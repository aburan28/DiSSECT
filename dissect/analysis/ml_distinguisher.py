#!/usr/bin/env python3
"""Supervised distinguisher: can a classifier tell standard curves from simulated ones?

Trains a model to separate the standard curves of one family from their
simulated pools, scores every standard curve out-of-fold (each is left out of
training once), and reports the leave-one-out AUC. On its own that number is
meaningless with a handful of positives, so the same procedure is run on a
null: a random subset of simulated curves is relabelled as "standard", the
true standards are dropped, and the AUC distribution that produces is the
reference. A real difference shows as an observed AUC the null rarely reaches.

Two things a distinguisher will find that are not properties of the curves:

* Generation conventions of the simulator. Every `x962_sim` curve has
  b < p/2: the X9.62 procedure fixes b only up to sign (b^2 r = a^3 mod p has
  two roots), and the simulator keeps the smaller one where the standards kept
  whichever the derivation produced. A tree model separated secp128r1 and
  secp224r1 perfectly on `weierstrass_b` alone. Such features are dropped by
  prefix (`--drop-prefix`).
* Duplicates across the positive class. Brainpool `t1` is the `r1` curve at
  the same bitlength, so leaving one out while its twin trains gives AUC 0.96
  by memorisation. Twins are dropped (`--drop-suffix`).
"""

import argparse
import sys

import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

from dissect.analysis.standard_vs_simulated import (build_features, matching_field_curves,
                                                    usable_features)

MODELS = {
    "logreg": lambda: LogisticRegression(C=0.1, max_iter=2000, class_weight="balanced"),
    "gbdt": lambda: GradientBoostingClassifier(n_estimators=100, max_depth=2, subsample=0.5,
                                               random_state=0),
}


def out_of_fold(X, y, model_fn, positives, folds=5, seed=0):
    """P(standard) for every row, never predicted by a model that saw it.

    Each positive is left out on its own; negatives are scored in `folds` folds
    with all positives in training.
    """
    n = len(y)
    scores = np.full(n, np.nan)
    for p in positives:
        train = np.ones(n, bool)
        train[p] = False
        scores[p] = model_fn().fit(X[train], y[train]).predict_proba(X[[p]])[0, 1]
    negatives = np.where(y == 0)[0]
    fold = np.random.default_rng(seed).permutation(len(negatives)) % folds
    for k in range(folds):
        test = negatives[fold == k]
        train = np.ones(n, bool)
        train[test] = False
        scores[test] = model_fn().fit(X[train], y[train]).predict_proba(X[test])[:, 1]
    return scores


def loo_auc(X, y, model_fn):
    positives = np.where(y == 1)[0]
    scores = out_of_fold(X, y, model_fn, positives)
    percentiles = [(scores[y == 0] < scores[p]).mean() for p in positives]
    return roc_auc_score(y, scores), scores, percentiles


def null_aucs(X_neg, n_positive, model_fn, reps, seed=1):
    """AUCs when `n_positive` random simulated curves play the standard ones."""
    rng = np.random.default_rng(seed)
    aucs = []
    for _ in range(reps):
        y = np.zeros(len(X_neg), int)
        y[rng.choice(len(X_neg), n_positive, replace=False)] = 1
        aucs.append(loo_auc(X_neg, y, model_fn)[0])
    return np.array(aucs)


def assemble(frames, drop_prefix, drop_suffix):
    """Stack per-bitlength frames on their shared features."""
    common = set.intersection(*[set(feats) for _, feats in frames])
    common = sorted(c for c in common if not any(c.startswith(p) for p in drop_prefix))
    df = pd.concat([df[["curve", "category"] + common] for df, _ in frames], ignore_index=True)
    for suffix in drop_suffix:
        df = df[~df.curve.str.endswith(suffix)]
    return df.reset_index(drop=True), common


def run(df, features, category, model_name, null_reps, subsample, seed=2):
    if subsample and len(df) > subsample:
        keep = (df.category == category) | (np.random.default_rng(seed).random(len(df)) < subsample / len(df))
        df = df[keep].reset_index(drop=True)
    X = df[features].values.astype(float)
    y = (df.category == category).astype(int).values
    if y.sum() < 2:
        print(f"  {model_name}: need at least two standard curves, have {y.sum()}")
        return
    print(f"\n  model={model_name}: {len(df)} curves, {y.sum()} standard, {len(features)} features")
    auc, scores, pct = loo_auc(X, y, MODELS[model_name])
    for p, q in zip(np.where(y == 1)[0], pct):
        print(f"    {df.curve[p]:22s} P(standard)={scores[p]:.3f}  above {100 * q:.1f}% of simulated")
    print(f"    leave-one-out AUC = {auc:.3f}")
    if null_reps:
        null = null_aucs(X[y == 0], int(y.sum()), MODELS[model_name], null_reps)
        print(f"    null ({null_reps} relabelled draws): AUC {null.mean():.3f} +/- {null.std():.3f}, "
              f"P(null >= observed) = {(null >= auc).mean():.2f}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--category", default="secg")
    parser.add_argument("--sim-category", default="x962_sim")
    parser.add_argument("--bits", type=int, nargs="+", default=[128, 160, 192, 224, 256])
    parser.add_argument("--source", default="https://dissect.crocs.fi.muni.cz/")
    parser.add_argument("--model", choices=[*MODELS, "both"], default="both")
    parser.add_argument("--null-reps", type=int, default=30)
    parser.add_argument("--subsample", type=int, default=20000,
                        help="cap on simulated curves per family, for speed (0 = all)")
    parser.add_argument("--drop-prefix", nargs="*", default=["weierstrass_b"],
                        help="features encoding a simulator convention rather than the curve")
    parser.add_argument("--drop-suffix", nargs="*", default=["t1"],
                        help="curve-name suffixes marking duplicates of another standard curve")
    args = parser.parse_args()

    frames = []
    for bits in args.bits:
        usable, skipped = matching_field_curves(args.source, args.category, args.sim_category, bits)
        if skipped:
            print(f"  {bits}: excluded (different field): {', '.join(skipped)}")
        if not usable:
            print(f"  {bits}: no standard curve shares the pool's field")
            continue
        df = build_features(args.source, args.category, bits, sim_category=args.sim_category,
                            standard_curves=usable)
        feats = usable_features(df, args.category, sim_category=args.sim_category)
        frames.append((df, feats))
        print(f"  {bits}: {len(df)} curves, {len(feats)} usable features", flush=True)
    if not frames:
        sys.exit("nothing to compare")

    df, features = assemble(frames, args.drop_prefix, args.drop_suffix)
    print(f"\n=== {args.category} vs {args.sim_category}: {len(features)} shared features ===")
    for name in (MODELS if args.model == "both" else [args.model]):
        run(df, features, args.category, name, args.null_reps, args.subsample)


if __name__ == "__main__":
    main()
