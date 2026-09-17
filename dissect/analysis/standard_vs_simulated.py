#!/usr/bin/env python3
"""Compare standard curves against their simulated counterparts.

Builds feature vectors for a (category, bitlength) pair, drops features whose
trait results are not actually present on *both* sides, and reports where the
standard curves land in the simulated LOF distribution.

The coverage filter is the point of this script. Trait results in the public
database are unevenly populated: `multiples_x`, for instance, is computed for
the standard X9.62 curves but for none of the simulated ones, while
`x962_invariant` is the other way round. `feature_builder` imputes such gaps
with -1.0, which falls outside the [0, 1] range that every real feature is
scaled into, so an uncomputed trait reads as a maximally extreme value. Running
outlier detection over the raw feature vectors therefore "discovers" whichever
curves happen to have the most unusual trait coverage, not the most unusual
mathematics.
"""

import argparse
import json
import sys
import urllib.request

import numpy as np
import pandas as pd
from sklearn.neighbors import LocalOutlierFactor

import dissect.analysis.data_processing as dp
from dissect.traits import TRAITS

IMPUTED = -1.0


def matching_field_curves(source, category, sim_category, bits):
    """Standard curves in `category` sharing the prime field of the simulated pool.

    Several standards name the same curve: P-192, secp192r1 and ansix9p192r1 are
    one curve. More usefully, the Koblitz `k1` curves sit at the same bitlength as
    their `r1` siblings but over a different prime, so they have no comparable
    simulated pool and must be excluded rather than compared across fields.
    """
    def fetch(cat):
        url = f"{source}db/curves?category={cat}&bits={bits}"
        with urllib.request.urlopen(url) as handle:
            return json.loads(handle.read())["data"]

    pool = fetch(sim_category)
    if not pool:
        return [], None
    prime = pool[0]["field"].get("p")

    usable, skipped = [], []
    for curve in fetch(category):
        if curve["field"].get("p") == prime:
            usable.append(curve["name"])
        else:
            skipped.append(curve["name"])
    return usable, skipped


def build_features(source, category, bits, traits=None, sim_category=None):
    """Join per-trait results for `category` and the simulated pool into one frame."""
    categories = [category, sim_category or f"{category}_sim"]
    query = {"category": categories, "bits": [str(bits)], "cofactors": "all", "example": None}

    curves = dp.get_curves(source, query)
    if curves.empty:
        # An unmatched query yields a columnless frame, so select nothing from it.
        return pd.DataFrame(columns=["curve", "category"])
    curves = curves[["curve", "category"]]
    for name in traits or sorted(TRAITS):
        try:
            trait_df = dp.get_trait(source, name, query, False)
        except Exception as error:
            print(f"  skipping {name}: {type(error).__name__}: {error}", file=sys.stderr)
            continue
        if trait_df.empty:
            continue

        for feature in TRAITS[name].numeric_outputs():
            dp.clean_feature(trait_df, feature)

        flat = dp.flatten_trait(name, trait_df)
        features = [c for c in flat.columns if c != "curve"]
        for feature in features:
            dp.scale_feature(flat, feature)

        curves = curves.merge(flat, "left", on="curve")
        for feature in features:
            dp.impute_feature(curves, feature, IMPUTED)

    return curves


def usable_features(df, category, max_missing=0.05, sim_category=None):
    """Features with real data on both sides, and some variation to speak of."""
    simulated = df[df.category == (sim_category or f"{category}_sim")]
    standard = df[df.category == category]

    keep = []
    for feature in (c for c in df.columns if c not in ("curve", "category")):
        if (simulated[feature] == IMPUTED).mean() > max_missing:
            continue
        if (standard[feature] == IMPUTED).any():
            continue
        if df[feature].nunique() <= 1:
            continue
        keep.append(feature)
    return keep


def report(df, category, features):
    lof = LocalOutlierFactor()
    lof.fit_predict(df[features].values)
    scores = -lof.negative_outlier_factor_

    df = df.assign(lof=scores)
    df["rank"] = df["lof"].rank(ascending=False).astype(int)

    total = len(df)
    print(f"  curves: {total} ({(df.category == category).sum()} standard)")
    print(f"  features used: {len(features)}")
    print(f"  flagged as outliers (LOF > 1.5): {(scores > 1.5).sum()}")
    for _, row in df[df.category == category].sort_values("rank").iterrows():
        percentile = 100 * (1 - row["rank"] / total)
        print(f"    {row['curve']:20s} LOF={row['lof']:.3f}  rank {row['rank']}/{total}"
              f"  (more outlying than {percentile:.1f}% of simulated)")
    return df


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--category", default="brainpool", help="standard category, e.g. brainpool")
    parser.add_argument("--sim-category", default=None,
                        help="simulated pool to compare against (default: <category>_sim). "
                             "NIST and SECG prime curves are generated by the X9.62 method, "
                             "so x962_sim is their counterpart.")
    parser.add_argument("--bits", type=int, nargs="+", default=[256])
    parser.add_argument("--source", default="https://dissect.crocs.fi.muni.cz/")
    parser.add_argument("--max-missing", type=float, default=0.05,
                        help="drop a feature if more than this fraction of simulated curves lack it")
    parser.add_argument("--no-coverage-filter", action="store_true",
                        help="keep coverage-broken features, reproducing the artefact")
    args = parser.parse_args()

    for bits in args.bits:
        print(f"\n=== {args.category} {bits}-bit ===")
        sim_category = args.sim_category or f"{args.category}_sim"
        if args.sim_category:
            usable, skipped = matching_field_curves(args.source, args.category, sim_category, bits)
            if skipped:
                print(f"  excluded (different prime field, no comparable pool): {', '.join(skipped)}")
            if not usable:
                print("  no standard curves share the simulated pool's field")
                continue
        df = build_features(args.source, args.category, bits, sim_category=sim_category)
        if args.sim_category:
            df = df[df.curve.isin(usable) | (df.category == sim_category)]
        if (df.category == args.category).sum() == 0:
            print("  no standard curves at this bitlength")
            continue

        every = [c for c in df.columns if c not in ("curve", "category")]
        features = every if args.no_coverage_filter else usable_features(
            df, args.category, args.max_missing, sim_category)
        if not args.no_coverage_filter:
            print(f"  dropped {len(every) - len(features)}/{len(every)} features lacking coverage on both sides")
        report(df, args.category, features)


if __name__ == "__main__":
    main()
