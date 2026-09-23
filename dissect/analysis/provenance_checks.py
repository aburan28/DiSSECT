"""Three provenance checks on the NIST prime curves that the database does not run.

1. Twist embedding degree. MOV on the twist needs p^k = 1 mod r' for small k,
   r' the largest prime factor of the twist order. Showing k > limit settles
   feasibility without factoring r' - 1.
2. Seed independence. If the ten NIST seeds came from a counter or a shared
   prefix, some pair would have a small difference or a Hamming distance far
   below the Binomial(160, 1/2) mean of 80.
3. Trace sign. All five NIST prime curves have t > 0. Tested against an
   out-of-sample set of other standard curves, with aliases removed.
   Brainpool is left out by default: its order < p requirement forces t > 0.

Curve data comes from J08nY/std-curves; every (p, t, n) is checked by
p + 1 - t == n * h before use.
"""

import argparse
import itertools
import json
import math
import statistics
import urllib.request

from dissect.analysis.isogeny_structure import factor, to_int

STD_CURVES = "https://raw.githubusercontent.com/J08nY/std-curves/master/"
SEED_BITS = 160
# other names for NIST curves; counting them again would double-count P-192/P-256
ALIASES = {"prime192v1", "prime256v1", "secp192r1", "secp224r1", "secp256r1",
           "secp384r1", "secp521r1", "ansip192r1", "ansip224r1", "ansip256r1",
           "ansip384r1", "ansip521r1"}


def embedding_degree(p, r, limit=1000):
    """Smallest k <= limit with p^k = 1 mod r, or None if there is none."""
    x = 1
    for k in range(1, limit + 1):
        x = x * p % r
        if x == 1:
            return k
    return None


def hamming(x, y):
    return bin(x ^ y).count("1")


def shared_leading_bits(x, y, bits=SEED_BITS):
    return bits - (x ^ y).bit_length()


def binomial_cdf(k, n=SEED_BITS):
    """P(Hamming <= k) for independent uniform n-bit strings."""
    return sum(math.comb(n, i) for i in range(k + 1)) / 2 ** n


def seed_pair_statistics(seeds, bits=SEED_BITS):
    """Pairwise statistics over a {name: int} dict of seeds."""
    pairs = []
    for a, b in itertools.combinations(sorted(seeds), 2):
        x, y = seeds[a], seeds[b]
        diff = abs(x - y)
        pairs.append({"a": a, "b": b, "hamming": hamming(x, y),
                      "log2_diff": math.log2(diff) if diff else 0.0,
                      "shared_leading": shared_leading_bits(x, y, bits)})
    if not pairs:
        return {"pairs": [], "n_pairs": 0}
    ham = [q["hamming"] for q in pairs]
    lo = min(ham)
    sd = math.sqrt(bits) / 2
    return {
        "pairs": pairs,
        "n_pairs": len(pairs),
        "hamming_mean": statistics.mean(ham),
        "z_mean": (statistics.mean(ham) - bits / 2) / (sd / math.sqrt(len(ham))),
        "hamming_min": lo,
        # pairs are not independent; this union-style bound is approximate
        "p_min": 1 - (1 - binomial_cdf(lo, bits)) ** len(ham),
        "min_log2_diff": min(q["log2_diff"] for q in pairs),
    }


def sign_test_p(positive, total):
    """One-sided P(at least `positive` of `total` fair coin flips are heads)."""
    return sum(math.comb(total, i) for i in range(positive, total + 1)) / 2 ** total


def _get(source, category):
    with urllib.request.urlopen(f"{source}{category}/curves.json") as handle:
        return json.loads(handle.read())["curves"]


def _checked(curve):
    """(p, t) for a prime-field curve whose data is self-consistent, else None."""
    if curve.get("field", {}).get("type") != "Prime":
        return None
    try:
        p = to_int(curve["field"]["p"])
        t = to_int(curve["characteristics"]["trace_of_frobenius"])
        n, h = to_int(curve["order"]), to_int(curve["cofactor"])
    except (KeyError, TypeError, ValueError):
        return None
    return (p, t) if p + 1 - t == n * h else None


def _seed(curve):
    s = curve.get("characteristics", {}).get("seed") or curve.get("seed")
    if not s:
        return None
    s = str(s)
    return to_int(s if s.lower().startswith("0x") else "0x" + s)


def report(source, limit, other_categories):
    nist = _get(source, "nist")

    print(f"=== Twist embedding degree (bound {limit}) ===")
    for curve in sorted(nist, key=lambda c: c["field"]["bits"]):
        pt = _checked(curve)
        if not pt:
            continue
        p, t = pt
        facs = factor(p + 1 + t)
        if not facs:
            print(f"  {curve['name']:6s} twist not factored")
            continue
        r = max(b for b, _ in facs)
        k = embedding_degree(p, r, limit)
        print(f"  {curve['name']:6s} largest twist prime {r.bit_length():3d}b  "
              f"embedding degree {'= ' + str(k) if k else '> ' + str(limit)}")

    print("\n=== Seed independence ===")
    seeds = {c["name"]: s for c in nist if (s := _seed(c)) is not None}
    stats = seed_pair_statistics(seeds)
    print(f"  {len(seeds)} seeds, {stats['n_pairs']} pairs")
    if stats["n_pairs"]:
        print(f"  Hamming mean {stats['hamming_mean']:.1f} (expect 80), "
              f"z {stats['z_mean']:+.2f}")
        print(f"  Hamming min {stats['hamming_min']}, "
              f"P(some pair this low) ~ {stats['p_min']:.2f}")
        print(f"  smallest |difference| 2^{stats['min_log2_diff']:.1f} "
              f"(a counter walk would be tiny)")
        print(f"  max shared leading bits "
              f"{max(q['shared_leading'] for q in stats['pairs'])}")

    print("\n=== Trace sign ===")
    nist_traces = [pt[1] for c in nist if (pt := _checked(c))]
    pos = sum(1 for t in nist_traces if t > 0)
    print(f"  NIST: {pos}/{len(nist_traces)} positive, "
          f"P = {sign_test_p(pos, len(nist_traces)):.3f} (chosen after looking)")
    seen, other = set(), []
    for category in other_categories:
        before = len(other)
        try:
            curves = _get(source, category)
        except Exception as error:
            print(f"  {category}: {type(error).__name__}")
            continue
        for c in curves:
            pt = _checked(c)
            if not pt or c["name"] in ALIASES or pt in seen:
                continue
            seen.add(pt)
            other.append(pt[1])
        mine = other[before:]
        print(f"  {category}: {sum(1 for t in mine if t > 0)}/{len(mine)} positive")
    if other:
        opos = sum(1 for t in other if t > 0)
        print(f"  out of sample: {opos}/{len(other)} positive, "
              f"P = {sign_test_p(opos, len(other)):.3f}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", default=STD_CURVES)
    parser.add_argument("--limit", type=int, default=1000,
                        help="largest embedding degree searched on the twist")
    parser.add_argument("--other", nargs="+", default=["secg", "x962"],
                        help="std-curves categories used as the out-of-sample set; "
                             "brainpool is excluded by default because RFC 5639 "
                             "requires order < p, which forces t > 0")
    args = parser.parse_args()
    report(args.source, args.limit, args.other)


if __name__ == "__main__":
    main()
