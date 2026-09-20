#!/usr/bin/env python3
"""Survey the *form* of the adopted parameters of standard curves.

Everything else in this analysis compares a standard curve against a pool of
curves generated the same way. That asks whether a curve is unusual *given* its
generation procedure, and is silent on the procedure itself. This script looks
at the published parameters directly -- the field, the coefficients, the
generator, the cofactor, the seed -- and reports the conventions and outliers
visible in them.

Nothing here is a security finding. A convention shared by a family of curves
is evidence about who wrote the standard and when, not about whether the curves
are sound. Its use is the other direction: a feature that encodes a convention
will separate standard from simulated curves perfectly while saying nothing
about the mathematics, and `ml_distinguisher` has to exclude such features to
mean anything. This script is how they are found.
"""

import argparse
import collections
import json
import sys
import urllib.request

# X9.62 A.3.3 fixes the seed length at 160 bits. Seeds are stored as integers,
# so a seed whose leading bits are zero reads as shorter than it is; see
# `leading_zero_profile`.
X962_SEED_BITS = 160

DEFAULT_CATEGORIES = ("secg", "nist", "x962", "brainpool", "nums", "anssi", "oscca",
                      "other", "bn", "mnt", "bls", "wtls", "gost")


def to_int(value):
    """Database integers arrive as ints or as strings in some base."""
    if isinstance(value, str):
        return int(value, 0)
    return int(value)


def factorize(n):
    """Trial-division factorization; extension degrees here are all < 1024."""
    factors = []
    remaining, divisor = n, 2
    while divisor * divisor <= remaining:
        while remaining % divisor == 0:
            factors.append(divisor)
            remaining //= divisor
        divisor += 1
    if remaining > 1:
        factors.append(remaining)
    return factors


def subfield_index(degree):
    """The largest power of two dividing a composite extension degree, else None.

    Weil descent needs a proper subfield, so a curve over F_(2^m) with prime m
    has none to offer. The X9.62 `c2pnb*w1` curves are the exception, and the
    index that makes them exceptional turns out to also be the size of their
    cofactor -- see `degree_cofactor_pairs`.
    """
    factors = factorize(degree)
    if len(factors) < 2:
        return None
    index = 1
    for factor in factors:
        if factor == 2:
            index *= 2
    return index if index > 1 else None


def degree_cofactor_pairs(curves):
    """Composite-extension-degree binary curves, with the ratio cofactor/index.

    A ratio near 1 means the cofactor is the size of the subfield the curve is
    defined over, which is the signature of a curve built over that subfield
    rather than a coincidence of the order.
    """
    rows = []
    for curve in curves:
        if curve["field"]["type"] == "Prime":
            continue
        degree = curve["field"]["bits"]
        index = subfield_index(degree)
        if index is None:
            continue
        cofactor = int(curve["cofactor"])
        subfield_size = 2 ** index
        rows.append({
            "curve": curve["name"],
            "category": curve["category"],
            "degree": degree,
            "factors": factorize(degree),
            "subfield_index": index,
            "subfield_size": subfield_size,
            "cofactor": cofactor,
            "cofactor_over_subfield": cofactor / subfield_size,
        })
    return sorted(rows, key=lambda row: row["degree"])


def leading_zero_profile(seed_bits, nominal=X962_SEED_BITS):
    """Compare observed seed bit-lengths against uniform `nominal`-bit seeds.

    A uniform n-bit value has bit_length n with probability 1/2, n-1 with
    probability 1/4, and so on, so a seed set that looks like it contains
    under-length seeds is usually just leading zeros. Returns observed and
    expected counts of seeds at or above each threshold, which is the
    comparison that is not sensitive to the tail.
    """
    total = len(seed_bits)
    rows = []
    for threshold in range(nominal, nominal - 4, -1):
        observed = sum(1 for bits in seed_bits if bits >= threshold)
        # P(bit_length >= threshold) = 1 - 2^-(nominal - threshold + 1)
        expected = total * (1 - 2.0 ** -(nominal - threshold + 1))
        rows.append({"threshold": threshold, "observed": observed, "expected": expected})
    return rows


def generator_scale(curve):
    """Bit length of the generator's x-coordinate, or None if absent.

    Standards split by era rather than by mathematics: the older ones publish a
    full-width x, the newer ones take the smallest x that lands on the curve.
    """
    try:
        x = to_int(curve["generator"]["x"]["raw"])
    except (KeyError, TypeError, ValueError):
        return None
    return x.bit_length()


def seed_of(curve):
    seed = curve.get("simulation", {}).get("seed")
    if seed in (None, ""):
        return None
    try:
        return to_int(seed)
    except (TypeError, ValueError):
        return None


def degenerate_parameters(curve, embedding_bound=20):
    """Textbook-weak parameter shapes. Pairing curves trip the last one by design."""
    flags = []
    properties = curve.get("properties", {})
    trace = properties.get("trace")
    embedding = properties.get("embedding_degree")
    if trace is not None:
        if to_int(trace) == 1:
            flags.append("anomalous (trace 1, #E = p)")
        if to_int(trace) == 0:
            flags.append("supersingular (trace 0)")
    if embedding is not None and to_int(embedding) < embedding_bound:
        flags.append(f"embedding degree {to_int(embedding)}")
    return flags


def fetch(source, category):
    if source.startswith("mongodb"):
        import dissect.utils.database_handler as database
        return list(database.get_curves(database.connect(source), {"category": [category]}))
    url = f"{source}db/curves?category={category}"
    with urllib.request.urlopen(url) as handle:
        return json.loads(handle.read())["data"]


def load(source, categories):
    """One record per curve name; a curve listed in two categories is one curve."""
    by_name = {}
    for category in categories:
        try:
            records = fetch(source, category)
        except Exception as error:  # a category absent from the database is not fatal
            print(f"  skipping {category}: {type(error).__name__}: {error}", file=sys.stderr)
            continue
        for record in records:
            by_name.setdefault(record["name"], record)
    return list(by_name.values())


def report(curves, embedding_bound=20):
    print(f"curves surveyed: {len(curves)}")
    binary = [c for c in curves if c["field"]["type"] != "Prime"]
    prime = [c for c in curves if c["field"]["type"] == "Prime"]
    print(f"  prime field: {len(prime)}, binary field: {len(binary)}")

    print("\n=== composite extension degree, against the cofactor ===")
    rows = degree_cofactor_pairs(curves)
    if not rows:
        print("  none")
    for row in rows:
        factors = "x".join(str(f) for f in row["factors"])
        ratio = row["cofactor_over_subfield"]
        shown = f"{ratio:.4f}" if ratio < 1000 else f"{ratio:.3e}"
        print(f"  {row['curve']:16s} {row['category']:7s} m={row['degree']:4d} = {factors:14s}"
              f" subfield=2^{row['subfield_index']:<3d} cofactor={row['cofactor']:<8d}"
              f" cofactor/|subfield|={shown}")
    print(f"  ({len(binary) - len(rows)} of {len(binary)} binary curves have prime m)")

    print("\n=== seed coverage ===")
    seeded = [(c, seed_of(c)) for c in curves]
    seeded = [(c, s) for c, s in seeded if s is not None]
    print(f"  curves carrying a seed: {len(seeded)}/{len(curves)}")
    missing = sorted(c["name"] for c in curves if seed_of(c) is None)
    print(f"  curves with no seed: {len(missing)}")
    by_category = collections.Counter(c["category"] for c in curves if seed_of(c) is None)
    for category, count in sorted(by_category.items()):
        print(f"    {category:12s} {count}")

    print("\n=== seed lengths, against uniform 160-bit seeds ===")
    lengths = [s.bit_length() for _, s in seeded if s.bit_length() <= 2 * X962_SEED_BITS]
    outsized = [(c["name"], s.bit_length()) for c, s in seeded
                if s.bit_length() > 2 * X962_SEED_BITS]
    if lengths:
        print(f"  {'>= bits':>8} {'observed':>9} {'expected':>9}")
        for row in leading_zero_profile(lengths):
            print(f"  {row['threshold']:8d} {row['observed']:9d} {row['expected']:9.1f}")
        print("  a seed shorter than 160 bits is a 160-bit seed with leading zeros")
    for name, bits in outsized:
        print(f"  OUTSIZED: {name} carries a {bits}-bit value where a {X962_SEED_BITS}-bit "
              f"seed belongs ({(bits + 7) // 8} bytes)")

    print("\n=== generator x-coordinate ===")
    small = [(c["name"], generator_scale(c)) for c in curves
             if generator_scale(c) is not None and generator_scale(c) <= 32]
    print(f"  curves whose generator x fits in 32 bits: {len(small)}")
    for name, bits in sorted(small, key=lambda pair: pair[1]):
        print(f"    {name:38s} x is {bits}-bit")

    print("\n=== one prime, several standards ===")
    by_prime = collections.defaultdict(list)
    for curve in prime:
        by_prime[to_int(curve["field"]["p"])].append((curve["category"], curve["name"]))
    for value, entries in sorted(by_prime.items(), key=lambda kv: kv[0].bit_length()):
        categories = sorted({category for category, _ in entries})
        if len(categories) > 1:
            names = ", ".join(name for _, name in entries)
            print(f"  {value.bit_length():4d}-bit prime, {categories}: {names}")

    print("\n=== degenerate parameter shapes ===")
    flagged = [(c["name"], flags) for c in curves
               if (flags := degenerate_parameters(c, embedding_bound))]
    print(f"  flagged: {len(flagged)} (pairing curves trip the embedding-degree bound by design)")
    for name, flags in sorted(flagged):
        print(f"    {name:22s} {', '.join(flags)}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", default="https://dissect.crocs.fi.muni.cz/")
    parser.add_argument("--categories", nargs="*", default=list(DEFAULT_CATEGORIES))
    parser.add_argument("--embedding-bound", type=int, default=20,
                        help="flag curves whose embedding degree is below this")
    args = parser.parse_args()

    curves = load(args.source, args.categories)
    if not curves:
        sys.exit("no curves returned")
    report(curves, args.embedding_bound)


if __name__ == "__main__":
    main()
