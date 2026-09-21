#!/usr/bin/env python3
"""Bound how many candidate curves the generator of a standard curve examined.

The seed-provenance question is usually posed as "was this curve selected for a
hidden property?", which the published parameters cannot answer. A narrower
question can be answered: *how much search would the visible criteria have
cost, and how much more could the era have afforded?* The gap between those two
is the room a hidden search could have occupied.

The floor is measured, not modelled. The simulated pools store the seed of every
curve they kept, and the pools were produced by walking seeds in order, so the
gap between consecutive accepted seeds is the number of seeds tried per curve
kept. Dividing by the fraction that reach cofactor 1 gives the seeds per
prime-order curve, which is what the standards actually required.

The ceiling is modelled and rests on an assumption stated in the open: the cost
of one point count on period hardware. This module measures that cost on the
machine it runs on and scales it back by a caller-supplied factor. Change the
factor and the ceiling moves; nothing else here depends on it.

What this does NOT support is an inference from the twist. A large search for
some hidden property X leaves the twist order distributed as it would be for any
curve, because the twist is determined by the curve order and its factorization
is uncorrelated with X. Weak twists therefore say that twist security was not a
selection criterion, which for curves generated before the invalid-curve
literature is unsurprising. They say nothing about the size of the search.
"""

import argparse
import json
import math
import statistics
import subprocess
import sys
import urllib.request

SECONDS_PER_YEAR = 365.25 * 24 * 3600

# Single-thread integer performance from the late 1990s to now spans roughly
# 100-300x, and period SEA implementations were weaker than a modern library by
# perhaps another 3-30x. The product is wide on purpose; it is the least
# defensible number here and the caller should be able to move it.
DEFAULT_SLOWDOWN = (2 ** 9, 2 ** 13)


def to_int(value):
    """Database values are hex with an 0x prefix OR plain decimal.

    Assuming base 16 for a bare decimal string does not raise, because decimal
    digits are all valid hex. It silently returns a different number.
    """
    if isinstance(value, int):
        return value
    text = str(value).strip()
    negative = text.startswith("-")
    if negative:
        text = text[1:]
    n = int(text, 16) if text.lower().startswith("0x") else int(text, 10)
    return -n if negative else n


def seed_gaps(seeds):
    """Differences between consecutive accepted seeds, sorted ascending."""
    ordered = sorted(seeds)
    return [b - a for a, b in zip(ordered, ordered[1:])]


def floor_from_gaps(gaps, kept, prime_order):
    """Seeds tried per prime-order curve.

    `gaps` measures seeds per curve *kept*; the pools keep curves of several
    cofactors, so scale by how many of those reach cofactor 1.
    """
    if not gaps or prime_order == 0:
        return None
    return statistics.mean(gaps) * kept / prime_order


def extrapolated_floor(bits, reference_bits, reference_trials):
    """Floor at a bitlength with no simulated pool.

    The acceptance rate is dominated by the chance that a number near p is
    prime, so trials scale with ln(p). Calibrated against a measured pool
    rather than an assumed constant.
    """
    scale = math.log(2 ** reference_bits)
    return reference_trials * math.log(2 ** bits) / scale


def ceiling(seconds_per_count, slowdown, cpu_years):
    """Candidates testable in `cpu_years` at `slowdown` times today's cost."""
    per_cpu_year_now = SECONDS_PER_YEAR / seconds_per_count
    return per_cpu_year_now / slowdown * cpu_years


def hiding_room_bits(floor_trials, ceiling_trials):
    """Negative log2 of the rarest property a search could have selected for.

    This bounds what *could* have been done. It is not evidence that anything
    was.
    """
    if floor_trials <= 0 or ceiling_trials <= 0:
        return None
    return math.log2(ceiling_trials) - math.log2(floor_trials)


def fetch(source, category, bits):
    url = f"{source}db/curves?category={category}&bits={bits}"
    with urllib.request.urlopen(url) as handle:
        return json.loads(handle.read())["data"]


def measure_point_count(p, a, reps=3, timeout=3600):
    """Median wall-clock seconds for one SEA point count over this field.

    Uses the real candidate shape: the standard's own p and a, a random b.
    Returns None if PARI is unavailable.
    """
    program = f"""
p={p}; a={a}; setrand(1);
for(i=1,{reps},
  b=random(p); E=ellinit([a,b],p);
  t0=getabstime(); ellcard(E); t1=getabstime(); print(t1-t0)
);
"""
    try:
        out = subprocess.run(["gp", "-q", "-f", "--stacksize", "800000000"],
                             input=program, capture_output=True, text=True,
                             timeout=timeout).stdout.split()
    except (FileNotFoundError, subprocess.TimeoutExpired) as error:
        print(f"  point count unavailable: {type(error).__name__}", file=sys.stderr)
        return None
    times = [float(x) / 1000.0 for x in out if x.strip()]
    return statistics.median(times) if times else None


def report(source, category, sim_category, bitlengths, cpu_years, slowdown, measure):
    print(f"floor: seeds per prime-order curve, measured from gaps between accepted seeds")
    print(f"ceiling: candidates testable in {cpu_years:g} CPU-years at "
          f"{slowdown[0]}-{slowdown[1]}x today's point-count cost\n")
    print(f"{'bits':>5} {'pool':>7} {'mean gap':>9} {'cof-1':>6} {'floor':>8} "
          f"{'s/count':>8} {'ceiling':>17} {'hiding room':>12}")

    measured = {}
    for bits in bitlengths:
        try:
            pool = fetch(source, sim_category, bits)
        except Exception as error:
            print(f"{bits:5d}  no pool: {type(error).__name__}")
            continue
        if not pool:
            print(f"{bits:5d}  no pool")
            continue
        seeds = [to_int(c["simulation"]["seed"]) for c in pool
                 if c.get("simulation", {}).get("seed")]
        gaps = seed_gaps(seeds)
        n1 = sum(1 for c in pool if int(c["cofactor"]) == 1)
        fl = floor_from_gaps(gaps, len(pool), n1)
        if fl is None:
            print(f"{bits:5d}  no usable seeds")
            continue
        measured[bits] = fl

        secs = None
        if measure:
            p = to_int(pool[0]["field"]["p"])
            secs = measure_point_count(p, p - 3)
        if secs:
            lo = ceiling(secs, slowdown[1], cpu_years)
            hi = ceiling(secs, slowdown[0], cpu_years)
            room = hiding_room_bits(fl, hi)
            cs = f"2^{math.log2(lo):.1f}-2^{math.log2(hi):.1f}"
            print(f"{bits:5d} {len(pool):7d} {statistics.mean(gaps):9.1f} "
                  f"{100*n1/len(pool):5.1f}% 2^{math.log2(fl):6.1f} {secs:8.2f} "
                  f"{cs:>17} {room:11.0f}b")
        else:
            print(f"{bits:5d} {len(pool):7d} {statistics.mean(gaps):9.1f} "
                  f"{100*n1/len(pool):5.1f}% 2^{math.log2(fl):6.1f} "
                  f"{'--':>8} {'(use --measure)':>17} {'--':>12}")
    return measured


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", default="https://dissect.crocs.fi.muni.cz/")
    parser.add_argument("--category", default="nist")
    parser.add_argument("--sim-category", default="x962_sim")
    parser.add_argument("--bits", type=int, nargs="+", default=[192, 224, 256])
    parser.add_argument("--cpu-years", type=float, default=1e4,
                        help="size of the hypothetical period campaign")
    parser.add_argument("--slowdown", type=float, nargs=2, default=list(DEFAULT_SLOWDOWN),
                        help="how much slower one point count was then than now; "
                             "this is the least defensible input and moves the ceiling")
    parser.add_argument("--measure", action="store_true",
                        help="time a real point count per field (needs PARI; slow)")
    args = parser.parse_args()

    report(args.source, args.category, args.sim_category, args.bits,
           args.cpu_years, tuple(args.slowdown), args.measure)


if __name__ == "__main__":
    main()
