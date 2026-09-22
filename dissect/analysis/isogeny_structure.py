#!/usr/bin/env python3
"""Isogeny structure of a standard curve, and how it compares to its pool.

Two quantities that matter for isogenies are not usable from the public
database, which is why they are computed here:

* `twist_order` stores the twist cardinality but its `factorization` field is
  `null` for the curves that matter, because factoring a 256-bit number
  exceeded the trait's 30-second budget.
* `class_number` returns "NO DATA (timed out)" on real curves, and even when it
  succeeds it reports bounds so wide (an upper bound of `sqrt(d) log(d) w / 2pi`
  against a lower bound carrying a 1/55 or 1/7000 factor) that two curves cannot
  be ranked against each other.

The class number is the size of the horizontal isogeny class, so its *relative*
value is the interesting one. `log_l_value` supplies that without computing the
class number at all: `h(D) = sqrt|D| L(1, chi_D) / pi` up to small factors, so
comparing `L(1, chi_D)` between curves over the same field compares their
isogeny class sizes directly. It is a proxy, and a truncated one.

One property here is forced rather than chosen, and saying so matters more than
measuring it. For any prime-order curve over an odd prime field the trace is
odd, so `D = t^2 - 4p` is 5 mod 8 and 2 is inert. Every such curve has no
rational 2-isogeny and odd order. Observing that on the NIST curves says
nothing about the NIST curves.
"""

import argparse
import json
import math
import statistics
import subprocess
import sys
import urllib.request

SMALL_PRIMES = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]


def to_int(value):
    """Hex with an 0x prefix, or plain decimal.

    Reading a bare decimal string as base 16 raises nothing, because decimal
    digits are all valid hex. It silently returns a different number, which is
    exactly how an earlier pass of this analysis went wrong.
    """
    if isinstance(value, int):
        return value
    text = str(value).strip()
    negative = text.startswith("-")
    if negative:
        text = text[1:]
    n = int(text, 16) if text.lower().startswith("0x") else int(text, 10)
    return -n if negative else n


def frobenius_discriminant(p, trace):
    """D = t^2 - 4p, negative for every ordinary curve."""
    return trace * trace - 4 * p


def two_is_forced_inert(p, trace):
    """True when D = 5 mod 8, which prime order over an odd prime field forces."""
    return frobenius_discriminant(p, trace) % 8 == 5


def legendre(D, l):
    """1 if l splits in the order, -1 if inert, 0 if ramified. Odd prime l."""
    if D % l == 0:
        return 0
    return 1 if pow(D % l, (l - 1) // 2, l) == 1 else -1


def kronecker_profile(D, primes=None):
    """Splitting behaviour at each small prime.

    For l not dividing the conductor this is the number of horizontal
    l-isogenies: 2 when it splits, 0 when inert, 1 when ramified.
    """
    return {l: legendre(D, l) for l in (primes or SMALL_PRIMES)}


def split_count(D, primes=None):
    return sum(1 for v in kronecker_profile(D, primes).values() if v == 1)


def conductor_small_part(D, primes=None):
    """{l: e} for each small l with l^e dividing D and e >= 2.

    l^2 | D means l divides the conductor, so the endomorphism ring is not
    maximal and a descending l-isogeny exists.
    """
    found = {}
    for l in (primes or SMALL_PRIMES):
        e, m = 0, abs(D)
        while m % l == 0:
            m //= l
            e += 1
        if e >= 2:
            found[l] = e
    return found


def log_l_value(D, bound=20000):
    """log of the truncated Euler product for L(1, chi_D).

    Proxy for class number relative to sqrt|D|, i.e. for the size of the
    horizontal isogeny class. Truncated, so only meaningful as a comparison
    between curves evaluated at the same bound.
    """
    total = 0.0
    sieve = _primes_upto(bound)
    for l in sieve:
        r = D % l
        if r == 0:
            continue
        chi = 1 if pow(r, (l - 1) // 2, l) == 1 else -1
        total -= math.log(1.0 - chi / l)
    return total


def _primes_upto(n):
    flags = bytearray([1]) * (n + 1)
    flags[0] = flags[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if flags[i]:
            flags[i * i::i] = bytearray(len(flags[i * i::i]))
    return [i for i in range(3, n + 1) if flags[i]]


def twist_order(p, trace):
    """#E' = p + 1 + t. Cross-checks against the database's twist_order trait."""
    return p + 1 + trace


def factor(n, timeout=900):
    """[(prime, exponent)] via PARI, or None if unavailable. Verified by product."""
    program = (f"f=factor({n}); for(i=1,matsize(f)[1], "
               f"print1(f[i,1],\":\",f[i,2],\" \"))")
    try:
        out = subprocess.run(["gp", "-q", "-f"], input=program, capture_output=True,
                             text=True, timeout=timeout).stdout.split()
    except (FileNotFoundError, subprocess.TimeoutExpired) as error:
        print(f"  factoring unavailable: {type(error).__name__}", file=sys.stderr)
        return None
    pairs = []
    for token in out:
        base, exponent = token.split(":")
        pairs.append((int(base), int(exponent)))
    product = 1
    for base, exponent in pairs:
        product *= base ** exponent
    if product != n:
        print(f"  factorization does not multiply back to n; discarding", file=sys.stderr)
        return None
    return pairs


def rho_security_bits(factors):
    """Half the bit length of the largest prime factor: Pollard rho on the
    largest subgroup, which is what an invalid-curve attack would face."""
    if not factors:
        return None
    return max(b for b, _ in factors).bit_length() // 2


def fetch(source, category, bits=None):
    url = f"{source}db/curves?category={category}"
    if bits:
        url += f"&bits={bits}"
    with urllib.request.urlopen(url) as handle:
        return json.loads(handle.read())["data"]


def report(source, category, sim_category, bitlengths, bound, do_factor):
    print("Forced, so not a distinguisher: prime order over an odd prime field makes")
    print("D = 5 mod 8, so 2 is inert and no rational 2-isogeny exists.\n")
    for bits in bitlengths:
        try:
            standard = [c for c in fetch(source, category, bits)
                        if c["field"]["type"] == "Prime"]
        except Exception as error:
            print(f"{bits}: {type(error).__name__}")
            continue
        if not standard:
            continue
        for curve in standard:
            p = to_int(curve["field"]["p"])
            trace = to_int(curve["properties"]["trace"])
            D = frobenius_discriminant(p, trace)
            print(f"=== {curve['name']} ({bits}-bit) ===")
            print(f"  2 inert (forced): {two_is_forced_inert(p, trace)}")
            cond = conductor_small_part(D)
            print(f"  conductor small part: "
                  f"{cond or 'none, endomorphism ring maximal below 50'}")
            print(f"  splits at {split_count(D)}/{len(SMALL_PRIMES)} of l = 3..47")

            try:
                pool = [c for c in fetch(source, sim_category, bits)
                        if to_int(c["field"]["p"]) == p]
            except Exception:
                pool = []
            if pool:
                obs = math.exp(log_l_value(D, bound))
                vals = []
                for s in pool[:1500]:
                    try:
                        st = to_int(s["properties"]["trace"])
                    except (KeyError, TypeError, ValueError):
                        continue
                    vals.append(math.exp(log_l_value(frobenius_discriminant(p, st), bound)))
                if vals:
                    pct = 100 * sum(1 for v in vals if v < obs) / len(vals)
                    print(f"  isogeny class size proxy L(1,chi) = {obs:.4f}; "
                          f"pool mean {statistics.mean(vals):.4f}, "
                          f"percentile {pct:.1f}% (n={len(vals)})")
            if do_factor:
                tw = twist_order(p, trace)
                facs = factor(tw)
                if facs:
                    shape = " * ".join(f"{b.bit_length()}b" + (f"^{e}" if e > 1 else "")
                                       for b, e in facs)
                    print(f"  twist #E' = {tw.bit_length()}b = {shape}")
                    print(f"  twist rho security ~{rho_security_bits(facs)} bits "
                          f"vs {p.bit_length()//2} on the curve")
            print()


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", default="https://dissect.crocs.fi.muni.cz/")
    parser.add_argument("--category", default="nist")
    parser.add_argument("--sim-category", default="x962_sim")
    parser.add_argument("--bits", type=int, nargs="+", default=[192, 224, 256])
    parser.add_argument("--l-bound", type=int, default=20000,
                        help="Euler product truncation; comparisons are only valid "
                             "between curves evaluated at the same bound")
    parser.add_argument("--factor-twist", action="store_true",
                        help="factor the twist order with PARI (slow, and the "
                             "database's own twist factorization is null because "
                             "it timed out)")
    args = parser.parse_args()
    report(args.source, args.category, args.sim_category, args.bits,
           args.l_bound, args.factor_twist)


if __name__ == "__main__":
    main()
