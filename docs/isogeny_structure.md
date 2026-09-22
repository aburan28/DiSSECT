---
layout: default
title: Isogeny structure of the NIST prime curves
---

# Isogeny structure of the NIST prime curves

Two isogeny quantities that matter are not usable from the public database, so
`dissect-isogeny_structure` computes them:

- **`twist_order`** stores the twist cardinality, but its `factorization` field
  is `null` for the curves of interest. Factoring a 256-bit number exceeded the
  trait's 30-second budget.
- **`class_number`** returns `NO DATA (timed out)` on real curves. Even when it
  succeeds, it reports an upper bound of `sqrt(d) log(d) w / 2pi` against a
  lower bound carrying a 1/55 or 1/7000 factor. Those are too far apart to rank
  two curves against each other.

Everything below is output of the tool. The twist cardinality it derives is
cross-checked against the database's own `twist_order` trait, and a test pins
that equality.

## One property is forced, so it distinguishes nothing

For any prime-order curve over an odd prime field the trace must be odd, so
`D = t^2 - 4p` is always 5 mod 8. That makes 2 inert: **no rational 2-isogeny
exists and the order is odd**, on every such curve ever generated.

All five NIST primes show it. So would any curve meeting the same requirement.
It is a consequence of demanding prime order, not a choice anyone made, and
reporting it as a property of the NIST curves would be misleading.

## Endomorphism rings, and what they do to the volcano

`l^2 | D` means `l` divides the conductor, the endomorphism ring is not
maximal, and a descending `l`-isogeny exists.

| curve | conductor small part | volcano shape |
|---|---|---|
| P-192 | none below 50 | bare crater, no descending isogeny |
| P-224 | 3^3 | depth-1 volcano at l = 3 |
| P-256 | none below 50 | bare crater |

This agrees with the stored `conductor` field where the database has one, 1 for
P-192 and P-256 and 3 for P-224. A non-maximal order is not rare: 18 to 35
percent of the simulated pool has one.

## Relative size of the horizontal isogeny class

The class number is how many curves share an endomorphism ring, and
`h(D) = sqrt|D| L(1, chi_D) / pi` up to small factors. Comparing `L(1, chi_D)`
between curves over the same field therefore compares isogeny class sizes
directly, without computing a class number at all.

| curve | L(1,chi) | pool mean | percentile |
|---|---|---|---|
| P-192 | 0.617 | 0.676 | 41 |
| P-224 | 0.539 | 0.884 | 6.8 |
| P-256 | 1.084 | 0.901 | 76 |

P-224's isogeny class is somewhat smaller than typical. Nothing here is
anomalous.

**A statistic that overstated it, recorded because the correction is the
point.** Counting how many primes under 50 split put P-224 at the 0.2nd
percentile, and it held there across four different prime ranges, which looked
robust. The L-value is the quantity the class number actually depends on, and
it moves P-224 to the 6.8th percentile. A crude count over a handful of tiny
primes was not measuring what it appeared to measure.

The proxy is truncated at a bound the caller sets, so comparisons are only
valid between curves evaluated at the same bound.

## Twist structure

Every factorization below is verified against its product before being
reported; the tool discards any that does not multiply back.

| curve | twist #E' | twist rho security | curve security |
|---|---|---|---|
| P-192 | 5b · 94b · 95b | ~47 bits | 96 bits |
| P-224 | 2b² · 4b · 6b · 22b · 26b · 48b · 118b | ~59 bits | 112 bits |
| P-256 | 2b · 3b · 4b · 8b · 241b | ~120 bits | 128 bits |
| P-384 | prime | 192 bits | 192 bits |
| P-521 | 3b · 3b · 27b · 30b · 461b | ~230 bits | 260 bits |

An implementation that skips point validation faces the twist, not the curve.
P-192 offers roughly 47 bits there.

This says nothing about how hard anyone searched. A search for a hidden
property leaves the twist distributed as it would be for any curve, since the
twist order is fixed by the curve order and its factorization is uncorrelated
with that property. See
[search effort]({{ '/search_effort.html' | relative_url }}) for that argument
and the retraction behind it.

## The headline is shared, not special

All five have `|D| ~ 4p`, so class numbers near 2^112 to 2^260. The isogeny
classes are astronomically large, which is precisely the anti-CM property these
curves were built for. The simulated curves have it too, by construction, which
is why it cannot separate them.

P-384 and P-521 have no simulated pool, so their percentiles cannot be computed
at all.

## Reproducing

```shell
dissect-isogeny_structure --bits 192 224 256
dissect-isogeny_structure --bits 224 256 --factor-twist   # needs PARI, slow
dissect-isogeny_structure --bits 256 --l-bound 50000      # tighter truncation
```
