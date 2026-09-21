---
layout: default
title: How much search went into a standard curve
---

# How much search went into a standard curve

The seed-provenance question is normally posed as *was this curve selected for
a hidden property?* The published parameters cannot answer that. A narrower
question can be answered, and it is the useful one:

> How much search would the **published** criteria have cost, and how much more
> could the era have afforded?

The gap between those two numbers is the room a hidden search could have
occupied. Everything below is output of `dissect-search_effort`.

## The floor is measured, not assumed

The simulated pools store the seed of every curve they kept, and they were
produced by walking seeds in order. So the gap between consecutive accepted
seeds *is* the number of seeds tried per curve kept. Scaling by the fraction
that reach cofactor 1 gives seeds per prime-order curve, which is what the
standards required.

| pool | curves kept | mean seed gap | cofactor 1 | seeds per prime-order curve |
|---|---|---|---|---|
| 192-bit | 18,836 | 265.4 | 44.3% | 2^9.2 |
| 224-bit | 22,211 | 225.1 | 42.6% | 2^9.0 |
| 256-bit | 18,502 | 270.2 | 44.1% | 2^9.3 |

About 600 seeds. On 1990s hardware that is a workstation afternoon. Nothing in
the published parameters requires more than this.

P-384 and P-521 have no simulated pool, so their floors are extrapolated
through `ln(p)`, calibrated against the 256-bit measurement rather than an
assumed constant: 2^9.9 and 2^10.3.

## The ceiling is modelled, and rests on one stated assumption

The binding cost is point counting, not hashing. Timing one SEA count per field
on the machine this was run on, using the standard's own `p` and `a = -3` with
a random `b`:

| curve | seconds per count | counts per CPU-year, today |
|---|---|---|
| P-192 | 1.44 | 2^24.4 |
| P-224 | 2.50 | 2^23.6 |
| P-256 | 5.05 | 2^22.6 |
| P-384 | 19.42 | 2^20.6 |
| P-521 | 125.29 | 2^17.9 |

Converting to the late 1990s needs a slowdown factor. The default is **2^9 to
2^13**, combining roughly 100-300x for single-thread integer performance with
another 3-30x for weaker period SEA implementations. That four-bit band is the
least defensible number on this page, which is why `--slowdown` exposes it. Move
it and the ceiling column moves with it; nothing else here depends on it.

## The bounds

At a generous 10,000 CPU-years, a dedicated campaign on one curve:

| curve | floor | 1997 ceiling | hiding room |
|---|---|---|---|
| P-192 | 2^9.2 | 2^24.7 - 2^28.7 | ~19 bits |
| P-224 | 2^9.0 | 2^23.9 - 2^27.9 | ~19 bits |
| P-256 | 2^9.3 | 2^22.9 - 2^26.9 | ~18 bits |
| P-384 | 2^9.9 | 2^20.9 - 2^24.9 | ~15 bits |
| P-521 | 2^10.3 | 2^18.2 - 2^22.2 | ~12 bits |

**Hiding room bounds what could have been done, not what was.** Nineteen bits
means a property occurring in roughly one curve in 500,000 was reachable for
P-192 under generous assumptions. It is not evidence that anything was
selected.

The structure worth noticing is that the ceiling falls with field size while the
floor barely moves. Point counting at 521 bits costs 87 times what it costs at
192 bits, while the acceptance rate drifts only from 2^9.2 to 2^10.3. If anyone
wanted to search widely, the small curves are where the room was.

Early abort, rejecting candidates by small-prime division polynomials before
paying for a full count, buys a constant factor of roughly 3-5x, about 2 bits.
It cannot do better here because the acceptance rate is already near 2^-9 and
the survivors still need the expensive count.

## What the twists do and do not show

The quadratic twist orders are mostly weak, and X9.62 never required otherwise.
All factorizations verify against their product.

| curve | twist factor sizes | twist security | curve security |
|---|---|---|---|
| P-192 | 5b · 94b · 95b | ~47 bits | 96 bits |
| P-224 | 2b² · 4b · 6b · 22b · 26b · 48b · 118b | ~59 bits | 112 bits |
| P-256 | 2b · 3b · 4b · 8b · 241b | ~120 bits | 128 bits |
| P-384 | prime | 192 bits | 192 bits |
| P-521 | 3b · 3b · 27b · 30b · 461b | ~230 bits | 260 bits |

**This is not evidence about search size, and an earlier draft of this analysis
wrongly said it was.** A large search for some hidden property X leaves the
twist distributed as it would be for any curve, because the twist order is
determined by the curve order and its factorization is uncorrelated with X.
Selecting hard on X produces exactly these weak twists.

What the twists show is narrower: twist security was not a selection criterion.
That is unsurprising for curves generated before the invalid-curve literature.
Biehl, Meyer and Müller published in 2000, key-validation work followed in 2003,
and twist security became a headline design goal with Curve25519 in 2006. These
curves predate all of it.

So the twists tell you what was on the checklist. They tell you nothing about
how many seeds were burned.

## The limit

The floor is a real lower bound. The ceiling is an order of magnitude with an
assumption attached. Between them sits a range that the parameters do not
narrow at all, because the seeds are unexplained 160-bit values with no
published derivation. That gap is the whole substance of the BADA55 argument,
and nothing in this analysis closes it.

## Reproducing

```shell
dissect-search_effort --bits 192 224 256
dissect-search_effort --bits 192 224 256 --measure          # times real point counts, slow
dissect-search_effort --bits 256 --measure --slowdown 64 256 # your own era assumption
```
