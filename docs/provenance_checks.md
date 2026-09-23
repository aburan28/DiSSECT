---
layout: default
title: Provenance checks on the NIST prime curves
---

# Provenance checks on the NIST prime curves

`dissect-provenance_checks` runs three checks the database does not. Curve data
comes from [J08nY/std-curves](https://github.com/J08nY/std-curves). Before any
curve is used, its `(p, t, n, h)` must satisfy `p + 1 - t = n·h`, so a trace
parsed in the wrong base cannot slip through. All numbers below are output of
the tool.

## 1. Twist embedding degree: MOV on the twist is infeasible

An invalid-point attack moves the discrete log to the quadratic twist. MOV
would reduce it further, to a finite field, if `p^k ≡ 1 mod r'` held for a
small `k`, where `r'` is the largest prime factor of the twist order. The tool
factors the twist with PARI and looks for such a `k` up to 1000. It does not
need to factor `r' - 1`, which is infeasible at these sizes.

| curve | largest twist prime | embedding degree |
|---|---|---|
| P-192 | 95 bits | > 1000 |
| P-224 | 118 bits | > 1000 |
| P-256 | 241 bits | > 1000 |
| P-384 | 385 bits | > 1000 |
| P-521 | 461 bits | > 1000 |

A pairing with embedding degree above 1000 lands in a field of more than
10^5 bits, so none of the twists is vulnerable to MOV. The weak twists of
P-192 and P-224 are weak against Pollard rho on the small factor (see the
[isogeny page]({{ '/isogeny_structure.html' | relative_url }})), not against
pairings.

## 2. The ten NIST seeds are independent

Suppose the seeds had come from a counter, or from values sharing a prefix.
Then some pair would have a small difference, or a Hamming distance far below
the Binomial(160, 1/2) mean of 80 (standard deviation 6.3).

| statistic | observed | expected if independent |
|---|---|---|
| pairs | 45 | |
| mean Hamming distance | 81.2 (z = +1.27) | 80 |
| minimum Hamming distance | 67 | P(some pair this low) ≈ 0.66 |
| smallest difference | 2^153.6 | near 2^160; a counter would give a few bits |
| most shared leading bits | 6 | about log2(45) ≈ 5.5 |

Nothing links the seeds. This does not show how the seeds were picked. It
shows only that they were not picked by a simple walk from one to the next.

## 3. All five NIST traces are positive, and that is chance

All five NIST prime curves have `t > 0`, so the order is below `p`. For fair
coin flips this happens with probability 1/32 = 0.031. The hypothesis was
chosen after looking at the data, though, so it needs an out-of-sample test.
That test uses other standard curves with NIST aliases removed (`prime192v1`,
`prime256v1`, `secp256r1`, and so on) and deduplicated by `(p, t)`:

| set | positive |
|---|---|
| SECG, excluding NIST aliases | 4/10 |
| X9.62, excluding NIST aliases | 4/5 |
| **out of sample** | **8/15, P = 0.50** |

Brainpool is excluded by default. RFC 5639 requires the group order to be
below `p`, which forces `t > 0` on all seven Brainpool curves, so including
them would count a rule as evidence. With Brainpool included the tool reports
15/22 (P = 0.067), and that figure is misleading for the same reason.

The sign pattern does not persist outside the five curves where it was found,
so it is not a signature of the NIST procedure.

## Reproduce

```
pip install -e .          # PARI/gp must be on PATH for the twist factorizations
dissect-provenance_checks
dissect-provenance_checks --other secg x962 brainpool   # shows the forced 7/7
```
