---
layout: default
title: Curve selection security review
description: >-
  Attack-oriented selection criteria, including the conditional one-third
  Brown-Gallant-Cheon profile and binary-field descent.
---

# Curve selection: attack-oriented review

## Conclusion first

The missing one-third criterion is real.  In the notation used here, let `r`
be the prime subgroup order.  If `r - 1` or `r + 1` has a suitable divisor `u`,
the Brown-Gallant-Cheon family can recover a repeatedly used static scalar with
roughly `u` oracle queries and `sqrt(r/u)` offline group operations.  A divisor
near `r^(1/3)` gives the familiar `O(r^(1/3))` balanced profile.

The existing `kn_factorization` trait does not test this.  It factors
`k * #E + 1` and `k * #E - 1`, where `#E` is the full curve cardinality, and
reports only the largest factor.  The relevant object is instead every divisor
of `r +/- 1`, especially a composite product near the attack optimum.

There is unusually strong documentary evidence for examining this property.
[RFC 5639, Section 2.1](https://www.rfc-editor.org/rfc/rfc5639.html#section-2.1)
explicitly says that the Brainpool subgroup orders may have a factor
`u` of order `q^(1/3)` to which the attack applies, then declines to make
resistance a curve requirement because the protocol precondition can usually
be avoided.  The measured parameters below confirm that property.

That is evidence of a deliberately omitted *selection safeguard*.  It is not,
without further provenance evidence, proof that the curves were selected
*because* they have the property.  The same RFC publishes a deterministic
seed-update procedure and says to accept a curve satisfying its stated
requirements; Cheon resistance is explicitly not one of them.  To establish
selection intent, the audit would need to replay that procedure, retain every
rejected candidate, and show an unexplained filter correlated with the
one-third score.

## Do not collapse three different “one-third” statements

| notation | attack family | applies to | status for ordinary ECDLP |
|---|---|---|---|
| `O(r^(1/3))` | Brown-Gallant-Cheon | Auxiliary powers, strong-DH assumptions, or an oracle applying one static scalar; divisor structure of `r +/- 1` matters | **Conditional**, not a generic ECDLP attack |
| `L_N(1/3, c)` | finite-field index calculus after MOV/Tate transfer, or after a successful descent | Small embedding degree or a concrete descent to a field/Jacobian where index calculus wins | Depends on the transfer; not implied by an extension degree alone |
| `2^(c sqrt(m log m))` | summation-polynomial work on curves over `F_(2^m)` | Binary-field ECDLP under algebraic-system and first-fall-degree heuristics | A research warning, not an established practical break of every binary curve |

The distinction changes the conclusion.  A normal ECDLP instance `(G, xG)`
does not give the auxiliary inputs needed by Cheon.  Conversely, a protocol
which publishes suitable powers of one secret or exposes chosen-input static
scalar multiplication cannot quote the ordinary Pollard-rho estimate as its
whole security argument.  The original
[static-DH analysis](https://eprint.iacr.org/2004/306) and the later
[elliptic-curve treatment](https://eprint.iacr.org/2010/177.pdf) make the oracle
precondition central.

## Measured subgroup-order profiles

The table uses the exact subgroup orders shipped by OpenSSL 3.0.13 and complete
factorizations of `r - 1` and `r + 1` verified with Sage.  The P-521 row is the
one stated exception: its complete `r - 1` factorization already demonstrates
the listed attack, while the available `r + 1` factorization is incomplete.
The row is therefore not a certificate that the listed profile is optimal.
The P-521 `r - 1` decomposition is also published in
[Appendix A.5 of ePrint 2017/055](https://eprint.iacr.org/2017/055.pdf); this
review uses its factor data, not the paper's parallel-work accounting.  For
every available divisor it computes

`C(u) = u + sqrt(r/u)`.

“Total” is `log2 C(u)` at the unconstrained minimum.  It is a conditional work
profile, not unconditional “security bits”; constants, memory, and the cost of
obtaining oracle answers are omitted.

| curve | generic rho | side | log2 queries | log2 offline work | log2 total | apparent reduction |
|---|---:|:---:|---:|---:|---:|---:|
| brainpoolP160r1/t1 | 79.933 | - | 52.616 | 53.625 | 54.207 | 25.726 |
| brainpoolP192r1/t1 | 95.804 | + | 62.898 | 64.355 | 64.803 | 31.001 |
| brainpoolP224r1/t1 | 111.877 | + | 73.736 | 75.008 | 75.508 | 36.368 |
| brainpoolP256r1/t1 | 127.705 | + | 81.750 | 86.829 | 86.871 | 40.833 |
| NIST P-192 | 96.000 | + | 76.627 | 57.686 | 76.627 | 19.373 |
| NIST P-224 | 112.000 | + | 47.278 | 88.361 | 88.361 | 23.639 |
| NIST P-256 | 128.000 | - | 84.664 | 85.668 | 86.252 | 41.748 |
| NIST P-384 | 192.000 | - | 103.347 | 140.327 | 140.327 | 51.673 |
| NIST P-521 [*] | 260.500 | - | 130.480 | 195.260 | 195.260 | 65.240 |
| secp160k1 | 80.000 | + | 47.331 | 56.335 | 56.337 | 23.663 |
| secp192k1 | 96.000 | + | 58.876 | 66.562 | 66.569 | 29.431 |
| secp224k1 | 112.000 | + | 78.170 | 72.915 | 78.207 | 33.793 |
| secp256k1 | 128.000 | - | 84.764 | 85.618 | 86.253 | 41.747 |

[*] P-521 uses the complete `r - 1` factorization only; searching a future
completion of `r + 1` can only improve, not invalidate, this demonstrated
profile.

The P-256 and secp256k1 rows are almost perfectly balanced around the
one-third exponent.  So are the Brainpool rows in the practical sense relevant
to this coarse asymptotic model.  P-224 and P-384 have valid but less balanced
profiles.  P-521 is **not immune** to the attack family: its best divisor from
the complete `r - 1` factorization has size about `2^130.480`, giving about
`2^195.260` offline work instead of the `2^260.500` generic baseline.  What is
true is narrower: this divisor is far from the balanced target
`r^(1/3) ~= 2^173.667`, and the attack needs an extraordinary number of oracle
answers.  The unresolved `r + 1` factorization also means that this audit cannot
certify P-521 as Cheon-resistant.

These attacks are not practical breaks of ordinary ECDSA or ephemeral ECDH,
nor of static ECDH that validates inputs and exposes only the output of a
one-way KDF.  They are serious parameter-selection facts for protocols whose
security model really supplies the auxiliary inputs.

For P-256, imposing a query cap changes the result as follows.  Each row picks
the best divisor that fits under the cap; the listed query count can be slightly
smaller than the cap because `u` must divide `r - 1` or `r + 1`.

| maximum queries | chosen log2 queries | log2 offline work | log2 total |
|---:|---:|---:|---:|
| `2^20` | 19.884 | 118.058 | 118.058 |
| `2^32` | 31.989 | 112.005 | 112.005 |
| `2^40` | 39.953 | 108.023 | 108.023 |
| `2^64` | 63.991 | 96.005 | 96.005 |

These remain enormous computations, but they show why a deployment-specific
query cap belongs in the criterion.  Calling the curve simply “86-bit” would be
wrong; calling its order fully Cheon-resistant would also be wrong.

This Cheon model is not an upper bound on every attack that can use a genuine
Diffie-Hellman oracle.  May and Schneider's 2023 implementation constructs a
smooth auxiliary curve over the subgroup-order field; with a precomputed
multi-terabyte codebook, it reports solving a P-256 discrete logarithm in about
30 seconds using fewer than 24,000 simulated DH-oracle calls.  That
[practical Maurer reduction](https://eprint.iacr.org/2023/539.pdf) is not driven
by factors of the base curve's `r +/- 1`, so it belongs in the deployment/oracle
model rather than the curve-characteristic score.  Its precondition is still a
raw DH oracle, not a lone ECDLP instance or a properly encapsulated KDF output.

### Matched-null check

The property is not unique to the published curve.  From DiSSECT's own
Brainpool simulations, a reproducible sample was selected by sorting curve
names by SHA-256:

| pool | sample / population | standard total | sample median | simulated scores no greater than standard |
|---|---:|---:|---:|---:|
| 160-bit | 200 / 3,184 | 54.2071 | 54.2291 | 7 / 200 |
| 192-bit | 25 / 2,640 | 64.8029 | 65.2085 | 9 / 25 |

At 160 bits the rank looks low, but the standard and median differ by only
0.022 bits because many curves saturate the continuous one-third optimum.  At
192 bits the standard is ordinary.  These bounded samples are enough to reject
the claim that merely observing an approximately one-third divisor is evidence
of intentional selection.  A full-family provenance replay is still required
before making an intent claim.

## Binary curves and ECC2K-130

ECC2K-130 is useful evidence, but for a different point.  It was a Certicom
challenge on the Koblitz curve over `F_(2^131)`, not a mainstream adopted
domain parameter.  The published solution used parallel Pollard rho with
negation and Frobenius orbits, estimating about `2^60.9` iterations; it did not
use a generic `O(r^(1/3))` algorithm.  See
[Breaking ECC2K-130](https://eprint.iacr.org/2009/541.pdf).

The binary-field screen in the earlier report also needed correction.  Every
`F_(2^m)` with `m > 1` contains `F_2`, even when `m` is prime.  Prime degree
means there is no *intermediate* subfield, not that there is no descent target.
The extended GHS analysis says that descent to `F_2` is infeasible for the
cryptographic prime degrees it studies; that attack conclusion cannot be
replaced by the syntactic test “m is prime.”  Composite degrees such as 155,
176, 185, 208, 272, 304, and 368 deserve priority because they expose more
descent choices, but composite degree alone is still not a vulnerability.

For each binary curve, a useful audit must record the chosen subfield, the
resulting Weil-restriction genus and splitting behavior, the best available
Jacobian/index-calculus cost, and whether the curve is in a special isogeny
class.  The primary references are the
[GHS construction](https://www.iacr.org/cryptodb/data/paper.php?pubkey=14115)
and its
[extended analysis](https://www.iacr.org/cryptodb/archive/2002/EUROCRYPT/2290/2290.pdf).
Semaev's
[summation-polynomial proposal](https://arxiv.org/abs/1504.01175) should be
tracked separately as heuristic; the follow-up
[limitations analysis](https://arxiv.org/abs/1503.08001) and
[last-fall-degree experiments](https://eprint.iacr.org/2015/573.pdf), which
raise doubt about the required first-fall-degree assumption, are part of the
same evidence package.

## Selection gate that should replace an outlier score

An LOF score answers whether measured traits are unusual.  It is not a security
acceptance test.  Curve selection should instead produce an attack ledger with
explicit preconditions and evidence:

| gate | required evidence | reject or flag condition |
|---|---|---|
| Parameter integrity | independently recompute field validity, nonsingularity, `#E`, trace, primality of `r`, cofactor, `G in E`, and `[r]G = O` | Any supplied/cached value disagrees; DiSSECT currently trusts several of these values |
| Generic group cost | largest prime factor of `#E`, rho cost after all usable automorphisms/endomorphisms | Below the target strength after speedups |
| Transfer attacks | embedding degree and an actual finite-field DLP estimate; anomalous and supersingular checks | Transfer cost below target, trace one, or supersingular outside an intentional pairing setting |
| Binary descent | all subfields including `F_2`, concrete GHS genus/splitting, best Jacobian and summation-polynomial estimates, with heuristic labels | A concrete transfer beats the target; never infer safety from prime `m` alone |
| Conditional one-third attacks | complete factorizations of `r +/- 1`, all composite divisors, query budget, protocol oracle/auxiliary-input model | Conditional total cost below target **and** the protocol supplies the precondition; otherwise retain as a documented hazard |
| Twist and invalid points | twist order, largest prime factor, twist embedding degree, input validation and formula completeness | Weak twist plus an attacker-reachable off-curve path, or inadequate subgroup handling |
| Provenance | executable generation transcript, seed/counter origin, every rejection reason, independent replay | Unexplained degrees of freedom or an undocumented rejection correlated with an attack metric |
| Deployment fit | static/ephemeral use, raw static-scalar or DH oracles, KDF behavior, attacker-chosen points, side channels, deprecation status | The implementation violates the assumptions under which the curve was accepted; include Maurer-style reductions when a raw DH oracle exists |

The pass/fail decision must be protocol-specific.  For the Cheon screen in
particular, first set a realistic maximum number of oracle answers and only
then minimize the offline cost over divisors within that budget.  Reporting the
unconstrained `r^(1/3)` optimum without its `r^(1/3)` query requirement is an
overstatement; omitting the criterion entirely is an understatement.

## Repository changes supporting this review

The new `cheon` trait:

- uses the prime subgroup order `r`, not the full cardinality;
- factors both `r - 1` and `r + 1` and searches composite divisors;
- reports queries, offline work, their sum, the rho baseline, and the apparent
  conditional reduction;
- accepts `max_query_bits = -1` for the unconstrained profile or a positive
  query budget for a deployment-specific result; and
- labels incomplete factorization or divisor search instead of silently
  treating partial data as a security verdict.

The `extension_degree` trait now includes degree 1 as the prime-field descent
target for every extension field.  Its description states that subfield
availability is only a screening condition, not a completed GHS analysis.
