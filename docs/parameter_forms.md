---
layout: default
title: The form of the adopted parameters
---

# The form of the adopted parameters

Every other page here compares a standard curve against a pool of curves
generated the same way. That asks whether a curve is unusual *given* its
generation procedure, and is silent about the procedure itself. This page looks
at the published parameters directly: the field, the coefficients, the
generator, the cofactor, the seed.

None of it is a security finding, and one of the most interesting entries is a
bug in the source dataset rather than anything about a curve. The reason to record it
is narrower and mostly methodological. A feature that encodes a *convention*
will separate standard from simulated curves perfectly while saying nothing
about the mathematics, and the
[supervised distinguisher]({{ '/ml_distinguisher_methodology.html' | relative_url }})
is worthless until such features are excluded. This page is how they get found.

Everything below is output of `dissect-parameter_forms` against the live
database — 202 distinct standard curves, 139 over prime fields and 63 over
binary fields.

## The five `c2pnb` curves encode the same 16 twice

Of the 63 binary-field standard curves, 57 have a prime extension degree. The
six that do not are five X9.62 curves and one pairing curve. In every one of the
five, the degree is 16 times an odd prime **and** the cofactor is within 0.7% of
2^16, the size of the subfield that composite degree exposes:

| curve | m | subfield | cofactor | cofactor / \|subfield\| |
|---|---|---|---|---|
| c2pnb176w1 | 16 × 11 | 2^16 | 65390 | 0.9978 |
| c2pnb208w1 | 16 × 13 | 2^16 | 65096 | 0.9933 |
| c2pnb272w1 | 16 × 17 | 2^16 | 65286 | 0.9962 |
| c2pnb304w1 | 16 × 19 | 2^16 | 65070 | 0.9929 |
| c2pnb368w1 | 16 × 23 | 2^16 | 65392 | 0.9978 |

Every other X9.62 binary curve has prime m and a cofactor between 2 and 10, the
two exceptions being `c2tnb359v1` at 76 and `c2tnb431r1` at 10080 — both with
prime m, so neither has a subfield to expose.

The 16 in the extension degree and the 16 in the cofactor are the same 16: these
curves were built over F(2^16) and the cofactor is the fingerprint of that
construction rather than an accident of the group order. This is the one place
in the standard set where the field structure and the group structure visibly
share a parameter, and a composite extension degree with a small-index subfield
is precisely the structure Weil descent targets. The curves are long deprecated
and the construction is documented; what is worth recording is that two
properties usually treated as independent are here the same fact.

## Most standard curves are not verifiably random at all

Only **57 of 202** curves carry a seed. The other 145 were not produced by a
public seed-and-hash procedure, so the entire seed-provenance question — the one
that motivates this analysis — simply does not apply to them:

| category | curves with no seed |
|---|---|
| nums | 24 |
| other | 21 |
| bn | 16 |
| secg | 16 |
| mnt | 14 |
| brainpool | 14 |
| x962 | 11 |
| gost | 9 |
| wtls | 8 |
| bls | 6 |
| nist | 5 |
| oscca | 1 |

Brainpool is the entry that deserves a footnote. Its selling point is seeds
derived from the digits of pi and e, and the database carries none of them for
its standard curves, only for its simulated ones. That is a gap in the record,
not in the standard. `secp256k1`, every Koblitz curve, SM2 and Ed448 are in this
list for the ordinary reason that they were constructed rather than searched.

## A seed-length alarm that checked out clean

Seeds are stored as integers, and read that way some are far under the 160-bit
length X9.62 mandates — `c2tnb239v3` reads as 128 bits, `secp128r1` as 148. That
is leading zeros, not short seeds. A uniform 160-bit value has a bit length of
160 only half the time:

| bit length at least | observed | expected |
|---|---|---|
| 160 | 24 | 28.0 |
| 159 | 37 | 42.0 |
| 158 | 49 | 49.0 |
| 157 | 50 | 52.5 |

The fit is close enough that no seed needs another explanation. Recorded because
it is exactly the shape of finding that looks alarming when tabulated and
evaporates when the null is written down.

## Every seed in the database derives its curve, except one

A seed is not decoration: a verifiably random curve commits to its coefficients
through *r*·*b*² ≡ *a*³ (mod *p*), where *r* comes from the seed by SHA-1 alone
(X9.62 A.3.3.1). That check is worth running over the whole corpus, and
`dissect-parameter_forms` now does.

**26 of the 27 prime-field seeds verify.** The single failure is `FRP256v1`,
the ANSSI curve.

That number is also the control. A verifier that reproduces P-192, P-224,
P-256, P-384, secp256r1 and twenty-one others from their published seeds is
working, so the one curve it rejects is being rejected on the evidence. (The
step that is easy to get wrong is clearing the leading bit of W0. Omit it and
*nothing* verifies — which is how the first version of this check failed, and
why the 26 matter as much as the 1.)

### What `FRP256v1` actually carries

Where a 160-bit seed belongs, the record has **1029 bytes** — 2058 hex
characters, against 40 for every other curve. Beside its siblings in the same
record it is fifty times the size of anything else:

| field | length (hex chars) |
|---|---|
| cm_disc | 79 |
| discriminant | 77 |
| embedding_degree | 77 |
| j_invariant | 77 |
| trace_of_frobenius | 39 |
| **seed** | **2058** |

Length alone would prove nothing, since X9.62 permits any seed of at least 160
bits. It is the verification that settles it, and truncating the value to its
first 20, first 32 or last 20 bytes does not verify either. The content is
high-entropy (7.82 bits per byte), is not DER, and contains none of the curve's
own *p*, *a*, *b* or order, so it is not another field misfiled into this one.

### Where the fault is, and where it is not

**Not in DiSSECT.** Its importer reads the field verbatim at
`database_handler.py:69-72`, taking the `characteristics.seed` branch because
the top-level `seed` is null for this curve, and that parsing is correct.

The value is in **`anssi/curves.json` in the `std-curves` dataset**, which
DiSSECT imports. The record's own cited source is the
[JORF notice](https://www.legifrance.gouv.fr/jorf/id/JORFTEXT000024668816),
which publishes the parameters and no seed — FRP256v1 arriving without a seed
or any selection justification is a documented criticism of the curve, so an
absent seed is the expected state. The fix belongs upstream: drop the field.

An earlier version of this page attributed the fault to DiSSECT's parsing. That
was wrong, and tracing it to the source dataset is what corrected it.

## Generator conventions split by era, not by mathematics

Eighteen curves publish a generator whose x-coordinate fits in 32 bits, against
a full-width coordinate everywhere else. They are not scattered: they are the
NUMS family, the GOST family, and a handful of others, taking the smallest x
that lands on the curve — 0, 1, 2, 3, 5, 8, 13, 18, 32.

This is cosmetic, and that is the point. Any trait reading the generator
separates curves by publication date and design philosophy rather than by
mathematics, so it belongs on the exclusion list of any distinguisher.

## One prime, several standards

Nine primes appear in more than one category. Most are the familiar aliasing —
P-256, secp256r1 and ansix9p256r1 are one curve under three names. Two are not:

- The **112-bit and 160-bit SECG primes** reappear in WTLS curves.
- **2^512 − 569** is used by both `numsp512d1`/`numsp512t1` and the GOST
  parameter sets `id-tc26-gost-3410-12-512-paramSetA` and
  `id-tc26-gost-3410-2012-512-paramSetC`.

The GOST case reads as convergence rather than copying: both rules select the
largest prime of that shape below 2^512, and there is only one answer.

## Nothing accidentally weak

Nineteen curves have an embedding degree below 20, the lowest being `mnt1` at 3.
Every one is an intentional pairing curve from the BN, BLS, MNT or `alnr`
families, where a small embedding degree is the entire design goal. No standard
curve in the database is anomalous (trace 1) or supersingular (trace 0).

## Reproducing

```shell
dissect-parameter_forms
dissect-parameter_forms --categories secg nist x962 --embedding-bound 40
```

The seed check is part of the default report: it recomputes X9.62 A.3.3.1 for
every prime-field curve carrying a seed and reports which ones fail to derive
their own `b`.
