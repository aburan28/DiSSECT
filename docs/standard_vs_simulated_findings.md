# Standard vs. simulated curves: analysis notes

Run against the public database at `https://dissect.crocs.fi.muni.cz/`
(334,858 curves), comparing standard curves with their simulated counterparts
at matching bitlengths. Reproduce with `dissect-standard_vs_simulated`.

Categories with enough simulated curves to compare against are Brainpool
(160/192/224/256-bit, 1,677–3,184 simulated each) and X9.62 (192/256-bit,
~18,500 each). Other standards either have no simulated counterpart at a
matching bitlength or too few curves.

## Result: no standard curve is an outlier

Using Local Outlier Factor over the traits that are actually populated for both
standard and simulated curves:

| curve | bits | LOF | rank | more outlying than |
|---|---|---|---|---|
| brainpoolP160r1 | 160 | 1.008 | 1614/3186 | 49.3% |
| brainpoolP160t1 | 160 | 1.006 | 1722/3186 | 46.0% |
| brainpoolP192r1 | 192 | 1.000 | 1844/2642 | 30.2% |
| brainpoolP192t1 | 192 | 1.001 | 1770/2642 | 33.0% |
| brainpoolP224r1 | 224 | 0.985 | 2325/2363 | 1.6% |
| brainpoolP224t1 | 224 | 0.998 | 1780/2363 | 24.7% |
| brainpoolP256r1 | 256 | 1.014 | 532/1679 | 68.3% |
| brainpoolP256t1 | 256 | 1.013 | 546/1679 | 67.5% |
| ansix9p192r1 | 192 | 1.002 | 12711/18837 | 32.5% |
| ansix9p256r1 | 256 | 1.020 | 7956/18503 | 57.0% |

Not one standard curve crosses the conventional LOF threshold of 1.5 — in fact
no curve in any of these populations does, standard or simulated. The standard
curves sit in the middle of the distribution, and `brainpoolP224r1` is among the
most ordinary curves in its entire population.

A per-feature scan agrees. Across 1,295 (curve, feature) comparisons, a standard
curve landed in the 5% tail of the simulated distribution 49 times where chance
alone predicts about 130, and in the 1% tail 4 times against about 26 expected.
Standard curves are, if anything, *more* typical than a randomly drawn simulated
curve. This is a negative result, and it is consistent with what DiSSECT's
authors report.

## The trap: missing trait data reads as a maximally extreme value

The headline finding of this run is methodological.

`feature_builder` scales every real feature into `[0, 1]` and then imputes
missing trait results with `-1.0`. An *uncomputed* trait is therefore not
neutral — it is further from every real value than any two real values are from
each other. Trait coverage in the published database is uneven, so this
manufactures outliers out of nothing:

| trait | X9.62 standard | X9.62 simulated |
|---|---|---|
| `multiples_x` | 100% | **0%** |
| `x962_invariant` | **absent** | 100% |
| `small_prime_order` | 100% | 73% |
| `class_number` | **0%** | 80% |
| `discriminant` | **0%** | 81% |
| `square_4p1` | 100% | 83% |

`multiples_x` is computed for the standard X9.62 curves and for none of the
18,836 simulated ones; `x962_invariant` is the exact reverse. Each such trait is
a *perfect* standard-vs-simulated separator that encodes nothing about the
curves at all.

The effect is not subtle. Run over the unfiltered feature vectors, the standard
X9.62 curves come out as the single most extreme curve in their population:

```
ansix9p192r1  LOF=4.289  rank 1/18837
ansix9p256r1  LOF=2.896  rank 1/18503
```

Their distance to the 20th nearest neighbour is 10.4 and 12.4 respectively,
against a simulated median of 2.5 and 3.8 — close to the 12.6 maximum attainable
in a 159-dimensional unit cube. Decomposing that distance shows it is built
almost entirely from features where the curve has a real value and its
neighbours have `-1.0`, or vice versa. Restricting to traits populated on both
sides collapses both curves to LOF ≈ 1.0 and mid-pack ranks.

This is worth stating plainly: the most striking "discovery" available from this
dataset is an artefact, and it points in the direction most likely to be
believed — that the standardised curves are special.

`dissect-standard_vs_simulated` applies the coverage filter by default and
reproduces the artefact under `--no-coverage-filter`.

## `brainpool_overlap` measures a construction artefact, faithfully simulated

The `brainpool_overlap` trait checks whether the low bits of `a` coincide with
the high bits of `b`, a consequence of Brainpool deriving both coefficients from
one seed stream. An exact overlap (`o == 0`) holds for:

- `brainpoolP192r1` and `brainpoolP256r1` — but not the `t1` twists, whose
  `a = -3` is not seed-derived
- roughly half of the simulated curves at every bitlength (1329/2640 at 192-bit,
  1171/2361 at 224-bit, 833/1677 at 256-bit)

The trait computes `a`'s low `cut` bits against `b`'s top `cut + 1` bits, so
equality turns on one extra bit behaving as a coin flip — hence the ~50% rate.
The property is real and it is by construction, but the simulation reproduces it
exactly, so it carries no distinguishing power. The trait is undefined below
193-bit fields (`cut <= 0`), which is why 160-bit Brainpool returns nothing.

## Bugs found and fixed

**`convert_dtypes()` overflows on arbitrary-precision trait outputs.**
`get_trait` called `convert_dtypes()` on the whole frame; pandas tries to fit
integer columns into a C long, and trait outputs such as curve orders and
x-coordinates routinely exceed 64 bits. This raised
`OverflowError: Python int too large to convert to C long` and took out
`multiples_x`, `pow_distance`, `small_prime_order`, `x962_invariant`,
`trace_factorization` and `twist_order` — 6 of 21 traits — at 256-bit and above.
Conversion is now per-column, leaving oversized columns as objects; the
downstream feature cleaning already routes through `Decimal`.

**A bare `except:` hid it.** `feature_builder` wrapped the fetch in
`for i in range(3): try: ... except: print("Reconnecting...")`, so a
deterministic type error was reported as a network problem, retried three times,
and then crashed with `UnboundLocalError` on the unbound `trait_df`. It now
reports the actual exception and exits non-zero.

Both bugs are silent in the sense that matters: the first removes traits from
the analysis, and the second disguises why.

## Caveats

- Standard curves per bitlength number 1–2. These ranks describe where those
  specific curves fall; they are not a powered hypothesis test, and no single
  curve's percentile should be read as meaningful on its own.
- LOF over 100+ correlated, partly discrete features is a blunt instrument.
  The negative result is robust (nothing is anywhere near the threshold); a
  positive one would have needed far more scrutiny.
- Only Brainpool and X9.62 have simulated counterparts large enough to compare.
  NIST, SECG and the rest are untested here.
