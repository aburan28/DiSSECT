---
layout: default
title: Detailed arithmetic and cover traits
---

# Detailed arithmetic and cover traits

Nine traits extend the original `canonical_geometry` and
`hyperelliptic_branch` traits. They compute arithmetic, moduli and divisor
geometry, automorphism structure, and local ramification for smooth
higher-genus curves and their explicit covers over finite fields.

## Run a complete report

```shell
dissect-higher-genus-report -i examples/detailed_covers.json -o report.json

# Equivalent module invocation in an uninstalled checkout:
sage -python -m dissect.utils.higher_genus_report \
  -i examples/detailed_covers.json -o report.json

# Select traits and adjust their limits:
dissect-higher-genus-report -i examples/detailed_covers.json \
  -t frobenius_polynomial hasse_witt cover_ramification \
  --params '{"timeout":60,"max_rational_fibers":128}'

# Individual traits use the existing JSON and MongoDB interfaces:
dissect-compute-json -t automorphism_structure \
  -i examples/higher_genus_curves.json \
  --params '{"extension_degree":2,"max_candidates":200000}'
dissect-compute-db -t cover_ramification --category detailed_cover_examples \
  --params '{"timeout":60}'
```

The combined report constructs each model once and shares successful
function-field and L-polynomial computations between traits. The individual
report includes the input model, field, effective parameters and Sage/Python
versions for reproducibility. The individual
trait registry now contains 44 traits. The nine detailed traits expose 132
top-level output fields, including status and scope fields; matrix entries,
local fibers, invariant coordinates and place records are structured data.
Use `TRAITS[name].OUTPUT` to inspect the full schema and descriptions.

| Trait | Main computations | Scope recorded explicitly |
| --- | --- | --- |
| `frobenius_polynomial` | L(T), characteristic polynomial, trace, curve counts and Jacobian orders over extensions | Curve order and Jacobian order are separate; optional metadata is checked |
| `jacobian_decomposition` | Honda-Tate simple factors, Schur indices, multiplicities, local invariants and extension decompositions | Base-field simplicity differs from geometric simplicity |
| `hasse_witt` | Cartier and Serre-dual Hasse-Witt matrices, semilinear ranks, p-rank and a-number | Matrix basis and extension-field coordinates are included |
| `newton_polygon` | Exact normalized slopes, vertices, multiplicities and ordinary/supersingular indicators | Supersingular is not superspecial |
| `weierstrass_gaps` | Generic gaps, all closed Weierstrass places, gap excess and numerical semigroups | Nonclassical generic gaps and output truncation are explicit |
| `automorphism_structure` | Finite-field group, order distribution, canonical actions and tame cyclic quotient genera | Full groups, cover stabilizers and proper search subgroups are distinguished |
| `moduli_invariants` | Igusa-Clebsch, Shioda and Dixmier-Ohno generators; invariant subsets in other cases | Normalization, characteristic restrictions and classification scope are explicit |
| `gonality_and_quotients` | Geometric gonality, base-field bounds, realizing projections and supplied quotient maps | Supplied maps do not exhaust all possible quotients |
| `cover_ramification` | All branch fibers, e/f degrees, different exponents, wild breaks, splitting fields and map checks | Cover degree differs from gonality; geometric and arithmetic Galois properties differ |

## Models and their normalization

The existing [model guide](higher_genus_traits.md) describes `Hyperelliptic`
and `SmoothPlane`. Geometric genus is computed or justified from the validated
model. The inputs never supply an unchecked family label.

- `Hyperelliptic`: `y^2+h(x)y=f(x)` with genus at least two, including
  characteristic two. Sage checks the smooth hyperelliptic model.
- `SmoothPlane`: a geometrically smooth homogeneous plane equation of degree
  d at least four; genus is `(d-1)(d-2)/2` and the curve is non-hyperelliptic.
- `CyclicTrigonal`: the smooth projective normalization of `y^3=f(x)`, with
  characteristic different from three and squarefree f of degree d at least
  four. The affine closure may have singularities at infinity; computations
  use the normalized function field, not its arithmetic plane genus.

For example:

```json
{
  "name": "nonplane-genus-four",
  "category": "detailed_cover_examples",
  "form": "CyclicTrigonal",
  "field": {"type": "Prime", "p": "0x7", "bits": 3},
  "params": {"f": [0, -1, 0, 0, 0, 1]},
  "example": true,
  "standard": false
}
```

This gives `y^3=x^5-x`, of genus four. A degree-three map in genus at least
three certifies non-hyperellipticity: an independent degree-two and
degree-three pair would bound the genus by two by Castelnuovo-Severi.
The cover is geometrically cyclic; it is Galois over the base field only when
that field contains the third roots of unity.

The trigonal genus is
`d-1-(gcd(3,d)-1)/2`. Singular or inseparable input forms, incorrect declared
genera and unmatched subgroup-order/cofactor metadata are rejected. Subgroup
metadata refers to the Jacobian and is never used as a curve point count.

General complete intersections and noncyclic trigonal equations need further
input adapters. This implementation does not claim to represent every
non-hyperelliptic curve or every possible cover.

## Status, bounds and exact serialization

Every detailed result has `status` and `reason`:

| Status | Meaning |
| --- | --- |
| `ok` | The declared computation finished within its stated mathematical scope |
| `partial` | A valid subset, subgroup or truncated record list is returned with a reason |
| `unsupported` | A required mathematical formula or Sage backend is unavailable |
| `resource_limit` | A configured size or search bound prevents the computation |
| `timeout` | The wall-clock alarm interrupted this trait |
| `not_applicable` | The model does not belong to a supported higher-genus family |

An uncomputed numeric field is null, never zero or minus one. A bounded
extension search does not silently turn an undecided property into a negative
answer. `ok` does not imply all conceivable geometric questions have been
decided: for example, a completed finite splitting search still leaves
`geometrically_simple` null when it has no witness.

Common defaults are 30 seconds, base-field size at most 4096 and genus at most
12. Exact point counting additionally requires `q^g <= 2,000,000`. Extension
counts are capped at 64; the decomposition extension bound is capped at 16.
Other defaults are 100,000 automorphism candidates, 128 transformation
records, 256 Weierstrass records, 16 rational place samples, 64 rational target
fibers, 16 declared maps and map degree/input exponent at most 256. Inspect each trait's `DEFAULT_PARAMS` and override
parameters with `--params`. Overrides replace a parameter choice rather than
creating duplicate jobs from the original parameter grid.

The wall-clock alarm is cancelled on every return path. Successful shared
computations are cached; failures are not. Some generic SageMath 9.8
function-field algorithms over nonprime constant fields are unavailable with
proof enabled; these return `unsupported` rather than globally disabling
Sage's proof settings. Hyperelliptic point counting and odd-characteristic
Cartier matrices have direct exact implementations that avoid that limitation.

Prime-field elements are JSON integers. Extension-field elements are lists
of prime-field coefficients in ascending powers of the defining generator,
padded to the absolute field degree. Automorphism results include the search
field modulus and the source-field generator's image in that basis. Rational slopes and local invariants are exact
`[numerator,denominator]` pairs. Lists and dictionaries survive the existing
MongoDB integer encoding and decoding, including negative polynomial
coefficients and null values.

## Frobenius and Jacobian decomposition

For genus g and field cardinality q, write

```text
Z(C,T) = L(T)/((1-T)(1-qT))
L(T) = product(1-alpha_i*T), i=1,...,2g
P(X) = X^(2g)*L(1/X)
#C(F_(q^n)) = q^n+1-sum(alpha_i^n)
#J(F_(q^n)) = product(1-alpha_i^n).
```

Coefficients are always ascending in the JSON output. Newton identities
compute extension characteristic polynomials with integral divisions checked.
The reciprocal symmetry `L[2g-i]=q^(g-i)*L[i]` is verified. Hyperelliptic
counts enumerate fibers over constant extensions, including the infinity
fiber. In characteristic two they use the absolute-trace criterion for the
Artin-Schreier quadratic. Other models use Sage's global function-field
L-polynomial algorithm. No elliptic trace recurrence is applied to genus g>1.

`metadata_order_matches` compares a supplied Jacobian subgroup order times
cofactor against `L(1)`; it is null when metadata is absent. A disagreement is
reported rather than changing the computed polynomial.

Factoring P over the rationals gives isotypic blocks, not automatically the
dimensions of simple factors. For each irreducible factor f, the implementation
constructs `Q(alpha)` and computes Honda-Tate local invariants. At a prime above
p, the invariant modulo one is

```text
v_P(alpha) * residue_degree(P) / [F_q:F_p].
```

Real places contribute denominator two. The least common multiple of all
denominators is the Schur index e. The corresponding simple factor has
dimension `degree(f)*e/2`; if f occurs m times, its multiplicity is `m/e`.
Integrality and divisibility are checked. This handles real linear factors and
quadratic Weil factors with Schur index two, both of which defeat the shortcut
“factor degree divided by two.”

Extension decompositions repeat this calculation over the reported field.
A split or repeated simple factor provides a witness against geometric
simplicity. A search without a witness leaves geometric simplicity undecided;
it is not a proof of absolute simplicity. An isogeny factor is not itself an
explicit morphism from the curve to a smaller curve.

## Cartier, Hasse-Witt and Newton data

The reported Cartier matrix C uses column images of a regular differential
basis. Over `F_(p^e)` the Cartier operator is inverse-Frobenius semilinear:

```text
Cartier(v) = C * v^(1/p)
C_(n+1) = C * C_n^(1/p)
H = transpose(C^p).
```

H is the Frobenius-semilinear matrix on Serre-dual `H^1(O)`. The p-rank is the
stable rank after g Cartier iterates, while the a-number is `g-rank(C)`.
These quantities can differ. For `y^2=x^5+1` over F_3 the first rank is one,
the second is zero, p-rank is zero and a-number is one.

For odd-characteristic hyperelliptic curves set `F=f+h^2/4`. If `c_k` are the
coefficients of `F^((p-1)/2)`, then in basis
`x^(j-1) dx/(y+h/2)` the Cartier entries are
`C[i,j]=c_(p*i-j)^(1/p)`, using one-based indices. The coefficient expansion
has an explicit term bound. Other supported cases use Cartier on Sage's
function-field regular differential space.

Newton slopes come from the lower convex hull of `(i,v_p(L_i))`, omitting
zero coefficients, then divide by e to normalize relative to q. Their exact
multiplicities and symmetry under `s -> 1-s` are reported. Ordinary means
p-rank g; supersingular means every slope is one half. Superspecial is detected
by a zero Cartier matrix and is a stronger property than supersingularity.

## Weierstrass places and numerical semigroups

The global function field supplies its actual generic gap sequence and the
complete set of closed Weierstrass places, including nonrational ones. Each
record contains its closed-point degree, gap sequence, gap excess over the
generic sequence, conductor and minimal numerical-semigroup generators.

Classical gap weight is `sum(gaps)-g*(g+1)/2`. Generic gap excess instead
subtracts the actual generic gap sum. The two need not coincide in positive
characteristic. These gap-excess sums are not advertised as a general formula
for a wild Wronskian divisor.

For example, the Hermitian quartic `x^4+y^4+z^4=0` in characteristic three has
generic gaps `[1,2,4]` and 28 geometric Weierstrass points. Replacing its generic
gaps by `[1,2,3]` or forcing a characteristic-zero total weight would be wrong.
The code's regression tests include this case.

The output cap truncates place records, not the computed locus totals.
`returned_locus_complete` records this difference. Additional rational-place
samples illustrate ordinary gaps but do not constitute an exhaustive rational
point listing.

## Automorphisms and quotient genera

The search field is `F_(q^r)`, where r is `extension_degree`. Its order and
modulus are reported. Counts in this field do not certify the automorphism
group over the algebraic closure.

For odd-characteristic hyperelliptic curves, every automorphism descends to
PGL2 because the degree-two pencil is unique. The algorithm enumerates PGL2
and tests the homogeneous branch form of degree `2g+2`. A transformation
preserving that form up to multiplier lambda lifts over the search field only
if lambda is a square. Both lifts are checked through their faithful action on
regular differentials. Branch stabilizers without a rational lift are counted
separately.

For characteristic-two curves with constant nonzero h, every automorphism
preserves the unique branch point at infinity and acts affinely on x. Polynomial
solutions of `s^2+s=A` produce all lifts. A regular permutation action certifies
the full group because the deck involution acts trivially on canonical
differentials in characteristic two. For nonconstant h, this backend certifies
the deck involution and labels the rest of the group search incomplete.

Smooth plane curves use projective linear transformations. Exhaustive PGL3
search is used when within budget. When a rational frame exists, every
automorphism must map that fixed frame into the complete rational point set;
enumerating such images can prove completeness much more efficiently.
Oversized searches return the exact monomial stabilizer as a subgroup with
`base_field_order=null`, not an asserted full group order.

Cyclic trigonal models enumerate the PGL2 stabilizer of the chosen cubic
pencil and all rational Kummer lifts, including the two possible nontrivial
deck characters. For genus greater than four the trigonal pencil is unique by
Castelnuovo-Severi, so a completed stabilizer search gives the full finite-field
automorphism group. In genera three and four it is reported as a complete
cover stabilizer, with full curve-group completeness left undecided.

Group results include a GAP structure description, center order, element-order
distribution and explicit transformations. For cyclic subgroups of order prime
to p, the invariant regular-differential dimension gives the quotient genus.
That formula is not applied to wild subgroups. The hyperelliptic deck quotient
has genus zero independently, including in characteristic two.

## Moduli invariants and normalization

The invariant values are exact field elements with their coefficient degrees
as weights. An explicitly chosen nonzero coordinate determines weight-zero
ratios; each ratio records both powers. These ratios can lose information and
are not by themselves a universal isomorphism test. Geometric moduli also do
not distinguish every twist over the base field.

- Genus two, characteristic greater than five: Sage's Igusa-Clebsch
  `(I2,I4,I6,I10)` on `f+h^2/4`.
- Hyperelliptic genus three, characteristic greater than seven: nine Shioda
  generators `J2,...,J10`, using the normalized factorial transvectant. The
  returned coordinates are not integrally rescaled.
- Plane quartics, characteristic greater than three: thirteen Dixmier-Ohno
  generators, constructed by binary restrictions, contravariants, Hessians,
  differential contractions and adjugates. The exact construction is in
  `curve_invariants.py`; it uses the displayed Girard-Kohel differential
  formulas, with `I27=Res(Fx,Fy,Fz)/2^54`. Coordinates must not be compared with
  another library without matching its normalization. Completeness is asserted
  only in characteristic greater than seven.
- Other odd-characteristic hyperelliptic cases: a binary discriminant and,
  when factorial denominators are invertible, a transvectant subset.
- Other plane degrees: a gradient-resultant discriminant when the
  characteristic does not divide the degree. This is an invariant subset,
  not a complete moduli system.
- Cyclic trigonal inputs: invariants of the unordered branch divisor of the
  chosen cover. The scope is that cover's branch configuration, not a complete
  invariant system for the underlying curve.

Characteristic-two hyperelliptic models with constant nonzero h use an exact
Artin-Schreier affine orbit instead of binary invariants: remove polynomial
`s^2+s` coboundaries, reduce constants by absolute trace, and minimize over all
base-field affine coordinate changes. This classifies base-field isomorphism
classes in the polynomial single-branch-point model class, including rational
twists. Its completeness is recorded in `classification_complete_in_stated_scope`;
weighted classical completeness is inapplicable. Nonconstant h still needs
another moduli backend and returns `unsupported`. Subsets and incomplete classical systems
return `partial`. Invariance tests include non-unimodular GL2/GL3 coordinate
changes, coefficient scaling and the Fermat quartic discriminant normalization.

## Gonality and explicit maps

Hyperelliptic gonality is two; cyclic trigonal gonality is three. A smooth
plane curve of degree d has geometric gonality d-1. When a base-field rational
plane point exists, the implementation constructs a separating projection
through it and checks the resulting function-field genus. Its equation and
coordinate-change matrix certify a base-field map of degree d-1. Without a
realizing map, base-field gonality remains bounded rather than asserted equal
to geometric gonality.

Declared maps live in a top-level `analysis.maps` list. The coefficient schema
is structured data; expressions are never evaluated as Python or Sage code.
Each source rational function uses two powers, referring to the source's
**reported function-field base and extension generators**, which need not be
the original plane coordinates. Each map result repeats its source
function-field equation and generator names, independently of any alternate
natural projection selected in the cover trait. For a numerator and denominator:

```json
{
  "numerator": [{"coeff": 1, "powers": [2, 0]}],
  "denominator": [{"coeff": 1, "powers": [0, 0]}]
}
```

This represents the square of the base generator. Omitting the denominator
uses one. A scalar is a constant field element. Missing `target` means P^1.

The full example `examples/detailed_covers.json` includes
`y^2=x^6+1 -> v^2=u^3+1` with `u=x^2,v=y`, an elliptic quotient of degree two:

```json
{
  "name": "elliptic-quotient",
  "x": {"numerator": [{"coeff": 1, "powers": [2, 0]}]},
  "y": {"numerator": [{"coeff": 1, "powers": [0, 1]}]},
  "target": {"form": "Hyperelliptic", "params": {"f": [1, 0, 0, 1]}}
}
```

Target hyperelliptic equations may have genus one even though source
`Hyperelliptic` input requires genus at least two. `SmoothPlane` targets use
their affine chart `F(u,v,1)=0`; `CyclicTrigonal` targets use `v^3=f(u)`.
Every target model is validated over the source constant field.

The images must satisfy the target relation. Map degree is the source pole
degree of u divided by the target pole degree of its separating coordinate.
Constants, vanishing denominators and inconsistent relations return
`invalid_map` per map. For inseparable maps, simultaneous p-th root extraction
of the generator images determines the inseparable degree. No separable
different formula is assigned to such a map.

For a verified separable map, the different is computed by

```text
div(phi^*(dt)) - phi^*(div(dt)).
```

Pullbacks of target prime ideals distinguish places above the same rational
coordinate. All local exponents must be nonnegative, all fibers satisfy their
degree identity, and the total different must satisfy Riemann-Hurwitz. This
works for rational, elliptic and higher-genus targets. Quotient discovery is
not exhaustive: the trait verifies supplied maps rather than claiming to have
found every morphism.

## Full local cover ramification

The natural function-field projection and verified separable declared maps
report every branch fiber. For each source place P above Q the record includes

- `degree(P)` and `degree(Q)` over the constant field;
- ramification index e and residue degree f;
- different exponent d, including wild contributions;
- whether p divides e, and the excess `d-(e-1)` over the tame exponent.

Each fiber checks `sum(e*f)=degree(map)`. The native different checks
`2g-2=-2*degree(map)+degree(Different)`. General target maps check
`2g_X-2=degree(map)*(2g_Y-2)+degree(Different)`.

For double covers, the lower ramification break at a ramified place is `d-1`.
For `y^2+y=x^5+x^3+1` over F_2 this gives one geometric branch point at
infinity, e=2, different exponent six and lower break five. Its different
degree is six; its number of distinct branch points is one. The old tame
`hyperelliptic_branch` trait remains inapplicable there, while the detailed
cover trait computes the wild data.

The complete infinity fiber is always included. Rational target fibers are
sampled up to the configured bound; only a complete set of q+1 fibers yields
`rational_curve_points_from_fibers`. The least common multiples of closed
branch and ramification-point degrees give the minimal constant extensions
making those respective point sets rational.

Separable double covers are recognized as Galois with group C2. Tame Kummer
equations `v^n=a(u)` have geometric monodromy Cn. Their base-field deck order
is `gcd(n,q-1)`; the minimal constant extension making the cover Galois has
degree the multiplicative order of q modulo n. For example, a quartic Kummer
projection over F_7 has geometric group C4 but rational deck order two. General
monodromy groups and Galois closures are not inferred from a branch partition.

## Validation and interpretation

Run `sage -python -m pytest -q` for the full suite. The regression tests cover
independent model point counts, Newton identities, Schur-index counterexamples,
semilinear rank drops, nonprime fields, nonclassical gaps, tame and wild
ramification, explicit positive-genus target maps, inseparability, group
completeness, coordinate invariance, output truncation, timeouts, CLI behavior
and nested database serialization. MongoDB tests use in-memory collections;
they do not require or modify a live server.

These traits are mathematical measurements, not automatic security verdicts.
An isogeny decomposition does not by itself give an efficient discrete-log
transfer. Gonality does not by itself specify an attack cost. Isogeny factors,
quotient maps and automorphisms over an extension must be interpreted over the
actual field and subgroup in use. Many outputs are correlated and should not
be treated as independent evidence in an outlier model.

Mathematical and implementation references:

- [SageMath global function fields](https://doc.sagemath.org/html/en/reference/function_fields/sage/rings/function_field/function_field.html)
- [SageMath differentials and Cartier operator](https://doc.sagemath.org/html/en/reference/function_fields/sage/rings/function_field/differential.html)
- [Girard and Kohel, Classification of genus 3 curves in special strata](https://arxiv.org/abs/math/0603555)
- [Lorenzo Garcia, On different expressions for invariants of hyperelliptic curves of genus 3](https://doi.org/10.2969/jmsj/83418341)
- [Magma invariant-system characteristic conventions](https://magma.maths.usyd.edu.au/magma/handbook/text/1602)
