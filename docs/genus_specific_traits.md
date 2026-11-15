---
layout: default
title: Genus-specific configurations
---

# Genus-specific configurations in genera two, three and four

Six traits extend the common higher-genus arithmetic with geometry specific
to each genus. They are registered in the existing trait registry and included
in the combined report. All computations use finite fields; geometric counts
refer to the algebraic closure, while explicit searches state their constant
field and completeness separately.

| Genus | Trait | Main results |
| --- | --- | --- |
| 2 | `genus2_rosenhain` | Six branch points, splitting field, Frobenius action, Rosenhain normalization orbit and rational two-torsion rank |
| 2 | `genus2_elliptic_quotients` | Exhaustive finite-field nonhyperelliptic involution search, degree-two elliptic maps, target equations and local ramification |
| 3 | `genus3_bitangents` | Bitangent lines, half-intersection divisors, Frobenius permutation and syzygetic contact conics |
| 3 | `genus3_flexes` | Flex points, tangent lines, contact multiplicities, actual generic behavior and incidence configurations |
| 4 | `genus4_canonical_quadric` | Canonical quadric/cubic, cone versus split/nonsplit quadric, explicit trigonal pencils, descent and automorphism action |
| 4 | `genus4_tritangents` | Odd-theta tritangent planes, contact divisors, Frobenius action, even-theta pencil planes and incidences |

The registry has 50 traits. These six expose 103 top-level fields, with contact
forms, matrices, permutations and maps nested in the records. Their schemas
and parameter bounds are available as `TRAITS[name].OUTPUT` and
`TRAITS[name].DEFAULT_PARAMS`. See the [common reference](detailed_cover_traits.md)
for status semantics, timeouts and exact finite-field serialization.

## Running the examples

```shell
# A mixed-genus report; inapplicable traits are identified automatically.
dissect-higher-genus-report -i examples/genus_specific_curves.json -o genera-report.json

# Full 28-line bitangent configuration over F_49, with Frobenius relative to F_7.
dissect-compute-json -t genus3_bitangents \
  -i examples/higher_genus_curves.json \
  --params '{"extension_degree":2,"timeout":120}'

# All 120 odd-theta planes over F_25, separately from the cone's 26 even-theta planes.
dissect-compute-json -t genus4_tritangents \
  -i examples/genus4_tritangent_curve.json \
  --params '{"extension_degree":2,"max_incidence_tests":10,"timeout":120}'

# Base-field descent of the two genus-four pencils.
dissect-compute-json -t genus4_canonical_quadric \
  -i examples/genus_specific_curves.json
```

The same traits work with `dissect-compute-db`. The full report retains the
input model, field, effective parameters and Sage/Python versions. An
uncomputed quantity remains null. A partial configuration is not presented
as a complete geometric configuration.

## Applicability and conventions

Genus checks run before the Sage computation. A trait for another genus
returns `not_applicable`. The genus-two traits require hyperelliptic models;
the genus-three traits require smooth plane quartics; the genus-four traits
require non-hyperelliptic curves with a recoverable canonical model.

Classical Rosenhain, bitangent and tritangent formulas require odd
characteristic. Flexes use actual function-field generic gaps and can also
handle nonclassical behavior in small characteristic when the Sage backend
supports the model. The current canonical genus-four input adapter and
quadratic-form classification require odd characteristic.

Projective coordinates are normalized by their first nonzero coordinate.
Permutation and incidence indices are zero-based. Extension-field elements
use power-basis coefficient vectors, with modulus coefficients, field size,
characteristic and source-generator embedding included in `search_field` or
the appropriate splitting-field record. Binary-form coefficients are ordered
by ascending power of the first variable: degree two means `[t^2,s*t,s^2]`;
degree three means `[t^3,s*t^2,s^2*t,s^3]`.

Output caps do not invalidate indices: full thin coordinate lists accompany
truncated flex or tritangent contact records. Incidence and record completeness
have separate flags. Optional incidence searches are bounded even when the
main geometric line or plane locus is complete.

## Genus two: Rosenhain configurations

Complete the square in `y^2+h(x)y=f(x)`. The binary branch form is the
homogenization of `h^2+4f` to degree six. An odd affine degree adds the point
at infinity. Irreducible-factor degrees, including infinity, give the
Frobenius orbit lengths and their least common multiple gives the minimal
constant extension splitting all six branch points.

When that field fits the size bound, the trait computes all branch coordinates
and their explicit Frobenius permutation. Each of the 120 ordered choices of
three distinct points is sent to `0,1,infinity` by its unique projectivity.
The other three coordinates form an unordered Rosenhain triple; duplicate
triples are removed. The canonical representative is the smallest such triple
in the declared field basis.

These configurations describe the geometric branch configuration, not a
quadratic-twist class. A Rosenhain equation over the original field, allowing
a scalar twist coefficient, exists exactly when all six branch points are
rational there. The trait does not claim the leading scalar can always be
normalized to one over that field.

Rational Jacobian two-torsion is also computed. In odd characteristic,

```text
J[2] = {even subsets of the six branch points}/<the full six-point set>.
```

The fixed quotient can contain vectors whose Frobenius difference is the full
set, so merely counting branch factors is insufficient. The code solves the
corresponding linear equations over F_2. Splitting-field limits can leave the
configuration orbit uncomputed while preserving exact orbit degrees and
two-torsion rank in a `partial` result.

The example `y^2=x^5-x` over F_7 has branch degrees `[1,1,1,1,2]`, splitting
degree two, rational two-torsion rank three, and five distinct unordered
Rosenhain triples in its full normalization orbit.

## Genus two: discovering elliptic quotients

This is a systematic search for **degree-two** elliptic quotient maps over the
reported finite field. It does not claim to exhaust elliptic maps of higher
degree or maps defined outside that field.

Every automorphism descends to the unique hyperelliptic P1 quotient. The trait
enumerates PGL2 stabilizers of the homogeneous branch form, tests which
transformations lift, and uses the faithful canonical action to select
nonhyperelliptic involutions. This excludes the hyperelliptic deck involution,
whose quotient has genus zero.

For each selected involution sigma, invariant coordinates are constructed from
traces or norms of x and traces of `x^k*(y+h/2)`. A resultant recovers a rational
equation `w^2=R(u)`; square factors are removed to give a squarefree genus-one
polynomial of degree three or four. The implementation verifies

```text
W^2 = F(U),
degree(divisor_of_poles(U)) = 4,
genus(target) = 1.
```

The map therefore has degree two, since U has pole degree two on the target.
Its local different is computed through the common cover engine and checked
against Riemann-Hurwitz. The returned map descriptors can be reused as
`analysis.maps` over the same constant field. Maps found after a constant
extension refer to the base-extended source and the declared embedding.

The record includes an elliptic j-invariant when its cubic conversion or
quartic formula is supported. A missing j-value is null and does not erase
the verified genus-one map. Distinct maps are retained even when their target
j-invariants coincide. `degree_two_search_complete` concerns this finite-field
involution search; `all_elliptic_quotients_complete` remains zero.

## Genus three: bitangents and syzygetic tetrads

A smooth plane quartic has 28 geometric bitangent lines in characteristic
different from two. The dual-plane search enumerates lines over the requested
constant field. On each line, the quartic restricts to a binary quartic.
The line is a bitangent precisely when its intersection divisor is even:

```text
F(P(s,t)) = scalar * q(s,t)^2.
```

The scalar need not be a square over the search field. Tests include the
point at infinity, so degree loss in a dehomogenized polynomial cannot hide an
odd intersection multiplicity. Coefficient recursion rejects nonsquares
without factoring every candidate; accepted half-divisors are factored for
their contact orbits and multiplicities.

Each line record includes its projective coordinates, parameterization,
section scalar, contact quadratic, closed contact factors and hyperflex flag.
A repeated root of the quadratic gives contact multiplicity four at one point.
The explicit Frobenius permutation acts on the line configuration relative to
the original field, not the larger search field.

Four bitangents form a syzygetic tetrad when a conic contains the sum of their
degree-two contact divisors. The code solves exact linear conditions for the
six conic coefficients, including scheme multiplicities at hyperflexes. A
successful tetrad has a unique contact conic. When all 28 lines and all
four-line tests are available, the count is checked to be 315.

For the Fermat quartic over F_7, four bitangents are rational. Over F_49 all 28
are found, including 12 hyperflex lines; all 315 syzygetic tetrads are recovered.
Finding only the base-field four produces `partial`, not a claim that the
curve has only four geometric bitangents.

## Genus three: flexes and nonclassical behavior

The flex trait enumerates points over the chosen field, forms the gradient
tangent line, and computes its exact intersection multiplicity by restriction
to a binary parameterization. Multiplicity at least three is a flex;
multiplicity four is a hyperflex. Coordinates, tangent equations and complete
line sections provide independently checkable certificates.

The actual generic tangent contact is obtained from the source function
field's generic canonical gaps. On a classical quartic this is two. The
trait compares special contact with that actual generic value and uses the
global closed Weierstrass locus to decide whether every finite special flex
has been found.

For `x^4+y^4+z^4=0` in characteristic three, the generic gaps are `[1,2,4]` and
generic tangent contact is three. Every geometric point is an ordinary flex,
so the ordinary flex locus is not finite. Its finite special locus consists
of 28 hyperflexes, found over F_9. The trait distinguishes these statements
rather than imposing a characteristic-zero count.

Configuration records include the Frobenius permutation, collinear triples of
found flex points and concurrent triples of their tangent lines. These
incidences have their own test and record bounds.

## Genus four: canonical input and reconstruction

`CanonicalGenus4` accepts a quadric Q and cubic C in four variables. Terms use
the existing structured coefficient schema with four nonnegative powers.
The model must be a geometrically smooth, one-dimensional intersection with
quadric rank three or four. Adjunction gives canonical bundle O(1), degree
six and genus four; this is a non-hyperelliptic canonical model.

```json
{
  "name": "canonical-genus-four",
  "category": "genus_specific_examples",
  "form": "CanonicalGenus4",
  "field": {"type": "Prime", "p": "0x7", "bits": 3},
  "params": {
    "quadric": [
      {"coeff": 1, "powers": [1,0,0,1]},
      {"coeff": -1, "powers": [0,1,1,0]}
    ],
    "cubic": [
      {"coeff": 1, "powers": [3,0,0,0]},
      {"coeff": 1, "powers": [0,3,0,0]},
      {"coeff": 1, "powers": [0,0,3,0]},
      {"coeff": 2, "powers": [0,0,0,3]}
    ]
  },
  "example": true,
  "standard": false
}
```

The shared arithmetic, Cartier, Weierstrass and cover traits use this model's
native Sage function field. This adds general smooth quadric/cubic input,
including nonsplit canonical quadrics, beyond the previous cyclic-trigonal
adapter.

For an existing non-hyperelliptic genus-four function field, the code instead
recovers the canonical model from a regular differential basis. Quadratic
relations have dimension one; cubic relations have dimension five. One cubic
is chosen outside the four-dimensional space of quadric multiples. Both
equations are verified in the source function field. The differential basis
and equation origin are retained in the result.

## Genus four: pencils, descent and automorphisms

In odd characteristic the symmetric quadric matrix determines three cases:

| Quadric | Geometric trigonal pencils | Pencils over the base field | Frobenius action |
| --- | --- | --- | --- |
| Rank-three cone | 1 | 1 over a finite field | Fixed |
| Smooth split quadric | 2 | 2 | Both fixed |
| Smooth nonsplit quadric | 2 | 0 | Exchanged |

For a smooth four-dimensional form, a square determinant gives the split
case. The nonsplit rulings split over the quadratic constant extension.
Explicit coordinate changes put the quadric into Segre form or cone form;
linear-form ratios then define the degree-three pencils. Their degrees are
certified by intersection with the canonical cubic, and pole degrees are
checked directly for pencils over the source field.

A cone also gives exactly one effective even theta characteristic, its ruling
degree-three line bundle. A smooth canonical quadric has none. The unique
pencil on a cone is preserved by every automorphism. This now upgrades a
completed genus-four cyclic-trigonal cover-stabilizer search to a full
finite-field automorphism search when the recovered quadric is a cone.

Canonical input models additionally support a bounded projective-frame search.
If a complete rational frame enumeration is unavailable, a monomial
stabilizer or the identity subgroup is reported. A witnessed automorphism
exchanging the two pencils proves the answer is yes, even in an incomplete
search. Without a witness, an incomplete search leaves the answer null.
Transformation-record caps do not erase a separately retained exchange witness.

Geometric gonality remains three on a nonsplit quadric, but the absence of a
base-field pencil raises the base-field lower bound to four. The code does not
assert equality with four without an explicit upper-bound certificate.

## Genus four: odd tritangents and even theta pencils

There are 120 geometric odd theta characteristics in odd characteristic. Their
unique effective divisors determine 120 tritangent planes, allowing coincident
contact points. The trait searches planes over a declared finite extension.

For a smooth plane section of the canonical quadric, the conic is parameterized
over the finite field. Restriction of the canonical cubic gives a binary sextic:

```text
Q(P(s,t)) = 0,
C(P(s,t)) = scalar * r(s,t)^2.
```

The cubic r is the degree-three half-intersection divisor. Its factors encode
contact fields and multiplicities. The returned scalar and parameterization
allow exact reconstruction of the intersection certificate.

Degenerate conic sections are handled separately. A pair of lines can still
give an odd tritangent when the cubic passes through their common point with
the appropriate multiplicities. Component restrictions and the shared-point
contribution are checked rather than discarding every singular conic.

On a quadric cone, a tangent plane cuts a double ruling line. Its half-section
is a divisor in the effective **even** theta pencil, with two sections. Such
planes form a positive-dimensional family and are not counted among the 120
odd tritangents. Over a search field of size q there are q+1 of these planes;
that count is checked independently.

For `y^3=x^5-x` over F_5, all 120 odd tritangent planes are found over F_25,
while the cone contributes 26 rational planes in the even-theta pencil family.
Frobenius over F_5 pairs the 120 odd planes into 60 degree-two orbits. The
implementation keeps these configurations distinct.

Plane configuration data include Frobenius cycles, contact schemes and bounded
four-plane concurrence tests. A finite search without all 120 odd planes
returns `partial`; a complete geometric plane set can still have bounded
optional incidence exploration, identified by its own completeness flags.

## Validation and limits

Run the full suite with `sage -python -m pytest -q`. Tests include independent
even-subset enumeration for rational two-torsion; coordinate changes and twists
of Rosenhain configurations; reusable elliptic quotient maps; exact bitangent
and tritangent contact identities; 28-line, 315-tetrad and 120-plane counts;
nonclassical flexes; cone/split/nonsplit descent; genuine automorphisms exchanging
pencils; degenerate-conic contacts; output caps; CLI and database serialization.

These are bounded finite-field configuration searches, not universal symbolic
classification algorithms. Bitangent searches scale with q_search squared and
tritangent searches with q_search cubed. Large splitting fields can exceed a
size, candidate or wall-clock limit. Complete searches over a finite field are
distinguished from complete geometric loci. Characteristic-two theta
configurations and elliptic-quotient discovery of degree greater than two need
additional algorithms.

References:

- [Shaska, Elliptic subfields and automorphisms of genus 2 function fields](https://arxiv.org/abs/math/0107142)
- [Girard and Kohel, Classification of genus 3 curves in special strata](https://arxiv.org/abs/math/0603555)
- [Kulkarni, Sayyary Namin and Sturmfels, Real space sextics and their tritangents](https://arxiv.org/abs/1712.06274)
- [SageMath projective curves](https://doc.sagemath.org/html/en/reference/curves/sage/schemes/curves/projective_curve.html)
