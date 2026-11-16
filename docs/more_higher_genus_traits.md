---
layout: default
title: Further traits for genus two and higher
---

# Further traits for genus two and higher

Six traits add new isogeny, torsion, divisor and canonical-ring data. The
registry contains 56 traits; these six expose 101 top-level output fields.
Nested records retain exact coefficients, kernels, matrices, equations and
verification results. All six work through `dissect-compute-json`,
`dissect-compute-db` and `dissect-higher-genus-report`.

| Trait | Applicability | New information |
| --- | --- | --- |
| `genus2_richelot` | Genus two, odd characteristic | Fifteen maximal isotropic kernels, Frobenius action, Richelot target equations and split products |
| `hyperelliptic_twists` | Hyperelliptic genus at least two | Explicit involution twist, twisted zeta polynomial, curve/Jacobian extension orders |
| `jacobian_torsion` | Genus at least two, prime ell different from p | Rational ell-torsion rank bounds; exact two-torsion modules for odd-characteristic hyperelliptic curves |
| `symmetric_power_divisors` | All supported higher-genus curves | Closed points, symmetric powers, effective degree-two Picard classes and Riemann-Roch identities |
| `frobenius_orders` | Simple ordinary Jacobians | CM field, Frobenius/Weil order discriminants and indices, maximal real subfield and endomorphism-order bounds |
| `canonical_relations` | All supported higher-genus curves within the monomial budget | Multiplication-map ranks, minimal canonical-image equations and linear quadric syzygies |

These complement the [arithmetic and cover traits](detailed_cover_traits.md)
and [genus-specific configurations](genus_specific_traits.md). They do not
apply elliptic formulas to higher-dimensional Jacobians.

## Running the examples

```shell
# Eight examples: ordinary genus two, nonmaximal Weil order, split branch
# Richelot kernels, characteristic two, plane quartic, hyperelliptic genus
# four/five, and a nonhyperelliptic genus-four cyclic trigonal curve.
dissect-higher-genus-report \
  -i examples/more_higher_genus_curves.json -o more-genera-report.json

dissect-compute-json -t genus2_richelot \
  -i examples/more_higher_genus_curves.json
dissect-compute-json -t jacobian_torsion \
  -i examples/more_higher_genus_curves.json \
  --params '{"ell":3,"extension_bound":12}'
dissect-compute-json -t canonical_relations \
  -i examples/more_higher_genus_curves.json \
  --params '{"relation_degree":4,"max_relation_records":128}'
```

Genus/model guards return `not_applicable`. Characteristic restrictions,
resource limits and timeouts have distinct statuses. Missing exact ranks,
unverified point counts and undetermined endomorphism orders remain null.
The report records effective parameters, model, field and software versions.
Use `TRAITS[name].OUTPUT`, `.INPUT` and `.DEFAULT_PARAMS` for the complete
machine-readable schema.

## Genus two: Richelot isogenies

Work in odd characteristic and complete the square in
`y^2+h(x)y=f(x)`. Write `B=f+h^2/4`. Its six geometric branch points,
including infinity when B has degree five, have fifteen partitions into
three unordered pairs. Each partition determines a maximal isotropic
subgroup of the four-dimensional geometric two-torsion module.

The kernel has four elements and the associated polarized abelian-surface
isogeny has degree four. This is an isogeny of Jacobians, rather than a
degree-two morphism from the source curve to an elliptic curve. The existing
`genus2_elliptic_quotients` trait searches those curve morphisms separately.

The implementation computes the branch splitting extension, all fifteen
partitions, and the Frobenius permutation on this finite set. A kernel is
defined over the base field exactly when its partition is Frobenius-stable;
its three pairs may still be permuted. Cycle lengths give each kernel's
minimal constant-field definition degree.

In the full branch field, form quadratic homogeneous factors, allowing an
affine linear factor for a pair involving infinity. Normalize them so that

```text
G1(x) G2(x) G3(x) = B(x).
```

Let Delta be the determinant of their three coefficient rows in the order
`1,x,x^2`. Use cyclic brackets

```text
H1 = G2' G3 - G2 G3',
H2 = G3' G1 - G3 G1',
H3 = G1' G2 - G1 G2'.
```

Every record is checked against the fundamental exact identity

```text
G1(u) H1(v) + G2(u) H2(v) + G3(u) H3(v) = -Delta (u-v)^2.
```

For nonzero Delta, the target is the smooth genus-two equation

```text
Y^2 = H1(X) H2(X) H3(X) / Delta.
```

The polynomial's degree and squarefreeness are checked. With this bracket
and coefficient convention, the scalar/sign is part of the equation:
discarding it can produce the wrong quadratic twist. Small targets are also
checked by direct counts over their field and its quadratic extension;
their Weil polynomial must equal the source's polynomial over the same field.

When Delta is zero, the quotient with its induced principal polarization
is a product of elliptic curves. The trait records this degeneration and
sets `target_components_computed=0`; it does not invent component equations.

`target_field` describes the field of the constructed equations and the
source embedding. A rational kernel does not by itself make the chosen
full-splitting-field equation a descended base-field equation. The trait
reports the kernel's smaller definition degree separately and does not
claim to construct that descent or an explicit rational isogeny map on
Mumford coordinates. Literal distinct target equations are counted without
claiming distinct isomorphism classes.

The curve `y^2=x^5-x` over F_7 has three rational kernels. Over F_49 the trait
finds nine genus-two targets and six split-product kernels. The split-branch
F_11 example has all fifteen rational kernels and nonzero odd Frobenius
coefficients, so its point-count tests detect a wrong twist sign.

If the branch field exceeds `max_field_size`, abstract pair partitions and
their Frobenius action remain exact in a `partial` result. Its branch labels
then refer to abstract factor cycles, not to uncomputed coordinates.
`max_target_point_work` bounds q_search squared per target (default 128).
Skipped verification remains null; it does not suppress the algebraic target.

## Hyperelliptic involution twists

The hyperelliptic involution acts as minus one on H1. Twisting by its
nontrivial constant cocycle therefore sends `L(T)` to `L(-T)` in every
genus. It is not valid to replace the entire degree-2g polynomial with an
elliptic trace.

In odd characteristic choose a nonsquare d and return

```text
y_new^2 = d * (f+h^2/4),     h_new=0.
```

In characteristic two, choose a constant a with absolute trace one and return

```text
y_new^2 + h*y_new = f + a*h^2.
```

The trace-one choice is made in the actual constant field. In an even-degree
extension of F_2, the constant one has trace zero and would give the wrong
twist. The coefficient field, its modulus and the selected parameter are
returned explicitly. Coefficients use the common power-basis encoding.

The trait gives both Jacobian orders, their gcd, both curve point counts,
and extension comparisons. In odd extensions the curve counts add to
`2(q^n+1)`; in even extensions both Frobenius polynomials agree. The field
`twist_parameter_trivialized` records the constant square/Artin-Schreier
class, not a complete isomorphism test for curves with extra automorphisms.
Such curves can have additional twists or already be isomorphic to this
chosen twist. Accordingly `all_twists_classified=0`.

## Jacobian torsion: exact modules versus polynomial constraints

For a prime ell different from the characteristic, geometric `J[ell]` has
dimension 2g over F_ell. Factor the Frobenius characteristic polynomial as

```text
P(T) mod ell = product f_i(T)^m_i.
```

This polynomial does not usually determine the Jordan blocks. In extension
degree n, write n as `ell^v * n0`, with n0 prime to ell. A primary factor
contributes only if its roots satisfy `lambda^n=1`. In that case the fixed
space dimension lies between

```text
degree(f_i) * min(m_i, ell^v)  and  degree(f_i) * m_i.
```

The lower bound places the entire multiplicity in one Jordan block; the
upper bound uses separate blocks. The upper bound is also capped by the
ell-adic valuation of the exact Jacobian order in that extension. Equal
bounds give an exact rank. A squarefree reduction gives exact ranks at
every tested degree. No unverified semisimplicity assumption is made.

The extension order is computed from the integral Frobenius polynomial,
not from the number of points of the curve. Divisibility by ell and the
fixed-space bounds are checked for consistency. The first searched
nonzero rank is the exact least torsion extension, because every smaller
degree was tested. Absence within the bound is not nonexistence.

When the semisimple order s is found, a certified full-torsion extension
degree is `s * ell^ceil(log_ell(max m_i))`. This is an upper bound. It is
reported as an exact full extension only when the unipotent part is known
to be trivial or an exact module is available.

### Exact hyperelliptic two-torsion

In odd characteristic an exact module is available in every hyperelliptic
genus:

```text
J[2] = {even subsets of the 2g+2 branch points}/<the full set>.
```

Factor degrees give a branch Frobenius permutation. In the even-subset
basis `e_i+e_last`, quotienting by the full set gives an explicit 2g by 2g
matrix over F_2. Its characteristic polynomial is checked against P mod 2.
Its fixed spaces and exact order give ranks and the minimal full-torsion
extension, including degrees beyond the finite rank-search bound. The least
nontrivial extension is found among divisors of that exact matrix order.
Basis labels follow abstract factor cycles; they are not claimed to be
individual branch coordinates.

The genus-two F_7 example has rational rank three and full two-torsion over
F_49. The hyperelliptic genus-four example has an eight-dimensional module
whose Frobenius has order nine. Characteristic-p torsion is left unsupported:
the full finite group scheme is not the same as an etale 2g-dimensional module.

## Symmetric powers and effective divisors

Let `N_n=#C(F_(q^n))`. The number B_n of closed points of degree n is

```text
B_n = (1/n) * sum_(d divides n) mu(d) N_(n/d).
```

The number A_n of rational effective divisors of degree n, equivalently the
number of F_q-points on `Sym^n(C)`, is the coefficient of T^n in

```text
L(T)/((1-T)(1-q*T)) = product_(d>=1) (1-T^d)^(-B_d).
```

The code computes the rational-function recurrence and independently
reconstructs the truncated closed-point Euler product. Integrality,
nonnegativity and equality are checked. Degree zero is included with A_0=1.

Over a finite field every Picard degree component has `#J(F_q)` rational
classes. For n at least 2g-1, Riemann-Roch makes every fiber a projective
space of constant dimension n-g, so

```text
A_n = #J(F_q) * (q^(n+1-g)-1)/(q-1).
```

For degree g, residual linear series give the useful identity

```text
A_g = #J(F_q) + q*A_(g-2).
```

Thus genus two gives `A_2=#J+q`, genus three gives `A_3=#J+q*N_1`, and
genus four gives `A_4=#J+q*A_2`. These are checked whenever the bound reaches
the required degree, rather than extrapolating constant high-degree fibers
into the special-divisor range.

In degree two, the hyperelliptic pencil is the unique class with two
sections. Its q+1 divisors must count as one Picard class. The number of
effective rational Picard degree-two classes is therefore `A_2-q` in the
hyperelliptic case, and A_2 in the nonhyperelliptic case. Remaining Picard
classes are reported as noneffective. This does not enumerate actual
Mumford representatives or provide general Jacobian group arithmetic.

## Frobenius orders and CM fields

For a simple ordinary Jacobian, the irreducible degree-2g Weil polynomial
defines the commutative base-field endomorphism algebra `K=Q(pi)`, a CM field.
Ordinarity is tested by the middle zeta coefficient being prime to p;
irreducibility is checked over Q. Other cases return `not_applicable` and
remain available to the Honda-Tate decomposition trait.

The code constructs

```text
Z[pi] <= Z[pi,q/pi] <= End_Fq(J) <= O_K.
```

It reports both generated-order discriminants and indices in O_K, the CM
field discriminant and signature, and the minimal polynomial and field
discriminant of the maximal real subfield generated by `pi+q/pi`.
The involution `pi -> q/pi` is constructed and checked exactly. Both orders
must satisfy `disc(order)=disc(K)*index^2`.

If the Weil order is already maximal, the containing-order bounds force
the actual base-field endomorphism order to be maximal. Otherwise the
actual order is undetermined; a nonmaximal generated order does not prove
the actual ring is nonmaximal. No algebraic-closure endomorphism-ring claim
is inferred from these base-field data.

For `y^2=x^5+x^3+2x` over F_5, K has discriminant 2304 and signature (0,2),
its real subfield has discriminant 24, and the Weil order is maximal.
The example `y^2=x^5+2x^3+x+1` has a Weil-order index of four and leaves the
actual intermediate order undetermined.

## Canonical equations and linear syzygies

The trait constructs canonical coordinate functions from a regular
differential basis divided by a common separating differential. Validated
canonical genus-four input uses its explicit coordinate functions. In each
degree n, monomials are evaluated in the source function field. Exact linear
relations give the kernel of the canonical multiplication map.

Lower-degree generator multiples are removed by linear algebra to find
new minimal generators. Every returned polynomial is checked to vanish in
the function field. A second kernel calculation finds linear combinations
of the quadratic generators with linear-form coefficients; each syzygy is
checked as a polynomial identity.

| Curve | Minimal canonical-image equations | Linear syzygies among quadrics |
| --- | --- | --- |
| Genus two | P1 has no image equations | 0 |
| Hyperelliptic genus three | One conic | 0 |
| Nonhyperelliptic genus three | One quartic | 0 |
| Hyperelliptic genus four | Three quadrics defining a twisted cubic | 2 |
| Nonhyperelliptic genus four | One quadric and one cubic | 0 |
| Hyperelliptic genus five | Six quadrics defining a rational normal quartic | 8 |

The degree profiles retain ambient monomial counts, multiplication ranks,
relation dimensions, new generator counts and the dimension of the full
pluricanonical space. For hyperelliptic curves the canonical image ring is
a proper subring of the full canonical ring in higher degrees. Its Hilbert
function is `n(g-1)+1`, whereas the full pluricanonical dimension for n>=2
is `(2n-1)(g-1)`. In the nonhyperelliptic case the degree-one generators
produce the full canonical ring. These dimensions are checked exactly.

`ideal_generation_complete` requires enough degrees: two for a rational
normal canonical image, four for a plane quartic, and three for other
nonhyperelliptic canonical curves. It uses the canonical-ideal generation
theorem; it is not a full minimal free resolution or a Green-conjecture
test. Only linear syzygies among quadrics are computed here.

`max_relation_records` caps verbose generators and syzygies independently.
`indexed_quadric_terms` retains every quadric's exact coefficients in the
syzygy ordering, even when verbose records are capped. Field metadata and
variable order make those certificates independently reconstructible.

## Bounds and validation

All six traits use the existing wall-clock, field and genus limits. Arithmetic
traits additionally bound point-count work. Torsion and twist extension
bounds cannot exceed 64; effective-divisor degree lies between 2 and 64.
The torsion prime is at most 65537 and must differ from the characteristic.
Canonical relations default to genus at most six, degree four and at most
1000 monomials; relation degree cannot exceed six. CM fields default to
degree at most twelve. Large computations return explicit resource/timeout
states rather than silently substituting zeros.

Run `sage -python -m pytest -q`. The tests check all fifteen Richelot kernels,
the bilinear identities and Weil-polynomial equality including nonzero h
and a nonsquare source twist; odd- and even-characteristic twist models;
trace-one selection over F_4; exact branch-subset torsion and ambiguous
Jordan cases; independent divisor multiset and Euler-product counts; both
maximal and nonmaximal CM orders; equations and syzygies through genus five;
record caps, bounds, CLI behavior and nested database encoding.

Mathematical background: Cassels and Flynn, *Prolegomena to a Middlebrow
Arithmetic of Curves of Genus 2* (Richelot quotients); Stichtenoth,
*Algebraic Function Fields and Codes* (zeta functions and Riemann-Roch);
Waterhouse, *Abelian varieties over finite fields* (ordinary endomorphism
algebras); Arbarello, Cornalba, Griffiths and Harris, *Geometry of Algebraic
Curves*, Volume I (canonical multiplication maps and Petri's theorem).
