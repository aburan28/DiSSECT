---
layout: default
title: Higher-genus curve traits
---

# Higher-genus curve traits

The fork supports two higher-genus models alongside the existing elliptic
models: `Hyperelliptic` and `SmoothPlane`. Input is validated with SageMath;
genus and curve family come from the model, rather than an unverified label.
Both `dissect-compute-json` and `dissect-compute-db` use this dispatch.

## Traits and mathematical conventions

Here **hyperelliptic** means genus at least two. Non-hyperelliptic curves in
this comparison have genus at least three. Every smooth genus-two curve is
geometrically hyperelliptic. The canonical geometry is over the algebraic
closure; rational branch counts and Frobenius orbits use the input base field.

| Property | Elliptic | Hyperelliptic, genus g >= 2 | Non-hyperelliptic, genus g >= 3 |
| --- | --- | --- | --- |
| Independent regular differentials | 1 | g | g |
| Jacobian dimension | 1 | g | g |
| Curve is its Jacobian after choosing an origin | Yes | No | No |
| Canonical target | P^0 | P^(g-1) | P^(g-1) |
| Canonical map | Constant | Degree 2 onto a rational normal curve | Embedding |
| Degree of canonical image curve | Not applicable | g-1 | 2g-2 |
| Degree-two map to P^1 exists geometrically | Yes | Yes | No |
| Degree-two map unique up to target automorphisms | No | Yes | Not applicable |

`canonical_geometry` reports this table as numeric features. The indicator
`unique_degree_two_map` is zero when no map exists. `plane_quartic` is one
exactly in non-hyperelliptic genus three, where the canonical model is a smooth
plane quartic. A degree-two covering, an involution or an equation `y^2=f(x)`
alone does not distinguish higher-genus curves from elliptic curves.

The points of a higher-genus curve do not carry the natural algebraic group
law of an elliptic curve. Group arithmetic instead uses degree-zero divisor
classes in its higher-dimensional Jacobian. These traits describe that
geometry; they do not implement Jacobian arithmetic or infer its group order.

`hyperelliptic_branch` factors `h(x)^2 + 4*f(x)` for the cover
`y^2 + h(x)*y = f(x)` in odd characteristic. A polynomial of odd degree also
branches at infinity. Its outputs are:

- `geometric_branch_points`: 2g+2 for a valid smooth model.
- `rational_branch_points`: linear factors plus a branch point at infinity.
- `branch_at_infinity`: zero or one.
- `branch_closed_points`: number of Frobenius orbits.
- `branch_orbit_degrees`: sorted irreducible-factor degrees, plus 1 for
  infinity when it branches; their sum is 2g+2.

The genus-one canonical map is constant, so its map degree and image-curve
degree are null. The branch trait returns null for other curve families and
for characteristic two, where wild ramification requires a different formula.
Existing elliptic traits default to supporting elliptic curves only and return
null outputs on higher-genus inputs. New traits declare their supported
families in `SUPPORTED_CURVE_TYPES`; direct calls to a trait's internal
`compute` method bypass this guard, so use the public `trait(curve, **params)`
interface.

These features describe family and genus, not security estimates. Several
are determined by the same inputs and should not be treated as independent
evidence in a statistical comparison. Compare populations within the same
genus, family and field regime before interpreting outlier scores.

## Input schema

See [the complete mixed-family example](../examples/higher_genus_curves.json).
Keep the usual `name`, `category`, `field` and `params` fields. `field` uses the
existing Prime, Binary or extension-field schema. Include `field.bits` when
importing records into the database. Set `example: true` and `standard: false`
for toy curves; the supplied examples are small mathematical examples.

For `form: "Hyperelliptic"`, `params.f` and optional `params.h` are coefficient
lists in **ascending powers**. Omitting `h` gives zero. For example, over F_7:

```json
{
  "name": "genus-two-example",
  "category": "geometry_examples",
  "form": "Hyperelliptic",
  "field": {"type": "Prime", "p": "0x7", "bits": 3},
  "params": {"f": [0, -1, 0, 0, 0, 1], "h": [0]},
  "example": true,
  "standard": false
}
```

This is `y^2 = x^5 - x`, a smooth genus-two curve. The model must be smooth,
including at infinity; genus-zero and genus-one models are rejected by this
input form. In characteristic two a suitable nonzero `h` is required for
the separable hyperelliptic equation. Canonical geometry still applies there.

For `form: "SmoothPlane"`, `params.polynomial` is a list of monomials with
`coeff` and three nonnegative `powers` for `x,y,z`. For example:

```json
{
  "polynomial": [
    {"coeff": 1, "powers": [4, 0, 0]},
    {"coeff": 1, "powers": [0, 4, 0]},
    {"coeff": 1, "powers": [0, 0, 4]}
  ]
}
```

Over F_7 this defines the smooth quartic `x^4+y^4+z^4=0`, of genus three.
Sage checks geometric smoothness. A smooth plane curve of degree d >= 4 is
non-hyperelliptic with genus `(d-1)*(d-2)/2`. Its given plane model coincides
with the canonical embedding only when d=4. For d>4 the canonical image lies
in a larger projective space.

Coefficients can be integers, hexadecimal strings, or existing field-element
objects with `raw` or `poly`. Use `poly` to specify arbitrary elements of an
odd-characteristic extension field. Polynomial input is structured data;
arbitrary expressions are not evaluated.

An optional top-level `genus` must agree with the computed genus. Higher-genus
records may omit `order` and `cofactor`; if supplied, provide both as positive
integers or hexadecimal strings, referring to a Jacobian subgroup and its
cofactor. They are metadata, not verified point counts. Missing values stay
null through database import and result storage. No elliptic Frobenius-trace
formula is applied to them.

This model set does not represent every non-hyperelliptic curve. For example,
general genus-four canonical curves need a quadric/cubic intersection in P^3.
Those require a future validated model adapter. The current traits compute
canonical-map invariants, not explicit canonical coordinates, Jacobian
decompositions, torsion groups or full Weierstrass gap sequences.

## Running and testing

```shell
dissect-compute-json -t canonical_geometry -i examples/higher_genus_curves.json
dissect-compute-json -t hyperelliptic_branch -i examples/higher_genus_curves.json

dissect-database mongodb://localhost:27017/ import -i examples/higher_genus_curves.json
dissect-compute-db -t canonical_geometry --category geometry_examples

# Formula/dispatch tests; Sage integration tests skip without Sage.
python -m pytest -q

# Full suite, including model validation, branch factorization and JSON CLI.
sage -pip install pytest pymongo==4.6.3
sage -python -m pytest -q
```

Database computation also fixes the producer's initial skip: `--skip 0`
processes the first curve, and `--skip N` skips exactly N combinations,
including when N exceeds the input size.
