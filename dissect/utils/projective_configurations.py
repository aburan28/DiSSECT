"""Exact projective normalization, Frobenius permutations and contact sections."""

from dissect.utils.curve_analysis import AnalysisUnavailable, field_value
from dissect.utils.curve_automorphisms import extension_field


def normalize(coordinates):
    coordinates = list(coordinates)
    pivot = next((c for c in coordinates if c), None)
    if pivot is None:
        raise ValueError("Zero projective coordinates")
    return tuple(c/pivot for c in coordinates)


def value_key(value):
    encoded = field_value(value)
    return tuple(encoded) if isinstance(encoded, list) else (encoded,)


def projective_key(coordinates):
    return tuple(value_key(c) for c in normalize(coordinates))


def encoded(coordinates):
    return [field_value(c) for c in coordinates]


def field_description(field, embedding, source):
    image = field_value(embedding(source.gen()))
    return {"size": int(field.cardinality()), "modulus": str(field.modulus()),
            "characteristic": int(field.characteristic()), "absolute_degree": int(field.degree()),
            "modulus_coefficients": [int(c) for c in field.modulus()],
            "source_generator_image": image if isinstance(image, list) else [image]}


def search_field(curve, params):
    extension = params.get("extension_degree", 1)
    if int(curve.q())**extension > params["max_field_size"]:
        raise AnalysisUnavailable("Search extension exceeds max_field_size", "resource_limit")
    return extension_field(curve.field(), extension)


def frobenius_configuration(points, q):
    normalized = [normalize(point) for point in points]
    indices = {point: i for i, point in enumerate(normalized)}
    permutation = [indices[normalize([c**q for c in point])] for point in normalized]
    seen, cycles = set(), []
    for i in range(len(points)):
        if i in seen:
            continue
        cycle, j = [], i
        while j not in seen:
            seen.add(j)
            cycle.append(j)
            j = permutation[j]
        cycles.append(cycle)
    return {"permutation": permutation, "cycles": cycles,
            "orbit_degrees": sorted(len(cycle) for cycle in cycles),
            "rational_over_base": sum(len(cycle) == 1 for cycle in cycles)}


def monomial_powers(variables, degree):
    if variables == 1:
        return [(degree,)]
    return [(i,)+rest for i in range(degree+1) for rest in monomial_powers(variables-1, degree-i)]


def substitute(polynomial, images, ring):
    result = ring.zero()
    for powers, coefficient in polynomial.dict().items():
        term = ring(coefficient)
        for value, exponent in zip(images, powers):
            term *= value**exponent
        result += term
    return result


def binary_even_section(section, degree):
    """Test even intersection divisor, including the point at infinity.

    A scalar times a square is sufficient; the scalar need not be a square
    over the search field. The returned half-divisor is a homogeneous form.
    """
    from sage.all import PolynomialRing

    if not section:
        raise ValueError("A component of the curve is contained in the tested section")
    ring, field = section.parent(), section.base_ring()
    s, t = ring.gens()
    unary_ring = PolynomialRing(field, "z")
    z = unary_ring.gen()
    unary = sum(coefficient*z**powers[0] for powers, coefficient in section.dict().items())
    infinity = degree-int(unary.degree())
    if infinity % 2:
        return None
    scalar = unary.leading_coefficient()
    normalized = unary/scalar
    half_degree_unary = int(unary.degree())//2
    if field.characteristic() != 2:
        coefficients = [field.zero()]*(half_degree_unary+1)
        coefficients[-1] = field.one()
        for i in range(half_degree_unary-1, -1, -1):
            power = half_degree_unary+i
            known = sum(coefficients[j]*coefficients[power-j] for j in range(i+1, half_degree_unary+1)
                        if i < power-j <= half_degree_unary)
            coefficients[i] = (normalized[power]-known)/2
        half = unary_ring(coefficients)
        if half*half != normalized:
            return None
    else:
        factors = unary.factor()
        if any(int(e) % 2 for _, e in factors):
            return None
        half = unary_ring.one()
        for factor, exponent in factors:
            half *= factor**(int(exponent)//2)
    records = []
    for factor, exponent in half.factor():
        records.append({"polynomial": [field_value(c) for c in factor], "degree": int(factor.degree()),
                        "half_divisor_multiplicity": int(exponent),
                        "intersection_multiplicity": 2*int(exponent)})
    half_degree = degree//2
    homogeneous = sum(coefficient*s**i*t**(half_degree-i) for i, coefficient in enumerate(half))
    if infinity:
        records.append({"at_infinity": True, "degree": 1,
                        "half_divisor_multiplicity": infinity//2, "intersection_multiplicity": infinity})
    return {"form": homogeneous, "scalar": scalar, "factors": records,
            "geometric_contact_points": sum(r["degree"] for r in records),
            "rational_contact_points": sum(r["degree"] == 1 for r in records)}


def form_coefficients(form, degree):
    return [field_value(form.monomial_coefficient(form.parent().monomial(*powers)))
            for powers in monomial_powers(form.parent().ngens(), degree)]


def line_basis(coefficients):
    from sage.all import matrix

    return list(matrix(coefficients[0].parent(), [coefficients]).right_kernel().basis())
