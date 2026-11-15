"""Bitangents, syzygetic conics and explicit flex incidences of smooth quartics."""

from itertools import combinations
from math import comb

from dissect.utils.curve_analysis import field_value, field_matrix, function_field, AnalysisUnavailable
from dissect.utils.projective_configurations import (
    normalize, projective_key, encoded, substitute, binary_even_section,
    form_coefficients, line_basis, frobenius_configuration,
)


def extended_quartic(curve, field, embedding):
    from sage.all import PolynomialRing

    ring = PolynomialRing(field, names=("x", "y", "z"))
    return ring({powers: embedding(c) for powers, c in curve.plane_polynomial().dict().items()})


def bitangents(curve, field, embedding, params):
    from sage.all import ProjectiveSpace, PolynomialRing, matrix

    q = int(field.cardinality())
    if q*q+q+1 > params["max_candidates"]:
        raise AnalysisUnavailable("Dual-plane enumeration exceeds max_candidates", "resource_limit")
    F = extended_quartic(curve, field, embedding)
    binary_ring = PolynomialRing(field, names=("s", "t"))
    s, t = binary_ring.gens()
    conic_powers = [(2, 0, 0), (1, 1, 0), (1, 0, 1), (0, 2, 0), (0, 1, 1), (0, 0, 2)]
    conics = [F.parent().monomial(*powers) for powers in conic_powers]
    records, lines, conditions = [], [], []
    for coefficients in ProjectiveSpace(field, 2).rational_points():
        line = normalize(coefficients)
        basis = line_basis(line)
        images = [basis[0][i]*s+basis[1][i]*t for i in range(3)]
        section = substitute(F, images, binary_ring)
        contact = binary_even_section(section, 4)
        if contact is None:
            continue
        Q = contact["form"]
        powers = [(2, 0), (1, 1), (0, 2)]
        coefficients_q = [Q.monomial_coefficient(binary_ring.monomial(*power)) for power in powers]
        pivot = next(i for i, c in enumerate(coefficients_q) if c)
        restrictions = [substitute(conic, images, binary_ring) for conic in conics]
        rows = []
        for i in range(3):
            if i != pivot:
                rows.append([coefficients_q[pivot]*R.monomial_coefficient(binary_ring.monomial(*powers[i]))
                             -coefficients_q[i]*R.monomial_coefficient(binary_ring.monomial(*powers[pivot])) for R in restrictions])
        lines.append(line)
        conditions.append(rows)
        records.append({"line": encoded(line), "line_parameterization": field_matrix(matrix(field, basis).transpose()),
                        "section_scalar": field_value(contact["scalar"]),
                        "contact_quadratic": form_coefficients(Q, 2),
                        "contact_factors": contact["factors"],
                        "geometric_contact_points": contact["geometric_contact_points"],
                        "rational_contact_points": contact["rational_contact_points"],
                        "hyperflex": int(contact["geometric_contact_points"] == 1), "theta_parity": "odd"})
    if len(lines) > 28:
        raise ValueError("A smooth odd-characteristic quartic has at most 28 geometric bitangents")
    tetrads = []
    tests = 0
    for indices in combinations(range(len(lines)), 4):
        if tests >= params["max_tetrad_tests"]:
            break
        tests += 1
        kernel = matrix(field, [row for i in indices for row in conditions[i]]).right_kernel()
        if kernel.dimension():
            if kernel.dimension() != 1:
                raise ValueError("A syzygetic tetrad should have a unique contact conic")
            tetrads.append({"bitangent_indices": list(indices), "conic_coefficients": encoded(normalize(kernel.basis()[0]))})
    full = len(lines) == 28
    configurations_complete = tests == comb(len(lines), 4)
    if full and configurations_complete and len(tetrads) != 315:
        raise ValueError("The complete bitangent configuration did not have 315 syzygetic tetrads")
    return {"bitangents": records, "bitangent_count": len(lines), "expected_geometric_bitangents": 28,
            "geometric_locus_complete": int(full),
            "frobenius": frobenius_configuration(lines, int(curve.q())),
            "hyperflex_bitangent_count": sum(record["hyperflex"] for record in records),
            "conic_coefficient_monomials": [list(power) for power in conic_powers],
            "syzygetic_tetrads": tetrads[:params["max_configuration_records"]],
            "syzygetic_tetrad_count": len(tetrads), "tetrad_tests": tests,
            "tetrad_search_complete_on_found_lines": int(configurations_complete),
            "configuration_records_complete": int(len(tetrads) <= params["max_configuration_records"]),
            "status": "ok" if full and configurations_complete and len(tetrads) <= params["max_configuration_records"] else "partial",
            "reason": None if full and configurations_complete and len(tetrads) <= params["max_configuration_records"] else "Missing geometric lines, untested tetrads or truncated configuration records remain"}


def flexes(curve, field, embedding, params):
    from sage.all import ProjectiveSpace, PolynomialRing, matrix, vector

    q = int(field.cardinality())
    if q*q+q+1 > params["max_candidates"]:
        raise AnalysisUnavailable("Projective-plane enumeration exceeds max_candidates", "resource_limit")
    F = extended_quartic(curve, field, embedding)
    L = function_field(curve)
    generic_gaps = [int(n) for n in L.gaps()]
    generic_contact = generic_gaps[-1]-1
    expected = sum(int(P.degree()) for P in L.weierstrass_places())
    ring = PolynomialRing(field, names=("s", "t"))
    s, t = ring.gens()
    records, points, tangents = [], [], []
    for point in ProjectiveSpace(field, 2).rational_points():
        P = normalize(point)
        if F(*P):
            continue
        tangent = normalize([F.derivative(variable)(*P) for variable in F.parent().gens()])
        basis = line_basis(tangent)
        R = next(v for v in basis if matrix(field, [list(P), list(v)]).rank() == 2)
        images = [P[i]*s+R[i]*t for i in range(3)]
        section = substitute(F, images, ring)
        multiplicity = min(power[1] for power, coefficient in section.dict().items() if coefficient)
        if multiplicity < 3:
            continue
        points.append(P)
        tangents.append(tangent)
        records.append({"point": encoded(P), "tangent_line": encoded(tangent),
                        "tangent_intersection_multiplicity": multiplicity,
                        "hyperflex": int(multiplicity == 4),
                        "special_above_generic": int(multiplicity > generic_contact),
                        "tangent_parameterization": field_matrix(matrix(field, [list(P), list(R)]).transpose()),
                        "tangent_section": form_coefficients(section, 4)})
    special_count = sum(record["special_above_generic"] for record in records)
    if special_count > expected:
        raise ValueError("Explicit special flexes exceed the exact Weierstrass locus")
    collinear, concurrent = [], []
    incidence_tests = 0
    for indices in combinations(range(len(points)), 3):
        if incidence_tests >= params["max_incidence_tests"]:
            break
        incidence_tests += 1
        if matrix(field, [points[i] for i in indices]).det() == 0:
            collinear.append(list(indices))
        if matrix(field, [tangents[i] for i in indices]).det() == 0:
            concurrent.append(list(indices))
    complete = special_count == expected
    return {"generic_gaps": generic_gaps, "generic_tangent_contact": generic_contact,
            "ordinary_flex_locus_finite": int(generic_contact == 2),
            "flexes": records[:params["max_point_records"]], "found_flex_count": len(records),
            "indexed_point_coordinates": [encoded(P) for P in points],
            "indexed_tangent_coordinates": [encoded(H) for H in tangents],
            "found_hyperflex_count": sum(record["hyperflex"] for record in records),
            "expected_geometric_special_flexes": expected, "found_special_flexes": special_count,
            "geometric_special_locus_complete": int(complete),
            "frobenius": frobenius_configuration(points, int(curve.q())),
            "collinear_point_triples": collinear[:params["max_configuration_records"]],
            "concurrent_tangent_triples": concurrent[:params["max_configuration_records"]],
            "incidence_tests": incidence_tests,
            "incidence_search_complete_on_found_points": int(incidence_tests == comb(len(points), 3)),
            "configuration_records_complete": int(max(len(collinear), len(concurrent)) <= params["max_configuration_records"]),
            "point_records_complete": int(len(records) <= params["max_point_records"]),
            "status": "ok" if complete and len(records) <= params["max_point_records"] else "partial",
            "reason": ("Point records exceed max_point_records" if len(records) > params["max_point_records"] else
                       None if complete else "The search field does not contain the entire finite special flex locus")}
