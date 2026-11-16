"""Odd-theta tritangents and the separate even-theta pencil on a quadric cone."""

from itertools import combinations
from math import comb

from dissect.utils.curve_analysis import AnalysisUnavailable, field_matrix, field_value
from dissect.utils.canonical_genus4 import canonical_data, quadric_matrix, isotropic_vector
from dissect.utils.projective_configurations import (
    normalize, encoded, substitute, line_basis, binary_even_section,
    form_coefficients, frobenius_configuration,
)


def conic_parameterization(conic):
    from sage.all import PolynomialRing, matrix, identity_matrix

    field = conic.base_ring()
    A = quadric_matrix(conic)
    P = isotropic_vector(A)
    basis = identity_matrix(field, 3).rows()
    change = next(matrix(field, [u, v, P]).transpose() for u in basis for v in basis
                  if matrix(field, [u, v, P]).det())
    binary = PolynomialRing(field, names=("s", "t"))
    s, t = binary.gens()
    # Q(change*(s,t,z)) = A(s,t)+z*B(s,t), since the last column is on Q.
    B = change.transpose()*A*change
    quadratic = B[0, 0]*s*s+2*B[0, 1]*s*t+B[1, 1]*t*t
    linear = 2*B[0, 2]*s+2*B[1, 2]*t
    coordinates = [s*linear, t*linear, -quadratic]
    images = [sum(change[i, j]*coordinates[j] for j in range(3)) for i in range(3)]
    if substitute(conic, images, binary):
        raise ValueError("The conic parameterization failed its equation")
    return images


def singular_conic_contact(conic, cubic, ambient_images, curve, search_embedding, params):
    from sage.all import PolynomialRing, matrix
    from dissect.utils.curve_automorphisms import extension_field

    field = conic.base_ring()
    original_search_field = field
    contact_embedding = lambda value: value
    A = quadric_matrix(conic)
    rank = int(A.rank())
    if rank == 1:
        row = next(row for row in A.rows() if any(row))
        basis = line_basis(row)
        ring = PolynomialRing(field, names=("s", "t"))
        s, t = ring.gens()
        images = [basis[0][i]*s+basis[1][i]*t for i in range(3)]
        section = substitute(cubic, images, ring)
        return {"kind": "even_theta_pencil", "half_divisor_dimension": 2,
                "contact_cubic": form_coefficients(section, 3),
                "contact_parameterization": [form_coefficients(substitute(image, images, ring), 1) for image in ambient_images]}
    if rank != 2:
        raise ValueError("Unexpected restricted-quadric rank")
    vertex = A.right_kernel().basis()[0]
    if cubic(*vertex):
        return None  # Two disjoint degree-three components cannot give an even divisor.
    factors = list(conic.factor())
    if len(factors) == 1:
        if int(field.cardinality())**2 > params["max_field_size"]:
            raise AnalysisUnavailable("A singular-conic contact needs a quadratic field beyond max_field_size", "resource_limit")
        contact_field, embedding = extension_field(field, 2)
        contact_embedding = embedding
        ring3 = PolynomialRing(contact_field, names=conic.parent().variable_names())
        conic = ring3({p: embedding(c) for p, c in conic.dict().items()})
        cubic = ring3({p: embedding(c) for p, c in cubic.dict().items()})
        ambient_images = [ring3({p: embedding(c) for p, c in image.dict().items()}) for image in ambient_images]
        vertex = [embedding(c) for c in vertex]
        field = contact_field
        factors = list(conic.factor())
    if len(factors) != 2 or any(f.degree() != 1 or e != 1 for f, e in factors):
        raise ValueError("A rank-two conic did not split into two distinct lines")
    binary = PolynomialRing(field, names=("s", "t"))
    s, t = binary.gens()
    pieces, multiplicities = [], []
    for factor, _ in factors:
        coefficients = [factor.monomial_coefficient(v) for v in factor.parent().gens()]
        basis = line_basis(coefficients)
        R = next(v for v in basis if matrix(field, [list(vertex), list(v)]).rank() == 2)
        images = [vertex[i]*s+R[i]*t for i in range(3)]
        section = substitute(cubic, images, binary)
        order = min(power[1] for power, coefficient in section.dict().items() if coefficient)
        if order % 2 == 0:
            return None
        residual = sum(c*s**p[0]*t**(p[1]-order) for p, c in section.dict().items())
        contact = binary_even_section(residual, 3-order)
        if contact is None:
            return None
        multiplicities.append(order)
        pieces.append({"line": encoded(coefficients), "vertex_intersection_multiplicity": order,
                       "residual_scalar": field_value(contact["scalar"]),
                       "residual_half_divisor": form_coefficients(contact["form"], (3-order)//2),
                       "residual_factors": contact["factors"],
                       "parameterization": [form_coefficients(substitute(image, images, binary), 1) for image in ambient_images]})
    from dissect.utils.projective_configurations import field_description

    return {"kind": "odd_theta", "conic_type": "two_lines", "component_contacts": pieces,
            "common_vertex": encoded([image(*vertex) for image in ambient_images]),
            "vertex_half_divisor_multiplicity": sum(multiplicities)//2,
            "contact_field": field_description(field, contact_embedding, original_search_field)}


def tritangent_configuration(curve, field, embedding, params):
    from sage.all import PolynomialRing, ProjectiveSpace, matrix

    q = int(field.cardinality())
    if q**3+q*q+q+1 > params["max_candidates"]:
        raise AnalysisUnavailable("Dual three-space enumeration exceeds max_candidates", "resource_limit")
    data = canonical_data(curve)
    ring4 = PolynomialRing(field, names=data["quadric"].parent().variable_names())
    Q, C = [ring4({power: embedding(c) for power, c in P.dict().items()}) for P in (data["quadric"], data["cubic"])]
    cone = quadric_matrix(Q).rank() == 3
    if cone:
        from dissect.utils.canonical_genus4 import ruling_transform
        change = ruling_transform(Q)
        variables = ring4.gens()
        C_standard = substitute(C, [sum(change[i, j]*variables[j] for j in range(4)) for i in range(4)], ring4)
    plane_ring = PolynomialRing(field, names=("u", "v", "w"))
    plane_variables = plane_ring.gens()
    planes, odd, even = [], [], []
    for coefficients in ProjectiveSpace(field, 3).rational_points():
        H = normalize(coefficients)
        if cone:
            standard_H = matrix(field, [H])*change
            values = list(standard_H.row(0))
            if values[3]:
                binary = PolynomialRing(field, names=("s", "t"))
                s, t = binary.gens()
                coordinates = [t*t, s*t, s*s, -(values[0]*t*t+values[1]*s*t+values[2]*s*s)/values[3]]
                section = substitute(C_standard, coordinates, binary)
                contact = binary_even_section(section, 6)
                if contact is None:
                    continue
                original = [sum(change[i, j]*coordinates[j] for j in range(4)) for i in range(4)]
                record = {"kind": "odd_theta", "conic_type": "smooth", "plane": encoded(H),
                          "section_scalar": field_value(contact["scalar"]),
                          "contact_cubic": form_coefficients(contact["form"], 3), "contact_factors": contact["factors"],
                          "geometric_contact_points": contact["geometric_contact_points"],
                          "parameterization": [form_coefficients(image, 2) for image in original]}
                planes.append(H)
                odd.append(record)
                continue
        basis = line_basis(H)
        images = [sum(basis[j][i]*plane_variables[j] for j in range(3)) for i in range(4)]
        conic = substitute(Q, images, plane_ring)
        cubic = substitute(C, images, plane_ring)
        rank = int(quadric_matrix(conic).rank())
        if rank == 3:
            parameterization = conic_parameterization(conic)
            binary = parameterization[0].parent()
            section = substitute(cubic, parameterization, binary)
            contact = binary_even_section(section, 6)
            if contact is None:
                continue
            record = {"kind": "odd_theta", "conic_type": "smooth",
                      "section_scalar": field_value(contact["scalar"]),
                      "contact_cubic": form_coefficients(contact["form"], 3),
                      "contact_factors": contact["factors"],
                      "geometric_contact_points": contact["geometric_contact_points"],
                      "parameterization": [form_coefficients(substitute(image, parameterization, binary), 2) for image in images]}
        else:
            record = singular_conic_contact(conic, cubic, images, curve, embedding, params)
            if record is None:
                continue
        record.update(plane=encoded(H), plane_parameterization=field_matrix(matrix(field, basis).transpose()),
                      restricted_quadric=str(conic), restricted_cubic=str(cubic))
        if record["kind"] == "even_theta_pencil":
            if not cone:
                raise ValueError("A smooth canonical quadric cannot have a double-line plane section")
            even.append(record)
        else:
            planes.append(H)
            odd.append(record)
    if len(odd) > 120:
        raise ValueError("Odd tritangent count exceeds the 120 odd theta characteristics")
    if cone and len(even) != q+1:
        raise ValueError("A quadric cone must have q+1 rational even-theta tangent planes")
    incidences, tests = [], 0
    for indices in combinations(range(len(planes)), 4):
        if tests >= params["max_incidence_tests"]:
            break
        tests += 1
        kernel = matrix(field, [planes[i] for i in indices]).right_kernel()
        if kernel.dimension():
            incidences.append({"plane_indices": list(indices), "intersection_basis": field_matrix(matrix(field, kernel.basis()))})
    full = len(odd) == 120
    return {"canonical_quadric": str(Q), "canonical_cubic": str(C),
            "tritangent_planes": odd[:params["max_plane_records"]],
            "indexed_plane_coordinates": [encoded(H) for H in planes],
            "tritangent_count": len(odd), "expected_odd_theta_characteristics": 120,
            "geometric_odd_locus_complete": int(full),
            "frobenius": frobenius_configuration(planes, int(curve.q())),
            "even_theta_pencil_planes": even[:params["max_plane_records"]],
            "even_theta_plane_count": len(even), "even_theta_locus_positive_dimensional": int(cone),
            "effective_even_theta_characteristics": int(cone),
            "concurrent_plane_quadruples": incidences[:params["max_configuration_records"]],
            "concurrent_plane_quadruple_count": len(incidences),
            "configuration_records_complete": int(len(incidences) <= params["max_configuration_records"]),
            "incidence_tests": tests,
            "incidence_search_complete_on_found_planes": int(tests == comb(len(planes), 4)),
            "plane_records_complete": int(max(len(odd), len(even)) <= params["max_plane_records"]),
            "status": "ok" if full and max(len(odd), len(even)) <= params["max_plane_records"] else "partial",
            "reason": ("Plane records exceed max_plane_records" if max(len(odd), len(even)) > params["max_plane_records"] else
                       None if full else "Odd tritangents outside the search field remain uncomputed; even-theta pencil planes are recorded separately")}
