"""Canonical-image equations and linear syzygies, checked in the source field."""

from math import comb, prod

from dissect.utils.curve_analysis import function_field, AnalysisUnavailable, field_value
from dissect.utils.canonical_genus4 import linear_relations, canonical_function_field
from dissect.utils.projective_configurations import monomial_powers, substitute, field_description


def canonical_relations(curve, params):
    from sage.all import PolynomialRing, matrix, VectorSpace

    g, degree = int(curve.genus()), params["relation_degree"]
    if degree > 6 or comb(g+degree-1, degree) > params["max_monomials"]:
        raise AnalysisUnavailable("Canonical monomial count or relation degree exceeds the bound", "resource_limit")
    L = function_field(curve)
    if curve.model_type() == "CanonicalGenus4":
        _, coordinates = canonical_function_field(curve)
        basis = [str(v) for v in coordinates]
        origin = "validated canonical input coordinates"
    else:
        dx = L(L.base_field().gen()).differential()
        differentials = L.basis_of_holomorphic_differentials()
        coordinates = [omega/dx for omega in differentials]
        basis, origin = [str(omega) for omega in differentials], "regular differentials divided by a common separating differential"
    if len(coordinates) != g:
        raise ValueError("Canonical coordinates must have dimension g")
    ring = PolynomialRing(curve.field(), names=tuple("X%d" % i for i in range(g)))
    generators, profiles = [], []
    for n in range(1, degree+1):
        powers = monomial_powers(g, n)
        values = [prod(v**e for v, e in zip(coordinates, power)) for power in powers]
        kernel = linear_relations(values, L)
        hilbert = len(powers)-kernel.dimension()
        expected = n*(g-1)+1 if curve.is_hyperelliptic() else (g if n == 1 else (2*n-1)*(g-1))
        if hilbert != expected:
            raise ValueError("The canonical multiplication map has an unexpected Hilbert function")
        vector_space = VectorSpace(curve.field(), len(powers))
        multiples = [P*ring.monomial(*power) for P in generators for power in monomial_powers(g, n-P.degree())]
        span = vector_space.subspace([[P.monomial_coefficient(ring.monomial(*power)) for power in powers] for P in multiples])
        new_generators = []
        for vector in kernel.basis():
            if vector not in span:
                polynomial = sum(c*ring.monomial(*power) for c, power in zip(vector, powers))
                if substitute(polynomial, coordinates, L):
                    raise ValueError("A recovered canonical equation is not zero in the source field")
                new_generators.append(polynomial)
                span = span+vector_space.subspace([vector])
        generators.extend(new_generators)
        profiles.append({"degree": n, "ambient_monomials": len(powers), "multiplication_rank": int(hilbert),
                         "relation_dimension": int(kernel.dimension()), "minimal_generator_count": len(new_generators),
                         "full_pluricanonical_dimension": g if n == 1 else (2*n-1)*(g-1),
                         "canonical_subring_deficiency": (g if n == 1 else (2*n-1)*(g-1))-int(hilbert)})
    quadrics = [P for P in generators if P.degree() == 2]
    syzygies = []
    if quadrics:
        powers = monomial_powers(g, 3)
        products = [P*variable for P in quadrics for variable in ring.gens()]
        rows = [[P.monomial_coefficient(ring.monomial(*power)) for power in powers] for P in products]
        for relation in matrix(curve.field(), rows).left_kernel().basis():
            forms = [sum(relation[i*g+j]*ring.gen(j) for j in range(g)) for i in range(len(quadrics))]
            if sum(a*b for a, b in zip(forms, quadrics)):
                raise ValueError("A linear syzygy failed its polynomial identity")
            syzygies.append([[field_value(relation[i*g+j]) for j in range(g)] for i in range(len(quadrics))])
    cap = params["max_relation_records"]
    records = [{"degree": int(P.degree()), "polynomial": str(P),
                "terms": [{"coefficient": field_value(c), "powers": list(power)} for power, c in P.dict().items()]} for P in generators]
    complete = len(records) <= cap and len(syzygies) <= cap
    return {"coefficient_field": field_description(curve.field(), lambda value: value, curve.field()),
            "coordinate_basis": basis, "coordinate_functions": [str(v) for v in coordinates], "coordinate_origin": origin,
            "variables": list(ring.variable_names()), "degree_profiles": profiles,
            "minimal_generators": records[:cap], "minimal_generator_count": len(records),
            "quadric_generator_indices": [i for i, P in enumerate(generators) if P.degree() == 2],
            "indexed_quadric_terms": [record["terms"] for record in records if record["degree"] == 2],
            "linear_quadric_syzygies": syzygies[:cap], "linear_quadric_syzygy_count": len(syzygies),
            "relation_degree_bound": degree, "equations_verified_in_source": 1, "records_complete": int(complete),
            "canonical_ring_generated_in_degree_one": int(not curve.is_hyperelliptic()),
            "ideal_generation_complete": int(not curve.is_hyperelliptic() and degree >= (4 if g == 3 else 3) or curve.is_hyperelliptic() and degree >= 2),
            "status": "ok" if complete else "partial", "reason": None if complete else "Equation or syzygy records exceed max_relation_records"}
