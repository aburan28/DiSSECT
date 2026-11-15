"""Exhaustive finite-field stabilizers and explicitly bounded subgroups.

Counts always refer to the reported finite field, not the algebraic closure.
"""

from itertools import product, permutations

from dissect.utils.curve_analysis import field_value, field_matrix, AnalysisUnavailable


def projective_matrices(field, dimension):
    from sage.all import matrix

    elements = list(field)
    size = dimension*dimension
    for pivot in range(size):
        for tail in product(elements, repeat=size-pivot-1):
            M = matrix(field, dimension, [field.zero()]*pivot + [field.one()] + list(tail))
            if M.det():
                yield M


def projective_candidate_count(q, dimension):
    return sum(q**i for i in range(dimension*dimension))


def extension_field(field, degree):
    if degree == 1:
        return field, lambda value: value
    extension = field.extension(degree, "a")
    if field.is_prime_field():
        return extension, field.hom(extension)
    root = field.modulus().change_ring(extension).roots(multiplicities=False)[0]
    return extension, field.hom([root], extension)


def group_summary(matrices, field, genus):
    from sage.all import MatrixGroup, identity_matrix

    group = MatrixGroup(matrices)
    order = int(group.order())
    if order != len(matrices):
        raise ValueError("Enumerated automorphisms are not a faithful closed group")
    distribution, quotients = {}, {}
    for M in matrices:
        n = int(M.multiplicative_order())
        distribution[n] = distribution.get(n, 0)+1
        if n % int(field.characteristic()):
            quotient_genus = genus-int((M-identity_matrix(field, genus)).rank())
            key = (n, quotient_genus)
            quotients[key] = quotients.get(key, 0)+1
    return {
        "computed_subgroup_order": order,
        "group_structure": group.structure_description(),
        "abelian": int(group.is_abelian()),
        "center_order": int(group.center().order()),
        "element_order_distribution": [{"order": n, "count": distribution[n]} for n in sorted(distribution)],
        "tame_cyclic_quotient_genera": [{"automorphism_order": n, "quotient_genus": g, "elements": count} for (n, g), count in sorted(quotients.items())],
        "canonical_action_faithful": 1,
    }


def hyperelliptic_odd(curve, field, embedding, params):
    from sage.all import PolynomialRing, matrix, identity_matrix

    q, genus = int(field.cardinality()), curve.genus()
    if projective_candidate_count(q, 2) > params["max_candidates"]:
        matrices = [identity_matrix(field, genus), -identity_matrix(field, genus)]
        result = group_summary(matrices, field, genus)
        return dict(result, search_complete=0, base_field_order=None, reduced_automorphism_order=1,
                    transformations=[], branch_stabilizer_order=None,
                    reason="PGL2 search exceeds max_candidates; only the hyperelliptic involution is reported")
    ring = PolynomialRing(field, "t")
    t = ring.gen()
    f, h = curve.hyperelliptic_polynomials()
    branch = ring([embedding(c) for c in h**2+4*f])
    degree = 2*genus+2
    pivot = next(i for i, c in enumerate(branch) if c)
    matrices, transformations = [], []
    branch_stabilizer, reduced = 0, 0
    for M in projective_matrices(field, 2):
        a, b, c, d = M.list()
        A, B = a*t+b, c*t+d
        transformed = sum(coefficient*A**i*B**(degree-i) for i, coefficient in enumerate(branch))
        multiplier = transformed[pivot]/branch[pivot]
        if not multiplier or transformed != multiplier*branch:
            continue
        branch_stabilizer += 1
        if not multiplier.is_square():
            continue
        reduced += 1
        square_root = multiplier.sqrt()
        for gamma in (square_root, -square_root):
            columns = [(M.det()/gamma)*A**i*B**(genus-1-i) for i in range(genus)]
            action = matrix(field, [[polynomial[j] for j in range(genus)] for polynomial in columns]).transpose()
            matrices.append(action)
            if len(transformations) < params["max_transformations"]:
                transformations.append({"pgl2_matrix": field_matrix(M), "square_root_multiplier": field_value(gamma),
                                        "differential_action": field_matrix(action)})
    result = group_summary(matrices, field, genus)
    return dict(result, search_complete=1, base_field_order=len(matrices),
                reduced_automorphism_order=reduced, branch_stabilizer_order=branch_stabilizer,
                transformations=transformations, reason=None)


def artin_schreier_solutions(polynomial):
    """Polynomial solutions of s^2+s=A, via removal of even leading powers."""
    ring, field = polynomial.parent(), polynomial.base_ring()
    x, remainder, solution = ring.gen(), polynomial, ring.zero()
    while remainder.degree() > 0:
        degree = int(remainder.degree())
        if degree % 2:
            return []
        term = remainder.leading_coefficient().sqrt()*x**(degree//2)
        solution += term
        remainder -= term**2+term
    constants = (x**2+x-ring(remainder[0])).roots(multiplicities=False)
    return [solution+constant for constant in constants]


def hyperelliptic_even(curve, field, embedding, params):
    from sage.all import PolynomialRing, matrix, identity_matrix, SymmetricGroup, PermutationGroup

    f, h = curve.hyperelliptic_polynomials()
    genus = curve.genus()
    if h.degree() != 0 or int(field.cardinality())*(int(field.cardinality())-1) > params["max_candidates"]:
        return {
            "computed_subgroup_order": 2, "group_structure": "C2", "abelian": 1, "center_order": 2,
            "element_order_distribution": [{"order": 1, "count": 1}, {"order": 2, "count": 1}],
            "tame_cyclic_quotient_genera": [{"automorphism_order": 1, "quotient_genus": genus, "elements": 1}],
            "canonical_action_faithful": 0, "search_complete": 0, "base_field_order": None,
            "reduced_automorphism_order": 1, "branch_stabilizer_order": None,
            "transformations": [{"x_image": "x", "y_image": "y+h(x)"}],
            "reason": "Only the deck involution is certified for nonconstant h or an oversized affine search",
        }
    ring = PolynomialRing(field, "t")
    t = ring.gen()
    normalized = ring([embedding(c) for c in f]) / embedding(h[0])**2
    automorphisms = []
    for a in field:
        if not a:
            continue
        for b in field:
            for s in artin_schreier_solutions(normalized(a*t+b)+normalized):
                automorphisms.append((a, b, s))
    lookup = {entry: i+1 for i, entry in enumerate(automorphisms)}
    symmetric = SymmetricGroup(len(automorphisms))
    permutations_list = []
    for a, b, s in automorphisms:
        images = []
        for c, d, r in automorphisms:
            composite = (a*c, c*b+d, s+r(a*t+b))
            images.append(lookup[composite])
        permutations_list.append(symmetric(images))
    group = PermutationGroup(permutations_list)
    if group.order() != len(automorphisms):
        raise ValueError("Artin-Schreier automorphism group is not closed")
    distribution, quotients, transformations = {}, {}, []
    for (a, b, s), permutation in zip(automorphisms, permutations_list):
        order = int(permutation.order())
        distribution[order] = distribution.get(order, 0)+1
        columns = [a*(a*t+b)**i for i in range(genus)]
        action = matrix(field, [[polynomial[j] for j in range(genus)] for polynomial in columns]).transpose()
        if order % 2:
            key = (order, genus-int((action-identity_matrix(field, genus)).rank()))
            quotients[key] = quotients.get(key, 0)+1
        if len(transformations) < params["max_transformations"]:
            transformations.append({"x_multiplier": field_value(a), "x_translation": field_value(b),
                                    "normalized_y_shift": [field_value(c) for c in s],
                                    "order": order, "differential_action": field_matrix(action)})
    return {
        "computed_subgroup_order": int(group.order()), "group_structure": group.structure_description(),
        "abelian": int(group.is_abelian()), "center_order": int(group.center().order()),
        "element_order_distribution": [{"order": n, "count": distribution[n]} for n in sorted(distribution)],
        "tame_cyclic_quotient_genera": [{"automorphism_order": n, "quotient_genus": g, "elements": count} for (n, g), count in sorted(quotients.items())],
        "canonical_action_faithful": 0, "search_complete": 1, "base_field_order": len(automorphisms),
        "reduced_automorphism_order": len(automorphisms)//2, "branch_stabilizer_order": None,
        "transformations": transformations, "reason": None,
    }


def monomial_matrices(field):
    from sage.all import matrix

    for permutation in permutations(range(3)):
        for a, b in product([c for c in field if c], repeat=2):
            M = matrix(field, 3)
            for i, (j, coefficient) in enumerate(zip(permutation, [field.one(), a, b])):
                M[i, j] = coefficient
            yield M


def frame_matrices(field, points):
    """All projectivities mapping one fixed rational frame into the given set."""
    from sage.all import matrix, diagonal_matrix, vector

    def frame(indices):
        basis = matrix(field, [points[i] for i in indices[:3]]).transpose()
        if not basis.det():
            return None
        coordinates = basis.inverse()*vector(field, points[indices[3]])
        if not all(coordinates):
            return None
        return basis*diagonal_matrix(field, coordinates)

    source = next((F for indices in permutations(range(len(points)), 4) if (F := frame(indices)) is not None), None)
    if source is None:
        return None
    inverse = source.inverse()

    def candidates():
        for indices in permutations(range(len(points)), 4):
            target = frame(indices)
            if target is not None:
                M = target*inverse
                pivot = next(c for c in M.list() if c)
                yield M/pivot
    return candidates()


def plane_automorphisms(curve, field, embedding, params):
    from sage.all import PolynomialRing, ProjectiveSpace, matrix, vector

    q = int(field.cardinality())
    complete = projective_candidate_count(q, 3) <= params["max_candidates"]
    ring = PolynomialRing(field, names=("x", "y", "z"))
    variables = ring.gens()
    F = ring({powers: embedding(coefficient) for powers, coefficient in curve.plane_polynomial().dict().items()})
    degree = int(F.degree())
    powers = [(i, j, degree-3-i-j) for i in range(degree-2) for j in range(degree-2-i)]
    forms = [variables[0]**i*variables[1]**j*variables[2]**k for i, j, k in powers]
    pivot = next(iter(F.dict()))
    # A few rational points reject most candidates before polynomial substitution.
    points = []
    if q*q+q+1 <= params["max_candidates"]:
        points = [list(point) for point in ProjectiveSpace(field, 2).rational_points() if F(*point) == 0]
    candidates = None
    frame_bound = len(points)*(len(points)-1)*(len(points)-2)*(len(points)-3)
    if not complete and len(points) >= 4 and frame_bound <= params["max_candidates"]:
        candidates = frame_matrices(field, points)
        complete = candidates is not None
    if candidates is None:
        if not complete and 6*(q-1)**2 > params["max_candidates"]:
            raise AnalysisUnavailable("Even the monomial subgroup search exceeds max_candidates", "resource_limit")
        candidates = projective_matrices(field, 3) if complete else monomial_matrices(field)
    matrices, transformations = [], []
    for M in candidates:
        if any(F(*list(M*vector(field, point))) for point in points[:8]):
            continue
        images = [sum(M[i, j]*variables[j] for j in range(3)) for i in range(3)]
        substitution = ring.hom(images, ring)
        transformed = substitution(F)
        multiplier = transformed.monomial_coefficient(ring.monomial(*pivot)) / F.dict()[pivot]
        if not multiplier or transformed != multiplier*F:
            continue
        canonical_forms = [M.det()/multiplier*substitution(form) for form in forms]
        action = matrix(field, [[form.monomial_coefficient(ring.monomial(*power)) for power in powers] for form in canonical_forms]).transpose()
        matrices.append(action)
        if len(transformations) < params["max_transformations"]:
            transformations.append({"projective_matrix": field_matrix(M), "equation_multiplier": field_value(multiplier),
                                    "differential_action": field_matrix(action)})
    result = group_summary(matrices, field, curve.genus())
    return dict(result, search_complete=int(complete), base_field_order=len(matrices) if complete else None,
                transformations=transformations, reduced_automorphism_order=None, branch_stabilizer_order=None,
                reason=None if complete else "Full PGL3 search exceeds max_candidates; the computed group is the monomial stabilizer only")


def rational_cube_roots(value):
    """All cube roots in a rational function field, using exact factorization."""
    numerator, denominator = value.numerator(), value.denominator()
    numerator_factors, denominator_factors = numerator.factor(), denominator.factor()
    if any(int(e) % 3 for _, e in list(numerator_factors)+list(denominator_factors)):
        return []
    scalar = numerator_factors.unit()/denominator_factors.unit()
    ring = numerator.parent()
    constants = (ring.gen()**3-scalar).roots(multiplicities=False)
    root = value.parent().one()
    for factor, exponent in numerator_factors:
        root *= value.parent()(factor)**(int(exponent)//3)
    for factor, exponent in denominator_factors:
        root /= value.parent()(factor)**(int(exponent)//3)
    return [constant*root for constant in constants]


def trigonal_automorphisms(curve, field, embedding, params):
    from sage.all import FunctionField, PolynomialRing, matrix

    K = FunctionField(field, "x")
    x = K.gen()
    ring = PolynomialRing(K, "y")
    f = PolynomialRing(field, "t")([embedding(c) for c in curve.cover_polynomial()])
    L = K.extension(ring.gen()**3-f(x), "y")
    V, to_differential, to_vector = L.space_of_holomorphic_differentials()
    basis = [to_differential(v) for v in V.basis()]
    complete_cover = projective_candidate_count(int(field.cardinality()), 2) <= params["max_candidates"]
    if complete_cover:
        candidates = projective_matrices(field, 2)
    else:
        from sage.all import identity_matrix
        candidates = [identity_matrix(field, 2)]
    matrices, transformations = [], []
    for M in candidates:
        a, b, c, d = M.list()
        image_x = (a*x+b)/(c*x+d)
        for exponent in (1, 2):
            ratio = f(image_x)/f(x)**exponent
            for root in rational_cube_roots(ratio):
                image_y = L(root)*L.gen()**exponent
                morphism = L.hom([image_y, L(image_x)])
                columns = [to_vector(morphism(w/L(x).differential())*L(image_x).differential()) for w in basis]
                action = matrix(field, columns).transpose()
                matrices.append(action)
                if len(transformations) < params["max_transformations"]:
                    transformations.append({"pgl2_matrix": field_matrix(M), "y_power": exponent,
                                            "y_multiplier": str(root), "differential_action": field_matrix(action)})
    result = group_summary(matrices, field, curve.genus())
    # Two distinct degree-three pencils would bound the genus by (3-1)^2=4.
    unique_pencil = curve.genus() > 4
    if curve.genus() == 4 and field.characteristic() != 2:
        from dissect.utils.canonical_genus4 import canonical_data, quadric_matrix

        unique_pencil = quadric_matrix(canonical_data(curve)["quadric"]).rank() == 3
    full = complete_cover and unique_pencil
    return dict(result, search_complete=int(full), cover_stabilizer_complete=int(complete_cover),
                base_field_order=len(matrices) if full else None,
                reduced_automorphism_order=None, branch_stabilizer_order=None,
                transformations=transformations,
                reason=None if full else "The reported group stabilizes the chosen trigonal pencil; pencil uniqueness is not certified for this model")
