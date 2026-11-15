"""Canonical complete intersections, quadric rulings and descent in genus four."""

from itertools import permutations
from math import prod

from dissect.utils.curve_analysis import cached, function_field, AnalysisUnavailable, field_value, field_matrix
from dissect.utils.projective_configurations import monomial_powers, substitute, normalize, encoded


def homogeneous_form(ring, terms, degree):
    from sage.all import ZZ
    from dissect.utils.higher_genus_curve import field_element

    polynomial = ring.zero()
    for term in terms:
        powers = tuple(ZZ(p) for p in term["powers"])
        if len(powers) != ring.ngens() or any(p < 0 for p in powers) or sum(powers) != degree:
            raise ValueError("Expected homogeneous monomials of degree %d in %d variables" % (degree, ring.ngens()))
        polynomial += field_element(term["coeff"], ring.base_ring())*ring.monomial(*powers)
    if not polynomial:
        raise ValueError("A defining homogeneous form cannot be zero")
    return polynomial


def linear_relations(values, L):
    from sage.all import matrix, lcm

    width = int(L.degree())
    coefficients = [list(v.list())+[L.base_field().zero()]*(width-len(v.list())) for v in values]
    denominator = lcm([c.denominator() for row in coefficients for c in row])
    polynomials = [[(c*denominator).numerator() for c in row] for row in coefficients]
    degree = max(int(p.degree()) for row in polynomials for p in row)
    rows = [[p[i] for p in row for i in range(degree+1)] for row in polynomials]
    return matrix(L.constant_base_field(), rows).left_kernel()


def canonical_data(curve):
    def compute():
        from sage.all import PolynomialRing, matrix

        if curve.genus() != 4 or curve.is_hyperelliptic():
            raise AnalysisUnavailable("Requires a non-hyperelliptic curve of genus four", "not_applicable")
        if curve.model_type() == "CanonicalGenus4":
            Q, C = curve.canonical_polynomials()
            _, coordinates = canonical_function_field(curve)
            return {"quadric": Q, "cubic": C, "coordinates": coordinates,
                    "basis": [str(v) for v in coordinates], "origin": "validated canonical input"}
        L = function_field(curve)
        dx = L(L.base_field().gen()).differential()
        differentials = L.basis_of_holomorphic_differentials()
        coordinates = [w/dx for w in differentials]
        ring = PolynomialRing(curve.field(), names=("X0", "X1", "X2", "X3"))

        def relations(degree):
            powers = monomial_powers(4, degree)
            monomials = [prod(v**e for v, e in zip(coordinates, power)) for power in powers]
            kernel = linear_relations(monomials, L)
            forms = [sum(c*ring.monomial(*power) for c, power in zip(vector, powers)) for vector in kernel.basis()]
            return powers, forms

        _, quadrics = relations(2)
        cubic_powers, cubics = relations(3)
        if len(quadrics) != 1 or len(cubics) != 5:
            raise ValueError("Canonical genus-four relation dimensions must be one and five")
        Q = quadrics[0]
        multiples = [Q*v for v in ring.gens()]
        span = matrix(curve.field(), [[P.monomial_coefficient(ring.monomial(*power)) for power in cubic_powers] for P in multiples]).row_space()
        C = next(P for P in cubics if [P.monomial_coefficient(ring.monomial(*power)) for power in cubic_powers] not in span)
        Q /= next(iter(Q.dict().values()))
        C /= next(iter(C.dict().values()))
        if substitute(Q, coordinates, L) or substitute(C, coordinates, L):
            raise ValueError("Recovered canonical equations failed in the source function field")
        return {"quadric": Q, "cubic": C, "coordinates": coordinates,
                "basis": [str(w) for w in differentials], "origin": "regular-differential relations"}
    return cached(curve, "canonical_genus4", compute)


def quadric_matrix(Q):
    from sage.all import matrix

    if Q.base_ring().characteristic() == 2:
        raise AnalysisUnavailable("This quadric rank/discriminant backend requires odd characteristic")
    variables = Q.parent().gens()
    return matrix(Q.base_ring(), [[Q.base_ring()(Q.derivative(a).derivative(b))/2 for b in variables] for a in variables])


def isotropic_vector(A):
    """An exact isotropic nonsingular vector of a rank >= 3 finite-field form."""
    from sage.all import matrix, vector, identity_matrix

    field = A.base_ring()
    remaining = list(identity_matrix(field, A.nrows()).rows())
    diagonal = []
    while remaining:
        pivot = next((v for v in remaining if v*A*v), None)
        if pivot is None:
            pivot = next((u+v for u, v in permutations(remaining, 2) if (u+v)*A*(u+v)), None)
        if pivot is None:
            break
        norm = pivot*A*pivot
        diagonal.append((pivot, norm))
        orthogonal = [v-(v*A*pivot)/norm*pivot for v in remaining]
        remaining = list(matrix(field, orthogonal).row_space().basis())
    if len(diagonal) < 3:
        raise ValueError("The canonical quadric must have rank at least three")
    (v0, a), (v1, b), (v2, c) = diagonal[:3]
    for x in field:
        square = -(a*x*x+b)/c
        if square.is_square():
            return x*v0+v1+square.sqrt()*v2
    raise ValueError("A finite-field ternary quadric failed to have an isotropic point")


def ruling_transform(Q):
    """Columns put Q into Segre form or a scalar multiple of X1^2-X0*X2."""
    from sage.all import matrix, vector, PolynomialRing, identity_matrix

    A, field = quadric_matrix(Q), Q.base_ring()
    rank = int(A.rank())
    v = isotropic_vector(A)
    standard = list(identity_matrix(field, 4).rows())
    w = next(u for u in standard if v*A*u)
    w = w-(w*A*w)/(2*(v*A*w))*v
    w /= 2*(v*A*w)
    complement = matrix(field, [list(v*A), list(w*A)]).right_kernel().basis()
    if rank == 3:
        radical = A.right_kernel().basis()[0]
        z = next((u for u in complement if u*A*u), None)
        if z is None:
            z = complement[0]+complement[1]
        M = matrix(field, [v, z, -(z*A*z)*w, radical]).transpose()
    elif rank == 4:
        B = matrix(field, complement).transpose()
        binary = B.transpose()*A*B
        ring = PolynomialRing(field, "t")
        t = ring.gen()
        isotropic = []
        if binary[0, 0] == 0:
            isotropic.append(B*vector(field, [1, 0]))
        for root in (binary[0, 0]*t*t+2*binary[0, 1]*t+binary[1, 1]).roots(multiplicities=False):
            isotropic.append(B*vector(field, [root, 1]))
        if len(isotropic) != 2:
            raise AnalysisUnavailable("The two rulings are not split over this field")
        u, z = isotropic
        z /= -2*(u*A*z)
        M = matrix(field, [v, u, z, w]).transpose()
    else:
        raise ValueError("Canonical quadric rank must be three or four")
    if not M.det():
        raise ValueError("Ruling coordinate transform is singular")
    return M


def canonical_function_field(curve):
    """A small-degree field model with explicit original canonical coordinates.

    Parameterizing Q first avoids the generic space-curve normalization. On a
    nonsplit quadric, projection from a rational tangent line gives degree four;
    a bounded search for its point falls back to a degree-at-most-six model.
    """
    def compute():
        from sage.all import FunctionField, PolynomialRing, matrix, vector, identity_matrix

        Q, C = curve.canonical_polynomials()
        field, A = curve.field(), quadric_matrix(Q)
        rank = int(A.rank())
        K = FunctionField(field, "x")
        ring = PolynomialRing(K, "y")
        x, y = ring(K.gen()), ring.gen()
        if rank == 3 or A.det().is_square():
            change = ruling_transform(Q)
            standard = [ring.one(), x, x*x, y] if rank == 3 else [ring.one(), x, y, x*y]
        else:
            basis = list(identity_matrix(field, 4).rows())

            def hyperbolic_frame(point):
                partner = next(v for v in basis if point*A*v)
                partner -= (partner*A*partner)/(2*(point*A*partner))*point
                partner /= 2*(point*A*partner)
                complement = matrix(field, [list(point*A), list(partner*A)]).right_kernel().basis()
                return partner, list(complement), point

            def affine_model(partner, complement, point):
                change = matrix(field, [partner]+complement+[point]).transpose()
                direction = vector(ring, complement[0])*x+vector(ring, complement[1])*y
                last = -sum(A[i,j]*direction[i]*direction[j] for i in range(4) for j in range(4))
                standard = [ring.one(), x, y, last]
                images = [sum(change[i,j]*standard[j] for j in range(4)) for i in range(4)]
                return change, standard, images

            point = isotropic_vector(A)
            partner, complement, point = hyperbolic_frame(point)
            change, standard, images = affine_model(partner, complement, point)
            # First look on Q \cap C, using at most 32 degree-six root searches.
            if C(*point):
                point = None
                polynomial = substitute(C, images, ring)
                constant_ring = PolynomialRing(field, "t")
                evaluate_x = lambda coefficient, value: coefficient.numerator()(value)/coefficient.denominator()(value)
                for index, value in enumerate(field):
                    if index == 32:
                        break
                    section = constant_ring([evaluate_x(coefficient, value) for coefficient in polynomial])
                    roots = section.roots(multiplicities=False)
                    if roots:
                        point = vector(field, [evaluate_x(image(roots[0]), value) for image in images])
                        break
            if point is not None:
                partner, complement, point = hyperbolic_frame(point)
                gradient = vector(field, [C.derivative(variable)(*point) for variable in C.parent().gens()])
                a, b = [gradient*v for v in complement]
                if not (a or b):
                    raise ValueError("A smooth canonical curve must have an independent cubic gradient")
                tangent = b*complement[0]-a*complement[1]
                transverse = next(v for v in complement if matrix(field, [v, tangent]).rank() == 2)
                change, standard, images = affine_model(partner, [transverse, tangent], point)
                if substitute(C, images, ring).degree() > 4:
                    raise ValueError("Projection from the canonical tangent line must have degree at most four")
        images = [sum(change[i,j]*standard[j] for j in range(4)) for i in range(4)]
        equation = substitute(C, images, ring)
        if substitute(Q, images, ring) or not equation or equation.degree() < 2:
            raise ValueError("The canonical quadric parameterization did not give a curve equation")
        L = K.extension(equation, "y")
        coordinates = [sum(L(coefficient)*L.gen()**i for i, coefficient in enumerate(image)) for image in images]
        if any(sum(L(c)*prod(v**e for v,e in zip(coordinates,p)) for p,c in polynomial.dict().items()) for polynomial in (Q,C)):
            raise ValueError("The function-field coordinates failed their canonical equations")
        return L, coordinates

    return cached(curve, "canonical_field_model", compute)


def quadric_profile(curve, params):
    from sage.all import PolynomialRing, matrix, vector
    from dissect.utils.curve_automorphisms import extension_field
    from dissect.utils.projective_configurations import field_description

    data = canonical_data(curve)
    Q, C = data["quadric"], data["cubic"]
    A = quadric_matrix(Q)
    rank = int(A.rank())
    if rank not in (3, 4):
        raise ValueError("A smooth canonical genus-four curve has quadric rank three or four")
    split = rank == 3 or A.det().is_square()
    extension_degree = 1 if split else 2
    result = {"canonical_quadric": str(Q), "canonical_cubic": str(C),
              "canonical_coordinate_basis": data["basis"], "canonical_model_origin": data["origin"],
              "quadric_matrix": field_matrix(A), "quadric_rank": rank,
              "quadric_type": "cone" if rank == 3 else ("smooth_split" if split else "smooth_nonsplit"),
              "quadric_determinant": field_value(A.det()),
              "geometric_trigonal_pencils": 1 if rank == 3 else 2,
              "base_field_trigonal_pencils": 1 if rank == 3 else (2 if split else 0),
              "pencil_splitting_extension_degree": extension_degree,
              "pencil_frobenius_permutation": [0] if rank == 3 else ([0, 1] if split else [1, 0]),
              "effective_even_theta_characteristics": int(rank == 3),
              "pencil_construction_complete": 0}
    if int(curve.q())**extension_degree > params["max_field_size"]:
        return dict(result, status="partial", reason="Ruling splitting field exceeds max_field_size")
    field, embedding = extension_field(curve.field(), extension_degree)
    ring = PolynomialRing(field, names=Q.parent().variable_names())
    extended_Q = ring({powers: embedding(c) for powers, c in Q.dict().items()})
    M = ruling_transform(extended_Q)
    inverse = M.inverse()
    variables = ring.gens()
    transformed = substitute(extended_Q, [sum(M[i, j]*variables[j] for j in range(4)) for i in range(4)], ring)
    expected = variables[1]**2-variables[0]*variables[2] if rank == 3 else variables[0]*variables[3]-variables[1]*variables[2]
    scalar = transformed.monomial_coefficient(next(iter(expected.monomials())))/expected.monomial_coefficient(next(iter(expected.monomials())))
    if transformed != scalar*expected:
        raise ValueError("Ruling transform failed its quadric identity")
    indices = [1] if rank == 3 else [2, 1]
    pencils = []
    L = function_field(curve)
    for i in indices:
        record = {"numerator_linear_form": encoded(inverse.row(i)),
                  "denominator_linear_form": encoded(inverse.row(0)), "degree": 3,
                  "degree_certificate": "A ruling line meets the canonical cubic in a degree-three divisor"}
        if extension_degree == 1:
            numerator = sum(inverse[i, j]*data["coordinates"][j] for j in range(4))
            denominator = sum(inverse[0, j]*data["coordinates"][j] for j in range(4))
            value = numerator/denominator
            if value.divisor_of_poles().degree() != 3:
                raise ValueError("The constructed trigonal function has incorrect pole degree")
            record.update(function=str(value), pole_degree_verified=1)
        pencils.append(record)
    result.update(pencil_field=field_description(field, embedding, curve.field()),
                  ruling_coordinate_change=field_matrix(M), standard_quadric=str(expected),
                  pencils=pencils, pencil_construction_complete=1)
    return result


def pencil_permutation(standard_transform, automorphism):
    """Detect whether a smooth-quadric automorphism exchanges its two rulings."""
    S = standard_transform.inverse()*automorphism*standard_transform
    points = [S.column(i) for i in (0, 2)]
    pairs = [list(P[:2]) if any(P[:2]) else list(P[2:]) for P in points]
    preserved = pairs[0][0]*pairs[1][1] == pairs[0][1]*pairs[1][0]
    return [0, 1] if preserved else [1, 0]


def canonical_automorphisms(curve, field, embedding, params):
    key = ("canonical_automorphisms", field, params["max_candidates"], params["max_transformations"])
    return dict(cached(curve, key, lambda: _canonical_automorphisms(curve, field, embedding, params)))


def _canonical_automorphisms(curve, field, embedding, params):
    """Exhaustive projective-frame search when the complete rational point set permits it."""
    from sage.all import PolynomialRing, ProjectiveSpace, matrix, vector, diagonal_matrix, identity_matrix
    from itertools import product
    from dissect.utils.curve_automorphisms import group_summary, extension_field

    data = canonical_data(curve)
    ring = PolynomialRing(field, names=data["quadric"].parent().variable_names())
    Q, C = [ring({power: embedding(c) for power, c in P.dict().items()}) for P in (data["quadric"], data["cubic"])]
    q = int(field.cardinality())
    points = []
    if q**3+q*q+q+1 <= params["max_candidates"]:
        points = [list(P) for P in ProjectiveSpace(field, 3).rational_points() if not Q(*P) and not C(*P)]
    count = prod(range(len(points)-4, len(points)+1)) if len(points) >= 5 else 0

    def frame(indices):
        A = matrix(field, [points[i] for i in indices[:4]]).transpose()
        if not A.det():
            return None
        coordinates = A.inverse()*vector(field, points[indices[4]])
        if not all(coordinates):
            return None
        return A*diagonal_matrix(field, coordinates)

    source = None
    if count and count <= params["max_candidates"]:
        source = next((M for indices in permutations(range(len(points)), 5) if (M := frame(indices)) is not None), None)
    complete = source is not None

    def candidates():
        if source is None:
            if 24*(q-1)**3 <= params["max_candidates"]:
                for permutation in permutations(range(4)):
                    for a, b, c in product([v for v in field if v], repeat=3):
                        M = matrix(field, 4)
                        for i, (j, value) in enumerate(zip(permutation, [field.one(), a, b, c])):
                            M[i, j] = value
                        yield M
            else:
                yield identity_matrix(field, 4)
            return
        inverse = source.inverse()
        for indices in permutations(range(len(points)), 5):
            target = frame(indices)
            if target is not None:
                M = target*inverse
                yield M/next(c for c in M.list() if c)

    quadric_ideal = ring.ideal(Q)
    reduced_cubic = C.reduce(quadric_ideal.groebner_basis())
    pivot_Q = next(iter(Q.dict()))
    pivot_C = next(iter(reduced_cubic.dict()))
    rank = int(quadric_matrix(Q).rank())
    if rank == 4:
        splitting_degree = 1 if quadric_matrix(Q).det().is_square() else 2
        split_field, split_embedding = extension_field(field, splitting_degree)
        split_ring = PolynomialRing(split_field, names=ring.variable_names())
        split_Q = split_ring({power: split_embedding(c) for power, c in Q.dict().items()})
        standard_transform = ruling_transform(split_Q)
    matrices, transformations, permutations_seen = [], [], set()
    exchange_witness = None
    for M in candidates():
        if any(Q(*list(M*vector(field, P))) or C(*list(M*vector(field, P))) for P in points[:8]):
            continue
        images = [sum(M[i, j]*ring.gen(j) for j in range(4)) for i in range(4)]
        transformed_Q = substitute(Q, images, ring)
        multiplier_Q = transformed_Q.dict().get(pivot_Q, field.zero())/Q.dict()[pivot_Q]
        if not multiplier_Q or transformed_Q != multiplier_Q*Q:
            continue
        transformed_C = substitute(C, images, ring).reduce(quadric_ideal.groebner_basis())
        multiplier_C = transformed_C.dict().get(pivot_C, field.zero())/reduced_cubic.dict()[pivot_C]
        if not multiplier_C or transformed_C != multiplier_C*reduced_cubic:
            continue
        action = M.det()/(multiplier_Q*multiplier_C)*M.transpose()
        matrices.append(action)
        permutation = [0] if rank == 3 else pencil_permutation(standard_transform, M.apply_map(split_embedding))
        permutations_seen.add(tuple(permutation))
        if permutation == [1, 0] and exchange_witness is None:
            exchange_witness = field_matrix(M)
        if len(transformations) < params["max_transformations"]:
            transformations.append({"projective_matrix": field_matrix(M), "differential_action": field_matrix(action),
                                    "quadric_multiplier": field_value(multiplier_Q), "cubic_multiplier": field_value(multiplier_C),
                                    "pencil_permutation": permutation})
    result = group_summary(matrices, field, 4)
    return dict(result, transformations=transformations, search_complete=int(complete),
                pencil_permutations_seen=[list(p) for p in sorted(permutations_seen)], pencil_exchange_witness=exchange_witness,
                base_field_order=len(matrices) if complete else None,
                reduced_automorphism_order=None, branch_stabilizer_order=None,
                reason=None if complete else "The rational point set did not permit an exhaustive frame search; a monomial stabilizer or the identity subgroup is reported")
