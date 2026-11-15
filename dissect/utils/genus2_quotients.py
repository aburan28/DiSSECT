"""Construct degree-two elliptic quotients from all rational nonhyperelliptic involutions."""

from dissect.utils.curve_analysis import AnalysisUnavailable, field_value, field_matrix
from dissect.utils.curve_automorphisms import projective_matrices, projective_candidate_count
from dissect.utils.projective_configurations import encoded


def field_input(value):
    result = field_value(value)
    if isinstance(result, int):
        return result
    return {"poly": [{"coeff": c, "power": i} for i, c in enumerate(result) if c]}


def invariant_rational_expression(u, value):
    """Recover R(U) with value=R(u), using a resultant and verify the identity."""
    from sage.all import PolynomialRing

    field = u.parent().constant_base_field()
    ring = PolynomialRing(field, names=("x", "U", "Z"))
    x, U, Z = ring.gens()

    def at(poly):
        return sum(c*x**i for i, c in enumerate(poly))

    first = at(u.numerator())-U*at(u.denominator())
    second = at(value.numerator())-Z*at(value.denominator())
    resultant = first.resultant(second, x)
    target_ring = PolynomialRing(field, "U")
    T = target_ring.gen()
    for factor, _ in resultant.factor():
        if factor.degree(Z) != 1:
            continue
        A, B = target_ring.zero(), target_ring.zero()
        for powers, coefficient in factor.dict().items():
            if powers[0]:
                raise ValueError("Resultant did not eliminate x")
            if powers[2] == 0:
                A += coefficient*T**powers[1]
            elif powers[2] == 1:
                B += coefficient*T**powers[1]
        rational = -A/B
        if rational.numerator()(u)/rational.denominator()(u) == value:
            return rational
    raise ValueError("Could not recover the quotient's invariant rational equation")


def rational_descriptor(value):
    """Structured expression in native [x,y] coordinates, over the current constant field."""
    from sage.all import lcm

    field = value.parent().constant_base_field()
    coefficients = value.list()
    denominator = lcm([c.denominator() for c in coefficients])

    def coefficient(c):
        result = field_value(field(c))
        if isinstance(result, int):
            return result
        return {"poly": [{"coeff": a, "power": i} for i, a in enumerate(result) if a]}

    numerator = []
    for j, entry in enumerate(coefficients):
        polynomial = (entry*denominator).numerator()
        for i, c in enumerate(polynomial):
            if c:
                numerator.append({"coeff": coefficient(c), "powers": [i, j]})
    return {"numerator": numerator,
            "denominator": [{"coeff": coefficient(c), "powers": [i, 0]} for i, c in enumerate(denominator) if c]}


def discover_elliptic_quotients(curve, field, embedding, params):
    from sage.all import PolynomialRing, FunctionField, matrix, identity_matrix
    from dissect.utils.cover_analysis import mapped_cover_ramification

    if projective_candidate_count(int(field.cardinality()), 2) > params["max_candidates"]:
        raise AnalysisUnavailable("PGL2 involution search exceeds max_candidates", "resource_limit")
    K = FunctionField(field, "x")
    x = K.gen()
    univariate = PolynomialRing(field, "t")
    f0, h0 = curve.hyperelliptic_polynomials()
    f, h = [univariate([embedding(c) for c in polynomial]) for polynomial in (f0, h0)]
    completed = f+h**2/4
    ring = PolynomialRing(K, "y")
    L = K.extension(ring.gen()**2+h(x)*ring.gen()-f(x), "y")
    v = L.gen()+h(x)/2
    results = []
    seen = set()
    t = univariate.gen()
    branch = 4*completed
    pivot = next(i for i, c in enumerate(branch) if c)
    I = identity_matrix(field, 2)
    for M in projective_matrices(field, 2):
        a, b, c, d = M.list()
        if M == I or not (M*M)[0, 1] == (M*M)[1, 0] == 0 or (M*M)[0, 0] != (M*M)[1, 1]:
            continue
        A, B = a*t+b, c*t+d
        transformed = sum(coefficient*A**i*B**(6-i) for i, coefficient in enumerate(branch))
        multiplier = transformed[pivot]/branch[pivot]
        if not multiplier or transformed != multiplier*branch or not multiplier.is_square():
            continue
        for gamma in (multiplier.sqrt(), -multiplier.sqrt()):
            # Canonical action is faithful in odd characteristic.
            columns = [(M.det()/gamma)*A**i*B**(1-i) for i in range(2)]
            action = matrix(field, [[p[j] for j in range(2)] for p in columns]).transpose()
            if action == I or action == -I or action*action != I:
                continue
            image_x = (a*x+b)/(c*x+d)
            image_v = gamma*v/(c*L(x)+d)**3
            u = x+image_x
            if u in field:
                u = x*image_x
            w = None
            for exponent in range(3):
                candidate = L(x)**exponent*v+L(image_x)**exponent*image_v
                if candidate:
                    w = candidate
                    break
            square = w*w
            if any(square.list()[1:]):
                raise ValueError("Quotient generator square is not rational in x")
            R = invariant_rational_expression(u, K(square.list()[0]))
            product = R.numerator()*R.denominator()
            square_part = product.parent().one()
            for factor, exponent in product.factor():
                square_part *= factor**(int(exponent)//2)
            quotient, remainder = product.quo_rem(square_part**2)
            if remainder:
                raise ValueError("Square removal in the quotient equation failed")
            target_f = univariate(list(quotient))
            scale = K(R.denominator()(u))/K(square_part(u))
            W = w*L(scale)
            U = L(u)
            if W*W != target_f(U) or target_f.degree() not in (3, 4):
                raise ValueError("The recovered elliptic quotient equation failed verification")
            if U.divisor_of_poles().degree() != 4:
                raise ValueError("Recovered quotient has unexpected map degree")
            key = tuple(M.list())+(gamma,)
            if key in seen:
                continue
            seen.add(key)
            target_K = FunctionField(field, "U")
            target_ring = PolynomialRing(target_K, "V")
            E = target_K.extension(target_ring.gen()**2-target_f(target_K.gen()), "V")
            if E.genus() != 1:
                raise ValueError("The involution quotient did not have genus one")
            ramification = mapped_cover_ramification(L, E, E.hom([W, U]), U, 2)
            if target_f.degree() == 3:
                from sage.all import EllipticCurve
                aa, bb, cc, dd = target_f[3], target_f[2], target_f[1], target_f[0]
                j_invariant = EllipticCurve(field, [0, bb, 0, aa*cc, aa**2*dd]).j_invariant()
            elif field.characteristic() != 3:
                aa, bb, cc, dd, ee = [target_f[i] for i in (4, 3, 2, 1, 0)]
                invariant_I = 12*aa*ee-3*bb*dd+cc**2
                invariant_J = 72*aa*cc*ee+9*bb*cc*dd-27*aa*dd**2-27*bb**2*ee-2*cc**3
                j_invariant = 6912*invariant_I**3/(4*invariant_I**3-invariant_J**2)
            else:
                j_invariant = None
            results.append({"involution_matrix": field_matrix(M), "lift_multiplier": field_value(gamma),
                            "canonical_action": field_matrix(action), "map_degree": 2,
                            "x_image": str(U), "y_image": str(W),
                            "x": rational_descriptor(U), "y": rational_descriptor(W),
                            "target": {"form": "Hyperelliptic", "params": {"f": [field_input(c) for c in target_f]}},
                            "target_polynomial": [field_value(c) for c in target_f],
                            "j_invariant": field_value(j_invariant) if j_invariant is not None else None,
                            "relation_verified": 1, "ramification": ramification})
    return results
