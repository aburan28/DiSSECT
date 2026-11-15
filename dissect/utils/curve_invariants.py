"""Classical invariants, with explicit transvectant and differential normalization.

Definitions: Girard--Kohel, arXiv:math/0603555, section 3; and
Lorenzo Garcia, J. Math. Soc. Japan 74 (2022), section 3.
This is an independent implementation of the mathematical constructions.
"""

from math import factorial, comb, gcd


def differentiate(polynomial, orders):
    for variable, order in zip(polynomial.parent().gens(), orders):
        if order:
            polynomial = polynomial.derivative(variable, order)
    return polynomial


def transvectant(f, g, order):
    m, n = int(f.degree()), int(g.degree())
    if not f or not g:
        return f.parent().zero()
    if order > min(m, n):
        return f.parent().zero()
    field = f.base_ring()
    coefficient = field(factorial(m-order)*factorial(n-order)) / field(factorial(m)*factorial(n))
    return coefficient * sum((-1)**i*comb(order, i)*differentiate(f, (order-i, i))*differentiate(g, (i, order-i)) for i in range(order+1))


def shioda_invariants(binary):
    g = transvectant(binary, binary, 4)
    k = transvectant(binary, binary, 6)
    h = transvectant(k, k, 2)
    m = transvectant(binary, k, 4)
    n = transvectant(binary, h, 4)
    p = transvectant(g, k, 4)
    q = transvectant(g, h, 4)
    pairs = [(binary, binary, 8), (binary, g, 8), (k, k, 4), (m, k, 4),
             (k, h, 4), (m, h, 4), (p, h, 4), (n, h, 4), (q, h, 4)]
    return [binary.base_ring()(transvectant(a, b, r)) for a, b, r in pairs]


def differential_pair(operator, polynomial):
    result = polynomial.parent().zero()
    for powers, coefficient in operator.dict().items():
        result += coefficient*differentiate(polynomial, powers)
    return result


def quartic_contravariants(F):
    from sage.all import PolynomialRing

    field = F.base_ring()
    ring = PolynomialRing(field, names=("a", "b", "u", "v"))
    a, b, u, v = ring.gens()
    restriction = ring.zero()
    for powers, coefficient in F.dict().items():
        restriction += coefficient*a**powers[0]*b**powers[1]*(-u*a-v*b)**powers[2]
    coefficients = [ring.zero() for _ in range(5)]
    for powers, coefficient in restriction.dict().items():
        i, j, r, s = powers
        coefficients[j] += coefficient*u**r*v**s / field(comb(4, j))
    c0, c1, c2, c3, c4 = coefficients
    sigma = c0*c4-4*c1*c3+3*c2**2
    psi = c0*c2*c4-c0*c3**2-c1**2*c4+2*c1*c2*c3-c2**3
    x, y, z = F.parent().gens()

    def homogenize(value, degree):
        result = F.parent().zero()
        for powers, coefficient in value.dict().items():
            if powers[0] or powers[1] or sum(powers[2:]) > degree:
                raise ValueError("Contravariant homogenization failed")
            r, s = powers[2:]
            result += coefficient*x**r*y**s*z**(degree-r-s)
        return result

    return homogenize(sigma, 4), homogenize(psi, 6)


def dixmier_ohno_invariants(F):
    from sage.all import matrix

    field, ring = F.base_ring(), F.parent()
    variables = ring.gens()
    sigma, psi = quartic_contravariants(F)
    H = matrix(ring, [[F.derivative(a).derivative(b) for b in variables] for a in variables]).det()
    rho = differential_pair(F, psi) / field(144)
    tau = differential_pair(rho, F) / field(12)
    xi = differential_pair(sigma, H) / field(72)
    eta = differential_pair(xi, sigma) / field(12)
    nu = differential_pair(eta, differential_pair(rho, H)) / field(8)

    def quadratic(form):
        return matrix(field, [[field(form.derivative(a).derivative(b)) / field(2) for b in variables] for a in variables])

    def dot(A, B):
        return sum(a*b for a, b in zip(A.list(), B.list()))

    R, T, X, E, N = [quadratic(f) for f in (rho, tau, xi, eta, nu)]
    I3 = field(differential_pair(sigma, F)) / field(144)
    I6 = (field(differential_pair(psi, H))-8*I3**2) / field(4608)
    discriminant_resultant = ring.macaulay_resultant([F.derivative(v) for v in variables])
    I27 = discriminant_resultant / field(2**54)
    values = [I3, I6, dot(T, R), dot(X, R), R.det(), dot(T, E), T.det(), X.det(),
              dot(T.adjugate(), R.adjugate()), dot(X.adjugate(), R.adjugate()),
              E.det(), dot(N, E), I27]
    names = ["I3", "I6", "I9", "J9", "I12", "J12", "I15", "J15", "I18", "J18", "I21", "J21", "I27"]
    weights = [3, 6, 9, 9, 12, 12, 15, 15, 18, 18, 21, 21, 27]
    return names, weights, values


def absolute_ratios(values, weights):
    """Weight-zero ratios on a declared nonzero-coordinate chart."""
    from dissect.utils.curve_analysis import field_value

    index = next((i for i, value in enumerate(values) if value), None)
    if index is None:
        return None, []
    pivot, pivot_weight = values[index], weights[index]
    ratios = []
    for i, (value, weight) in enumerate(zip(values, weights)):
        common = gcd(weight, pivot_weight)
        ratios.append({"index": i, "numerator_power": pivot_weight//common,
                       "pivot_power": weight//common,
                       "value": field_value(value**(pivot_weight//common) / pivot**(weight//common))})
    return index, ratios


def binary_discriminant(polynomial, homogeneous_degree):
    value = polynomial.discriminant()
    if polynomial.degree() == homogeneous_degree-1:
        value *= polynomial.leading_coefficient()**2
    elif polynomial.degree() != homogeneous_degree:
        raise ValueError("The binary form has a multiple root at infinity")
    return value


def higher_binary_invariants(polynomial, degree):
    from sage.all import PolynomialRing

    field = polynomial.base_ring()
    names, weights = ["Delta"], [2*degree-2]
    values = [binary_discriminant(polynomial, degree)]
    if field.characteristic() > degree:
        ring = PolynomialRing(field, names=("x", "z"))
        x, z = ring.gens()
        binary = sum(c*x**i*z**(degree-i) for i, c in enumerate(polynomial))
        k = transvectant(binary, binary, degree-2)
        h = transvectant(k, k, 2)
        names = ["T2", "T4", "T6"]+names
        weights = [2, 4, 6]+weights
        values = [field(transvectant(binary, binary, degree)),
                  field(transvectant(k, k, 4)), field(transvectant(k, h, 4))]+values
    return names, weights, values


def plane_discriminant(F):
    degree, ring, field = int(F.degree()), F.parent(), F.base_ring()
    if degree % int(field.characteristic()) == 0:
        from dissect.utils.curve_analysis import AnalysisUnavailable

        raise AnalysisUnavailable("The gradient resultant discriminant formula needs characteristic not dividing the plane degree")
    value = ring.macaulay_resultant([F.derivative(v) for v in ring.gens()])
    return value / field(degree)**(degree*degree-3*degree+3)


def artin_schreier_orbit(f, h, max_candidates):
    """Exact base-field isomorphism class for polynomial, single-pole AS models."""
    from dissect.utils.curve_analysis import AnalysisUnavailable, field_value

    if h.degree() != 0:
        raise AnalysisUnavailable("Artin-Schreier orbit classification currently requires constant nonzero h")
    ring, field = f.parent(), f.base_ring()
    q = int(field.cardinality())
    if q*(q-1) > max_candidates:
        raise AnalysisUnavailable("Artin-Schreier affine orbit exceeds max_candidates", "resource_limit")
    x = ring.gen()
    trace_one = next(c for c in field if c.trace() == 1)

    def reduced(polynomial):
        for half in range(int(polynomial.degree())//2, 0, -1):
            coefficient = polynomial[2*half]
            if coefficient:
                term = coefficient.sqrt()*x**half
                polynomial -= term**2+term
        polynomial += polynomial[0]
        return polynomial

    def normalize(polynomial):
        constant_trace = polynomial[0].trace()
        # Removing even terms can change the constant only through the
        # subsequent constant AS equation, whose quotient is F_2 via trace.
        polynomial = reduced(polynomial)
        if constant_trace:
            polynomial += trace_one
        return polynomial

    normalized = f/h[0]**2
    degree = int(normalize(normalized).degree())
    orbit = {}
    for a in field:
        if not a:
            continue
        for b in field:
            polynomial = normalize(normalized(a*x+b))
            coordinates = [field_value(polynomial[i]) for i in range(degree+1)]
            key = tuple(tuple(value) if isinstance(value, list) else (value,) for value in coordinates)
            orbit[key] = polynomial
    canonical = orbit[min(orbit)]
    return {"polynomial": canonical, "orbit_size": len(orbit),
            "affine_stabilizer_order": q*(q-1)//len(orbit)}
