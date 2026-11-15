"""Exact shared arithmetic for finite-field curves; Sage imports stay lazy."""

from fractions import Fraction


class AnalysisUnavailable(Exception):
    def __init__(self, reason, status="unsupported"):
        super().__init__(reason)
        self.status = status


def check_budget(curve, params, point_count=False):
    q, g = int(curve.q()), int(curve.genus())
    if q > params.get("max_field_size", 4096) or g > params.get("max_genus", 12):
        raise AnalysisUnavailable("Field size or genus exceeds the configured bound", "resource_limit")
    if point_count and q**g > params.get("max_point_work", 2000000):
        raise AnalysisUnavailable("q^genus exceeds max_point_work for exact point counting", "resource_limit")


def cached(curve, key, computation):
    cache = getattr(curve, "_analysis_cache", None)
    if cache is None:
        cache = curve._analysis_cache = {}
    if key not in cache:
        cache[key] = computation()
    return cache[key]


def function_field(curve):
    def build():
        if curve.model_type() == "CyclicTrigonal":
            return curve.cover_function_field()
        if curve.is_hyperelliptic():
            from sage.all import FunctionField, PolynomialRing

            K = FunctionField(curve.field(), "x")
            x = K.gen()
            ring = PolynomialRing(K, "y")
            y = ring.gen()
            f, h = curve.hyperelliptic_polynomials()
            return K.extension(y**2 + h(x)*y - f(x), "y")
        return curve.curve().function_field()
    return cached(curve, "function_field", build)


def l_polynomial(curve, params):
    check_budget(curve, params, point_count=True)
    def compute():
        if not curve.is_hyperelliptic():
            return function_field(curve).L_polynomial()
        from sage.all import PolynomialRing, ZZ
        from dissect.utils.curve_automorphisms import extension_field

        q, genus = int(curve.q()), curve.genus()
        f, h = curve.hyperelliptic_polynomials()
        traces = []
        for degree in range(1, genus+1):
            field, embedding = extension_field(curve.field(), degree)
            ring = PolynomialRing(field, "x")
            F, H = [ring([embedding(c) for c in polynomial]) for polynomial in (f, h)]

            def quadratic_roots(a, b):
                if field.characteristic() != 2:
                    discriminant = a*a+4*b
                    return 1 if not discriminant else (2 if discriminant.is_square() else 0)
                return 1 if not a else (2 if (b/a**2).trace() == 0 else 0)

            count = sum(quadratic_roots(H(x), F(x)) for x in field)
            count += quadratic_roots(H[genus+1], F[2*genus+2])
            traces.append(q**degree+1-count)
        coefficients = [1]
        for n in range(1, genus+1):
            numerator = -sum(coefficients[n-i]*traces[i-1] for i in range(1, n+1))
            coefficient, remainder = divmod(numerator, n)
            if remainder:
                raise ValueError("Point counts produced nonintegral zeta coefficients")
            coefficients.append(coefficient)
        coefficients.extend([0]*genus)
        for i in range(genus):
            coefficients[2*genus-i] = q**(genus-i)*coefficients[i]
        return PolynomialRing(ZZ, "t")(coefficients)
    return cached(curve, "l_polynomial", compute)


def field_value(value):
    """JSON-safe exact value; extension-field elements use power-basis coordinates."""
    field = value.parent()
    if field.is_prime_field():
        return int(value)
    coefficients = list(value.polynomial())
    return [int(c) for c in coefficients] + [0] * (int(field.degree()) - len(coefficients))


def field_matrix(matrix):
    return [[field_value(c) for c in row] for row in matrix.rows()]


def power_sums(characteristic_coefficients, length):
    """Newton sums of roots of a monic polynomial, coefficients ascending."""
    coefficients = list(reversed([int(c) for c in characteristic_coefficients]))
    degree = len(coefficients) - 1
    if coefficients[0] != 1:
        raise ValueError("A monic characteristic polynomial is required")
    sums = [degree]
    for n in range(1, length + 1):
        value = sum(coefficients[i]*sums[n-i] for i in range(1, min(n, degree + 1)))
        if n <= degree:
            value += n * coefficients[n]
        sums.append(-value)
    return sums


def extension_polynomial(coefficients, extension):
    """Characteristic polynomial of the extension-th powers of all roots."""
    degree = len(coefficients) - 1
    sums = power_sums(coefficients, degree * extension)
    descending = [1]
    for n in range(1, degree + 1):
        numerator = -sum(descending[n-i] * sums[i*extension] for i in range(1, n+1))
        quotient, remainder = divmod(numerator, n)
        if remainder:
            raise ValueError("Nonintegral extension polynomial")
        descending.append(quotient)
    return list(reversed(descending))


def newton_slopes(coefficients, prime, field_degree):
    """Exact normalized slopes of the lower Newton polygon of L(T)."""
    points = []
    for i, coefficient in enumerate(coefficients):
        coefficient = abs(int(coefficient))
        if not coefficient:
            continue
        valuation = 0
        while coefficient % prime == 0:
            coefficient //= prime
            valuation += 1
        point = (i, valuation)
        while len(points) >= 2:
            a, b = points[-2:]
            if (b[1]-a[1])*(i-b[0]) < (valuation-b[1])*(b[0]-a[0]):
                break
            points.pop()
        points.append(point)
    slopes = []
    for a, b in zip(points, points[1:]):
        slope = Fraction(b[1]-a[1], (b[0]-a[0])*field_degree)
        slopes.extend([slope] * (b[0]-a[0]))
    return slopes, points


def semigroup_profile(gaps):
    gaps = sorted(int(n) for n in gaps)
    genus = len(gaps)
    nongaps = [n for n in range(1, 2*genus+2) if n not in gaps]
    generators = [n for n in nongaps if not any(n-a in nongaps for a in nongaps if a < n)]
    return {
        "gaps": gaps,
        "classical_gap_weight": sum(gaps) - genus*(genus+1)//2,
        "semigroup_conductor": max(gaps, default=-1) + 1,
        "semigroup_generators": generators,
    }
