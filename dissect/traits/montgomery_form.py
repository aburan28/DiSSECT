from dissect.traits import Trait


def is_quadratic_residue(value, p):
    """Euler's criterion, for odd `p`. Zero is not counted as a residue here."""
    value %= p
    if value == 0:
        return False
    return pow(value, (p - 1) // 2, p) == 1


def admissible_roots(roots, a, p):
    """Roots of $x^3 + ax + b$ that yield a Montgomery form.

    A short Weierstrass curve over $\\mathbb{F}_p$, $p > 3$, is birationally
    equivalent to a Montgomery curve exactly when its cubic has a root
    $\\alpha$ for which $3\\alpha^2 + a$ is a nonzero square: that square root
    is the scaling $B$ taking the curve to $By^2 = x^3 + Ax^2 + x$.
    """
    return sum(1 for alpha in roots if is_quadratic_residue(3 * alpha * alpha + a, p))


class MontgomeryFormTrait(Trait):
    NAME = "montgomery_form"
    DESCRIPTION = (
        "Whether the curve is birationally equivalent to a Montgomery curve, and so also to a"
        " twisted Edwards curve. This is what admits an x-only Montgomery ladder, and it forces"
        " the order to be divisible by 4."
    )
    INPUT = {}
    OUTPUT = {
        "montgomery": (int, "1 when a Montgomery form exists, 0 otherwise"),
        "admissible_roots": (int, "Roots of the cubic giving such a form"),
    }
    DEFAULT_PARAMS = {}

    def compute(self, curve, params):
        """Montgomery representability over a prime field.

        Defined for $p > 3$ only: the criterion is stated for short Weierstrass
        form, which binary and small-characteristic fields do not share. A
        twisted Edwards form exists under exactly the same condition, being
        birationally equivalent, so it is not reported separately.
        """
        from sage.all import PolynomialRing

        field = curve.field()
        p = int(field.characteristic())
        if not curve.is_over_prime() or p <= 3:
            return {"montgomery": None, "admissible_roots": None}

        a, b = int(curve.a()) % p, int(curve.b()) % p
        x = PolynomialRing(field, "x").gen()
        roots = [int(root) for root, _ in (x**3 + a * x + b).roots()]
        count = admissible_roots(roots, a, p)
        return {"montgomery": 1 if count else 0, "admissible_roots": count}


def _brute_force_roots(a, b, p):
    return [x for x in range(p) if (x * x * x + a * x + b) % p == 0]


def _brute_force_cardinality(a, b, p):
    """Point count by trying every affine point, plus the point at infinity."""
    squares = {}
    for y in range(p):
        squares[(y * y) % p] = squares.get((y * y) % p, 0) + 1
    return 1 + sum(squares.get((x * x * x + a * x + b) % p, 0) for x in range(p))


def test_is_quadratic_residue():
    # Squares mod 11 are {1, 3, 4, 5, 9}.
    assert {v for v in range(1, 11) if is_quadratic_residue(v, 11)} == {1, 3, 4, 5, 9}
    assert not is_quadratic_residue(0, 11)


def test_admissibility_implies_order_divisible_by_four():
    """The condition is only worth trusting if it agrees with an independent count."""
    p = 101
    checked = 0
    for a in range(p):
        for b in range(p):
            if (4 * a**3 + 27 * b**2) % p == 0:
                continue  # singular
            roots = _brute_force_roots(a, b, p)
            if admissible_roots(roots, a, p):
                assert _brute_force_cardinality(a, b, p) % 4 == 0
                checked += 1
    assert checked > 100, checked


def test_criterion_does_not_accept_everything():
    p = 101
    negatives = sum(
        1
        for a in range(p)
        for b in range(p)
        if (4 * a**3 + 27 * b**2) % p
        and not admissible_roots(_brute_force_roots(a, b, p), a, p)
    )
    assert negatives > 100, negatives


def test_prime_order_curve_is_never_montgomery():
    p = 101
    for a in range(p):
        for b in range(p):
            if (4 * a**3 + 27 * b**2) % p == 0:
                continue
            if _brute_force_cardinality(a, b, p) % 4 != 0:
                assert admissible_roots(_brute_force_roots(a, b, p), a, p) == 0
