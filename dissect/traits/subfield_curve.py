import math

from dissect.traits import Trait


def subfield_degree(degree, in_subfield):
    """Smallest $d \\mid m$ whose subfield contains both curve coefficients.

    `in_subfield(d)` reports whether both coefficients lie in the subfield of
    degree `d`. Returns `degree` itself when no proper subfield qualifies.
    """
    for d in range(1, degree + 1):
        if degree % d == 0 and in_subfield(d):
            return d
    return degree


def rho_speedup_bits(relative_degree):
    """Security margin lost to the Frobenius endomorphism, in bits.

    Frobenius has order `relative_degree` on a subfield curve, enlarging the
    Pollard rho equivalence classes by that factor and cutting the expected
    number of steps by its square root.
    """
    return 0.5 * math.log2(relative_degree)


class SubfieldCurveTrait(Trait):
    NAME = "subfield_curve"
    DESCRIPTION = (
        "Smallest subfield containing both curve coefficients. When that subfield is proper,"
        " Frobenius acts as an endomorphism and Pollard rho gains a factor of $\\sqrt{m/d}$;"
        " the Koblitz curves are the case $d = 1$."
    )
    INPUT = {}
    OUTPUT = {
        "subfield_degree": (int, "Degree $d$ of the smallest subfield containing $a$ and $b$"),
        "relative_degree": (int, "$m/d$: 1 for a generic curve, $m$ for a Koblitz curve"),
        "rho_speedup_bits": (float, "$\\log_2\\sqrt{m/d}$, the margin lost to Frobenius"),
    }
    DEFAULT_PARAMS = {}

    def compute(self, curve, params):
        """Locate the smallest subfield of definition for the curve coefficients."""
        field = curve.field()
        degree = int(field.degree())
        if degree == 1:
            # A prime field has no proper subfield to descend to.
            return {"subfield_degree": 1, "relative_degree": 1, "rho_speedup_bits": 0.0}

        characteristic = int(field.characteristic())
        a, b = field(curve.a()), field(curve.b())

        def in_subfield(d):
            power = characteristic**d
            return a**power == a and b**power == b

        d = subfield_degree(degree, in_subfield)
        relative = degree // d
        return {
            "subfield_degree": d,
            "relative_degree": relative,
            "rho_speedup_bits": rho_speedup_bits(relative),
        }


def test_subfield_degree_generic_curve():
    # Coefficients in no proper subfield: the smallest is the whole field.
    assert subfield_degree(163, lambda d: d == 163) == 163


def test_subfield_degree_koblitz():
    # Coefficients in the prime field, as for the K- and sect*k1 curves.
    assert subfield_degree(163, lambda d: True) == 1


def test_subfield_degree_intermediate():
    # A curve defined over F_(2^4) but presented over F_(2^12).
    assert subfield_degree(12, lambda d: d in (4, 6, 12)) == 4


def test_subfield_degree_only_considers_divisors():
    seen = []

    def in_subfield(d):
        seen.append(d)
        return d == 12

    assert subfield_degree(12, in_subfield) == 12
    assert seen == [1, 2, 3, 4, 6, 12]


def test_rho_speedup_bits():
    assert rho_speedup_bits(1) == 0.0
    assert math.isclose(rho_speedup_bits(163), 3.6743, abs_tol=1e-4)
    assert math.isclose(rho_speedup_bits(571), 4.5787, abs_tol=1e-4)
