import math

from dissect.traits import Trait


def normalized_trace(trace, q):
    """Frobenius trace scaled into $[-1, 1]$ by the Hasse bound.

    Hasse gives $|t| \\le 2\\sqrt{q}$, so $\\tau = t / 2\\sqrt{q}$ is comparable
    across field sizes. Integer square root keeps the scaling exact for the
    several-hundred-bit values involved; the result is clamped because the floor
    can put a curve sitting exactly on the bound a hair outside it.
    """
    tau = trace / (2 * math.isqrt(q))
    return max(-1.0, min(1.0, tau))


class SatoTateTrait(Trait):
    NAME = "sato_tate"
    DESCRIPTION = (
        "Frobenius trace normalized by the Hasse bound, $\\tau = t / 2\\sqrt{q}$. Over a fixed"
        " field with the coefficients varying, Birch's theorem makes $\\tau$ semicircular on"
        " $[-1, 1]$, so unlike the raw trace this has a null distribution to be read against."
    )
    INPUT = {}
    OUTPUT = {
        "normalized_trace": (float, "$t / 2\\sqrt{q}$, in $[-1, 1]$"),
    }
    DEFAULT_PARAMS = {}

    def compute(self, curve, params):
        """Normalized trace.

        Reported alone rather than alongside the Sato-Tate angle
        $\\theta = \\arccos\\tau$: the angle is a monotone transform carrying no
        extra information, and a second copy of one quantity would simply weight
        the trace twice in any feature vector built from these traits.
        """
        return {"normalized_trace": normalized_trace(int(curve.trace()), int(curve.q()))}


def test_zero_trace():
    assert normalized_trace(0, 2**256) == 0.0


def test_hasse_bound_endpoints():
    q = 2**256
    assert normalized_trace(2 * math.isqrt(q), q) == 1.0
    assert normalized_trace(-2 * math.isqrt(q), q) == -1.0
    # Beyond the bound is impossible for a real curve, but the clamp must hold.
    assert normalized_trace(3 * math.isqrt(q), q) == 1.0


def test_scale_invariance():
    # The same relative position in the Hasse interval scores the same at any size.
    small = normalized_trace(math.isqrt(2**64), 2**64)
    large = normalized_trace(math.isqrt(2**512), 2**512)
    assert math.isclose(small, 0.5, abs_tol=1e-9)
    assert math.isclose(large, 0.5, abs_tol=1e-9)


def test_brainpoolp256r1():
    trace = 300418416528525664980082381967979838673
    p = 76884956397045344220809746629001649093037950200943055203735601445031516197751
    assert math.isclose(normalized_trace(trace, p), 0.5417212486209377, abs_tol=1e-12)


def test_sign_tracks_cardinality_against_field():
    # #E = q + 1 - t, so a positive trace is exactly the curve whose cardinality
    # falls below q. Brainpool requires that; X9.62 does not.
    assert normalized_trace(1000, 2**128) > 0
    assert normalized_trace(-1000, 2**128) < 0
