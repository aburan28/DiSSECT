from typing import List

from dissect.traits.detailed_trait import DetailedTrait
from dissect.utils.curve_analysis import l_polynomial, power_sums, extension_polynomial


class FrobeniusPolynomialTrait(DetailedTrait):
    NAME = "frobenius_polynomial"
    DESCRIPTION = "Exact Weil polynomial, zeta numerator, curve counts and Jacobian orders over extensions."
    INPUT = dict(DetailedTrait.INPUT, max_point_work=(int, "Bound on q^genus"), extension_degree=(int, "Largest reported constant extension"))
    DEFAULT_PARAMS = dict(DetailedTrait.DEFAULT_PARAMS, max_point_work=[2000000], extension_degree=[3])
    OUTPUT = dict(DetailedTrait.OUTPUT, **{
        "l_polynomial": (List[int], "Zeta numerator coefficients in ascending powers"),
        "frobenius_polynomial": (List[int], "Monic Frobenius polynomial coefficients in ascending powers"),
        "frobenius_trace": (int, "Trace on H^1"),
        "curve_cardinality": (int, "Number of curve points over the base field"),
        "jacobian_cardinality": (int, "L(1), not the curve point count"),
        "curve_extension_counts": (List[int], "Curve counts for extensions 1 through extension_degree"),
        "jacobian_extension_orders": (List[int], "Jacobian orders for the same extensions"),
        "functional_equation_verified": (int, "Verification of the reciprocal coefficient symmetry"),
        "metadata_order_matches": (int, "Whether supplied order*cofactor matches L(1); null if absent"),
    })

    def compute(self, curve, params):
        from dissect.utils.curve_analysis import AnalysisUnavailable

        if params["extension_degree"] > 64:
            raise AnalysisUnavailable("extension_degree exceeds 64", "resource_limit")
        L = [int(c) for c in l_polynomial(curve, params)]
        P = list(reversed(L))
        q, g = int(curve.q()), curve.genus()
        count = params["extension_degree"]
        traces = power_sums(P, count)
        jacobian = sum(L)
        metadata = curve.order()
        return {
            "l_polynomial": L, "frobenius_polynomial": P,
            "frobenius_trace": -L[1], "curve_cardinality": q+1+L[1],
            "jacobian_cardinality": jacobian,
            "curve_extension_counts": [q**n+1-traces[n] for n in range(1, count+1)],
            "jacobian_extension_orders": [sum(extension_polynomial(P, n)) for n in range(1, count+1)],
            "functional_equation_verified": int(len(L) == 2*g+1 and all(L[2*g-i] == q**(g-i)*L[i] for i in range(g+1))),
            "metadata_order_matches": None if metadata is None else int(metadata*curve.cofactor() == jacobian),
        }
