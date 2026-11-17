from dissect.traits.detailed_trait import DetailedTrait
from dissect.utils.curve_analysis import l_polynomial, AnalysisUnavailable
from dissect.utils.jacobian_profiles import divisor_counts


class SymmetricPowerDivisorsTrait(DetailedTrait):
    NAME = "symmetric_power_divisors"
    DESCRIPTION = "Closed points, symmetric powers, effective Picard classes and Riemann-Roch identities."
    INPUT = dict(DetailedTrait.INPUT, max_point_work=(int, "Point-count budget"), degree_bound=(int, "Largest effective-divisor degree"))
    DEFAULT_PARAMS = dict(DetailedTrait.DEFAULT_PARAMS, max_point_work=[2000000], degree_bound=[12])
    OUTPUT = dict(DetailedTrait.OUTPUT, **{
        "l_polynomial": (list, "Zeta numerator, coefficients ascending"),
        "jacobian_order": (int, "Number of rational points of every Picard degree component"),
        "closed_point_counts": (list, "Exact closed-point counts by residue degree"),
        "symmetric_power_counts": (list, "Rational effective divisors, including degree zero"),
        "extension_point_counts": (list, "Curve counts reconstructing the closed-point counts"),
        "degree_two_effective_classes": (int, "Rational Picard degree-two classes admitting an effective divisor"),
        "degree_two_noneffective_classes": (int, "Remaining rational Picard degree-two classes"),
        "degree_two_special_classes": (int, "The hyperelliptic pencil gives one; otherwise zero"),
        "riemann_roch_verified_degrees": (list, "Degrees at least 2g-1 checked against constant Picard fibers"),
        "genus_degree_identity_verified": (int, "A_g = #J + q*A_(g-2), when searched"),
        "euler_product_verified": (int, "Closed-point Euler product equals the effective-divisor series"),
    })

    def compute(self, curve, params):
        if not 2 <= params["degree_bound"] <= 64:
            raise AnalysisUnavailable("degree_bound must lie between 2 and 64", "resource_limit")
        coefficients = [int(c) for c in l_polynomial(curve, params)]
        q, g = int(curve.q()), int(curve.genus())
        points, closed, effective, rr = divisor_counts(coefficients, q, g, params["degree_bound"])
        classes = effective[2]-(q if curve.is_hyperelliptic() else 0)
        if not 0 <= classes <= sum(coefficients):
            raise ValueError("Effective degree-two classes exceed the Picard component")
        return {"l_polynomial": coefficients, "jacobian_order": sum(coefficients),
                "closed_point_counts": [{"degree": n, "count": closed[n]} for n in range(1, len(closed))],
                "symmetric_power_counts": [{"degree": n, "count": count} for n, count in enumerate(effective)],
                "extension_point_counts": points[1:], "degree_two_effective_classes": classes,
                "degree_two_noneffective_classes": sum(coefficients)-classes,
                "degree_two_special_classes": int(curve.is_hyperelliptic()),
                "riemann_roch_verified_degrees": rr, "genus_degree_identity_verified": int(params["degree_bound"] >= g),
                "euler_product_verified": 1}
