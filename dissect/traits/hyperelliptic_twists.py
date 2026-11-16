from math import gcd

from dissect.traits.detailed_trait import DetailedTrait
from dissect.utils.curve_analysis import l_polynomial, extension_polynomial, field_value, AnalysisUnavailable
from dissect.utils.projective_configurations import field_description


class HyperellipticTwistsTrait(DetailedTrait):
    NAME = "hyperelliptic_twists"
    DESCRIPTION = "Explicit involution twists, curve and Jacobian orders, and extension-parity identities."
    SUPPORTED_CURVE_TYPES = ("hyperelliptic",)
    INPUT = dict(DetailedTrait.INPUT, max_point_work=(int, "Point-count budget"), extension_bound=(int, "Largest twist extension"))
    DEFAULT_PARAMS = dict(DetailedTrait.DEFAULT_PARAMS, max_point_work=[2000000], extension_bound=[6])
    OUTPUT = dict(DetailedTrait.OUTPUT, **{
        "field": (dict, "Coefficient field and basis"), "twist_kind": (str, "Quadratic or Artin-Schreier involution twist"),
        "twist_parameter": (object, "Nonsquare or absolute-trace-one constant"),
        "twisted_f": (list, "Ascending coefficients of the explicit twisted hyperelliptic equation"),
        "twisted_h": (list, "Ascending coefficients of its h polynomial"),
        "source_l_polynomial": (list, "Original zeta numerator"), "twist_l_polynomial": (list, "L(-T), not a genus-one trace formula"),
        "source_jacobian_order": (int, "Original L(1)"), "twist_jacobian_order": (int, "Twisted L(-1)"),
        "jacobian_order_gcd": (int, "Common divisor of the two base-field Jacobian orders"),
        "extension_comparisons": (list, "Curve/Jacobian orders and equality over even extensions"),
        "same_geometric_branch_divisor": (int, "The deck-involution twist preserves branch geometry"),
        "all_twists_classified": (int, "0: extra automorphisms may give other twists"),
    })

    def compute(self, curve, params):
        if params["extension_bound"] > 64:
            raise AnalysisUnavailable("extension_bound exceeds 64", "resource_limit")
        field, q = curve.field(), int(curve.q())
        coefficients = [int(c) for c in l_polynomial(curve, params)]
        twisted = [(-1)**i*c for i, c in enumerate(coefficients)]
        f, h = curve.hyperelliptic_polynomials()
        if field.characteristic() == 2:
            parameter = next(field.gen()**i for i in range(int(field.degree())) if (field.gen()**i).trace() == 1)
            F, H, kind = f+parameter*h*h, h, "Artin-Schreier constant twist"
        else:
            parameter = next(value for value in field if value and not value.is_square())
            F, H, kind = parameter*(f+h*h/4), h.parent().zero(), "quadratic deck-involution twist"
        comparisons = []
        for n in range(1, params["extension_bound"]+1):
            original = extension_polynomial(list(reversed(coefficients)), n)
            other = extension_polynomial(list(reversed(twisted)), n)
            source_points, twist_points = q**n+1+original[-2], q**n+1+other[-2]
            if source_points+twist_points != 2*(q**n+1) and n % 2:
                raise ValueError("The odd-extension twist trace identity failed")
            if n % 2 == 0 and original != other:
                raise ValueError("The quadratic twist did not split over an even extension")
            comparisons.append({"extension_degree": n, "source_curve_points": source_points, "twist_curve_points": twist_points,
                                "source_jacobian_order": sum(original), "twist_jacobian_order": sum(other),
                                "twist_parameter_trivialized": int(n % 2 == 0), "frobenius_polynomials_equal": int(original == other)})
        return {"field": field_description(field, lambda value: value, field), "twist_kind": kind,
                "twist_parameter": field_value(parameter), "twisted_f": [field_value(c) for c in F], "twisted_h": [field_value(c) for c in H],
                "source_l_polynomial": coefficients, "twist_l_polynomial": twisted,
                "source_jacobian_order": sum(coefficients), "twist_jacobian_order": sum(twisted),
                "jacobian_order_gcd": gcd(sum(coefficients), sum(twisted)), "extension_comparisons": comparisons,
                "same_geometric_branch_divisor": 1, "all_twists_classified": 0}
