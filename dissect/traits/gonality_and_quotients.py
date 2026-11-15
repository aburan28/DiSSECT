from dissect.traits.detailed_trait import DetailedTrait
from dissect.utils.curve_analysis import check_budget, function_field
from dissect.utils.cover_analysis import declared_maps, plane_gonality_cover


class GonalityAndQuotientsTrait(DetailedTrait):
    NAME = "gonality_and_quotients"
    DESCRIPTION = "Proved geometric gonality, base-field bounds and verified low-genus quotient maps."
    INPUT = dict(DetailedTrait.INPUT, max_maps=(int, "Maximum declared maps"), max_map_degree=(int, "Maximum declared-map degree and input exponent"))
    DEFAULT_PARAMS = dict(DetailedTrait.DEFAULT_PARAMS, max_maps=[16], max_map_degree=[256])
    OUTPUT = dict(DetailedTrait.OUTPUT, **{
        "geometric_gonality": (int, "Exact gonality over the algebraic closure"),
        "base_field_gonality": (int, "Exact when a map realizes the geometric lower bound"),
        "base_field_lower_bound": (int, "Geometric gonality is a lower bound over the base field"),
        "base_field_upper_bound": (int, "Degree of a verified base-field rational function"),
        "gonality_certificate": (dict, "Theorem and an explicit function-field model of a realizing cover"),
        "rational_point_found": (int, "Whether a rational plane point was found"),
        "declared_maps": (list, "Map verifications, including rational and elliptic quotients"),
        "verified_elliptic_quotients": (int, "Number of supplied, verified separable genus-one quotient maps"),
        "quotient_search_complete": (int, "0: supplied maps do not exhaust all quotients"),
    })

    def compute(self, curve, params):
        check_budget(curve, params)
        L = function_field(curve)
        upper = int(L.degree())
        if curve.is_hyperelliptic():
            gonality, rational = 2, None
            certificate = {"theorem": "Unique hyperelliptic double cover", "equation": str(L.polynomial())}
        elif curve.model_type() == "CyclicTrigonal":
            gonality, rational = 3, None
            certificate = {"theorem": "A degree-three cover in genus >= 3 cannot be hyperelliptic by Castelnuovo-Severi",
                           "equation": str(L.polynomial())}
        else:
            degree = int(curve.plane_polynomial().degree())
            gonality = degree-1
            projection = plane_gonality_cover(curve)
            rational = int(projection is not None)
            certificate = {"theorem": "Smooth plane curve of degree d >= 4 has geometric gonality d-1"}
            if projection:
                upper = gonality
                certificate.update(center=projection["center"], coordinate_change=projection["coordinate_change"],
                                   equation=str(projection["field"].polynomial()))
        maps = declared_maps(curve, params)
        for morphism in maps:
            if morphism["status"] == "ok" and morphism["target_genus"] == 0:
                upper = min(upper, morphism["degree"])
        if upper < gonality:
            raise ValueError("A supplied map contradicts the proved geometric gonality")
        return {
            "geometric_gonality": gonality,
            "base_field_gonality": gonality if upper == gonality else None,
            "base_field_lower_bound": gonality, "base_field_upper_bound": upper,
            "gonality_certificate": certificate, "rational_point_found": rational,
            "declared_maps": maps,
            "verified_elliptic_quotients": sum(m["status"] == "ok" and m.get("target_genus") == 1 and m.get("separable") for m in maps),
            "quotient_search_complete": 0,
        }
