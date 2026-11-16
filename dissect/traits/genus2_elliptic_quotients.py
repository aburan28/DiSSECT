from dissect.traits.detailed_trait import DetailedTrait
from dissect.utils.curve_analysis import check_budget, AnalysisUnavailable
from dissect.utils.projective_configurations import search_field, field_description
from dissect.utils.genus2_quotients import discover_elliptic_quotients


class Genus2EllipticQuotientsTrait(DetailedTrait):
    NAME = "genus2_elliptic_quotients"
    DESCRIPTION = "Systematic discovery of all degree-two elliptic quotient maps over a specified finite field."
    GENUS = 2
    SUPPORTED_CURVE_TYPES = ("hyperelliptic",)
    INPUT = dict(DetailedTrait.INPUT, extension_degree=(int, "Constant-field search extension"), max_candidates=(int, "PGL2 search budget"))
    DEFAULT_PARAMS = dict(DetailedTrait.DEFAULT_PARAMS, extension_degree=[1], max_candidates=[100000])
    OUTPUT = dict(DetailedTrait.OUTPUT, **{
        "search_field": (dict, "Constant field and embedding in which maps are defined"),
        "quotients": (list, "Verified involutions, explicit maps, target models and local different divisors"),
        "degree_two_quotient_count": (int, "Number of nonhyperelliptic rational involution quotients"),
        "distinct_j_invariant_count": (int, "Distinct geometric elliptic j-invariants among computed maps"),
        "degree_two_search_complete": (int, "Whether all degree-two maps over this search field were considered"),
        "all_elliptic_quotients_complete": (int, "0: degree greater than two is not searched"),
    })

    def compute(self, curve, params):
        check_budget(curve, params)
        if curve.field().characteristic() == 2:
            raise AnalysisUnavailable("This involution quotient construction requires odd characteristic")
        field, embedding = search_field(curve, params)
        quotients = discover_elliptic_quotients(curve, field, embedding, params)
        keys = {str(record["j_invariant"]) for record in quotients if record["j_invariant"] is not None}
        return {"search_field": field_description(field, embedding, curve.field()),
                "quotients": quotients, "degree_two_quotient_count": len(quotients),
                "distinct_j_invariant_count": len(keys) if all(Q["j_invariant"] is not None for Q in quotients) else None, "degree_two_search_complete": 1,
                "all_elliptic_quotients_complete": 0}
