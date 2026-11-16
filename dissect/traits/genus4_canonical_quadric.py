from dissect.traits.detailed_trait import DetailedTrait
from dissect.utils.curve_analysis import check_budget
from dissect.utils.canonical_genus4 import quadric_profile, canonical_automorphisms


class Genus4CanonicalQuadricTrait(DetailedTrait):
    NAME = "genus4_canonical_quadric"
    DESCRIPTION = "Canonical quadric type, explicit trigonal pencils, Frobenius descent and even theta-null."
    GENUS = 4
    SUPPORTED_CURVE_TYPES = ("non_hyperelliptic",)
    INPUT = dict(DetailedTrait.INPUT, max_candidates=(int, "Canonical automorphism frame budget"), max_transformations=(int, "Maximum automorphism action records"))
    DEFAULT_PARAMS = dict(DetailedTrait.DEFAULT_PARAMS, max_candidates=[100000], max_transformations=[128])
    OUTPUT = dict(DetailedTrait.OUTPUT, **{
        "canonical_quadric": (str, "Unique quadric containing the canonical curve"),
        "canonical_cubic": (str, "Cubic relation independent of the quadric multiples"),
        "canonical_coordinate_basis": (list, "Input canonical coordinates or source regular differential basis"),
        "canonical_model_origin": (str, "Validated input or recovered regular-differential relations"),
        "quadric_matrix": (list, "Symmetric quadratic-form matrix in odd characteristic"),
        "quadric_rank": (int, "Rank three for a cone and four for a smooth quadric"),
        "quadric_type": (str, "cone, smooth_split or smooth_nonsplit"),
        "quadric_determinant": (object, "Field value determining ruling splitting in rank four"),
        "geometric_trigonal_pencils": (int, "One for a cone; two for a smooth quadric"),
        "base_field_trigonal_pencils": (int, "Number of pencils defined over the original field"),
        "pencil_splitting_extension_degree": (int, "One or two according to ruling descent"),
        "pencil_frobenius_permutation": (list, "Frobenius preserves or exchanges the two geometric pencils"),
        "effective_even_theta_characteristics": (int, "One precisely when the canonical quadric is a cone"),
        "pencil_field": (dict, "Field and embedding of the explicit ruling maps"),
        "ruling_coordinate_change": (list, "Columns give standard ruling coordinates"),
        "standard_quadric": (str, "Standard Segre or cone equation after the coordinate change"),
        "pencils": (list, "Linear-form ratios defining degree-three maps"),
        "pencil_construction_complete": (int, "Whether all geometric pencils have explicit certificates"),
        "automorphisms_exchange_pencils": (int, "1 if witnessed; 0 if excluded; null if the search cannot decide"),
        "pencil_automorphism_actions": (list, "Verified projective automorphisms and their pencil permutations"),
        "automorphism_search_complete": (int, "Completeness of the finite-field automorphism search used here"),
        "pencil_exchange_witness": (list, "A verified projective automorphism exchanging the two pencils"),
    })

    def compute(self, curve, params):
        check_budget(curve, params)
        result = quadric_profile(curve, params)
        if result["quadric_rank"] == 3:
            result.update(automorphisms_exchange_pencils=0, pencil_automorphism_actions=[], automorphism_search_complete=None)
        else:
            group = canonical_automorphisms(curve, curve.field(), lambda c: c, params)
            exchange = [1, 0] in group["pencil_permutations_seen"]
            result.update(automorphisms_exchange_pencils=1 if exchange else (0 if group["search_complete"] else None),
                          pencil_automorphism_actions=group["transformations"], automorphism_search_complete=group["search_complete"],
                          pencil_exchange_witness=group["pencil_exchange_witness"])
        return result
