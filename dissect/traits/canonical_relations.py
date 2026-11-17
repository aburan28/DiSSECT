from dissect.traits.detailed_trait import DetailedTrait
from dissect.utils.curve_analysis import check_budget
from dissect.utils.canonical_relations import canonical_relations


class CanonicalRelationsTrait(DetailedTrait):
    NAME = "canonical_relations"
    DESCRIPTION = "Canonical multiplication maps, minimal image equations and linear syzygies across genera."
    INPUT = dict(DetailedTrait.INPUT, relation_degree=(int, "Largest relation degree"), max_monomials=(int, "Largest monomial basis"),
                 max_relation_records=(int, "Maximum equation and syzygy records"))
    DEFAULT_PARAMS = dict(DetailedTrait.DEFAULT_PARAMS, max_genus=[6], relation_degree=[4], max_monomials=[1000], max_relation_records=[128])
    OUTPUT = dict(DetailedTrait.OUTPUT, **{
        "coefficient_field": (dict, "Coefficient field and exact power-basis convention"),
        "coordinate_basis": (list, "Regular differential basis or validated input canonical coordinates"),
        "coordinate_functions": (list, "Canonical coordinate ratios in the source function field"),
        "coordinate_origin": (str, "Construction of the canonical coordinates"), "variables": (list, "Polynomial variable order"),
        "degree_profiles": (list, "Exact multiplication ranks, relation dimensions and new generators by degree"),
        "minimal_generators": (list, "Minimal canonical-image equations with exact structured coefficients"),
        "minimal_generator_count": (int, "Number of generators through the searched degree"),
        "quadric_generator_indices": (list, "Indices of quadrics referred to by syzygies"),
        "indexed_quadric_terms": (list, "Full ordered quadric equations retained when verbose generator records are capped"),
        "linear_quadric_syzygies": (list, "Linear-form coefficients giving zero combinations of the quadrics"),
        "linear_quadric_syzygy_count": (int, "Dimension of the linear syzygy space"),
        "relation_degree_bound": (int, "Largest checked degree"), "equations_verified_in_source": (int, "All equations and syzygies verified exactly"),
        "records_complete": (int, "Whether all computed equation and syzygy records fit the output cap"),
        "canonical_ring_generated_in_degree_one": (int, "Canonical image ring is the full canonical ring only in the nonhyperelliptic case"),
        "ideal_generation_complete": (int, "All degrees required by canonical-ideal generation theorems were searched"),
    })

    def compute(self, curve, params):
        check_budget(curve, params)
        return canonical_relations(curve, params)
