from dissect.traits.detailed_trait import DetailedTrait
from dissect.utils.curve_analysis import check_budget
from dissect.utils.projective_configurations import search_field, field_description
from dissect.utils.plane_quartic_configurations import flexes


class Genus3FlexesTrait(DetailedTrait):
    NAME = "genus3_flexes"
    DESCRIPTION = "Explicit flexes, hyperflexes, tangent contact orders, Frobenius action and incidences."
    GENUS = 3
    MODEL_TYPES = ("SmoothPlane",)
    SUPPORTED_CURVE_TYPES = ("non_hyperelliptic",)
    INPUT = dict(DetailedTrait.INPUT, extension_degree=(int, "Constant-field search extension"), max_candidates=(int, "Projective-plane budget"),
                 max_incidence_tests=(int, "Maximum point/tangent triple tests"), max_configuration_records=(int, "Maximum incidence records"),
                 max_point_records=(int, "Maximum explicit flex records"))
    DEFAULT_PARAMS = dict(DetailedTrait.DEFAULT_PARAMS, extension_degree=[1], max_candidates=[100000],
                          max_incidence_tests=[50000], max_configuration_records=[1000], max_point_records=[256])
    OUTPUT = dict(DetailedTrait.OUTPUT, **{
        "search_field": (dict, "Search field basis and source embedding"),
        "generic_gaps": (list, "Actual generic canonical gap sequence"),
        "generic_tangent_contact": (int, "Generic tangent intersection multiplicity"),
        "ordinary_flex_locus_finite": (int, "0 when every geometric point is a flex on a nonclassical quartic"),
        "flexes": (list, "Points, tangent lines, exact sections and contact multiplicities"),
        "indexed_point_coordinates": (list, "Full found point index for Frobenius and incidence records"),
        "indexed_tangent_coordinates": (list, "Tangent coordinates in the same index"),
        "found_flex_count": (int, "Flexes over the search field"),
        "found_hyperflex_count": (int, "Tangent multiplicity four"),
        "expected_geometric_special_flexes": (int, "Exact geometric Weierstrass-place count"),
        "found_special_flexes": (int, "Found points with contact exceeding the actual generic contact"),
        "geometric_special_locus_complete": (int, "All finite special flexes found"),
        "frobenius": (dict, "Frobenius permutation on the found flex points"),
        "collinear_point_triples": (list, "Collinear triples of recorded flex points"),
        "concurrent_tangent_triples": (list, "Concurrent triples of tangent lines"),
        "incidence_tests": (int, "Number of triples examined"),
        "incidence_search_complete_on_found_points": (int, "Whether all found-point triples were examined"),
        "point_records_complete": (int, "Whether the full found point set fits the output bound"),
        "configuration_records_complete": (int, "Whether all found incidence records fit the output bound"),
    })

    def compute(self, curve, params):
        check_budget(curve, params)
        field, embedding = search_field(curve, params)
        return dict(flexes(curve, field, embedding, params), search_field=field_description(field, embedding, curve.field()))
