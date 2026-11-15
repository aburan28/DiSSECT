from dissect.traits.detailed_trait import DetailedTrait
from dissect.utils.curve_analysis import check_budget, AnalysisUnavailable
from dissect.utils.projective_configurations import search_field, field_description
from dissect.utils.genus4_tritangents import tritangent_configuration


class Genus4TritangentsTrait(DetailedTrait):
    NAME = "genus4_tritangents"
    DESCRIPTION = "Odd-theta tritangent planes, contact divisors and Frobenius configurations, separating even-theta pencils."
    GENUS = 4
    SUPPORTED_CURVE_TYPES = ("non_hyperelliptic",)
    INPUT = dict(DetailedTrait.INPUT, extension_degree=(int, "Constant-field search extension"), max_candidates=(int, "Dual-space enumeration budget"),
                 max_incidence_tests=(int, "Maximum four-plane incidence tests"), max_configuration_records=(int, "Maximum incidence records"),
                 max_plane_records=(int, "Maximum full plane/contact records"))
    DEFAULT_PARAMS = dict(DetailedTrait.DEFAULT_PARAMS, extension_degree=[1], max_candidates=[250000],
                          max_incidence_tests=[50000], max_configuration_records=[1000], max_plane_records=[256])
    OUTPUT = dict(DetailedTrait.OUTPUT, **{
        "search_field": (dict, "Search constant field and source embedding"),
        "canonical_quadric": (str, "Quadric used for conic plane sections"),
        "canonical_cubic": (str, "Cubic used for degree-six intersection divisors"),
        "tritangent_planes": (list, "Odd-theta planes with exact contact half-divisors"),
        "indexed_plane_coordinates": (list, "Full coordinate index for Frobenius and incidence records"),
        "tritangent_count": (int, "Odd-theta planes over the search field"),
        "expected_odd_theta_characteristics": (int, "120 in characteristic different from two"),
        "geometric_odd_locus_complete": (int, "Whether all 120 geometric odd-theta planes were found"),
        "frobenius": (dict, "Frobenius action on the found odd tritangents"),
        "even_theta_pencil_planes": (list, "Double-line cone sections belonging to the effective even theta pencil"),
        "even_theta_plane_count": (int, "q_search+1 for a cone; zero for a smooth quadric"),
        "even_theta_locus_positive_dimensional": (int, "1 for the even-theta pencil family on a cone"),
        "effective_even_theta_characteristics": (int, "One for a cone, zero for a smooth canonical quadric"),
        "concurrent_plane_quadruples": (list, "Four-plane incidence records in the found configuration"),
        "concurrent_plane_quadruple_count": (int, "Number of incidences found before the record cap"),
        "configuration_records_complete": (int, "Whether all found incidence records fit the output bound"),
        "incidence_tests": (int, "Number of quadruples examined"),
        "incidence_search_complete_on_found_planes": (int, "Completeness within the finite found set"),
        "plane_records_complete": (int, "Whether all full contact records fit the output bound"),
    })

    def compute(self, curve, params):
        check_budget(curve, params)
        if curve.field().characteristic() == 2:
            raise AnalysisUnavailable("The classical 120 odd-theta tritangents require characteristic different from two")
        field, embedding = search_field(curve, params)
        return dict(tritangent_configuration(curve, field, embedding, params), search_field=field_description(field, embedding, curve.field()))
