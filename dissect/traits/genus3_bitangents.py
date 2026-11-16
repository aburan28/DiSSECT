from dissect.traits.detailed_trait import DetailedTrait
from dissect.utils.curve_analysis import check_budget, AnalysisUnavailable
from dissect.utils.projective_configurations import search_field, field_description
from dissect.utils.plane_quartic_configurations import bitangents


class Genus3BitangentsTrait(DetailedTrait):
    NAME = "genus3_bitangents"
    DESCRIPTION = "Explicit bitangents, contact divisors, Frobenius action and syzygetic contact conics."
    GENUS = 3
    MODEL_TYPES = ("SmoothPlane",)
    SUPPORTED_CURVE_TYPES = ("non_hyperelliptic",)
    INPUT = dict(DetailedTrait.INPUT, extension_degree=(int, "Constant-field search extension"), max_candidates=(int, "Dual-plane search budget"),
                 max_tetrad_tests=(int, "Maximum four-line conic tests"), max_configuration_records=(int, "Maximum tetrad records"))
    DEFAULT_PARAMS = dict(DetailedTrait.DEFAULT_PARAMS, extension_degree=[1], max_candidates=[100000],
                          max_tetrad_tests=[25000], max_configuration_records=[1000])
    OUTPUT = dict(DetailedTrait.OUTPUT, **{
        "search_field": (dict, "Search field basis and source embedding"),
        "bitangents": (list, "Line coordinates, parameterizations and exact half-intersection divisors"),
        "bitangent_count": (int, "Number of bitangents over the search field"),
        "expected_geometric_bitangents": (int, "28 odd theta characteristics in characteristic different from two"),
        "geometric_locus_complete": (int, "Whether all 28 geometric bitangents were found"),
        "frobenius": (dict, "Frobenius permutation on the found lines"),
        "hyperflex_bitangent_count": (int, "Bitangents with one contact point of multiplicity four"),
        "conic_coefficient_monomials": (list, "Basis for contact conic equations"),
        "syzygetic_tetrads": (list, "Four bitangents whose contact divisors lie on the returned conic"),
        "syzygetic_tetrad_count": (int, "Number found; 315 when the full configuration is searched"),
        "tetrad_tests": (int, "Number of four-line tests performed"),
        "tetrad_search_complete_on_found_lines": (int, "Whether all tetrads of the found set were tested"),
        "configuration_records_complete": (int, "Whether all conic records fit the output bound"),
    })

    def compute(self, curve, params):
        check_budget(curve, params)
        if curve.field().characteristic() == 2:
            raise AnalysisUnavailable("The classical 28-bitangent configuration requires characteristic different from two")
        field, embedding = search_field(curve, params)
        return dict(bitangents(curve, field, embedding, params), search_field=field_description(field, embedding, curve.field()))
