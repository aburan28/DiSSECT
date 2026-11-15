"""Algebraic identities, independent point counts, and actual Sage backends."""

import copy
import json
import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("sage.all")
from sage.all import GF, QQ, PolynomialRing, matrix

from dissect.traits import TRAITS
from dissect.traits.detailed_trait import DetailedTrait
from dissect.traits.jacobian_decomposition import honda_tate_factors
from dissect.utils.curve_factory import curve_from_dict
from dissect.utils.curve_invariants import dixmier_ohno_invariants, shioda_invariants
from dissect.utils.higher_genus_report import REPORT_TRAITS, curve_report
from dissect.utils.database_handler import _format_curve, _decode_ints, store_trait_result


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = json.loads((ROOT/"examples/higher_genus_curves.json").read_text())["curves"]


def hyper(f=None, h=None, prime=7):
    record = {
        "name": "hyper", "category": "test", "form": "Hyperelliptic",
        "field": {"type": "Prime", "p": prime, "bits": int(prime).bit_length()},
        "params": {"f": [0, -1, 0, 0, 0, 1] if f is None else f},
        "example": True, "standard": False,
    }
    if h is not None:
        record["params"]["h"] = h
    return record


@pytest.fixture(scope="module")
def hyper_curve():
    return curve_from_dict(hyper())


@pytest.fixture(scope="module")
def plane_curve():
    return curve_from_dict(EXAMPLES[2])


def monomial(a, b, coefficient=1):
    return {"numerator": [{"coeff": coefficient, "powers": [a, b]}]}


def assert_ok(result):
    assert result["status"] == "ok", result.get("reason")


@pytest.mark.parametrize("fixture_name,genus,order", [("hyper_curve", 2, 64), ("plane_curve", 3, 512)])
def test_frobenius_counts_match_independent_sage_curve_enumeration(request, fixture_name, genus, order):
    curve = request.getfixturevalue(fixture_name)
    result = TRAITS["frobenius_polynomial"](curve)
    assert_ok(result)
    assert result["jacobian_cardinality"] == order
    assert result["functional_equation_verified"] == 1
    assert len(result["frobenius_polynomial"]) == 2*genus+1
    # Independent model-level enumeration, rather than cached function-field places.
    field = curve.field()
    if curve.is_hyperelliptic():
        f, _ = curve.hyperelliptic_polynomials()
        independent = 1+sum(y*y == f(x) for x in field for y in field)
    else:
        independent = len(curve.curve().rational_points())
    assert result["curve_cardinality"] == independent
    assert result["curve_cardinality"] != result["jacobian_cardinality"]


def test_supplied_jacobian_order_is_checked_not_assumed():
    curve = curve_from_dict(dict(hyper(), order=3, cofactor=1))
    result = TRAITS["frobenius_polynomial"](curve)
    assert_ok(result)
    assert result["metadata_order_matches"] == 0


@pytest.mark.parametrize("fixture_name", ["hyper_curve", "plane_curve"])
def test_cartier_and_newton_arithmetic_agree(request, fixture_name):
    curve = request.getfixturevalue(fixture_name)
    cartier = TRAITS["hasse_witt"](curve)
    newton = TRAITS["newton_polygon"](curve)
    assert_ok(cartier)
    assert_ok(newton)
    assert cartier["p_rank"] == newton["p_rank"]
    assert cartier["a_number"] == curve.genus()
    assert cartier["superspecial"] == newton["supersingular"] == 1
    assert newton["slope_symmetry_verified"] == 1


def test_cartier_rank_is_not_mislabeled_p_rank():
    curve = curve_from_dict(hyper([1, 0, 0, 0, 0, 1], prime=3))
    result = TRAITS["hasse_witt"](curve)
    assert_ok(result)
    assert result["cartier_rank"] == 1
    assert result["stable_cartier_ranks"] == [1, 0]
    assert result["p_rank"] == 0
    assert result["a_number"] == 1
    assert result["superspecial"] == 0
    assert result["nilpotence_index"] == 2
    assert TRAITS["newton_polygon"](curve)["p_rank"] == 0


def test_semilinear_cartier_over_a_nonprime_field():
    record = hyper([1, 0, 0, 0, 0, 1], prime=3)
    record["field"] = {"type": "Extension", "base": "0x3", "degree": 2,
                       "poly": [{"coeff": "0x1", "power": "0x2"}, {"coeff": "0x1", "power": "0x0"}]}
    curve = curve_from_dict(record)
    result = TRAITS["hasse_witt"](curve)
    assert_ok(result)
    assert result["stable_cartier_ranks"] == [1, 0]
    assert all(isinstance(c, list) and len(c) == 2 for row in result["cartier_matrix"] for c in row)
    assert result["p_rank"] == TRAITS["newton_polygon"](curve)["p_rank"]
    from dissect.utils.curve_analysis import extension_polynomial
    base = TRAITS["frobenius_polynomial"](curve_from_dict(hyper([1, 0, 0, 0, 0, 1], prime=3)))
    extended = TRAITS["frobenius_polynomial"](curve)
    assert_ok(extended)
    assert extended["frobenius_polynomial"] == extension_polynomial(base["frobenius_polynomial"], 2)


def test_honda_tate_schur_index_prevents_false_elliptic_decompositions():
    R = PolynomialRing(QQ, "T")
    T = R.gen()
    # Over F_49, 7*zeta_3 is not an elliptic Weil number: 7 splits in Q(zeta_3).
    blocks = honda_tate_factors((T*T+7*T+49)**2, 7, 2)
    assert blocks[0]["schur_index"] == 2
    assert blocks[0]["simple_factor_dimension"] == 2
    assert blocks[0]["simple_factor_multiplicity"] == 1
    # A real rational Weil number also needs Schur index two.
    real = honda_tate_factors((T+7)**4, 7, 2)[0]
    assert real["simple_factor_dimension"] == 1
    assert real["simple_factor_multiplicity"] == 2


def test_jacobian_decompositions_are_over_the_reported_field(hyper_curve):
    result = TRAITS["jacobian_decomposition"](hyper_curve)
    assert_ok(result)
    assert result["base_elliptic_factor_dimension"] == 2
    assert result["base_simple"] == 0
    assert result["geometrically_simple"] == 0
    assert result["extension_decompositions"][0]["factors"][0]["schur_index"] == 2
    assert result["extension_search_complete"] == 0


def test_all_closed_weierstrass_places_include_nonrational_points(hyper_curve):
    result = TRAITS["weierstrass_gaps"](hyper_curve)
    assert_ok(result)
    assert result["weierstrass_closed_points"] == 5
    assert result["weierstrass_geometric_points"] == 6
    assert result["rational_weierstrass_points"] == 4
    assert [P["degree"] for P in result["weierstrass_places"]] == [1, 1, 1, 1, 2]
    assert all(P["gaps"] == [1, 3] for P in result["weierstrass_places"])
    assert all(P["semigroup_generators"] == [2, 5] for P in result["weierstrass_places"])


def test_plane_quartic_hyperflex_weights(plane_curve):
    result = TRAITS["weierstrass_gaps"](plane_curve)
    assert_ok(result)
    assert result["weierstrass_geometric_points"] == 12
    assert result["weighted_gap_excess"] == 24
    assert all(P["gaps"] == [1, 2, 5] for P in result["weierstrass_places"])


def test_nonclassical_curves_use_their_actual_generic_gaps():
    record = copy.deepcopy(EXAMPLES[2])
    record["field"] = {"type": "Prime", "p": 3}
    result = TRAITS["weierstrass_gaps"](curve_from_dict(record))
    assert_ok(result)
    assert result["generic_gaps"] == [1, 2, 4]
    assert result["classical"] == 0
    assert result["weierstrass_geometric_points"] == 28
    assert result["weighted_gap_excess"] == 28


def test_weierstrass_output_truncation_does_not_change_locus_totals(hyper_curve):
    result = TRAITS["weierstrass_gaps"](hyper_curve, max_places=2)
    assert result["status"] == "partial"
    assert len(result["weierstrass_places"]) == 2
    assert result["weierstrass_geometric_points"] == 6
    assert result["returned_locus_complete"] == 0


@pytest.mark.parametrize("fixture_name,degree,different,branch", [("hyper_curve", 2, 6, 6), ("plane_curve", 4, 12, 4)])
def test_cover_fibers_verify_local_and_global_degree_identities(request, fixture_name, degree, different, branch):
    curve = request.getfixturevalue(fixture_name)
    result = TRAITS["cover_ramification"](curve)
    assert_ok(result)
    assert result["cover_degree"] == degree
    assert result["different_degree"] == different
    assert result["branch_geometric_points"] == branch
    assert result["riemann_hurwitz_verified"] == 1
    assert result["discriminant_degree"] == different
    assert all(F["degree_sum"] == degree for F in result["branch_fibers"]+result["rational_fibers"])
    assert result["rational_curve_points_from_fibers"] == TRAITS["frobenius_polynomial"](curve)["curve_cardinality"]
    assert result["tame"] == 1


def test_tame_kummer_monodromy_is_distinct_from_rational_deck_group(plane_curve):
    result = TRAITS["cover_ramification"](plane_curve)
    assert_ok(result)
    assert result["geometric_monodromy_group"] == "C4"
    assert result["galois"] == 0
    assert result["deck_group_order"] == 2
    assert result["galois_constant_extension_degree"] == 2
    assert result["branch_splitting_extension_degree"] == 2


def test_wild_double_cover_reports_different_exponents_and_lower_breaks():
    curve = curve_from_dict(hyper([1, 0, 0, 1, 0, 1], [1], prime=2))
    result = TRAITS["cover_ramification"](curve)
    assert_ok(result)
    assert result["branch_geometric_points"] == 1
    assert result["different_degree"] == 6
    assert result["tame"] == 0
    point = result["branch_fibers"][0]["points"][0]
    assert point["ramification_index"] == 2
    assert point["different_exponent"] == 6
    assert point["lower_ramification_break"] == 5
    assert point["wild_different_excess"] == 5


def test_rational_fiber_truncation_does_not_invent_a_point_count(hyper_curve):
    result = TRAITS["cover_ramification"](hyper_curve, max_rational_fibers=2)
    assert_ok(result)
    assert len(result["rational_fibers"]) == 2
    assert result["rational_fibers_complete"] == 0
    assert result["rational_curve_points_from_fibers"] is None
    assert result["infinity_fiber"]["degree_sum"] == 2


def quotient_record():
    record = hyper([1, 0, 0, 0, 0, 0, 1])
    record["analysis"] = {"maps": [
        {"name": "elliptic", "x": monomial(2, 0), "y": monomial(0, 1),
         "target": {"form": "Hyperelliptic", "params": {"f": [1, 0, 0, 1]}}},
        {"name": "degree-four", "x": monomial(2, 0)},
        {"name": "inseparable", "x": monomial(7, 0)},
        {"name": "invalid", "x": monomial(2, 0), "y": monomial(0, 1, 2),
         "target": {"form": "Hyperelliptic", "params": {"f": [1, 0, 0, 1]}}},
    ]}
    return record


def test_explicit_elliptic_quotient_has_verified_local_ramification():
    curve = curve_from_dict(quotient_record())
    result = TRAITS["gonality_and_quotients"](curve)
    assert_ok(result)
    assert result["verified_elliptic_quotients"] == 1
    elliptic, rational, inseparable, invalid = result["declared_maps"]
    assert elliptic["degree"] == 2
    assert elliptic["target_genus"] == 1
    assert elliptic["ramification"]["different_degree"] == 2
    assert elliptic["ramification"]["riemann_hurwitz_verified"] == 1
    assert elliptic["ramification"]["branch_geometric_points"] == 2
    assert rational["degree"] == 4
    assert rational["ramification"]["different_degree"] == 10
    assert inseparable["inseparable_degree"] == 7
    assert inseparable["separable_degree"] == 2
    assert inseparable["ramification"] is None
    assert invalid["status"] == "invalid_map"


def test_verified_map_to_a_nonhyperelliptic_target(plane_curve):
    record = copy.deepcopy(EXAMPLES[2])
    x = {"numerator": [{"coeff": 1, "powers": [0, 0]}], "denominator": [{"coeff": 1, "powers": [0, 1]}]}
    y = {"numerator": [{"coeff": 1, "powers": [1, 0]}], "denominator": [{"coeff": 1, "powers": [0, 1]}]}
    record["analysis"] = {"maps": [{"name": "identity", "x": x, "y": y,
                                     "target": {"form": "SmoothPlane", "params": record["params"]}}]}
    result = TRAITS["cover_ramification"](curve_from_dict(record))
    assert_ok(result)
    identity = result["declared_maps"][0]
    assert identity["status"] == "ok"
    assert identity["degree"] == 1
    assert identity["target_genus"] == 3
    assert identity["ramification"]["different_degree"] == 0


def test_plane_gonality_has_an_explicit_realizing_cover(plane_curve):
    result = TRAITS["gonality_and_quotients"](plane_curve)
    assert_ok(result)
    assert result["geometric_gonality"] == result["base_field_gonality"] == 3
    assert result["gonality_certificate"]["center"]
    cover = TRAITS["cover_ramification"](plane_curve, projection_mode=2)
    assert_ok(cover)
    assert cover["cover_degree"] == 3
    assert cover["different_degree"] == 10
    assert cover["rational_curve_points_from_fibers"] == 8


def test_hyperelliptic_branch_stabilizers_must_lift_over_the_field(hyper_curve):
    result = TRAITS["automorphism_structure"](hyper_curve)
    assert_ok(result)
    assert result["branch_stabilizer_order"] == 8
    assert result["reduced_automorphism_order"] == 4
    assert result["base_field_order"] == 8
    assert result["center_order"] == 2
    assert result["geometric_automorphism_order"] is None
    assert any(Q["quotient_genus"] == 1 for Q in result["tame_cyclic_quotient_genera"])


def test_plane_automorphisms_are_exhaustive_via_rational_frames(plane_curve):
    result = TRAITS["automorphism_structure"](plane_curve)
    assert_ok(result)
    assert result["base_field_order"] == 24
    assert result["search_complete"] == 1
    assert sum(row["count"] for row in result["element_order_distribution"]) == 24
    partial = TRAITS["automorphism_structure"](plane_curve, max_candidates=400)
    assert partial["status"] == "partial"
    assert partial["base_field_order"] is None
    assert partial["computed_subgroup_order"] == 24


def test_characteristic_two_automorphisms_use_a_faithful_regular_action():
    curve = curve_from_dict(hyper([1, 0, 0, 1, 0, 1], [1], prime=2))
    result = TRAITS["automorphism_structure"](curve)
    assert_ok(result)
    assert result["base_field_order"] == 4
    assert result["group_structure"] == "C4"
    assert result["canonical_action_faithful"] == 0
    assert result["hyperelliptic_quotient_genus"] == 0


def test_automorphism_fields_of_definition_are_explicit():
    curve = curve_from_dict(hyper(prime=3))
    base = TRAITS["automorphism_structure"](curve)
    extension = TRAITS["automorphism_structure"](curve, extension_degree=2)
    assert_ok(base)
    assert_ok(extension)
    assert extension["search_field_size"] == 9
    assert extension["source_field_generator_image"]
    assert extension["base_field_order"] % base["base_field_order"] == 0
    assert extension["geometric_automorphism_order"] is None


def test_dixmier_ohno_generators_transform_with_their_actual_gl3_weights():
    field = GF(101)
    ring = PolynomialRing(field, names=("x", "y", "z"))
    x, y, z = ring.gens()
    F = x**4+2*y**4+3*z**4+x*x*y*z+5*x*y**3+7*x**3*z+y*z**3
    names, weights, values = dixmier_ohno_invariants(F)
    M = matrix(field, [[1, 2, 1], [0, 3, 2], [1, 0, 1]])
    substitution = ring.hom([sum(M[i, j]*ring.gen(j) for j in range(3)) for i in range(3)], ring)
    _, _, changed = dixmier_ohno_invariants(substitution(F))
    assert len(names) == 13
    assert all(a == b*M.det()**(4*weight//3) for a, b, weight in zip(changed, values, weights))
    _, _, scaled = dixmier_ohno_invariants(3*F)
    assert all(a == 3**weight*b for a, b, weight in zip(scaled, values, weights))
    assert values[-1] != 0


def test_dixmier_discriminant_normalization_on_fermat_quartic():
    ring = PolynomialRing(QQ, names=("x", "y", "z"))
    x, y, z = ring.gens()
    _, _, values = dixmier_ohno_invariants(x**4+y**4+z**4)
    assert values[0] == QQ(1)/2
    assert values[-1] == 1


def test_shioda_transvectants_transform_with_gl2_weights():
    field = GF(101)
    ring = PolynomialRing(field, names=("x", "z"))
    x, z = ring.gens()
    F = x**8+x**7*z+3*x**4*z**4+7*x*z**7+z**8
    values = shioda_invariants(F)
    substitution = ring.hom([2*x+z, 3*x+4*z], ring)
    changed = shioda_invariants(substitution(F))
    assert all(a == b*field(5)**(4*weight) for a, b, weight in zip(changed, values, range(2, 11)))


def test_moduli_reports_characteristic_limits_without_false_completeness(plane_curve):
    result = TRAITS["moduli_invariants"](plane_curve)
    assert result["status"] == "partial"
    assert result["invariant_count"] == 13
    assert result["complete_classical_system"] == 0
    record = copy.deepcopy(EXAMPLES[2])
    record["field"] = {"type": "Prime", "p": 11}
    full = TRAITS["moduli_invariants"](curve_from_dict(record))
    assert_ok(full)
    assert full["complete_classical_system"] == 1


def test_artin_schreier_moduli_remove_coboundaries_and_preserve_rational_twists():
    original = curve_from_dict(hyper([1, 0, 0, 1, 0, 1], [1], prime=2))
    equivalent = curve_from_dict(hyper([1, 1, 0, 1, 1, 1], [1], prime=2))
    twist = curve_from_dict(hyper([0, 0, 0, 1, 0, 1], [1], prime=2))
    result = TRAITS["moduli_invariants"](original)
    same = TRAITS["moduli_invariants"](equivalent)
    other = TRAITS["moduli_invariants"](twist)
    assert_ok(result)
    assert_ok(same)
    assert_ok(other)
    assert result["values"] == same["values"]
    assert result["values"] != other["values"]
    assert result["classification_complete_in_stated_scope"] == 1
    assert result["complete_classical_system"] is None
    assert 2*result["artin_schreier_affine_stabilizer_order"] == TRAITS["automorphism_structure"](original)["base_field_order"]
    extension = {"type": "Binary", "degree": 2,
                 "poly": [{"coeff": "0x1", "power": "0x2"}, {"coeff": "0x1", "power": "0x1"},
                          {"coeff": "0x1", "power": "0x0"}]}
    first = hyper([1, 0, 0, 1, 0, 1], [1], prime=2)
    second = hyper([0, 0, 0, 1, 0, 1], [1], prime=2)
    first["field"] = second["field"] = extension
    assert TRAITS["moduli_invariants"](curve_from_dict(first))["values"] == TRAITS["moduli_invariants"](curve_from_dict(second))["values"]


def test_higher_genus_moduli_return_valid_subsets_with_partial_status():
    record = hyper([1, 1, 0, 0, 0, 0, 0, 0, 0, 1], prime=11)
    result = TRAITS["moduli_invariants"](curve_from_dict(record))
    assert result["status"] == "partial"
    assert result["invariant_count"] == 4
    assert result["weights"] == [2, 4, 6, 18]


@pytest.mark.parametrize("coefficients,prime", [([1, 0, 0, 0, 1], 3), ([1, 0, 0, 0, 0, 0, 1], 7)])
def test_cyclic_trigonal_models_reject_inseparable_or_nonsquarefree_input(coefficients, prime):
    record = hyper(coefficients, prime=prime)
    record["form"] = "CyclicTrigonal"
    # Second polynomial is squarefree in characteristic seven; make a repeated-root polynomial.
    if prime == 7:
        record["params"]["f"] = [0, 0, 0, 0, 1]
    with pytest.raises(ValueError):
        curve_from_dict(record)


def trigonal_record(degree=5):
    record = hyper([0, -1]+[0]*(degree-2)+[1])
    record["form"] = "CyclicTrigonal"
    return record


def test_nonplane_genus_four_trigonal_cover_has_exact_arithmetic_and_gonality():
    record = trigonal_record()
    report = curve_report(record)
    assert report["genus"] == 4
    assert report["family"] == "non_hyperelliptic"
    cover = report["traits"]["cover_ramification"]
    gonality = report["traits"]["gonality_and_quotients"]
    assert_ok(cover)
    assert_ok(gonality)
    assert cover["cover_degree"] == gonality["geometric_gonality"] == 3
    assert cover["different_degree"] == 12
    assert cover["geometric_monodromy_group"] == "C3"
    assert cover["galois"] == 1
    assert cover["rational_curve_points_from_fibers"] == report["traits"]["frobenius_polynomial"]["curve_cardinality"]
    assert report["traits"]["automorphism_structure"]["cover_stabilizer_complete"] == 1
    # The genus-four canonical-quadric trait now proves that this is a cone,
    # so its unique ruling pencil makes the cover stabilizer the full group.
    assert report["traits"]["automorphism_structure"]["search_complete"] == 1


def test_high_genus_trigonal_pencil_uniqueness_certifies_full_automorphism_search():
    curve = curve_from_dict(trigonal_record(7))
    result = TRAITS["automorphism_structure"](curve)
    assert curve.genus() == 6
    assert_ok(result)
    assert result["cover_stabilizer_complete"] == result["search_complete"] == 1
    assert result["base_field_order"] == 42


def test_explicit_map_to_a_trigonal_target_is_checked():
    record = trigonal_record()
    record["analysis"] = {"maps": [{"name": "identity", "x": monomial(1, 0), "y": monomial(0, 1),
                                     "target": {"form": "CyclicTrigonal", "params": record["params"]}}]}
    result = TRAITS["cover_ramification"](curve_from_dict(record))
    assert_ok(result)
    identity = result["declared_maps"][0]
    assert identity["status"] == "ok"
    assert identity["target_genus"] == 4
    assert identity["degree"] == 1
    assert identity["ramification"]["different_degree"] == 0


def test_tame_trigonal_covers_are_supported_in_characteristic_two():
    record = hyper([1, 0, 1, 0, 0, 1], prime=2)
    record["form"] = "CyclicTrigonal"
    curve = curve_from_dict(record)
    result = TRAITS["cover_ramification"](curve)
    assert_ok(result)
    assert curve.genus() == 4
    assert result["tame"] == 1
    assert result["cover_degree"] == 3
    assert result["different_degree"] == 12
    assert result["galois"] == 0
    assert result["galois_constant_extension_degree"] == 2


@pytest.mark.parametrize("analysis", [[], {"maps": {}}, {"maps": ["bad"]}])
def test_malformed_map_containers_are_rejected(analysis):
    with pytest.raises(ValueError):
        curve_from_dict(dict(hyper(), analysis=analysis))


def test_large_map_exponents_have_a_separate_resource_status():
    record = dict(hyper(), analysis={"maps": [{"name": "oversized", "x": monomial(257, 0)}]})
    result = TRAITS["cover_ramification"](curve_from_dict(record))
    assert_ok(result)
    assert result["declared_maps"][0]["status"] == "resource_limit"


def test_reports_keep_all_trait_statuses_and_database_maps_round_trip():
    report = curve_report(quotient_record())
    assert set(report["traits"]) == set(REPORT_TRAITS)
    assert report["trait_parameters"]["frobenius_polynomial"]["timeout"] == 30
    assert sum(report["coverage"].values()) == len(REPORT_TRAITS)
    record = quotient_record()
    restored = _decode_ints(_format_curve(record))
    assert restored["analysis"] == record["analysis"]
    assert curve_from_dict(restored).analysis()["maps"][0]["name"] == "elliptic"


def test_nested_trait_results_round_trip_database_encoding(hyper_curve):
    result = TRAITS["jacobian_decomposition"](hyper_curve)

    class Collection:
        def insert_one(self, record):
            self.result = record
            return type("Ack", (), {"acknowledged": True})()

    collection = Collection()
    assert store_trait_result({"trait_jacobian_decomposition": collection}, hyper_curve, "jacobian_decomposition", {}, result)
    decoded = _decode_ints(collection.result)["result"]
    assert decoded["base_factors"][0]["simple_factor_dimension"] == 1
    assert decoded["status"] == "ok"


def test_timeouts_are_explicit_and_alarm_is_cancelled(hyper_curve):
    class Slow(DetailedTrait):
        def compute(self, curve, params):
            while True:
                pass

    result = Slow()(hyper_curve, timeout=1)
    assert result["status"] == "timeout"
    assert_ok(TRAITS["hasse_witt"](hyper_curve))


def test_resource_limits_are_distinct_from_mathematical_zero(hyper_curve):
    result = TRAITS["frobenius_polynomial"](hyper_curve, max_point_work=1)
    assert result["status"] == "resource_limit"
    assert result["jacobian_cardinality"] is None
    assert result["curve_cardinality"] is None


@pytest.mark.parametrize("module,arguments", [
    ("dissect.utils.compute_json", ["-t", "cover_ramification", "--params", '{"max_rational_fibers":2}']),
    ("dissect.utils.higher_genus_report", ["-t", "frobenius_polynomial", "hasse_witt", "cover_ramification"]),
])
def test_json_entrypoints_accept_mixed_families_and_budget_overrides(module, arguments):
    completed = subprocess.run([sys.executable, "-m", module, "-i", str(ROOT/"examples/higher_genus_curves.json")]+arguments,
                               check=True, capture_output=True, text=True)
    result = json.loads(completed.stdout)
    assert len(result["data"]) == 3
    if "report" in module:
        assert result["schema_version"] == 1
        assert result["software"]["sage"]
        assert result["data"][0]["traits"]["hasse_witt"]["status"] == "not_applicable"
        assert result["data"][1]["traits"]["hasse_witt"]["status"] == "ok"
    else:
        assert len(result["data"][1]["result"]["rational_fibers"]) == 2
