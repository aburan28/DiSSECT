import copy
import json
import subprocess
import sys
from itertools import product
from pathlib import Path

import pytest

pytest.importorskip("sage.all")
from sage.all import GF, PolynomialRing, matrix, vector

from dissect.traits import TRAITS
from dissect.utils.curve_factory import curve_from_dict
from dissect.utils.curve_analysis import function_field
from dissect.utils.canonical_genus4 import canonical_data
from dissect.utils.curve_automorphisms import extension_field
from dissect.utils.projective_configurations import monomial_powers, substitute, normalize
from dissect.utils.higher_genus_report import curve_report
from dissect.utils.database_handler import _format_curve, _decode_ints, store_trait_result


ROOT = Path(__file__).resolve().parents[1]
BASE = json.loads((ROOT/"examples/higher_genus_curves.json").read_text())["curves"]


def assert_ok(result):
    assert result["status"] == "ok", result.get("reason")


def decode(value, field):
    if isinstance(value, list):
        return sum(field(c)*field.gen()**i for i, c in enumerate(value))
    return field(value)


def binary_form(values, ring, degree):
    return sum(decode(c, ring.base_ring())*ring.monomial(*powers)
               for c, powers in zip(values, monomial_powers(2, degree)))


def hyper(f=None, prime=7):
    result = copy.deepcopy(BASE[1])
    result["field"] = {"type": "Prime", "p": prime, "bits": int(prime).bit_length()}
    if f is not None:
        result["params"]["f"] = f
    return result


def trigonal(prime=7):
    result = hyper(prime=prime)
    result["form"] = "CyclicTrigonal"
    return result


def canonical_record(nonsplit=False, degenerate_contact=False):
    ring = PolynomialRing(GF(7), names=("X0", "X1", "X2", "X3"))
    a, b, c, d = ring.gens()
    if nonsplit:
        Q = a*d-b*b-c*c
        C = a**3+b**3+2*c**3+3*d**3+a*b*c+b*c*d
    elif degenerate_contact:
        Q = a*d-b*c
        C = a*a*(b+c)-2*a*(b*b+c*c)+b**3+c**3+a*b*c+d*(a*a+b*b+c*c)+2*d**3
    else:
        Q = a*d-b*c
        C = a**3+b**3+c**3+2*d**3
    terms = lambda P: [{"coeff": int(coefficient), "powers": list(power)} for power, coefficient in P.dict().items()]
    return {"name": "canonical-example", "category": "genus_specific_examples", "example": True, "standard": False,
            "form": "CanonicalGenus4", "field": {"type": "Prime", "p": 7, "bits": 3},
            "params": {"quadric": terms(Q), "cubic": terms(C)}}


@pytest.fixture(scope="module")
def g2_curve():
    return curve_from_dict(hyper())


@pytest.fixture(scope="module")
def quartic_curve():
    return curve_from_dict(BASE[2])


@pytest.fixture(scope="module")
def bitangent_result(quartic_curve):
    return TRAITS["genus3_bitangents"](quartic_curve, extension_degree=2, timeout=120)


@pytest.fixture(scope="module")
def flex_result(quartic_curve):
    return TRAITS["genus3_flexes"](quartic_curve, extension_degree=2)


@pytest.fixture(scope="module")
def cone_curve():
    return curve_from_dict(trigonal())


@pytest.fixture(scope="module")
def split_curve():
    return curve_from_dict(canonical_record())


@pytest.fixture(scope="module")
def nonsplit_curve():
    return curve_from_dict(canonical_record(nonsplit=True))


@pytest.fixture(scope="module")
def complete_tritangents():
    return TRAITS["genus4_tritangents"](curve_from_dict(trigonal(prime=5)), extension_degree=2,
                                       max_incidence_tests=10, timeout=120)


def test_rosenhain_full_orbit_and_frobenius(g2_curve):
    result = TRAITS["genus2_rosenhain"](g2_curve)
    assert_ok(result)
    assert result["branch_orbit_degrees"] == [1, 1, 1, 1, 2]
    assert result["branch_splitting_extension_degree"] == 2
    assert result["branch_frobenius"]["orbit_degrees"] == [1, 1, 1, 1, 2]
    assert result["rational_two_torsion_rank"] == 3
    assert result["ordered_normalizations_tested"] == 120
    assert result["configuration_count"] == 5
    assert result["configuration_complete"] == 1


def test_rosenhain_two_torsion_matches_independent_subset_enumeration(g2_curve):
    result = TRAITS["genus2_rosenhain"](g2_curve)
    permutation = result["branch_frobenius"]["permutation"]
    classes = set()
    fixed = set()
    for bits in product((0, 1), repeat=6):
        if sum(bits) % 2:
            continue
        key = min(bits, tuple(1-b for b in bits))
        classes.add(key)
        image = [0]*6
        for i, j in enumerate(permutation):
            image[j] = bits[i]
        if min(tuple(image), tuple(1-b for b in image)) == key:
            fixed.add(key)
    assert len(classes) == 16
    assert len(fixed) == 2**result["rational_two_torsion_rank"]


def test_rosenhain_configuration_is_invariant_under_translation_and_twist(g2_curve):
    field = GF(7)
    ring = PolynomialRing(field, "x")
    x = ring.gen()
    original = x**5-x
    transformed = 3*original(x+2)
    other = curve_from_dict(hyper([int(c) for c in transformed]))
    first = TRAITS["genus2_rosenhain"](g2_curve)
    second = TRAITS["genus2_rosenhain"](other)
    assert_ok(second)
    assert first["rosenhain_configurations"] == second["rosenhain_configurations"]


def test_rosenhain_partial_result_preserves_exact_base_arithmetic(g2_curve):
    result = TRAITS["genus2_rosenhain"](g2_curve, max_field_size=7)
    assert result["status"] == "partial"
    assert result["rational_two_torsion_rank"] == 3
    assert result["branch_splitting_extension_degree"] == 2
    assert result["canonical_configuration"] is None


def test_rosenhain_rational_branch_points_give_full_rational_two_torsion():
    ring = PolynomialRing(GF(7), "x")
    x = ring.gen()
    f = ring.one()
    for i in range(5):
        f *= x-i
    result = TRAITS["genus2_rosenhain"](curve_from_dict(hyper([int(c) for c in f])))
    assert_ok(result)
    assert result["rosenhain_model_over_base"] == 1
    assert result["rational_two_torsion_rank"] == 4
    assert result["branch_splitting_extension_degree"] == 1


def test_classical_rosenhain_excludes_characteristic_two():
    record = hyper([1, 0, 0, 1, 0, 1], prime=2)
    record["params"]["h"] = [1]
    assert TRAITS["genus2_rosenhain"](curve_from_dict(record))["status"] == "unsupported"


def test_elliptic_quotients_exhaust_nonhyperelliptic_involutions_and_are_reusable(g2_curve):
    result = TRAITS["genus2_elliptic_quotients"](g2_curve)
    assert_ok(result)
    group = TRAITS["automorphism_structure"](g2_curve)
    order_two = next(row["count"] for row in group["element_order_distribution"] if row["order"] == 2)
    assert result["degree_two_quotient_count"] == order_two-1 == 4
    assert result["degree_two_search_complete"] == 1
    assert result["all_elliptic_quotients_complete"] == 0
    record = hyper()
    record["analysis"] = {"maps": [{"name": "found-%d" % i, "x": Q["x"], "y": Q["y"], "target": Q["target"]}
                                    for i, Q in enumerate(result["quotients"])]}
    checked = TRAITS["cover_ramification"](curve_from_dict(record))
    assert_ok(checked)
    assert all(Q["status"] == "ok" and Q["target_genus"] == 1 and Q["degree"] == 2 for Q in checked["declared_maps"])
    assert all(Q["ramification"]["different_degree"] == 2 for Q in result["quotients"])


def test_bitangent_complete_configuration_has_classical_counts(bitangent_result):
    assert_ok(bitangent_result)
    assert bitangent_result["bitangent_count"] == 28
    assert bitangent_result["hyperflex_bitangent_count"] == 12
    assert bitangent_result["syzygetic_tetrad_count"] == 315
    assert bitangent_result["geometric_locus_complete"] == 1
    assert bitangent_result["frobenius"]["rational_over_base"] == 4
    assert sum(bitangent_result["frobenius"]["orbit_degrees"]) == 28


def test_bitangent_contact_certificates_reconstruct_the_exact_sections(bitangent_result):
    field, embedding = extension_field(GF(7), 2)
    R = PolynomialRing(field, names=("x", "y", "z"))
    x, y, z = R.gens()
    F = x**4+y**4+z**4
    B = PolynomialRing(field, names=("s", "t"))
    s, t = B.gens()
    for record in bitangent_result["bitangents"]:
        parameterization = [[decode(c, field) for c in row] for row in record["line_parameterization"]]
        images = [row[0]*s+row[1]*t for row in parameterization]
        contact = binary_form(record["contact_quadratic"], B, 2)
        assert substitute(F, images, B) == decode(record["section_scalar"], field)*contact**2
        line = [decode(c, field) for c in record["line"]]
        assert sum(c*value for c, value in zip(line, images)) == 0


def test_syzygetic_tetrad_conics_contain_the_four_half_divisors(bitangent_result):
    field, _ = extension_field(GF(7), 2)
    ring = PolynomialRing(field, names=("x", "y", "z"))
    binary = PolynomialRing(field, names=("s", "t"))
    s, t = binary.gens()
    powers = bitangent_result["conic_coefficient_monomials"]
    for tetrad in bitangent_result["syzygetic_tetrads"][:10]:
        conic = sum(decode(c, field)*ring.monomial(*power) for c, power in zip(tetrad["conic_coefficients"], powers))
        for index in tetrad["bitangent_indices"]:
            record = bitangent_result["bitangents"][index]
            images = [decode(row[0], field)*s+decode(row[1], field)*t for row in record["line_parameterization"]]
            restricted = substitute(conic, images, binary)
            contact = binary_form(record["contact_quadratic"], binary, 2)
            rows = [[P.monomial_coefficient(binary.monomial(*power)) for power in monomial_powers(2, 2)] for P in (restricted, contact)]
            assert matrix(field, rows).rank() <= 1


def test_bitangent_base_field_search_is_labeled_partial(quartic_curve):
    result = TRAITS["genus3_bitangents"](quartic_curve)
    assert result["status"] == "partial"
    assert result["bitangent_count"] == 4
    assert result["geometric_locus_complete"] == 0


def test_flexes_include_explicit_contact_and_frobenius(flex_result):
    assert_ok(flex_result)
    assert flex_result["generic_tangent_contact"] == 2
    assert flex_result["found_flex_count"] == flex_result["found_hyperflex_count"] == 12
    assert flex_result["geometric_special_locus_complete"] == 1
    assert flex_result["frobenius"]["orbit_degrees"] == [2]*6
    assert all(P["tangent_intersection_multiplicity"] == 4 for P in flex_result["flexes"])
    assert flex_result["collinear_point_triples"]


def test_nonclassical_flexes_do_not_force_characteristic_zero_counts():
    record = copy.deepcopy(BASE[2])
    record["field"] = {"type": "Prime", "p": 3}
    result = TRAITS["genus3_flexes"](curve_from_dict(record), extension_degree=2)
    assert_ok(result)
    assert result["generic_gaps"] == [1, 2, 4]
    assert result["generic_tangent_contact"] == 3
    assert result["ordinary_flex_locus_finite"] == 0
    assert result["found_hyperflex_count"] == result["found_special_flexes"] == 28


def test_flex_record_cap_keeps_a_full_index(quartic_curve):
    result = TRAITS["genus3_flexes"](quartic_curve, extension_degree=2, max_point_records=2)
    assert result["status"] == "partial"
    assert len(result["flexes"]) == 2
    assert len(result["indexed_point_coordinates"]) == len(result["frobenius"]["permutation"]) == 12


@pytest.mark.parametrize("fixture_name,kind,rank,geometric,rational", [
    ("cone_curve", "cone", 3, 1, 1), ("split_curve", "smooth_split", 4, 2, 2),
    ("nonsplit_curve", "smooth_nonsplit", 4, 2, 0),
])
def test_genus_four_quadric_types_and_pencil_descent(request, fixture_name, kind, rank, geometric, rational):
    curve = request.getfixturevalue(fixture_name)
    result = TRAITS["genus4_canonical_quadric"](curve)
    assert_ok(result)
    assert result["quadric_type"] == kind
    assert result["quadric_rank"] == rank
    assert result["geometric_trigonal_pencils"] == geometric
    assert result["base_field_trigonal_pencils"] == rational
    assert result["effective_even_theta_characteristics"] == int(rank == 3)
    assert result["pencil_construction_complete"] == 1
    if rational:
        assert all(P["pole_degree_verified"] == 1 for P in result["pencils"])
    else:
        assert result["pencil_frobenius_permutation"] == [1, 0]
        assert result["pencil_field"]["size"] == 49


def test_canonical_relations_are_actual_source_field_identities(cone_curve):
    data = canonical_data(cone_curve)
    L = function_field(cone_curve)
    assert substitute(data["quadric"], data["coordinates"], L) == 0
    assert substitute(data["cubic"], data["coordinates"], L) == 0
    assert len(data["coordinates"]) == 4


def test_automorphism_exchanging_pencils_has_a_verifiable_witness(split_curve):
    result = TRAITS["genus4_canonical_quadric"](split_curve)
    assert_ok(result)
    assert result["automorphisms_exchange_pencils"] == 1
    witness = matrix(GF(7), result["pencil_exchange_witness"])
    assert witness.det()
    Q, C = split_curve.canonical_polynomials()
    ring = Q.parent()
    images = [sum(witness[i, j]*ring.gen(j) for j in range(4)) for i in range(4)]
    transformed_Q = substitute(Q, images, ring)
    assert ring.ideal(Q).reduce(transformed_Q) == 0
    assert ring.ideal(Q, C).reduce(substitute(C, images, ring)) == 0


def test_pencil_exchange_decision_does_not_depend_on_record_truncation(split_curve):
    result = TRAITS["genus4_canonical_quadric"](split_curve, max_transformations=1)
    assert_ok(result)
    assert len(result["pencil_automorphism_actions"]) == 1
    assert result["automorphisms_exchange_pencils"] == 1
    assert result["pencil_exchange_witness"]


def test_unexhausted_automorphism_search_does_not_claim_pencils_are_never_exchanged(split_curve):
    result = TRAITS["genus4_canonical_quadric"](split_curve, max_candidates=1)
    assert_ok(result)
    assert result["automorphism_search_complete"] == 0
    assert result["automorphisms_exchange_pencils"] is None


def test_no_rational_trigonal_pencil_does_not_imply_base_field_gonality_three(nonsplit_curve):
    result = TRAITS["gonality_and_quotients"](nonsplit_curve)
    assert_ok(result)
    assert result["geometric_gonality"] == 3
    assert result["base_field_lower_bound"] == 4
    assert result["base_field_gonality"] != 3


def test_complete_tritangent_configuration_excludes_even_theta_pencil(complete_tritangents):
    assert_ok(complete_tritangents)
    assert complete_tritangents["tritangent_count"] == 120
    assert complete_tritangents["geometric_odd_locus_complete"] == 1
    assert complete_tritangents["even_theta_plane_count"] == 26
    assert complete_tritangents["even_theta_locus_positive_dimensional"] == 1
    assert complete_tritangents["effective_even_theta_characteristics"] == 1
    assert complete_tritangents["frobenius"]["orbit_degrees"] == [2]*60
    assert all(P["kind"] == "odd_theta" for P in complete_tritangents["tritangent_planes"])


def test_tritangent_certificates_have_even_intersection_and_lie_on_the_plane(complete_tritangents):
    curve = curve_from_dict(trigonal(prime=5))
    data = canonical_data(curve)
    field, _ = extension_field(GF(5), 2)
    ring = PolynomialRing(field, names=data["quadric"].parent().variable_names())
    Q, C = [ring({p: field(c) for p, c in F.dict().items()}) for F in (data["quadric"], data["cubic"])]
    binary = PolynomialRing(field, names=("s", "t"))
    for record in complete_tritangents["tritangent_planes"]:
        images = [binary_form(row, binary, 2) for row in record["parameterization"]]
        contact = binary_form(record["contact_cubic"], binary, 3)
        assert substitute(Q, images, binary) == 0
        assert substitute(C, images, binary) == decode(record["section_scalar"], field)*contact**2
        H = [decode(c, field) for c in record["plane"]]
        assert sum(c*value for c, value in zip(H, images)) == 0


def test_smooth_quadric_tritangents_do_not_inherit_an_even_pencil(nonsplit_curve):
    result = TRAITS["genus4_tritangents"](nonsplit_curve)
    assert result["status"] == "partial"
    assert result["tritangent_count"] == 2
    assert result["even_theta_plane_count"] == 0
    assert result["even_theta_locus_positive_dimensional"] == 0


def test_degenerate_conic_sections_can_still_be_odd_tritangents():
    curve = curve_from_dict(canonical_record(degenerate_contact=True))
    result = TRAITS["genus4_tritangents"](curve)
    assert result["status"] == "partial"
    contact = next(record for record in result["tritangent_planes"] if record["plane"] == [0, 0, 0, 1])
    assert contact["kind"] == "odd_theta"
    assert contact["conic_type"] == "two_lines"
    assert contact["vertex_half_divisor_multiplicity"] == 1
    assert [piece["vertex_intersection_multiplicity"] for piece in contact["component_contacts"]] == [1, 1]


def test_canonical_input_works_with_shared_arithmetic(split_curve):
    result = TRAITS["frobenius_polynomial"](split_curve)
    assert_ok(result)
    assert result["curve_cardinality"] == len(split_curve.curve().rational_points()) == 6
    assert function_field(split_curve).genus() == 4
    assert_ok(TRAITS["hasse_witt"](split_curve))


def test_invalid_canonical_quadric_cubic_intersections_are_rejected():
    record = canonical_record()
    record["params"]["cubic"] = [{"coeff": 1, "powers": [3, 0, 0, 0]}]
    with pytest.raises(ValueError, match="smooth"):
        curve_from_dict(record)
    record = canonical_record()
    record["params"]["quadric"] = [{"coeff": 1, "powers": [2, 0, 0, 0]}]
    with pytest.raises(ValueError, match="rank"):
        curve_from_dict(record)


def test_new_configuration_results_survive_database_encoding(quartic_curve):
    result = TRAITS["genus3_bitangents"](quartic_curve)

    class Collection:
        def insert_one(self, record):
            self.record = record
            return type("Ack", (), {"acknowledged": True})()

    collection = Collection()
    assert store_trait_result({"trait_genus3_bitangents": collection}, quartic_curve, "genus3_bitangents", {}, result)
    decoded = _decode_ints(collection.record)["result"]
    assert decoded["bitangent_count"] == 4
    assert decoded["frobenius"]["orbit_degrees"] == [1]*4


def test_combined_report_selects_genus_specific_traits():
    report = curve_report(hyper(), traits=("genus2_rosenhain", "genus3_bitangents", "genus4_canonical_quadric"))
    assert report["traits"]["genus2_rosenhain"]["status"] == "ok"
    assert report["traits"]["genus3_bitangents"]["status"] == "not_applicable"
    assert report["traits"]["genus4_canonical_quadric"]["status"] == "not_applicable"


def test_genus_specific_json_cli_and_parameter_override():
    completed = subprocess.run([sys.executable, "-m", "dissect.utils.compute_json", "-t", "genus2_rosenhain",
                                "-i", str(ROOT/"examples/higher_genus_curves.json"), "--params", '{"max_field_size":7}'],
                               check=True, capture_output=True, text=True)
    rows = json.loads(completed.stdout)["data"]
    assert rows[0]["result"]["status"] == "not_applicable"
    assert rows[1]["result"]["status"] == "partial"
