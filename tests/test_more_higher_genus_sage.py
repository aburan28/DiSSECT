import copy
import json
import subprocess
import sys
from itertools import combinations, product
from math import prod
from pathlib import Path

import pytest

pytest.importorskip("sage.all")
from sage.all import GF, PolynomialRing, matrix, identity_matrix, companion_matrix

from dissect.traits import TRAITS
from dissect.utils.curve_factory import curve_from_dict
from dissect.utils.curve_analysis import function_field, l_polynomial
from dissect.utils.curve_automorphisms import extension_field
from dissect.utils.higher_genus_report import curve_report
from dissect.utils.jacobian_profiles import torsion_constraints, polynomial_power_mod
from dissect.utils.projective_configurations import substitute
from dissect.utils.richelot import pair_partitions
from dissect.utils.database_handler import store_trait_result, _decode_ints


ROOT = Path(__file__).resolve().parents[1]
BASE = json.loads((ROOT/"examples/higher_genus_curves.json").read_text())["curves"]
MORE = ("genus2_richelot", "hyperelliptic_twists", "jacobian_torsion", "symmetric_power_divisors", "frobenius_orders", "canonical_relations")


def hyper(coefficients=None, prime=7, h=None):
    record = copy.deepcopy(BASE[1])
    record["field"] = {"type": "Prime", "p": prime, "bits": prime.bit_length()}
    if coefficients is not None:
        record["params"] = {"f": coefficients, "h": h or []}
    return record


def decode(value, field):
    return sum(field(c)*field.gen()**i for i,c in enumerate(value)) if isinstance(value, list) else field(value)


def assert_ok(result):
    assert result["status"] == "ok", result.get("reason")


@pytest.fixture(scope="module")
def g2():
    return curve_from_dict(hyper())


@pytest.fixture(scope="module")
def ordinary():
    return curve_from_dict(hyper([0,2,0,1,0,1], prime=5))


@pytest.fixture(scope="module")
def quartic():
    return curve_from_dict(BASE[2])


@pytest.fixture(scope="module")
def split_richelot():
    ring = PolynomialRing(GF(11), "x")
    x = ring.gen()
    polynomial = prod(x-a for a in (0,1,2,3,5))
    return TRAITS["genus2_richelot"](curve_from_dict(hyper([int(c) for c in polynomial], prime=11)))


def test_richelot_partitions_exhaust_the_fifteen_lagrangians():
    partitions = list(pair_partitions(tuple(range(6))))
    assert len(partitions) == len(set(partitions)) == 15
    for partition in partitions:
        assert sorted(i for pair in partition for i in pair) == list(range(6))
        # Disjoint even branch pairs are isotropic in the even-subset model.
        assert all(len(set(a) & set(b)) % 2 == 0 for a,b in combinations(partition, 2))


def test_richelot_all_split_branch_kernels_are_rational_and_targets_isogenous(split_richelot):
    result = split_richelot
    assert_ok(result)
    assert result["base_field_kernel_count"] == 15
    assert result["genus_two_target_count"]+result["split_product_kernel_count"] == 15
    assert result["point_count_verified_targets"] == result["genus_two_target_count"] > 0
    assert all(row.get("source_target_weil_polynomial_verified") == 1 for row in result["kernel_partitions"] if row["determinant"])
    assert result["targets_complete"] == 1


def test_richelot_bilinear_certificates_and_source_factorization(split_richelot):
    field = GF(11)
    ring = PolynomialRing(field, names=("u", "v"))
    u,v = ring.gens()
    source = prod(u-a for a in (0,1,2,3,5))
    for row in split_richelot["kernel_partitions"]:
        factors = [sum(field(c)*u**i for i,c in enumerate(P)) for P in row["quadratic_factors"]]
        brackets = [sum(field(c)*v**i for i,c in enumerate(P)) for P in row["brackets"]]
        assert prod(factors) == source
        assert sum(a*b for a,b in zip(factors,brackets)) == -field(row["determinant"])*(u-v)**2


def test_richelot_splits_products_and_tracks_nonrational_kernels(g2):
    result = TRAITS["genus2_richelot"](g2)
    assert_ok(result)
    assert result["base_field_kernel_count"] == 3
    assert result["genus_two_target_count"] == 9
    assert result["split_product_kernel_count"] == 6
    assert sum(result["kernel_frobenius"]["orbit_degrees"]) == 15
    assert all(row["target_components_computed"] == 0 for row in result["kernel_partitions"] if not row["determinant"])


def test_richelot_field_limit_preserves_exact_kernel_frobenius(g2):
    result = TRAITS["genus2_richelot"](g2, max_field_size=7)
    assert result["status"] == "partial"
    assert result["base_field_kernel_count"] == 3
    assert len(result["kernel_partitions"]) == 15
    assert result["target_field"] is None


def test_richelot_handles_nonzero_h_and_a_nonsquare_source_twist():
    ring = PolynomialRing(GF(11),"x")
    x = ring.gen()
    h = x*x+1
    f = 2*prod(x-a for a in (0,1,2,3,5))-h*h/4
    record = hyper([int(c) for c in f],prime=11,h=[int(c) for c in h])
    result = TRAITS["genus2_richelot"](curve_from_dict(record))
    assert_ok(result)
    assert result["point_count_verified_targets"] == result["genus_two_target_count"]


@pytest.mark.parametrize("prime,coefficients,h", [(7,[0,-1,0,0,0,1],[]), (2,[1,0,0,1,0,1],[1])])
def test_explicit_involution_twists_match_independent_weil_polynomials(prime,coefficients,h):
    curve = curve_from_dict(hyper(coefficients, prime, h))
    result = TRAITS["hyperelliptic_twists"](curve)
    assert_ok(result)
    twisted = curve_from_dict(hyper(result["twisted_f"], prime, result["twisted_h"]))
    assert [int(c) for c in l_polynomial(twisted, {"max_point_work":2000000})] == result["twist_l_polynomial"]
    assert result["all_twists_classified"] == 0
    for row in result["extension_comparisons"]:
        if row["extension_degree"] % 2 == 0:
            assert row["source_jacobian_order"] == row["twist_jacobian_order"]
            assert row["twist_parameter_trivialized"] == 1
        else:
            assert row["source_curve_points"]+row["twist_curve_points"] == 2*(prime**row["extension_degree"]+1)


def test_characteristic_two_richelot_and_etale_two_torsion_are_not_misapplied():
    curve = curve_from_dict(hyper([1,0,0,1,0,1], 2, [1]))
    assert TRAITS["genus2_richelot"](curve)["status"] == "unsupported"
    assert TRAITS["jacobian_torsion"](curve)["status"] == "unsupported"


def test_artin_schreier_twist_uses_trace_one_in_an_even_degree_constant_field():
    record = hyper([1,0,0,1,0,1],2,[1])
    record["field"] = {"type":"Binary","degree":2,"bits":2,
                       "poly":[{"coeff":1,"power":2},{"coeff":1,"power":1},{"coeff":1,"power":0}]}
    curve = curve_from_dict(record)
    result = TRAITS["hyperelliptic_twists"](curve)
    assert_ok(result)
    assert curve.field().one().trace() == 0
    assert decode(result["twist_parameter"],curve.field()).trace() == 1
    encode = lambda value: {"poly":[{"coeff":c,"power":i} for i,c in enumerate(value)]} if isinstance(value,list) else value
    other = copy.deepcopy(record)
    other["params"] = {"f":[encode(c) for c in result["twisted_f"]],"h":[encode(c) for c in result["twisted_h"]]}
    assert [int(c) for c in l_polynomial(curve_from_dict(other),{})] == result["twist_l_polynomial"]


def test_two_torsion_matrix_certifies_rational_ranks_and_minimal_full_field(g2):
    result = TRAITS["jacobian_torsion"](g2)
    assert_ok(result)
    assert result["base_rational_rank"] == 3
    assert result["full_torsion_extension_degree"] == 2
    assert result["least_nontrivial_extension_degree"] == 1
    assert result["extension_ranks"][1]["rank"] == 4
    F = matrix(GF(2), result["frobenius_matrix"])
    assert F != identity_matrix(GF(2),4)
    assert F*F == identity_matrix(GF(2),4)


def test_two_torsion_quotient_agrees_with_independent_even_subset_orbits(g2):
    result = TRAITS["jacobian_torsion"](g2)
    permutation = TRAITS["genus2_rosenhain"](g2)["branch_frobenius"]["permutation"]
    for n in (1,2,3):
        fixed = set()
        for bits in product((0,1),repeat=6):
            if sum(bits) % 2:
                continue
            image = list(bits)
            for _ in range(n):
                new = [0]*6
                for i,j in enumerate(permutation):
                    new[j] = image[i]
                image = new
            key = min(bits,tuple(1-b for b in bits))
            if min(tuple(image),tuple(1-b for b in image)) == key:
                fixed.add(key)
        assert len(fixed) == 2**result["extension_ranks"][n-1]["rank"]


def test_exact_two_torsion_extends_beyond_genus_two():
    curve = curve_from_dict(hyper([1,1]+[0]*7+[1]))
    result = TRAITS["jacobian_torsion"](curve)
    assert_ok(result)
    assert result["module_dimension"] == 8
    F = matrix(GF(2),result["frobenius_matrix"])
    assert result["full_torsion_extension_degree"] == F.multiplicative_order() == 9


def test_repeated_frobenius_factors_do_not_assume_semisimplicity():
    # The Weil polynomial of E^2 with q=7 and trace(E)=2 reduces to (T-1)^4 mod 3.
    rows, semisimple, unipotent, _ = torsion_constraints([49,-28,18,-4,1],3,9)
    assert rows[0]["rank_lower_bound"] == 1
    assert rows[0]["rank_upper_bound"] == 2
    assert rows[0]["rank"] is None
    assert semisimple == 1
    assert unipotent == 9
    assert rows[8]["rank"] == 4


def test_linear_factor_modular_exponentiation_avoids_ntl_unreduced_base():
    ring = PolynomialRing(GF(2), "T")
    T = ring.gen()
    assert polynomial_power_mod(T,17,T+1) == 1
    assert polynomial_power_mod(T,17,T) == 0


def test_squarefree_torsion_characteristic_polynomial_determines_ranks(ordinary):
    result = TRAITS["jacobian_torsion"](ordinary,ell=7)
    assert_ok(result)
    assert all(row["rank"] is not None for row in result["extension_ranks"])
    assert result["unipotent_order_upper_bound"] == 1


def test_torsion_prime_and_extension_constraints(g2):
    with pytest.raises(ValueError,match="prime"):
        TRAITS["jacobian_torsion"](g2,ell=4)
    assert TRAITS["jacobian_torsion"](g2,ell=7)["status"] == "unsupported"
    assert TRAITS["jacobian_torsion"](g2,extension_bound=65)["status"] == "resource_limit"


def test_genus_two_symmetric_square_is_jacobian_blowup_count(g2):
    result = TRAITS["symmetric_power_divisors"](g2)
    assert_ok(result)
    assert result["symmetric_power_counts"][2]["count"] == result["jacobian_order"]+7
    assert result["degree_two_effective_classes"] == result["jacobian_order"]
    assert result["degree_two_noneffective_classes"] == 0
    assert result["degree_two_special_classes"] == 1
    assert result["riemann_roch_verified_degrees"] == list(range(3,13))


def test_genus_three_degree_two_classes_do_not_have_hyperelliptic_pencil(quartic):
    result = TRAITS["symmetric_power_divisors"](quartic)
    assert_ok(result)
    assert result["degree_two_effective_classes"] == result["symmetric_power_counts"][2]["count"]
    assert result["degree_two_special_classes"] == 0
    assert result["symmetric_power_counts"][3]["count"] == result["jacobian_order"]+7*result["extension_point_counts"][0]


def test_ordinary_cm_orders_satisfy_index_identities_and_maximality(ordinary):
    result = TRAITS["frobenius_orders"](ordinary)
    assert_ok(result)
    assert result["frobenius_polynomial"] == [25,20,8,4,1]
    assert result["cm_field_discriminant"] == 2304
    assert result["real_subfield_discriminant"] == 24
    assert result["cm_field_signature"] == [0,2]
    for name in ("frobenius","weil"):
        assert result[name+"_order_discriminant"] == result["cm_field_discriminant"]*result[name+"_order_index"]**2
    assert result["actual_endomorphism_order_determined"] == 1
    assert result["actual_endomorphism_order_maximal"] == 1


def test_cm_orders_do_not_apply_ordinary_formula_to_supersingular_jacobian(g2):
    assert TRAITS["frobenius_orders"](g2)["status"] == "not_applicable"


def test_nonmaximal_weil_order_does_not_determine_the_actual_endomorphism_order():
    result = TRAITS["frobenius_orders"](curve_from_dict(hyper([1,1,0,2,0,1],prime=5)))
    assert_ok(result)
    assert result["weil_order_index"] == 4
    assert result["actual_endomorphism_order_determined"] == 0
    assert result["actual_endomorphism_order_maximal"] is None


def test_cm_number_field_degree_is_bounded(ordinary):
    assert TRAITS["frobenius_orders"](ordinary,max_order_degree=2)["status"] == "resource_limit"


@pytest.mark.parametrize("genus,quadrics,syzygies", [(2,0,0),(3,1,0),(4,3,2),(5,6,8)])
def test_hyperelliptic_canonical_images_are_rational_normal_curves(genus,quadrics,syzygies):
    coefficients = [0,-1,0,0,0,1] if genus == 2 else [1,1]+[0]*(2*genus-1)+[1]
    curve = curve_from_dict(hyper(coefficients))
    result = TRAITS["canonical_relations"](curve)
    assert_ok(result)
    assert curve.genus() == genus
    assert result["minimal_generator_count"] == quadrics
    assert result["linear_quadric_syzygy_count"] == syzygies
    assert result["canonical_ring_generated_in_degree_one"] == 0
    assert [row["multiplication_rank"] for row in result["degree_profiles"]] == [n*(genus-1)+1 for n in range(1,5)]


def test_nonhyperelliptic_genus_three_has_one_quartic_equation(quartic):
    result = TRAITS["canonical_relations"](quartic)
    assert_ok(result)
    assert [row["minimal_generator_count"] for row in result["degree_profiles"]] == [0,0,0,1]
    assert result["canonical_ring_generated_in_degree_one"] == 1
    assert all(row["canonical_subring_deficiency"] == 0 for row in result["degree_profiles"])


def test_genus_four_canonical_image_has_one_quadric_and_one_cubic():
    record = hyper()
    record["form"] = "CyclicTrigonal"
    curve = curve_from_dict(record)
    result = TRAITS["canonical_relations"](curve)
    assert_ok(result)
    assert [row["minimal_generator_count"] for row in result["degree_profiles"]] == [0,1,1,0]


def test_canonical_syzygy_certificates_survive_record_caps():
    curve = curve_from_dict(hyper([1,1]+[0]*7+[1]))
    result = TRAITS["canonical_relations"](curve,max_relation_records=1)
    assert result["status"] == "partial"
    assert result["minimal_generator_count"] == 3
    assert len(result["indexed_quadric_terms"]) == 3
    ring = PolynomialRing(curve.field(),names=result["variables"])
    quadrics = [sum(decode(term["coefficient"],curve.field())*ring.monomial(*term["powers"]) for term in terms)
                for terms in result["indexed_quadric_terms"]]
    for row in result["linear_quadric_syzygies"]:
        forms = [sum(decode(c,curve.field())*v for c,v in zip(coefficients,ring.gens())) for coefficients in row]
        assert sum(a*b for a,b in zip(forms,quadrics)) == 0


def test_canonical_monomial_bound_and_incomplete_degree_scope(quartic):
    assert TRAITS["canonical_relations"](quartic,max_monomials=2)["status"] == "resource_limit"
    result = TRAITS["canonical_relations"](quartic,relation_degree=3)
    assert_ok(result)
    assert result["ideal_generation_complete"] == 0


def test_new_traits_round_trip_nested_database_results(g2):
    class Collection:
        def insert_one(self,record):
            self.record=record
            return type("Ack",(),{"acknowledged":True})()
    collection = Collection()
    result = TRAITS["genus2_richelot"](g2)
    assert store_trait_result({"trait_genus2_richelot":collection},g2,"genus2_richelot",{},result)
    assert _decode_ints(collection.record)["result"] == result


def test_new_traits_run_in_combined_report_and_json_cli():
    report = curve_report(hyper(),traits=MORE)
    assert set(report["traits"]) == set(MORE)
    assert all(row["status"] in ("ok","not_applicable") for row in report["traits"].values())
    completed = subprocess.run([sys.executable,"-m","dissect.utils.compute_json","-t","jacobian_torsion",
                                "-i",str(ROOT/"examples/higher_genus_curves.json"),"--params",'{"ell":3,"extension_bound":3}'],
                               capture_output=True,text=True,check=True)
    rows = json.loads(completed.stdout)["data"]
    assert rows[0]["result"]["status"] == "not_applicable"
    assert rows[1]["result"]["module_dimension"] == 4
