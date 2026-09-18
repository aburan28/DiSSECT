import math

from dissect.traits.cheon import (
    best_oracle_attack,
    divisors_from_factorization,
    oracle_attack_profile,
)
from dissect.traits.extension_degree import proper_divisors


def test_balanced_oracle_profile_counts_both_phases():
    profile = oracle_attack_profile(2**90 + 1, 2**30)
    assert math.isclose(profile["query_bits"], 30.0)
    assert math.isclose(profile["offline_bits"], 30.0, abs_tol=1e-12)
    assert math.isclose(profile["total_bits"], 31.0, abs_tol=1e-12)


def test_composite_divisors_are_searched():
    result = best_oracle_attack(1009, [(2, 4), (3, 2), (7, 1)])
    assert result["divisor"] == 6


def test_query_budget_changes_selected_divisor():
    result = best_oracle_attack(
        1009, [(2, 4), (3, 2), (7, 1)], max_query_bits=2
    )
    assert result["divisor"] == 4


def test_r_plus_1_itself_is_not_an_oracle_divisor():
    result = best_oracle_attack(101, [(2, 1), (3, 1), (17, 1)])
    assert result["divisor"] == 3


def test_p521_complete_r_minus_1_is_not_mislabeled_cheon_immune():
    order = int(
        "1ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
        "fa51868783bf2f966b7fcc0148f709a5d03bb5c9b8899c47aebb6fb71e91386409",
        16,
    )
    factors = [
        (2, 3),
        (7, 1),
        (11, 1),
        (1283, 1),
        (1458105463, 1),
        (1647781915921980690468599, 1),
        (
            int(
                "3615194794881930010216942559103847593050265703173292383701371712"
                "350878926821661243755933835426896058418509759880171943"
            ),
            1,
        ),
    ]

    result = best_oracle_attack(order, factors)

    assert result["divisor"] == 1898873518475180724503002533770555108536
    assert math.isclose(result["query_bits"], 130.48033951323217)
    assert math.isclose(result["total_bits"], 195.25983024338393)
    assert result["total_bits"] < 0.5 * math.log2(order)


def test_divisor_enumeration_has_a_resource_bound():
    try:
        divisors_from_factorization([(2, 20)], limit=20)
    except ValueError as error:
        assert "21 divisors" in str(error)
    else:
        raise AssertionError("oversized divisor search was not rejected")


def test_prime_extension_degree_still_has_the_prime_subfield():
    assert proper_divisors(163) == [1]


def test_composite_extension_degree_includes_all_subfield_degrees():
    assert proper_divisors(155) == [1, 5, 31]
