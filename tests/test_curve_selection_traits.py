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
