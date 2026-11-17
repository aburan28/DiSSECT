import pytest

from dissect.traits import TRAITS
from dissect.utils.higher_genus_report import REPORT_TRAITS
from dissect.utils.jacobian_profiles import divisor_counts, mobius


NAMES = ("genus2_richelot", "hyperelliptic_twists", "jacobian_torsion",
         "symmetric_power_divisors", "frobenius_orders", "canonical_relations")


def test_more_traits_are_registered_for_both_reporting_commands():
    assert len(TRAITS) == 56
    for name in NAMES:
        assert name in REPORT_TRAITS
        assert {"status", "reason"} <= set(TRAITS[name].OUTPUT)
        assert {"timeout", "max_genus", "max_field_size"} <= set(TRAITS[name].DEFAULT_PARAMS)


@pytest.mark.parametrize("name", NAMES)
def test_more_traits_do_not_apply_elliptic_formulas(name):
    assert TRAITS[name](None)["status"] == "not_applicable"


@pytest.mark.parametrize("name", NAMES)
def test_more_traits_validate_limits_before_loading_sage(name):
    with pytest.raises(ValueError, match="positive integer"):
        TRAITS[name](None, timeout=0)
    with pytest.raises(ValueError, match="Unknown trait parameters"):
        TRAITS[name](None, invented_parameter=1)


def test_mobius_counts_include_square_factors_and_empty_product():
    assert [mobius(n) for n in range(1, 13)] == [1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0]


def test_genus_two_effective_divisors_match_independent_multiset_count():
    points, closed, effective, rr = divisor_counts([1, 0, -2, 0, 49], 7, 2, 8)
    assert points[1:3] == [8, 46]
    assert closed[1:3] == [8, 19]
    assert effective[2] == 8*9//2+19 == 48+7
    assert effective[3] == 48*8
    assert rr == list(range(3, 9))
    for n in range(1, 9):
        assert points[n] == sum(d*closed[d] for d in range(1, n+1) if n % d == 0)


def test_genus_three_picard_fiber_identity_and_degree_zero():
    points, _, effective, rr = divisor_counts([1, 0, 21, 0, 147, 0, 343], 7, 3, 8)
    assert effective[0] == 1
    assert effective[3] == 512+7*points[1]
    assert rr == list(range(5, 9))
