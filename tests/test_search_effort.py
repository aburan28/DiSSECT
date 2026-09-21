import math

import pytest

from dissect.analysis.search_effort import (
    SECONDS_PER_YEAR,
    ceiling,
    extrapolated_floor,
    floor_from_gaps,
    hiding_room_bits,
    seed_gaps,
    to_int,
)


def test_to_int_does_not_silently_misread_a_decimal_string_as_hex():
    # The bug this guards: int("31", 16) == 49 and raises nothing.
    assert to_int("31") == 31
    assert to_int("0x31") == 49
    assert to_int(31) == 31
    assert to_int("-0x10") == -16
    assert to_int("-31") == -31


def test_seed_gaps_sorts_before_differencing():
    assert seed_gaps([30, 10, 20]) == [10, 10]
    assert seed_gaps([5]) == []
    assert seed_gaps([]) == []


def test_floor_scales_by_the_prime_order_fraction():
    gaps = [10] * 9  # mean gap 10, so 10 seeds per curve kept
    # half the kept curves have cofactor 1, so 20 seeds per prime-order curve
    assert floor_from_gaps(gaps, kept=100, prime_order=50) == pytest.approx(20.0)
    assert floor_from_gaps(gaps, kept=100, prime_order=100) == pytest.approx(10.0)


def test_floor_is_none_when_it_cannot_be_measured():
    assert floor_from_gaps([], kept=10, prime_order=5) is None
    assert floor_from_gaps([1, 2], kept=10, prime_order=0) is None


def test_extrapolated_floor_grows_with_ln_p():
    # calibrated at 256 bits, a 512-bit field should cost about twice as many trials
    got = extrapolated_floor(512, reference_bits=256, reference_trials=600)
    assert got == pytest.approx(1200.0)
    # and the calibration point reproduces itself
    assert extrapolated_floor(256, 256, 600) == pytest.approx(600.0)


def test_ceiling_is_linear_in_cpu_years_and_inverse_in_slowdown():
    one = ceiling(seconds_per_count=1.0, slowdown=1.0, cpu_years=1.0)
    assert one == pytest.approx(SECONDS_PER_YEAR)
    assert ceiling(1.0, 1.0, 10.0) == pytest.approx(10 * one)
    assert ceiling(1.0, 100.0, 1.0) == pytest.approx(one / 100)
    assert ceiling(2.0, 1.0, 1.0) == pytest.approx(one / 2)


def test_hiding_room_is_the_log_gap_between_floor_and_ceiling():
    assert hiding_room_bits(2**9, 2**27) == pytest.approx(18.0)
    assert hiding_room_bits(100, 100) == pytest.approx(0.0)
    assert hiding_room_bits(0, 10) is None
    assert hiding_room_bits(10, 0) is None


def test_a_bigger_curve_has_less_hiding_room_when_counting_costs_more():
    """The structural point: the ceiling falls with field size faster than the
    floor rises, so larger curves leave less room for a hidden search."""
    small = hiding_room_bits(2**9.2, ceiling(1.44, 2**9, 1e4))
    large = hiding_room_bits(2**10.3, ceiling(125.29, 2**9, 1e4))
    assert large < small
    assert small - large == pytest.approx(math.log2(125.29 / 1.44) + (10.3 - 9.2), abs=0.1)
