"""Tests use P-256's published constants, cross-checked by p + 1 - t == n."""

import math

import pytest

from dissect.analysis.isogeny_structure import (
    SMALL_PRIMES,
    conductor_small_part,
    frobenius_discriminant,
    kronecker_profile,
    legendre,
    log_l_value,
    rho_security_bits,
    split_count,
    to_int,
    twist_order,
    two_is_forced_inert,
)

P256_P = 0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF
P256_TRACE = 89188191154553853111372247798585809583
P256_ORDER = 0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551
# the database's own twist_order trait reports this cardinality for P-256
P256_TWIST = 115792089210356248762697446949407573530175331606444868048645003556665683663535


def test_the_test_vector_is_self_consistent():
    assert P256_P + 1 - P256_TRACE == P256_ORDER


def test_to_int_does_not_read_a_decimal_string_as_hex():
    assert to_int("31") == 31
    assert to_int("0x31") == 49
    assert to_int("-0x10") == -16
    assert to_int(7) == 7


def test_frobenius_discriminant_is_negative_and_about_4p():
    D = frobenius_discriminant(P256_P, P256_TRACE)
    assert D < 0
    assert abs(D).bit_length() in (P256_P.bit_length(), P256_P.bit_length() + 2)


def test_two_is_inert_and_that_is_forced_not_chosen():
    assert two_is_forced_inert(P256_P, P256_TRACE)
    D = frobenius_discriminant(P256_P, P256_TRACE)
    assert D % 8 == 5
    # forced for ANY odd trace over an odd prime field, not special to P-256
    for fake_trace in (1, 3, 5, 12345):
        assert frobenius_discriminant(P256_P, fake_trace) % 8 == 5
    # an even trace breaks it, which is why prime order is the precondition
    assert frobenius_discriminant(P256_P, 2) % 8 != 5


def test_legendre_classifies_split_inert_and_ramified():
    assert legendre(9, 3) == 0          # 3 divides
    assert legendre(1, 3) == 1          # 1 is a square mod 3
    assert legendre(2, 3) == -1         # 2 is not
    assert legendre(4, 5) == 1
    assert legendre(2, 5) == -1


def test_kronecker_profile_covers_every_requested_prime():
    D = frobenius_discriminant(P256_P, P256_TRACE)
    profile = kronecker_profile(D)
    assert set(profile) == set(SMALL_PRIMES)
    assert set(profile.values()) <= {-1, 0, 1}
    assert split_count(D) == sum(1 for v in profile.values() if v == 1)


def test_p256_splits_at_nine_of_the_fourteen_small_primes():
    assert split_count(frobenius_discriminant(P256_P, P256_TRACE)) == 9


def test_conductor_small_part_finds_repeated_factors_only():
    assert conductor_small_part(-(3 ** 3) * 7) == {3: 3}
    assert conductor_small_part(-3 * 5 * 7) == {}
    # P-256 has a maximal endomorphism ring as far as 47, matching stored conductor 1
    assert conductor_small_part(frobenius_discriminant(P256_P, P256_TRACE)) == {}


def test_l_value_is_positive_and_monotone_in_the_truncation_bound():
    D = frobenius_discriminant(P256_P, P256_TRACE)
    small = log_l_value(D, 500)
    large = log_l_value(D, 5000)
    assert math.exp(small) > 0
    # more Euler factors changes the estimate; they must not be equal
    assert small != large


def test_l_value_separates_two_discriminants():
    a = log_l_value(frobenius_discriminant(P256_P, P256_TRACE), 2000)
    b = log_l_value(frobenius_discriminant(P256_P, P256_TRACE - 2), 2000)
    assert a != b


def test_twist_order_matches_the_database_trait():
    assert twist_order(P256_P, P256_TRACE) == P256_TWIST


def test_twist_order_and_curve_order_sum_to_twice_p_plus_two():
    assert twist_order(P256_P, P256_TRACE) + P256_ORDER == 2 * P256_P + 2


def test_rho_security_is_half_the_largest_prime_factor():
    assert rho_security_bits([(2, 3), (7, 1), (2 ** 127 - 1, 1)]) == 63
    assert rho_security_bits([(3, 1)]) == 1
    assert rho_security_bits([]) is None
