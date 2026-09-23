"""Tests for the provenance checks, with P-256's published constants."""

import math

from dissect.analysis.provenance_checks import (
    binomial_cdf,
    embedding_degree,
    hamming,
    seed_pair_statistics,
    shared_leading_bits,
    sign_test_p,
)

P256_P = 0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF
P256_ORDER = 0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551


def test_embedding_degree_finds_small_k():
    # 2 has order 10 mod 11
    assert embedding_degree(2, 11) == 10
    assert embedding_degree(12, 11) == 1


def test_embedding_degree_respects_limit():
    assert embedding_degree(2, 11, limit=9) is None


def test_p256_curve_embedding_degree_is_large():
    assert embedding_degree(P256_P, P256_ORDER, limit=200) is None


def test_hamming_and_shared_prefix():
    assert hamming(0b1010, 0b0110) == 2
    assert shared_leading_bits(1 << 159, 1 << 159) == 160
    assert shared_leading_bits(1 << 159, 0) == 0


def test_binomial_cdf_is_symmetric_at_the_mean():
    assert math.isclose(binomial_cdf(79) + binomial_cdf(80) , 1.0)
    assert binomial_cdf(160) == 1.0


def test_counter_seeds_are_flagged():
    stats = seed_pair_statistics({"a": 1 << 159, "b": (1 << 159) + 1})
    assert stats["hamming_min"] == 1
    assert stats["min_log2_diff"] == 0.0
    assert stats["p_min"] < 1e-40
    assert stats["pairs"][0]["shared_leading"] == 159


def test_complementary_seeds_have_maximal_distance():
    x = 0x123456789ABCDEF0123456789ABCDEF012345678
    stats = seed_pair_statistics({"a": x, "b": x ^ ((1 << 160) - 1)})
    assert stats["hamming_mean"] == 160


def test_sign_test():
    assert sign_test_p(5, 5) == 1 / 32
    assert sign_test_p(0, 7) == 1.0
