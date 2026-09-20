import math

from dissect.analysis.parameter_forms import (
    degenerate_parameters,
    degree_cofactor_pairs,
    factorize,
    generator_scale,
    leading_zero_profile,
    seed_of,
    subfield_index,
    to_int,
)


def test_to_int_accepts_ints_and_based_strings():
    assert to_int(255) == 255
    assert to_int("0xff") == 255
    assert to_int("255") == 255


def test_factorize_handles_primes_and_prime_powers():
    assert factorize(163) == [163]
    assert factorize(176) == [2, 2, 2, 2, 11]
    assert factorize(16) == [2, 2, 2, 2]


def test_prime_degree_offers_no_subfield():
    assert subfield_index(163) is None
    assert subfield_index(571) is None


def test_composite_degree_reports_the_power_of_two_index():
    assert subfield_index(176) == 16
    assert subfield_index(208) == 16
    assert subfield_index(368) == 16
    # an odd composite has no power-of-two subfield index
    assert subfield_index(9) is None


def _binary(name, degree, cofactor):
    return {"name": name, "category": "x962", "cofactor": cofactor,
            "field": {"type": "Binary", "bits": degree}}


def test_degree_cofactor_pairs_skips_prime_degree_and_prime_fields():
    curves = [
        _binary("c2pnb176w1", 176, 65390),
        _binary("ansit163r2", 163, 2),
        {"name": "p256", "category": "nist", "cofactor": 1,
         "field": {"type": "Prime", "bits": 256, "p": 7}},
    ]
    rows = degree_cofactor_pairs(curves)
    assert [row["curve"] for row in rows] == ["c2pnb176w1"]
    # the cofactor is almost exactly the size of the subfield, 2^16
    assert rows[0]["subfield_index"] == 16
    assert rows[0]["subfield_size"] == 65536
    assert math.isclose(rows[0]["cofactor_over_subfield"], 65390 / 65536)
    assert 0.99 < rows[0]["cofactor_over_subfield"] < 1.0


def test_leading_zero_profile_matches_the_uniform_tail():
    # 80 uniform 160-bit seeds: half reach 160 bits, three quarters reach 159.
    rows = {row["threshold"]: row for row in leading_zero_profile([160] * 80)}
    assert math.isclose(rows[160]["expected"], 40.0)
    assert math.isclose(rows[159]["expected"], 60.0)
    assert math.isclose(rows[158]["expected"], 70.0)
    assert rows[160]["observed"] == 80


def test_seed_of_treats_missing_and_empty_alike():
    assert seed_of({"simulation": {"seed": 12}}) == 12
    assert seed_of({"simulation": {"seed": ""}}) is None
    assert seed_of({"simulation": {}}) is None
    assert seed_of({}) is None


def test_generator_scale_reports_bit_length_and_tolerates_absence():
    assert generator_scale({"generator": {"x": {"raw": 13}}}) == 4
    assert generator_scale({"generator": {"x": {"raw": 0}}}) == 0
    assert generator_scale({}) is None


def test_degenerate_parameters_names_each_shape():
    assert "anomalous (trace 1, #E = p)" in degenerate_parameters(
        {"properties": {"trace": 1}})
    assert "supersingular (trace 0)" in degenerate_parameters(
        {"properties": {"trace": 0}})
    assert degenerate_parameters({"properties": {"trace": 5, "embedding_degree": 12}}) == [
        "embedding degree 12"]
    assert degenerate_parameters({"properties": {"trace": 5, "embedding_degree": 10**9}}) == []
    assert degenerate_parameters({}) == []
