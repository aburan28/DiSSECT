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


# secp256r1, whose seed is published and verifiable. Fixed vectors, no network.
SECP256R1 = {
    "name": "secp256r1",
    "category": "secg",
    "cofactor": 1,
    "field": {
        "type": "Prime",
        "bits": 256,
        "p": 0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF,
    },
    "params": {
        "a": {"raw": 0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFC},
        "b": {"raw": 0x5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B},
    },
    "simulation": {"seed": 0xC49D360886E704936A6678E1139D26B7819F7E90},
}


def test_x962_derivation_reproduces_a_published_curve():
    from dissect.analysis.parameter_forms import verifies_x962_seed

    p = SECP256R1["field"]["p"]
    a = SECP256R1["params"]["a"]["raw"]
    b = SECP256R1["params"]["b"]["raw"]
    seed = SECP256R1["simulation"]["seed"].to_bytes(20, "big")
    assert verifies_x962_seed(seed, p, a, b)


def test_a_wrong_seed_does_not_verify():
    from dissect.analysis.parameter_forms import verifies_x962_seed

    p = SECP256R1["field"]["p"]
    a = SECP256R1["params"]["a"]["raw"]
    b = SECP256R1["params"]["b"]["raw"]
    wrong = (SECP256R1["simulation"]["seed"] ^ 1).to_bytes(20, "big")
    assert not verifies_x962_seed(wrong, p, a, b)


def test_seed_verification_reports_success_on_a_real_curve():
    from dissect.analysis.parameter_forms import seed_verification

    verified, note = seed_verification(SECP256R1)
    assert verified is True
    assert "20-byte" in note


def test_seed_verification_skips_what_it_cannot_check():
    from dissect.analysis.parameter_forms import seed_verification

    # binary field: the derivation is not defined the same way
    assert seed_verification({"field": {"type": "Binary", "bits": 163},
                              "simulation": {"seed": 5}}) is None
    # no seed to check
    assert seed_verification({"field": {"type": "Prime", "bits": 256, "p": 7},
                              "simulation": {}}) is None


def test_unreadable_parameters_are_not_reported_as_a_failed_derivation():
    """A curve we cannot check must not be accused of a bad seed."""
    from dissect.analysis.parameter_forms import seed_verification

    broken = {
        "name": "unreadable",
        "field": {"type": "Prime", "bits": 256, "p": "not-a-number"},
        "params": {"a": {"raw": 1}, "b": {"raw": 2}},
        "simulation": {"seed": 0xC49D360886E704936A6678E1139D26B7819F7E90},
    }
    verified, note = seed_verification(broken)
    assert verified is None, "unreadable parameters must be 'unknown', not False"
    assert "unavailable" in note


def test_report_separates_unknown_from_failed(capsys):
    from dissect.analysis.parameter_forms import report

    curves = [
        SECP256R1,
        {"name": "unreadable", "category": "x", "cofactor": 1,
         "field": {"type": "Prime", "bits": 256, "p": "not-a-number"},
         "params": {"a": {"raw": 1}, "b": {"raw": 2}},
         "simulation": {"seed": 0xC49D360886E704936A6678E1139D26B7819F7E90}},
    ]
    report(curves)
    out = capsys.readouterr().out
    assert "1/1 prime-field seeds derive their own coefficients" in out
    assert "FAILS: unreadable" not in out
    assert "could not be checked" in out
