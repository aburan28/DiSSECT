"""Independent arithmetic identities and resource/parameter contracts; no Sage required."""

from fractions import Fraction

import pytest

from dissect.traits import TRAITS
from dissect.traits.detailed_trait import DetailedTrait
from dissect.utils.curve_analysis import extension_polynomial, newton_slopes, power_sums, semigroup_profile


@pytest.mark.parametrize("extension", range(1, 9))
def test_frobenius_extensions_agree_with_explicit_supersingular_eigenvalues(extension):
    polynomial = [49, 0, 14, 0, 1]  # (T^2+7)^2
    trace = power_sums(polynomial, extension)[extension]
    expected_trace = 0 if extension % 2 else 4*(-7)**(extension//2)
    assert trace == expected_trace
    order = sum(extension_polynomial(polynomial, extension))
    expected_order = (1+7**extension)**2 if extension % 2 else (1-(-7)**(extension//2))**4
    assert order == expected_order


def test_newton_polygon_distinguishes_intermediate_slopes():
    slopes, _ = newton_slopes([1, 1, 14, 7, 49], 7, 1)
    assert slopes == [Fraction(0), Fraction(1, 2), Fraction(1, 2), Fraction(1)]


def test_newton_slopes_are_normalized_by_the_field_degree():
    slopes, _ = newton_slopes([1, 0, 49], 7, 2)
    assert slopes == [Fraction(1, 2)]*2


def test_newton_polygon_handles_zero_coefficients_and_nonbinary_denominators():
    slopes, _ = newton_slopes([1, 0, 0, 7, 0, 0, 343], 7, 1)
    assert slopes == [Fraction(1, 3)]*3+[Fraction(2, 3)]*3


@pytest.mark.parametrize("gaps,generators,conductor", [
    ([1, 3], [2, 5], 4),
    ([1, 2], [3, 4, 5], 3),
    ([1, 2, 5], [3, 4], 6),
    ([1, 3, 5, 7], [2, 9], 8),
])
def test_known_weierstrass_numerical_semigroups(gaps, generators, conductor):
    profile = semigroup_profile(gaps)
    assert profile["semigroup_generators"] == generators
    assert profile["semigroup_conductor"] == conductor


def test_overrides_replace_parameter_choices_without_duplicate_jobs():
    trait = TRAITS["twist_order"]
    assert list(trait.params_iter({"deg": 2})) == [{"deg": 2}]


@pytest.mark.parametrize("value", [0, -1, 0.5, "30", True, None])
def test_invalid_budget_parameters_are_rejected_before_sage_import(value):
    with pytest.raises(ValueError, match="positive integer"):
        TRAITS["hasse_witt"](None, timeout=value)


def test_unknown_parameters_are_rejected():
    with pytest.raises(ValueError, match="Unknown"):
        TRAITS["hasse_witt"](None, timout=1)


def test_all_detailed_traits_are_registered_with_limits_and_statuses():
    names = ["frobenius_polynomial", "jacobian_decomposition", "hasse_witt", "newton_polygon",
             "weierstrass_gaps", "automorphism_structure", "moduli_invariants",
             "gonality_and_quotients", "cover_ramification"]
    for name in names:
        trait = TRAITS[name]
        assert isinstance(trait, DetailedTrait)
        assert "timeout" in trait.DEFAULT_PARAMS
        assert "status" in trait.OUTPUT
        assert trait(None)["status"] == "not_applicable"
