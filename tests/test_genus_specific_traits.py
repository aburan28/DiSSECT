import pytest

from dissect.traits import TRAITS
from dissect.utils.curve_factory import HIGHER_GENUS_FORMS
from dissect.utils.higher_genus_report import REPORT_TRAITS


NAMES = ("genus2_rosenhain", "genus2_elliptic_quotients", "genus3_bitangents",
         "genus3_flexes", "genus4_canonical_quadric", "genus4_tritangents")


def test_genus_specific_traits_are_registered_with_report_and_limits():
    assert len(TRAITS) == 50
    assert "CanonicalGenus4" in HIGHER_GENUS_FORMS
    for name in NAMES:
        trait = TRAITS[name]
        assert name in REPORT_TRAITS
        assert trait.GENUS in (2, 3, 4)
        assert "timeout" in trait.DEFAULT_PARAMS
        assert "status" in trait.OUTPUT


@pytest.mark.parametrize("name", NAMES)
def test_wrong_genus_is_not_applicable_before_sage_import(name):
    trait = TRAITS[name]

    class OtherGenus:
        def genus(self):
            return 5

        def curve_type(self):
            return trait.SUPPORTED_CURVE_TYPES[0]

    result = trait(OtherGenus())
    assert result["status"] == "not_applicable"
    assert "genus" in result["reason"]


@pytest.mark.parametrize("name", NAMES)
def test_elliptic_curves_are_excluded_from_genus_specific_traits(name):
    assert TRAITS[name](None)["status"] == "not_applicable"
