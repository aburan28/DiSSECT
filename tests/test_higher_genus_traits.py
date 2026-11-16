import pytest
from contextlib import nullcontext
from types import SimpleNamespace

from dissect.traits import TRAITS
from dissect.traits.canonical_geometry import canonical_profile
from dissect.utils.database_handler import _format_curve


@pytest.mark.parametrize("genus,hyperelliptic,map_degree,image_degree", [
    (2, True, 2, 1),
    (3, True, 2, 2),
    (3, False, 1, 4),
    (4, False, 1, 6),
])
def test_known_canonical_models(genus, hyperelliptic, map_degree, image_degree):
    profile = canonical_profile(genus, hyperelliptic)
    assert profile["canonical_map_degree"] == map_degree
    assert profile["canonical_image_degree"] == image_degree
    assert profile["jacobian_dimension"] == genus
    assert profile["holomorphic_differentials"] == genus
    assert profile["canonical_embedding"] == int(not hyperelliptic)
    assert profile["degree_two_map"] == int(hyperelliptic)
    assert profile["unique_degree_two_map"] == int(hyperelliptic)
    assert profile["jacobian_is_curve"] == 0


def test_elliptic_double_cover_is_not_an_exclusive_higher_genus_trait():
    profile = canonical_profile(1, False)
    assert profile["degree_two_map"] == 1
    assert profile["unique_degree_two_map"] == 0
    assert profile["canonical_map_degree"] is None
    assert profile["canonical_image_degree"] is None
    assert profile["canonical_embedding"] == 0
    assert profile["jacobian_is_curve"] == 1


@pytest.mark.parametrize("genus,hyperelliptic", [(0, False), (1, True), (2, False)])
def test_impossible_families_are_rejected(genus, hyperelliptic):
    with pytest.raises(ValueError):
        canonical_profile(genus, hyperelliptic)


def test_genus_three_nonhyperelliptic_means_plane_quartic():
    assert canonical_profile(3, False)["plane_quartic"] == 1
    assert canonical_profile(3, True)["plane_quartic"] == 0
    assert canonical_profile(4, False)["plane_quartic"] == 0


@pytest.mark.parametrize("family", ["hyperelliptic", "non_hyperelliptic"])
def test_existing_traits_do_not_apply_elliptic_formulas_to_higher_genus(family):
    class CurveWithoutEllipticMethods:
        def curve_type(self):
            return family

    curve = CurveWithoutEllipticMethods()
    for name, trait in TRAITS.items():
        if trait.SUPPORTED_CURVE_TYPES == ("elliptic",):
            assert trait(curve) == dict.fromkeys(trait.OUTPUT)


def test_higher_genus_database_records_do_not_need_invented_group_orders():
    record = {
        "name": "example", "category": "examples", "form": "Hyperelliptic",
        "field": {"type": "Prime", "p": "0x7", "bits": 3},
        "params": {"f": [0, -1, 0, 0, 0, 1]}, "genus": 2,
    }
    formatted = _format_curve(record)
    assert formatted["order"] is None
    assert formatted["cofactor"] is None
    assert formatted["genus"] == 2
    assert formatted["params"] == record["params"]


@pytest.mark.parametrize("skip,expected", [(0, ["first", "second"]), (1, ["second"]), (10, [])])
def test_database_producer_skips_exactly_the_requested_number(monkeypatch, skip, expected):
    from dissect.utils import compute_db

    monkeypatch.setattr(compute_db, "connect", lambda database: None)
    monkeypatch.setattr(compute_db, "create_trait_index", lambda *args: None)
    monkeypatch.setattr(compute_db, "get_curves_count", lambda *args, **kwargs: 2)
    monkeypatch.setattr(compute_db, "get_trait_results_count", lambda *args, **kwargs: 0)
    monkeypatch.setattr(compute_db, "get_curves", lambda *args, **kwargs: ["first", "second"])
    monkeypatch.setattr(compute_db, "curve_from_dict", lambda record: record)
    queued = []
    queue = SimpleNamespace(put=lambda value, **kwargs: queued.append(value))
    args = SimpleNamespace(skip=skip, jobs=1, verbose=False)
    compute_db.producer(None, "canonical_geometry", args, queue, nullcontext())
    assert queued == [(name, {}) for name in expected] + [(None, None)]
