"""Run with sage -python -m pytest; ordinary Python skips this module."""

import json
import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("sage.all")

from dissect.traits import TRAITS
from dissect.utils.curve_factory import curve_from_dict
from dissect.utils.database_handler import _decode_ints, _format_curve, store_trait_result
from dissect.utils.custom_curve import CustomCurve


EXAMPLES = Path(__file__).resolve().parents[1] / "examples" / "higher_genus_curves.json"


def hyperelliptic(f=None, h=None):
    record = {
        "name": "genus2", "category": "examples", "form": "Hyperelliptic",
        "field": {"type": "Prime", "p": "0x7", "bits": 3},
        "params": {"f": [0, -1, 0, 0, 0, 1] if f is None else f},
    }
    if h is not None:
        record["params"]["h"] = h
    return record


def plane(degree=4):
    return {
        "name": "fermat", "category": "examples", "form": "SmoothPlane",
        "field": {"type": "Prime", "p": "0x7", "bits": 3},
        "params": {"polynomial": [
            {"coeff": 1, "powers": [degree, 0, 0]},
            {"coeff": 1, "powers": [0, degree, 0]},
            {"coeff": 1, "powers": [0, 0, degree]},
        ]},
    }


def test_hyperelliptic_cover_counts_infinity_and_frobenius_orbits():
    curve = curve_from_dict(hyperelliptic())
    assert curve.genus() == 2
    branch = TRAITS["hyperelliptic_branch"](curve)
    assert branch == {
        "geometric_branch_points": 6, "rational_branch_points": 4,
        "branch_at_infinity": 1, "branch_closed_points": 5,
        "branch_orbit_degrees": [1, 1, 1, 1, 2],
    }


def test_even_degree_cover_is_unramified_at_infinity():
    curve = curve_from_dict(hyperelliptic([1, 0, 0, 0, 0, 0, 1]))
    branch = TRAITS["hyperelliptic_branch"](curve)
    assert branch["branch_at_infinity"] == 0
    assert branch["geometric_branch_points"] == 6
    assert branch["branch_orbit_degrees"] == [2, 2, 2]


def test_nonzero_h_uses_the_completed_square_discriminant():
    curve = curve_from_dict(hyperelliptic([-1, -1, 0, 0, 0, 1], [2]))
    assert TRAITS["hyperelliptic_branch"](curve) == TRAITS["hyperelliptic_branch"](
        curve_from_dict(hyperelliptic())
    )


def test_characteristic_two_geometry_does_not_assume_tame_branching():
    record = hyperelliptic([1, 0, 0, 1, 0, 1], [1])
    record["field"] = {"type": "Prime", "p": "0x2", "bits": 2}
    curve = curve_from_dict(record)
    assert TRAITS["canonical_geometry"](curve)["genus"] == 2
    branch = TRAITS["hyperelliptic_branch"]
    assert branch(curve) == dict.fromkeys(branch.OUTPUT)


def test_branch_factorization_is_over_the_actual_extension_field():
    record = hyperelliptic()
    record["field"] = {
        "type": "Extension", "base": "0x3", "degree": 2, "bits": 4,
        "poly": [{"coeff": "0x1", "power": "0x2"}, {"coeff": "0x1", "power": "0x0"}],
    }
    curve = curve_from_dict(record)
    assert TRAITS["hyperelliptic_branch"](curve)["branch_orbit_degrees"] == [1] * 6


@pytest.mark.parametrize("degree,genus", [(4, 3), (5, 6)])
def test_smooth_plane_models_have_canonical_embeddings(degree, genus):
    curve = curve_from_dict(plane(degree))
    geometry = TRAITS["canonical_geometry"](curve)
    assert geometry["genus"] == genus
    assert geometry["canonical_embedding"] == 1
    assert geometry["canonical_image_degree"] == 2 * genus - 2
    assert geometry["degree_two_map"] == 0
    assert geometry["plane_quartic"] == int(degree == 4)
    branch = TRAITS["hyperelliptic_branch"]
    assert branch(curve) == dict.fromkeys(branch.OUTPUT)


@pytest.mark.parametrize("record", [
    hyperelliptic([0, 0, 0, 0, 0, 1]),
    hyperelliptic([1, 1, 0, 1]),
    dict(hyperelliptic(), genus=3),
    dict(hyperelliptic(), order=101),
    dict(hyperelliptic(), order=0, cofactor=1),
    plane(3),
])
def test_invalid_models_are_rejected(record):
    with pytest.raises(ValueError):
        curve_from_dict(record)


def test_singular_plane_curve_is_rejected():
    record = plane()
    record["params"]["polynomial"].pop()
    with pytest.raises(ValueError, match="smooth"):
        curve_from_dict(record)


def test_nonhomogeneous_plane_curve_is_rejected():
    record = plane()
    record["params"]["polynomial"].append({"coeff": 1, "powers": [0, 0, 0]})
    with pytest.raises(ValueError, match="homogeneous"):
        curve_from_dict(record)


def test_database_round_trip_preserves_model_and_optional_metadata():
    record = dict(hyperelliptic(), genus=2)
    curve = curve_from_dict(_decode_ints(_format_curve(record)))
    assert curve.genus() == 2
    assert curve.order() is None
    assert curve.cofactor() is None

    class Collection:
        def insert_one(self, result):
            self.result = result
            return type("Acknowledged", (), {"acknowledged": True})()

    collection = Collection()
    assert store_trait_result({"trait_canonical_geometry": collection}, curve,
                              "canonical_geometry", {}, TRAITS["canonical_geometry"](curve))
    assert collection.result["curve"]["cofactor"] is None
    assert _decode_ints(collection.result)["result"]["genus"] == 2


def test_existing_elliptic_models_still_use_custom_curve():
    record = json.loads(EXAMPLES.read_text())["curves"][0]
    curve = curve_from_dict(record)
    assert isinstance(curve, CustomCurve)
    assert curve.ec().cardinality() == 5
    assert TRAITS["canonical_geometry"](curve)["jacobian_is_curve"] == 1
    assert TRAITS["cofactor"](curve) == {"order": 5, "cofactor": 1}


@pytest.mark.parametrize("trait_name", ["canonical_geometry", "hyperelliptic_branch", "cofactor"])
def test_compute_json_on_mixed_curve_families(trait_name):
    result = subprocess.run([
        sys.executable, "-m", "dissect.utils.compute_json", "-t", trait_name,
        "-i", str(EXAMPLES),
    ], check=True, capture_output=True, text=True)
    records = json.loads(result.stdout)["data"]
    assert len(records) == 3
    if trait_name == "canonical_geometry":
        assert [row["result"]["genus"] for row in records] == [1, 2, 3]
    elif trait_name == "hyperelliptic_branch":
        assert records[1]["result"]["geometric_branch_points"] == 6
        assert all(value is None for value in records[0]["result"].values())
    else:
        assert records[0]["result"] == {"order": 5, "cofactor": 1}
        assert all(value is None for value in records[1]["result"].values())
