"""Compute a complete higher-genus trait report with shared exact caches."""

import argparse
import json
import sys
from collections import Counter

from dissect.traits import TRAITS
from dissect.utils.curve_factory import curve_from_dict
from dissect.utils.database_handler import _cast_sage_types


REPORT_TRAITS = (
    "canonical_geometry", "hyperelliptic_branch", "frobenius_polynomial",
    "jacobian_decomposition", "hasse_witt", "newton_polygon", "weierstrass_gaps",
    "automorphism_structure", "moduli_invariants", "gonality_and_quotients", "cover_ramification",
    "genus2_rosenhain", "genus2_elliptic_quotients", "genus3_bitangents", "genus3_flexes",
    "genus4_canonical_quadric", "genus4_tritangents",
    "genus2_richelot", "hyperelliptic_twists", "jacobian_torsion",
    "symmetric_power_divisors", "frobenius_orders", "canonical_relations",
)


def curve_report(record, traits=REPORT_TRAITS, overrides=None):
    curve = curve_from_dict(record)
    results, parameters = {}, {}
    for name in traits:
        trait = TRAITS[name]
        options = {key: value for key, value in (overrides or {}).items() if key in trait.DEFAULT_PARAMS}
        parameters[name] = {key: values[0] for key, values in trait.DEFAULT_PARAMS.items()}
        parameters[name].update(options)
        value = trait(curve, **options)
        if "status" not in value:
            applicable = any(item is not None for item in value.values())
            value = dict(status="ok" if applicable else "not_applicable", reason=None, **value)
        results[name] = value
    statuses = Counter(value.get("status", "ok") for value in results.values())
    return {"curve": curve.name(), "family": curve.curve_type(), "genus": curve.genus(),
            "field_size": int(curve.q()), "field": record["field"],
            "model": {"form": record["form"], "params": record["params"]},
            "trait_parameters": parameters, "coverage": dict(statuses), "traits": results}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", help="Curve JSON file; defaults to stdin")
    parser.add_argument("-o", "--output", help="Report JSON file; defaults to stdout")
    parser.add_argument("-t", "--traits", nargs="+", choices=REPORT_TRAITS, default=REPORT_TRAITS)
    parser.add_argument("--params", type=json.loads, default={}, help="Global JSON parameter overrides, filtered by trait")
    args = parser.parse_args()
    if not isinstance(args.params, dict):
        parser.error("--params must be a JSON object")
    known = set().union(*(TRAITS[name].DEFAULT_PARAMS for name in args.traits))
    if set(args.params)-known:
        parser.error("Unknown report parameters: " + ", ".join(sorted(set(args.params)-known)))
    if args.input:
        with open(args.input) as stream:
            records = json.load(stream)
    else:
        records = json.load(sys.stdin)
    if not isinstance(records, list):
        records = records["curves"]
    from sage.env import SAGE_VERSION
    import platform

    result = {"schema_version": 1, "software": {"sage": SAGE_VERSION, "python": platform.python_version()},
              "data": [curve_report(record, args.traits, args.params) for record in records]}
    result = _cast_sage_types(result)
    if args.output:
        with open(args.output, "w") as stream:
            json.dump(result, stream, indent=2)
    else:
        json.dump(result, sys.stdout, indent=2)


if __name__ == "__main__":
    main()
