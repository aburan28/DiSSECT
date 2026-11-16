import argparse
import json
import sys

from dissect.traits import TRAITS
from dissect.utils.curve_factory import curve_from_dict
from dissect.utils.database_handler import _cast_sage_types


def main():
    parser = argparse.ArgumentParser(description="DiSSECT trait computation script.")
    parser.add_argument(
        "-t",
        "--trait",
        type=str,
        choices=list(TRAITS.keys()),
        help="Trait identifier",
        required=True,
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        help="Input curves file (stdin by default)",
        default=None,
    )
    parser.add_argument(
        "-o", "--output", type=str, help="Output file (stdout by default)", default=None
    )
    parser.add_argument("--params", type=json.loads, default={},
                        help='JSON object overriding trait parameters, e.g. {"timeout":60}')

    args = parser.parse_args()
    if not isinstance(args.params, dict):
        parser.error("--params must be a JSON object")
    unknown = set(args.params)-set(TRAITS[args.trait].DEFAULT_PARAMS)
    if unknown:
        parser.error("Unknown trait parameters: " + ", ".join(sorted(unknown)))
    if args.trait not in TRAITS:
        print(f"Trait {args.trait} is not implemented", file=sys.stderr)
        sys.exit(1)

    results = {"data": []}
    if args.input:
        with open(args.input, "r") as f:
            curves = json.load(f)
    else:
        curves = json.load(sys.stdin)

    if not isinstance(curves, list):
        curves = curves["curves"]

    for curve in curves:
        curve = curve_from_dict(curve)
        for params in TRAITS[args.trait].params_iter(args.params):
            result = {"curve": curve.name()}
            result["params"] = params
            result["result"] = TRAITS[args.trait](curve, **params)
            results["data"].append(result)

    json.dump(
        _cast_sage_types(results),
        sys.stdout if not args.output else open(args.output, "w"),
        indent=2,
    )


if __name__ == "__main__":
    main()
