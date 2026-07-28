from __future__ import annotations

import argparse
import json
from pathlib import Path

from nfc_cryst.canonical import write_canonical_json
from nfc_cryst.conventional_gate import GateThresholds, qualify_experiments
from nfc_cryst.evidence import load_evidence, markdown_table
from nfc_cryst.methods import verify_method_sources
from nfc_cryst.paths import release_root


def _qualify(args: argparse.Namespace) -> int:
    result = qualify_experiments(
        {
            "full": args.full,
            "half_a": args.half_a,
            "half_b": args.half_b,
        },
        {
            "full": args.assigned_full,
            "half_a": args.assigned_half_a,
            "half_b": args.assigned_half_b,
        },
        GateThresholds(
            minimum_assigned_fraction=args.minimum_assigned,
            maximum_relative_metric_difference=args.maximum_metric,
            maximum_orientation_difference_degrees=args.maximum_orientation,
            basis_search_bound=args.search_bound,
        ),
    )
    write_canonical_json(args.output, result)
    print(json.dumps(result, sort_keys=True))
    return 0 if result["principal_outcome"].endswith("_STABLE") else 2


def _evidence(args: argparse.Namespace) -> int:
    document = load_evidence()
    if args.format == "markdown":
        print(markdown_table(document))
    else:
        print(json.dumps(document, sort_keys=True, indent=2))
    return 0


def _verify_methods(_: argparse.Namespace) -> int:
    result = verify_method_sources()
    print(json.dumps(result, sort_keys=True, indent=2))
    return 0 if result["passed"] else 1


def _release_root(_: argparse.Namespace) -> int:
    print(release_root())
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="nfc-cryst")
    subparsers = parser.add_subparsers(dest="command", required=True)

    evidence = subparsers.add_parser("evidence", help="print the evidence table")
    evidence.add_argument("--format", choices=("markdown", "json"), default="markdown")
    evidence.set_defaults(handler=_evidence)

    verify = subparsers.add_parser(
        "verify-methods", help="verify bundled frozen-source identities"
    )
    verify.set_defaults(handler=_verify_methods)

    root = subparsers.add_parser(
        "release-root",
        help="print the installed or source-checkout release-payload root",
    )
    root.set_defaults(handler=_release_root)

    qualify = subparsers.add_parser(
        "qualify-conventional",
        help="compare full/split DIALS crystal models modulo basis changes",
    )
    qualify.add_argument("--full", type=Path, required=True)
    qualify.add_argument("--half-a", type=Path, required=True)
    qualify.add_argument("--half-b", type=Path, required=True)
    qualify.add_argument("--assigned-full", type=float, required=True)
    qualify.add_argument("--assigned-half-a", type=float, required=True)
    qualify.add_argument("--assigned-half-b", type=float, required=True)
    qualify.add_argument("--minimum-assigned", type=float, default=0.5)
    qualify.add_argument("--maximum-metric", type=float, default=0.05)
    qualify.add_argument("--maximum-orientation", type=float, default=2.0)
    qualify.add_argument("--search-bound", type=int, default=1)
    qualify.add_argument("--output", type=Path, required=True)
    qualify.set_defaults(handler=_qualify)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    return int(args.handler(args))


if __name__ == "__main__":
    raise SystemExit(main())
