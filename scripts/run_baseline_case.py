#!/usr/bin/env python3
"""Run exact baseline-B evaluator source on a compact truth-free D4.5 case."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
RUNTIME = ROOT / "methods/baseline_B/runtime"


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode(
        "utf-8"
    )


def load_json(path: Path) -> Any:
    if path.suffix == ".gz":
        with gzip.open(path, "rb") as handle:
            return json.loads(handle.read())
    return json.loads(path.read_text(encoding="utf-8"))


def load_module(path: Path) -> Any:
    spec = importlib.util.spec_from_file_location("nfc_baseline_b_evaluator", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load evaluator: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    case = load_json(args.case)
    if case.get("truth", {}).get("basis") is not None:
        raise SystemExit("refusing case with a populated truth basis")
    evaluator = load_module(RUNTIME / "rebuild/evaluate.py")
    scripts = RUNTIME / "restored/candidate/scripts"
    frozen = load_json(scripts / "frozen_core/d5_reciprocal_rule_freeze_0_1_0.json")
    successor = load_json(scripts / "d5_successor_rule_0_1_0.json")
    rule = load_json(RUNTIME / "rebuild/d5_multiscale_bidirectional_rule_0_1_0.json")
    portability = evaluator.old.load_portability(scripts / "pilatus_portability.py")
    result = evaluator.evaluate_case(case, None, portability, frozen, successor, rule)
    if result["decision"].get("truth_consulted"):
        raise SystemExit("evaluator unexpectedly reports truth access")
    record = {
        "method": "BASELINE_B",
        "release_zip_sha256": (
            "ba18310a04a45c13f1fdf100599c77f2da9fa8ba2f43ec9e942719871b6edf48"
        ),
        "case_file_sha256": hashlib.sha256(args.case.read_bytes()).hexdigest(),
        "result": result,
    }
    record["semantic_sha256"] = hashlib.sha256(canonical_bytes(record)).hexdigest()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(canonical_bytes(record))
    print(
        json.dumps(
            {
                "decision": result["decision"],
                "semantic_sha256": record["semantic_sha256"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
