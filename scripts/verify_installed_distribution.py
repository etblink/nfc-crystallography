#!/usr/bin/env python3
"""Verify the installed wheel carries the complete frozen release payload."""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

from nfc_cryst.evidence import load_evidence
from nfc_cryst.methods import verify_method_sources
from nfc_cryst.paths import release_root


def load_module(path: Path) -> object:
    spec = importlib.util.spec_from_file_location(
        "nfc_raw_builder_installed_distribution_smoke",
        path,
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load module from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> int:
    root = release_root()
    if "_release" not in root.parts:
        raise RuntimeError(
            f"resolved checkout payload instead of wheel payload: {root}"
        )

    methods = verify_method_sources()
    if not methods["passed"]:
        raise RuntimeError("installed frozen-source verification failed")
    evidence = load_evidence()

    builder_path = (
        root / "methods/baseline_B/raw_pipeline/build_pilatus_d45_feed.py"
    )
    candidate_scripts = (
        root / "methods/baseline_B/runtime/restored/candidate/scripts"
    ).resolve()
    builder = load_module(builder_path)
    portability, d45, feed_builder, imported_directory = builder.import_sources(
        argparse.Namespace(candidate_scripts=candidate_scripts)
    )

    expected = {
        "build_fixed_control_feeds": "build_fixed_control_feeds.py",
        "d45_scan_direction_diagnostic": "d45_scan_direction_diagnostic.py",
        "d45_successor": "d45_successor.py",
        "pilatus_portability": "pilatus_portability.py",
        "prior_run_generalization": "prior_run_generalization.py",
    }
    if imported_directory != candidate_scripts:
        raise RuntimeError("raw builder imported from an unexpected directory")
    if Path(portability.__file__).resolve() != (
        candidate_scripts / expected["pilatus_portability"]
    ):
        raise RuntimeError("PILATUS portability module identity mismatch")
    if Path(d45.__file__).resolve() != (
        candidate_scripts / expected["d45_successor"]
    ):
        raise RuntimeError("D4.5 successor module identity mismatch")
    if Path(feed_builder.__file__).resolve() != (
        candidate_scripts / expected["build_fixed_control_feeds"]
    ):
        raise RuntimeError("fixed-feed builder module identity mismatch")
    for module_name, filename in expected.items():
        module_path = Path(sys.modules[module_name].__file__).resolve()
        if module_path != candidate_scripts / filename:
            raise RuntimeError(f"{module_name} imported from {module_path}")

    print(
        json.dumps(
            {
                "evidence_rows": len(evidence["rows"]),
                "method_sources_verified": sum(
                    len(method["sources"]) for method in methods["methods"]
                ),
                "principal_outcome": "INSTALLED_DISTRIBUTION_VERIFIED",
                "release_root": str(root),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
