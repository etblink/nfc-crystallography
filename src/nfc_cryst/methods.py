from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from nfc_cryst.canonical import sha256_file

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def verify_method_sources() -> dict[str, Any]:
    records: list[dict[str, Any]] = []
    passed = True
    for manifest_path in (
        REPOSITORY_ROOT / "methods/baseline_B/method.json",
        REPOSITORY_ROOT / "methods/candidate_C_experimental/method.json",
    ):
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        source_results = []
        for item in manifest.get("bundled_sources", []):
            path = REPOSITORY_ROOT / item["path"]
            actual = sha256_file(path) if path.is_file() else None
            matches = actual == item["sha256"]
            source_results.append(
                {
                    "path": item["path"],
                    "expected_sha256": item["sha256"],
                    "actual_sha256": actual,
                    "matches": matches,
                }
            )
            passed = passed and matches
        records.append(
            {
                "method_id": manifest["method_id"],
                "release_zip_sha256": manifest["release_zip_sha256"],
                "sources": source_results,
            }
        )
    return {"passed": passed, "methods": records}
