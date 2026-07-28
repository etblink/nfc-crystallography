from __future__ import annotations

import json
from typing import Any

from nfc_cryst.canonical import sha256_file
from nfc_cryst.paths import release_root


def verify_method_sources() -> dict[str, Any]:
    root = release_root()
    records: list[dict[str, Any]] = []
    passed = True
    for manifest_path in (
        root / "methods/baseline_B/method.json",
        root / "methods/candidate_C_experimental/method.json",
    ):
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        source_results = []
        for item in manifest.get("bundled_sources", []):
            path = root / item["path"]
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
