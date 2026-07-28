from __future__ import annotations

import json
from pathlib import Path
from typing import Any

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def load_evidence() -> dict[str, Any]:
    return json.loads(
        (REPOSITORY_ROOT / "results/evidence_table.json").read_text(encoding="utf-8")
    )


def markdown_table(document: dict[str, Any]) -> str:
    lines = [
        "| Corpus | Role | Bound method | Outcome | Interpretation |",
        "| --- | --- | --- | --- | --- |",
    ]
    for row in document["rows"]:
        lines.append(
            "| {corpus} | {role} | {method} | {outcome} | {interpretation} |".format(
                corpus=row["corpus"],
                role=row["role"].replace("_", " "),
                method=row["method"],
                outcome=row["outcome"].replace("_", " "),
                interpretation=row["interpretation"],
            )
        )
    return "\n".join(lines)
