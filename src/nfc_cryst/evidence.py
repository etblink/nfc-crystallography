from __future__ import annotations

import json
from typing import Any

from nfc_cryst.paths import release_root


def load_evidence() -> dict[str, Any]:
    return json.loads(
        (release_root() / "results/evidence_table.json").read_text(encoding="utf-8")
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
