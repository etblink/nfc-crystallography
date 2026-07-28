from __future__ import annotations

import json
from pathlib import Path

from nfc_cryst.evidence import load_evidence, markdown_table
from nfc_cryst.methods import verify_method_sources

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def test_method_source_bindings() -> None:
    assert verify_method_sources()["passed"]


def test_evidence_table_preserves_claim_boundaries() -> None:
    evidence = load_evidence()
    boundaries = set(evidence["claim_boundaries"])
    assert "D6C1_THROUGH_D7_TRANSFER_NOT_ESTABLISHED" in boundaries
    assert "NFC_AS_COSMOLOGY_OR_THEORY_OF_EVERYTHING_NOT_ESTABLISHED" in boundaries
    table = markdown_table(evidence)
    assert "8VTD" in table
    assert "9JZO" in table
    assert "6CKT" in table


def test_historical_gate_development_cases_are_resolved() -> None:
    for corpus in ("6ckt", "6tpi"):
        result = json.loads(
            (
                REPOSITORY_ROOT / f"results/{corpus}_invariant_gate_development.json"
            ).read_text(encoding="utf-8")
        )
        assert (
            result["principal_outcome"]
            == "CONVENTIONAL_FULL_AND_SPLIT_INDEXABILITY_STABLE"
        )
        assert (
            result["maximum_pairwise_invariant_orientation_difference_degrees"] < 0.04
        )


def test_repository_layout_replayed_8vtd_decision() -> None:
    result = json.loads(
        (REPOSITORY_ROOT / "results/8vtd_repository_baseline_replay.json").read_text(
            encoding="utf-8"
        )
    )
    assert result["decision"]["decision"] == "LATTICE_RECOVERED"
    assert result["decision"]["truth_consulted"] is False
    comparison = result["comparison_to_archived_result"]
    assert comparison["nonfloating_value_mismatches"] == 0
    assert comparison["maximum_absolute_floating_difference"] < 3e-15
