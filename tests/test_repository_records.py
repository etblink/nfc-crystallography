from __future__ import annotations

import hashlib
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
    assert "4JX2" in table
    assert "4G2A" in table
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


def test_4jx2_manifest_preserves_prospective_commitment() -> None:
    manifest = json.loads(
        (REPOSITORY_ROOT / "data_manifests/4jx2.json").read_text(
            encoding="utf-8"
        )
    )
    assert (
        manifest["download"]["sha256"]
        == "1d833cd501fd86dd846f3a7eeb22062eb9917012a1c1cc37fc6800596c3c504a"
    )
    commitment = manifest["truth_free_commitment"]
    commitment_path = REPOSITORY_ROOT / commitment["file"]
    commitment_record = json.loads(commitment_path.read_text(encoding="utf-8"))
    assert hashlib.sha256(commitment_path.read_bytes()).hexdigest() == commitment[
        "file_sha256"
    ]
    assert commitment["decision"] == "INSUFFICIENT_SIGNAL"
    assert commitment["reason"] == "NO_PERSISTENT_FAMILY"
    assert (
        commitment["semantic_sha256"]
        == "2aac64f387f036f887146dc0fb1b5603dba680c1a18fbedf76b87be720fd6609"
    )
    assert commitment_record["semantic_sha256"] == commitment["semantic_sha256"]
    assert commitment_record["conventional_truth_accessed"] is False
    assert manifest["method"]["candidate_c_executed"] is False
    scoring = manifest["postcommitment_scoring"]
    assert scoring["classification"] == "CONSERVATIVE_MISSED_RECOVERY"
    assert scoring["baseline_b_decision_promoted"] is False


def test_4g2a_remains_a_preexecution_exclusion() -> None:
    exclusions = json.loads(
        (REPOSITORY_ROOT / "data_manifests/exclusions.json").read_text(
            encoding="utf-8"
        )
    )["exclusions"]
    exclusion = next(
        item for item in exclusions if item["identifier"] == "4G2A"
    )
    assert (
        exclusion["reason"]
        == "PREEXECUTION_UNSUPPORTED_FIXED_RAW_INTERFACE_NO_NFC_RUN"
    )
    assert exclusion["frame_count"] == 578
