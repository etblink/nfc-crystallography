#!/usr/bin/env python3
"""Build fixed D4.5-successor inputs for exploratory D5 development.

Every real feed is constructed from a raw-derived D4 interface without using
the conventional cell or orientation.  Synthetic truth is carried only in
the scoring metadata.  Full and split D4.5 representations are constructed
independently.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

from d45_aggregation import canonical_bytes, semantic_sha256
from d45_successor import aggregate_d45_successor
from prior_run_generalization import (
    frame_local_synthetic_case,
    load_gzip_json,
    primary_9jq9,
    primary_9z6f,
)


VALIDITY = "CERTIFIED_POSITIVE_SEPARATION_INTERIOR"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def load_portability(path: Path) -> Any:
    source = path.resolve()
    sys.path.insert(0, str(source.parent))
    spec = importlib.util.spec_from_file_location(
        "nfc_portability_fixed_feed_builder", source
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load frozen portability source")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def basis_from_dials_expt(path: Path, crystal_index: int) -> np.ndarray:
    record = json.loads(path.read_text(encoding="utf-8"))
    crystal = record["crystal"][crystal_index]
    direct = np.column_stack(
        [
            crystal["real_space_a"],
            crystal["real_space_b"],
            crystal["real_space_c"],
        ]
    )
    return np.linalg.inv(direct).T


def primary_public(
    corpus: dict[str, Any], portability: Any
) -> tuple[
    np.ndarray,
    np.ndarray,
    list[dict[str, Any]],
    dict[str, dict[str, Any]],
]:
    q, radius, ordered, _ = portability.compute_all_q(
        corpus["frames"], corpus["peaks"]
    )
    mask = np.asarray(
        [
            portability.FROZEN_D5.tier_for_peak(peak)
            == "PRIMARY_LATTICE_TIER"
            for peak in ordered
        ],
        dtype=bool,
    )
    peaks = [
        peak
        for peak, selected in zip(ordered, mask, strict=True)
        if selected
    ]
    frames = {frame["frame_id"]: frame for frame in corpus["frames"]}
    return q[mask], radius[mask], peaks, frames


def subset_interface(
    q: np.ndarray,
    radius: np.ndarray,
    peaks: list[dict[str, Any]],
    frames: dict[str, dict[str, Any]],
    allowed_frame_ids: set[str],
) -> tuple[
    np.ndarray,
    np.ndarray,
    list[dict[str, Any]],
    dict[str, dict[str, Any]],
]:
    mask = np.asarray(
        [peak["f"] in allowed_frame_ids for peak in peaks], dtype=bool
    )
    return (
        q[mask],
        radius[mask],
        [
            peak
            for peak, selected in zip(peaks, mask, strict=True)
            if selected
        ],
        {
            frame_id: frames[frame_id]
            for frame_id in sorted(allowed_frame_ids)
        },
    )


def ordered_frame_ids(
    frames: dict[str, dict[str, Any]]
) -> list[str]:
    by_scan: dict[str, list[dict[str, Any]]] = {}
    for frame in frames.values():
        by_scan.setdefault(str(frame["scan_configuration_id"]), []).append(
            frame
        )
    ordered: list[str] = []
    for scan_id in sorted(by_scan):
        local = by_scan[scan_id]
        increments = {
            float(frame["geometry"]["angle_increment_degrees"])
            for frame in local
        }
        signs = {1 if value > 0 else -1 for value in increments}
        if 0.0 in increments or len(signs) != 1:
            raise ValueError(f"invalid increment convention in {scan_id}")
        sign = signs.pop()
        ordered.extend(
            str(frame["frame_id"])
            for frame in sorted(
                local,
                key=lambda frame: (
                    sign
                    * float(frame["geometry"]["start_angle_degrees"]),
                    str(frame["frame_id"]),
                ),
            )
        )
    return ordered


def certified_feed(
    q: np.ndarray,
    radius: np.ndarray,
    peaks: list[dict[str, Any]],
    frames: dict[str, dict[str, Any]],
) -> tuple[dict[str, Any], dict[str, Any]]:
    aggregate_q, aggregate_peaks, construction = aggregate_d45_successor(
        q, radius, peaks, frames
    )
    mask = np.asarray(
        [
            peak["d45"]["disposition"]
            == "UNIQUE_FORMAL_CONSENSUS_AGGREGATE"
            and int(peak["d45"]["member_count"]) >= 2
            for peak in aggregate_peaks
        ],
        dtype=bool,
    )
    selected_q = aggregate_q[mask]
    selected_peaks = [
        peak
        for peak, keep in zip(aggregate_peaks, mask, strict=True)
        if keep
    ]
    used_frame_ids = sorted({str(peak["f"]) for peak in selected_peaks})
    feed = {
        "q_cycles_per_angstrom": selected_q.tolist(),
        "objects": [
            {
                "object_id": str(peak["peak_id"]),
                "representative_frame_id": str(peak["f"]),
                "scan_configuration_id": str(
                    frames[peak["f"]]["scan_configuration_id"]
                ),
                "signal_to_noise": float(
                    peak["integrated_signal_to_formal_noise"]
                ),
                "member_count": int(peak["d45"]["member_count"]),
                "member_peak_ids": list(peak["d45"]["member_peak_ids"]),
            }
            for peak in selected_peaks
        ],
        "frames": {
            frame_id: {
                "frame_id": frame_id,
                "scan_configuration_id": str(
                    frames[frame_id]["scan_configuration_id"]
                ),
            }
            for frame_id in used_frame_ids
        },
    }
    feed["q_sha256"] = hashlib.sha256(
        np.asarray(selected_q, dtype="<f8", order="C").tobytes(order="C")
    ).hexdigest()
    feed["membership_sha256"] = semantic_sha256(
        [
            {
                "aggregate_peak_id": item["object_id"],
                "member_peak_ids": item["member_peak_ids"],
            }
            for item in feed["objects"]
        ]
    )
    return feed, construction


def independently_constructed_splits(
    q: np.ndarray,
    radius: np.ndarray,
    peaks: list[dict[str, Any]],
    frames: dict[str, dict[str, Any]],
    groups: dict[str, set[str]],
) -> tuple[dict[str, Any], dict[str, Any]]:
    feeds: dict[str, Any] = {}
    constructions: dict[str, Any] = {}
    full_feed, full_construction = certified_feed(q, radius, peaks, frames)
    feeds["FULL"] = full_feed
    constructions["FULL"] = full_construction
    for label, frame_ids in groups.items():
        local = subset_interface(q, radius, peaks, frames, frame_ids)
        feeds[label], constructions[label] = certified_feed(*local)
    return feeds, constructions


def chronological_halves(
    frames: dict[str, dict[str, Any]]
) -> dict[str, set[str]]:
    ordered = ordered_frame_ids(frames)
    midpoint = len(ordered) // 2
    return {
        "HALF_A": set(ordered[:midpoint]),
        "HALF_B": set(ordered[midpoint:]),
    }


def scan_hash_halves(
    frames: dict[str, dict[str, Any]]
) -> dict[str, set[str]]:
    scans = sorted(
        {str(frame["scan_configuration_id"]) for frame in frames.values()}
    )
    scan_groups = {
        "HALF_A": {
            scan
            for scan in scans
            if hashlib.sha256(scan.encode("utf-8")).digest()[0] % 2 == 0
        },
        "HALF_B": {
            scan
            for scan in scans
            if hashlib.sha256(scan.encode("utf-8")).digest()[0] % 2 == 1
        },
    }
    return {
        label: {
            frame_id
            for frame_id, frame in frames.items()
            if str(frame["scan_configuration_id"]) in selected
        }
        for label, selected in scan_groups.items()
    }


def truth_metadata(
    basis: np.ndarray | None,
    role: str,
) -> dict[str, Any]:
    return {
        "basis": basis.tolist() if basis is not None else None,
        "role": role,
        "used_by_construction": False,
        "used_by_candidate_generation_or_decision": False,
    }


def load_prior_object_identity(path: Path) -> dict[str, Any]:
    record = load_gzip_json(path)
    objects = [
        item
        for item in record["objects"]
        if item["aggregate_peak"]["d45"]["disposition"]
        == "UNIQUE_FORMAL_CONSENSUS_AGGREGATE"
        and int(item["aggregate_peak"]["d45"]["member_count"]) >= 2
    ]
    q = np.asarray(
        [
            item["aggregate_q_cycles_per_angstrom"]
            for item in objects
        ],
        dtype=float,
    )
    return {
        "repeat_certified_object_count": len(objects),
        "q_sha256": hashlib.sha256(
            np.asarray(q, dtype="<f8", order="C").tobytes(order="C")
        ).hexdigest(),
        "membership_sha256": semantic_sha256(
            [
                {
                    "aggregate_peak_id": item["aggregate_peak"]["peak_id"],
                    "member_peak_ids": item["aggregate_peak"]["d45"][
                        "member_peak_ids"
                    ],
                }
                for item in objects
            ]
        ),
    }


def build_real_case(
    case_id: str,
    q: np.ndarray,
    radius: np.ndarray,
    peaks: list[dict[str, Any]],
    frames: dict[str, dict[str, Any]],
    groups: dict[str, set[str]],
    truth_basis: np.ndarray | None,
    required_outcome: str,
    role: str,
    prior_identity: dict[str, Any] | None = None,
) -> dict[str, Any]:
    print(f"BUILD {case_id}", flush=True)
    feeds, constructions = independently_constructed_splits(
        q, radius, peaks, frames, groups
    )
    result = {
        "case_id": case_id,
        "case_kind": "REAL_PUBLIC_OR_HISTORICAL_CORPUS",
        "role": role,
        "required_outcome": required_outcome,
        "truth": truth_metadata(truth_basis, role),
        "feeds": feeds,
        "constructions": constructions,
    }
    if prior_identity is not None:
        observed = {
            "repeat_certified_object_count": len(
                feeds["FULL"]["q_cycles_per_angstrom"]
            ),
            "q_sha256": feeds["FULL"]["q_sha256"],
            "membership_sha256": feeds["FULL"]["membership_sha256"],
        }
        result["prior_fixed_interface_identity_check"] = {
            "expected": prior_identity,
            "observed": observed,
            "exact_equal": observed == prior_identity,
        }
        if observed != prior_identity:
            raise RuntimeError(
                f"{case_id}: successor changed established positive-scan feed"
            )
    return result


def build_synthetic_case(
    case: dict[str, Any],
    persistence_fraction: float,
    truth_basis: np.ndarray,
) -> dict[str, Any]:
    (
        q,
        radius,
        peaks,
        frames,
        _latent_q,
        scan_by_observation,
    ) = frame_local_synthetic_case(case, persistence_fraction)
    member_to_latent = {
        peak["peak_id"]: peak["synthetic_latent_observation_id"]
        for peak in peaks
    }
    groups: dict[str, set[str]] = {}
    for label, scans in {
        "HALF_A": {0, 1},
        "HALF_B": {2, 3},
    }.items():
        latent = {
            observation_id
            for observation_id, scan_index in scan_by_observation.items()
            if scan_index in scans
        }
        groups[label] = {
            peak["f"]
            for peak in peaks
            if member_to_latent[peak["peak_id"]] in latent
        }
    feeds, constructions = independently_constructed_splits(
        q, radius, peaks, frames, groups
    )
    return {
        "case_id": str(case["case_id"]),
        "case_kind": "SYNTHETIC_D45_NATIVE_CONTROL",
        "persistence_fraction": persistence_fraction,
        "required_outcome": str(case["required_outcome"]),
        "truth": truth_metadata(
            truth_basis, "POST_CONSTRUCTION_SYNTHETIC_SCORING_ONLY"
        ),
        "feeds": feeds,
        "constructions": constructions,
    }


def write_gzip_json(path: Path, value: Any) -> None:
    raw = canonical_bytes(value)
    with path.open("wb") as handle:
        with gzip.GzipFile(
            filename="", mode="wb", fileobj=handle, mtime=0
        ) as compressed:
            compressed.write(raw)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    general = (
        root
        / "NFC_CRYST_D45_GENERALIZATION_SCAN_DIRECTION_AND_D5_ABSTENTION_EVALUATION_0_1_0"
    )
    causal = (
        root
        / "release_stage"
        / "NFC_CRYST_D45_REPEAT_CERTIFICATION_AND_UNCHANGED_FROZEN_D5_CAUSAL_EVALUATION_0_1_0"
    )
    portability = load_portability(
        root
        / "d45_d5_successor_work"
        / "scripts"
        / "pilatus_portability.py"
    )

    real_cases: list[dict[str, Any]] = []
    public_specs = [
        (
            "6GN2",
            causal / "inputs/6GN2_FULL_FROZEN_D4_D5_INPUT_CORPUS_0_1_0.json.gz",
            causal / "inputs/6GN2_INDEXED_SCORING_REFERENCE.expt",
            0,
            causal / "results/D45_6GN2_FULL_OBJECTS_0_1_0.json.gz",
        ),
        (
            "6GN3_SWEEP1",
            causal / "inputs/6GN3_SWEEP1_FROZEN_D4_D5_INPUT_CORPUS_0_1_0.json.gz",
            causal / "inputs/6GN3_INDEPENDENT_SWEEPS_INDEXED_SCORING_REFERENCE.expt",
            0,
            causal / "results/D45_6GN3_SWEEP1_OBJECTS_0_1_0.json.gz",
        ),
        (
            "6GN3_SWEEP2",
            causal / "inputs/6GN3_SWEEP2_FROZEN_D4_D5_INPUT_CORPUS_0_1_0.json.gz",
            causal / "inputs/6GN3_INDEPENDENT_SWEEPS_INDEXED_SCORING_REFERENCE.expt",
            1,
            causal / "results/D45_6GN3_SWEEP2_OBJECTS_0_1_0.json.gz",
        ),
    ]
    for case_id, corpus_path, expt_path, crystal_index, prior_path in public_specs:
        corpus = load_gzip_json(corpus_path)
        q, radius, peaks, frames = primary_public(corpus, portability)
        real_cases.append(
            build_real_case(
                case_id,
                q,
                radius,
                peaks,
                frames,
                chronological_halves(frames),
                basis_from_dials_expt(expt_path, crystal_index),
                "RECOVER_PRIMITIVE_LATTICE",
                "PUBLIC_GEOMETRIC_POSITIVE_CONTROL",
                load_prior_object_identity(prior_path),
            )
        )

    q, radius, peaks, frames, _ = primary_9z6f(
        general / "inputs/NFC_CRYST_9Z6F_D45_PRIMARY_INTERFACE_0_1_0.json.gz",
        general / "inputs/NFC_CRYST_9Z6F_D45_PRIMARY_Q_INTERFACE_0_1_0.json.gz",
    )
    historical = json.loads(
        (
            general
            / "inputs/NFC_CRYST_9Z6F_D5_RAW_RECIPROCAL_LATTICE_RESULT_0_1_0.json"
        ).read_text(encoding="utf-8")
    )
    basis_9z6f = np.asarray(
        historical["lattice_candidate"]["metric"][
            "reciprocal_basis_columns_cycles_per_angstrom"
        ],
        dtype=float,
    )
    prior_generalization = json.loads(
        (
            general
            / "results"
            / "NFC_CRYST_D45_SYNTHETIC_9Z6F_9JQ9_GENERALIZATION_AND_SCAN_DIRECTION_DIAGNOSTIC_RESULT_0_1_0.json"
        ).read_text(encoding="utf-8")
    )

    def prior_sign_aware_identity(corpus_key: str) -> dict[str, Any]:
        prior = prior_generalization[corpus_key][
            "scan_local_increment_sign_aware_diagnostic"
        ]
        return {
            "repeat_certified_object_count": prior[
                "repeat_certified_object_count"
            ],
            "q_sha256": prior["repeat_certified_q_sha256"],
            "membership_sha256": prior[
                "repeat_certified_membership_sha256"
            ],
        }

    real_cases.append(
        build_real_case(
            "9Z6F",
            q,
            radius,
            peaks,
            frames,
            scan_hash_halves(frames),
            basis_9z6f,
            "RECOVER_PRIMITIVE_LATTICE",
            "HISTORICAL_DEVELOPMENT_POSITIVE_CONTROL",
            prior_sign_aware_identity("9z6f"),
        )
    )

    q, radius, peaks, frames, deposited_record = primary_9jq9(
        general / "inputs/9JQ9_NFC_D4_D5_NUMERICAL_CORPUS_0_1_0.npz",
        general / "inputs/9JQ9_NFC_D4_D5_NUMERICAL_CORPUS_0_1_0.json",
    )
    del deposited_record
    real_cases.append(
        build_real_case(
            "9JQ9",
            q,
            radius,
            peaks,
            frames,
            chronological_halves(frames),
            None,
            "INSUFFICIENT_SIGNAL_OR_AMBIGUOUS_LATTICE",
            "UNRESOLVED_PUBLIC_ADVERSE_DIAGNOSTIC",
            prior_sign_aware_identity("9jq9"),
        )
    )

    synthetic = load_gzip_json(
        general / "inputs/SYNTHETIC_D5_CONTROL_CORPUS_0_1_0.json.gz"
    )
    synthetic_truth = np.asarray(
        synthetic["truth"][
            "reciprocal_basis_columns_cycles_per_angstrom"
        ],
        dtype=float,
    )
    synthetic_models: dict[str, list[dict[str, Any]]] = {}
    for label, fraction in {
        "ALL_OBJECTS_PERSIST": 1.0,
        "HASHED_FIFTY_PERCENT_PERSISTENCE": 0.5,
    }.items():
        cases = []
        for case in synthetic["cases"]:
            print(f"BUILD SYNTHETIC {label} {case['case_id']}", flush=True)
            cases.append(
                build_synthetic_case(case, fraction, synthetic_truth)
            )
        synthetic_models[label] = cases

    source_paths = [
        item[1] for item in public_specs
    ] + [
        general / "inputs/NFC_CRYST_9Z6F_D45_PRIMARY_INTERFACE_0_1_0.json.gz",
        general / "inputs/NFC_CRYST_9Z6F_D45_PRIMARY_Q_INTERFACE_0_1_0.json.gz",
        general / "inputs/9JQ9_NFC_D4_D5_NUMERICAL_CORPUS_0_1_0.npz",
        general / "inputs/9JQ9_NFC_D4_D5_NUMERICAL_CORPUS_0_1_0.json",
        general / "inputs/SYNTHETIC_D5_CONTROL_CORPUS_0_1_0.json.gz",
        general
        / "results"
        / "NFC_CRYST_D45_SYNTHETIC_9Z6F_9JQ9_GENERALIZATION_AND_SCAN_DIRECTION_DIAGNOSTIC_RESULT_0_1_0.json",
    ]
    result: dict[str, Any] = {
        "artifact_id": "NFC_CRYST_D45_SUCCESSOR_FIXED_CONTROL_FEEDS_0_1_0",
        "scientific_scope": "OPEN_EXPLORATORY_D45_D5_SUCCESSOR_DEVELOPMENT",
        "d45_successor_contract": {
            "scan_local": True,
            "increment_sign_aware": True,
            "cross_scan_association": False,
            "association": "FORMAL_REGION_OVERLAP_ONLY",
            "ambiguity_policy": "PRESERVE_UNRESOLVED",
            "centroid": "SIGNAL_WEIGHTED",
            "certification": (
                "UNIQUE_FORMAL_CONSENSUS_AGGREGATE_MEMBER_COUNT_GEQ_2"
            ),
            "cell_or_orientation_input": False,
            "processed_reflection_input": False,
        },
        "real_cases": real_cases,
        "synthetic_models": synthetic_models,
        "input_bindings": {
            path.name: {
                "byte_count": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path in source_paths
        },
        "truth_firewall": {
            "truth_used_by_d45": False,
            "truth_used_by_future_d5_candidate_generation": False,
            "truth_used_by_future_d5_decision": False,
            "truth_role": "POST_DECISION_SCORING_ONLY",
        },
    }
    body = dict(result)
    result["semantic_sha256"] = semantic_sha256(body)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_gzip_json(args.output, result)
    print(
        json.dumps(
            {
                "output": str(args.output),
                "sha256": sha256_file(args.output),
                "semantic_sha256": result["semantic_sha256"],
                "real_counts": {
                    case["case_id"]: {
                        label: len(feed["q_cycles_per_angstrom"])
                        for label, feed in case["feeds"].items()
                    }
                    for case in real_cases
                },
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
