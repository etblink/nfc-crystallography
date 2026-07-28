#!/usr/bin/env python3
"""Evaluate exploratory D4.5 with unchanged frozen D5 beyond public controls.

The campaign has three components:

1. Project the existing D5-boundary synthetic suite into deterministic,
   frame-local repeat observations and rerun the published D4.5 construction.
2. Apply the published D4.5 construction to the compact 9Z6F and 9JQ9 D4
   corpora.
3. Run a separately labelled scan-local, increment-sign-aware diagnostic to
   isolate the frame-order failure found on negative-increment 9Z6F scans.

No D5 source or rule is modified.  No D6 or D7 operation is run.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import math
import platform
import signal
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Callable

import numpy as np
import scipy

from d45_aggregation import aggregate_d4, canonical_bytes, semantic_sha256
from d45_scan_direction_diagnostic import (
    aggregate_d4_scan_direction_aware,
)


FROZEN_BINDINGS = {
    "published_d45_source_sha256": (
        "856f21d2b5f3fe5f7c14e6ec93c5a7702424ccef892c34d9e0e50621a8d33f4c"
    ),
    "portability_source_sha256": (
        "a04aadc1b4cf2c7c24b45bedb8d8cd8c36ed84ed9477ab0806c4a601c2c5b585"
    ),
    "frozen_d4_source_sha256": (
        "909ce94e7519f3c10ebe4d874a8ca03697d1d4bee99a0dc54f22ef9a0508be12"
    ),
    "frozen_d4_rule_sha256": (
        "991ab9ec59fc325a177b18397380a1e87b31ea19e04ff76a07585a22b6e02708"
    ),
    "frozen_d5_source_sha256": (
        "740a5fb62694dd08ab9499ddec7fe834b07f955e61a389d703972f5a52a684fd"
    ),
    "frozen_d5_rule_sha256": (
        "16497549cc55344c9a31e9b0e65919a2186f4c1ff87e08b73987bada9343555a"
    ),
}

VALIDITY = "CERTIFIED_POSITIVE_SEPARATION_INTERIOR"
SYNTHETIC_FORMAL_RADIUS = 0.003
SYNTHETIC_DETECTION_OFFSET = 0.0006
SYNTHETIC_OBJECTS_PER_FRAME_BLOCK = 16
DEFAULT_D5_FIT_BUDGET_SECONDS = 60


class D5FitBudgetExceeded(TimeoutError):
    """Raised by the evaluation wrapper, not by the frozen D5 method."""


def _budget_alarm_handler(signum: int, frame: Any) -> None:
    del signum, frame
    raise D5FitBudgetExceeded("declared per-fit execution budget exceeded")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(canonical_bytes(value))


def seal(record: dict[str, Any]) -> dict[str, Any]:
    body = dict(record)
    body.pop("semantic_sha256", None)
    record["semantic_sha256"] = semantic_sha256(body)
    return record


def load_gzip_json(path: Path) -> Any:
    with gzip.open(path, "rb") as handle:
        return json.loads(handle.read())


def load_portability(path: Path) -> Any:
    source = path.resolve()
    sys.path.insert(0, str(source.parent))
    spec = importlib.util.spec_from_file_location(
        "nfc_frozen_portability_generalization", source
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load frozen portability source")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def verify_runtime(
    portability_source: Path,
    published_d45_source: Path,
    d4_source: Path,
    d4_rule: Path,
    d5_source: Path,
    d5_rule: Path,
) -> dict[str, str]:
    paths = {
        "published_d45_source_sha256": published_d45_source,
        "portability_source_sha256": portability_source,
        "frozen_d4_source_sha256": d4_source,
        "frozen_d4_rule_sha256": d4_rule,
        "frozen_d5_source_sha256": d5_source,
        "frozen_d5_rule_sha256": d5_rule,
    }
    observed = {key: sha256_file(path) for key, path in paths.items()}
    if observed != FROZEN_BINDINGS:
        raise RuntimeError(f"scientific source identity mismatch: {observed}")
    return observed


def lattice_equivalence(
    truth_basis: np.ndarray,
    candidate_basis: np.ndarray,
) -> dict[str, Any]:
    transform = np.linalg.inv(truth_basis) @ candidate_basis
    nearest = np.rint(transform).astype(np.int64)
    deviation = float(np.max(np.abs(transform - nearest)))
    determinant = int(round(float(np.linalg.det(nearest))))
    return {
        "nearest_integer_basis_transform": nearest.tolist(),
        "maximum_transform_integer_deviation": deviation,
        "integer_transform_determinant": determinant,
        "primitive_lattice_equivalent": bool(
            deviation <= 0.03 and abs(determinant) == 1
        ),
        "criterion": (
            "MAX_TRANSFORM_INTEGER_DEVIATION_LEQ_0.03_AND_"
            "ABS_DETERMINANT_EQ_1"
        ),
    }


def score_basis(q: np.ndarray, basis: np.ndarray) -> dict[str, Any]:
    q = np.asarray(q, dtype=float)
    if not len(q):
        return {"count": 0, "status": "NO_POINTS"}
    h = np.rint(q @ np.linalg.inv(basis).T).astype(np.int64)
    residual = np.linalg.norm(q - h @ basis.T, axis=1)
    return {
        "count": len(q),
        "residual_quantiles_cycles_per_angstrom": {
            "p50": float(np.quantile(residual, 0.50)),
            "p90": float(np.quantile(residual, 0.90)),
            "p99": float(np.quantile(residual, 0.99)),
        },
        "support_fraction": {
            "0.002": float(np.mean(residual <= 0.002)),
            "0.003": float(np.mean(residual <= 0.003)),
            "0.006": float(np.mean(residual <= 0.006)),
            "0.012": float(np.mean(residual <= 0.012)),
        },
        "maximum_abs_index": int(np.max(np.abs(h))),
    }


def descriptive_cell_comparison(
    metric: dict[str, Any],
    reference_cell: dict[str, Any],
) -> dict[str, Any]:
    candidate_lengths = np.sort(
        np.asarray(metric["direct_lengths_angstrom_length_sorted"], dtype=float)
    )
    reference_lengths = np.sort(
        np.asarray(reference_cell["lengths_angstrom"], dtype=float)
    )
    candidate_angles = np.sort(
        np.asarray(
            metric["direct_angles_degrees_alpha_beta_gamma"], dtype=float
        )
    )
    reference_angles = np.sort(
        np.asarray(reference_cell["angles_degrees"], dtype=float)
    )
    relative_length_error = (
        candidate_lengths - reference_lengths
    ) / reference_lengths
    volume = float(metric["direct_cell_volume_angstrom_cubed"])
    reference_volume = float(reference_cell["volume_angstrom_cubed"])
    return {
        "comparison_status": (
            "DESCRIPTIVE_SORTED_CELL_COMPARISON_NOT_A_BASIS_EQUIVALENCE_TEST"
        ),
        "candidate_lengths_angstrom_sorted": candidate_lengths.tolist(),
        "reference_lengths_angstrom_sorted": reference_lengths.tolist(),
        "rms_relative_length_error": float(
            np.sqrt(np.mean(relative_length_error**2))
        ),
        "candidate_angles_degrees_sorted": candidate_angles.tolist(),
        "reference_angles_degrees_sorted": reference_angles.tolist(),
        "rms_angle_error_degrees": float(
            np.sqrt(np.mean((candidate_angles - reference_angles) ** 2))
        ),
        "relative_volume_error": float(
            (volume - reference_volume) / reference_volume
        ),
    }


def compact_seed(record: dict[str, Any]) -> dict[str, Any]:
    return {
        "balanced_peak_count": record["balanced_peak_count"],
        "difference_vector_count": record["difference_vector_count"],
        "voxel_width_cycles_per_angstrom": record[
            "voxel_width_cycles_per_angstrom"
        ],
        "maximum_voxel_count": record["maximum_voxel_count"],
        "mode_minimum_bin_count": record["mode_minimum_bin_count"],
        "deduplicated_mode_count": record["deduplicated_mode_count"],
        "selected_seed_vectors": record["selected_seed_vectors"],
        "selected_seed_determinant": record["selected_seed_determinant"],
        "first_50_modes_by_norm": record["deduplicated_modes"][:50],
    }


def run_frozen_d5(
    q: np.ndarray,
    peaks: list[dict[str, Any]],
    frames: dict[str, dict[str, Any]],
    portability: Any,
    d5_rule: dict[str, Any],
    truth_basis: np.ndarray | None = None,
    reference_cell: dict[str, Any] | None = None,
    fit_budget_seconds: int = DEFAULT_D5_FIT_BUDGET_SECONDS,
) -> dict[str, Any]:
    q = np.asarray(q, dtype=float)
    if len(q) == 0:
        return {
            "execution_status": "D45_NO_REPEAT_CERTIFIED_OBJECTS",
            "candidate_returned": False,
            "structured_frozen_d5_abstention": False,
            "input_object_count": 0,
            "interpretation": (
                "UPSTREAM_D45_EMPTY_FEED_NOT_A_FROZEN_D5_ABSTENTION_RULE"
            ),
        }
    if fit_budget_seconds <= 0:
        raise ValueError("fit budget must be positive")
    prior_handler = signal.getsignal(signal.SIGALRM)
    signal.signal(signal.SIGALRM, _budget_alarm_handler)
    signal.setitimer(signal.ITIMER_REAL, float(fit_budget_seconds))
    try:
        seed, seed_record = portability.FROZEN_D5.seed_lattice(
            q, peaks, frames, d5_rule
        )
        basis, refinement = portability.FROZEN_D5.refine_lattice(
            q, seed, d5_rule
        )
        direct = np.linalg.inv(basis).T
        metric = portability.FROZEN_D5.cell_metric(direct, basis)
        result: dict[str, Any] = {
            "execution_status": "CANDIDATE_RETURNED",
            "candidate_returned": True,
            "structured_frozen_d5_abstention": False,
            "input_object_count": len(q),
            "basis": basis.tolist(),
            "metric": metric,
            "score_on_input_q": score_basis(q, basis),
            "seed_search": compact_seed(seed_record),
            "refinement": refinement,
        }
        if truth_basis is not None:
            result["truth_equivalence"] = lattice_equivalence(
                np.asarray(truth_basis, dtype=float), basis
            )
            result["truth_score_on_input_q"] = score_basis(q, truth_basis)
        if reference_cell is not None:
            result["reference_cell_comparison"] = (
                descriptive_cell_comparison(metric, reference_cell)
            )
        return result
    except D5FitBudgetExceeded:
        return {
            "execution_status": "EVALUATION_BUDGET_EXCEEDED",
            "candidate_returned": False,
            "structured_frozen_d5_abstention": False,
            "input_object_count": len(q),
            "evaluation_budget_seconds": fit_budget_seconds,
            "interpretation": (
                "OPERATIONAL_NONCOMPLETION_NOT_A_FROZEN_D5_ABSTENTION_"
                "OR_SCIENTIFIC_FAILURE"
            ),
        }
    except Exception as error:
        return {
            "execution_status": "UNSTRUCTURED_EXCEPTION",
            "candidate_returned": False,
            "structured_frozen_d5_abstention": False,
            "input_object_count": len(q),
            "exception_type": type(error).__name__,
            "exception_message": str(error),
        }
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0.0)
        signal.signal(signal.SIGALRM, prior_handler)


def progress(stage: str, **details: Any) -> None:
    print(
        json.dumps({"stage": stage, **details}, sort_keys=True),
        flush=True,
    )


def repeat_certified(
    aggregate_q: np.ndarray,
    aggregate_peaks: list[dict[str, Any]],
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    mask = np.asarray(
        [
            peak["d45"]["disposition"]
            == "UNIQUE_FORMAL_CONSENSUS_AGGREGATE"
            and int(peak["d45"]["member_count"]) >= 2
            for peak in aggregate_peaks
        ],
        dtype=bool,
    )
    return (
        aggregate_q[mask],
        [
            peak
            for peak, selected in zip(
                aggregate_peaks, mask, strict=True
            )
            if selected
        ],
    )


Aggregator = Callable[
    [
        np.ndarray,
        np.ndarray,
        list[dict[str, Any]],
        dict[str, dict[str, Any]],
    ],
    tuple[np.ndarray, list[dict[str, Any]], dict[str, Any]],
]


def evaluate_d45_feed(
    q: np.ndarray,
    radius: np.ndarray,
    peaks: list[dict[str, Any]],
    frames: dict[str, dict[str, Any]],
    aggregator: Aggregator,
    portability: Any,
    d5_rule: dict[str, Any],
    truth_basis: np.ndarray | None = None,
    reference_cell: dict[str, Any] | None = None,
) -> dict[str, Any]:
    aggregate_q, aggregate_peaks, construction = aggregator(
        q, radius, peaks, frames
    )
    certified_q, certified_peaks = repeat_certified(
        aggregate_q, aggregate_peaks
    )
    certified_member_count = sum(
        int(peak["d45"]["member_count"]) for peak in certified_peaks
    )
    return {
        "construction": construction,
        "all_d45_object_count": len(aggregate_q),
        "repeat_certified_object_count": len(certified_q),
        "repeat_certified_input_member_count": certified_member_count,
        "repeat_certified_q_sha256": hashlib.sha256(
            np.asarray(certified_q, dtype="<f8", order="C").tobytes(
                order="C"
            )
        ).hexdigest(),
        "repeat_certified_membership_sha256": semantic_sha256(
            [
                {
                    "aggregate_peak_id": peak["peak_id"],
                    "member_peak_ids": peak["d45"]["member_peak_ids"],
                }
                for peak in certified_peaks
            ]
        ),
        "frozen_d5": run_frozen_d5(
            certified_q,
            certified_peaks,
            frames,
            portability,
            d5_rule,
            truth_basis,
            reference_cell,
        ),
    }


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
    local_peaks = [
        peak
        for peak, selected in zip(peaks, mask, strict=True)
        if selected
    ]
    return (
        q[mask],
        radius[mask],
        local_peaks,
        {
            frame_id: frames[frame_id]
            for frame_id in sorted(allowed_frame_ids)
        },
    )


def primary_9z6f(
    d4_corpus_path: Path,
    q_corpus_path: Path,
) -> tuple[
    np.ndarray,
    np.ndarray,
    list[dict[str, Any]],
    dict[str, dict[str, Any]],
    dict[str, Any],
]:
    d4 = load_gzip_json(d4_corpus_path)
    q_corpus = load_gzip_json(q_corpus_path)
    if q_corpus["record_count"] != len(d4["peaks"]):
        raise RuntimeError("9Z6F D4/Q record-count mismatch")
    peak_by_id = {peak["peak_id"]: peak for peak in d4["peaks"]}
    records = sorted(
        (
            record
            for record in q_corpus["records"]
            if record["tier"] == "PRIMARY_LATTICE_TIER"
        ),
        key=lambda item: item["peak_id"],
    )
    if any(record["peak_id"] not in peak_by_id for record in records):
        raise RuntimeError("9Z6F Q record has no D4 peak")
    peaks = [peak_by_id[record["peak_id"]] for record in records]
    frames = {frame["frame_id"]: frame for frame in d4["frames"]}
    return (
        np.asarray(
            [record["q_cycles_per_angstrom"] for record in records],
            dtype=float,
        ),
        np.asarray(
            [
                record["formal_q_region_corner_radius"]
                for record in records
            ],
            dtype=float,
        ),
        peaks,
        frames,
        q_corpus.get("source_full_tier_counts", q_corpus["tier_counts"]),
    )


def primary_9jq9(
    array_path: Path,
    record_path: Path,
) -> tuple[
    np.ndarray,
    np.ndarray,
    list[dict[str, Any]],
    dict[str, dict[str, Any]],
    dict[str, Any],
]:
    record = json.loads(record_path.read_text(encoding="utf-8"))
    if sha256_file(array_path) != record["array_object"]["sha256"]:
        raise RuntimeError("9JQ9 numerical array hash mismatch")
    with np.load(array_path, allow_pickle=False) as archive:
        arrays = {key: archive[key] for key in archive.files}
    mask = arrays["tier_code"] == 0
    frame_by_ordinal = {
        int(frame["ordinal"]): frame for frame in record["frame_records"]
    }
    frames = {
        str(frame["frame_id"]): {
            "frame_id": str(frame["frame_id"]),
            "source_archive_path": frame["source_archive_path"],
            "source_file_sha256": frame["source_file_sha256"],
            "geometry": frame["geometry"],
            "scan_configuration_id": frame["scan_configuration_id"],
            "detection_summary": frame["detection_summary"],
        }
        for frame in record["frame_records"]
    }
    indices = np.flatnonzero(mask)
    peaks: list[dict[str, Any]] = []
    for index in indices:
        ordinal = int(arrays["frame_ordinal"][index])
        frame = frame_by_ordinal[ordinal]
        peaks.append(
            {
                "peak_id": arrays["peak_id"][index].decode("ascii"),
                "f": str(frame["frame_id"]),
                "x": float(arrays["x"][index]),
                "y": float(arrays["y"][index]),
                "Bhat": float(arrays["background"][index]),
                "Shat": float(arrays["signal"][index]),
                "sigma_S": float(arrays["sigma"][index]),
                "integrated_signal_to_formal_noise": float(
                    arrays["snr"][index]
                ),
                "q": {"validity": VALIDITY},
            }
        )
    if len(peaks) != int(
        record["replay_checks"]["tier_counts"]["PRIMARY_LATTICE_TIER"]
    ):
        raise RuntimeError("9JQ9 primary-tier count mismatch")
    return (
        arrays["q"][mask],
        arrays["radius"][mask],
        peaks,
        frames,
        record,
    )


def compare_feed_identities(
    left: dict[str, Any],
    right: dict[str, Any],
) -> dict[str, Any]:
    return {
        "repeat_certified_object_count_equal": (
            left["repeat_certified_object_count"]
            == right["repeat_certified_object_count"]
        ),
        "repeat_certified_q_sha256_equal": (
            left["repeat_certified_q_sha256"]
            == right["repeat_certified_q_sha256"]
        ),
        "repeat_certified_membership_sha256_equal": (
            left["repeat_certified_membership_sha256"]
            == right["repeat_certified_membership_sha256"]
        ),
        "frozen_d5_basis_equal": (
            left["frozen_d5"].get("basis")
            == right["frozen_d5"].get("basis")
        ),
    }


def evaluate_9z6f(
    q: np.ndarray,
    radius: np.ndarray,
    peaks: list[dict[str, Any]],
    frames: dict[str, dict[str, Any]],
    historical_result: dict[str, Any],
    portability: Any,
    d5_rule: dict[str, Any],
) -> dict[str, Any]:
    truth_basis = np.asarray(
        historical_result["lattice_candidate"]["metric"][
            "reciprocal_basis_columns_cycles_per_angstrom"
        ],
        dtype=float,
    )
    reference_cell = {
        "lengths_angstrom": [18.4619, 22.3151, 40.7662],
        "angles_degrees": [90.5380, 90.0, 90.0],
        "volume_angstrom_cubed": 16794.08,
        "role": "RETROSPECTIVE_REFERENCE_ONLY",
    }
    baseline = run_frozen_d5(
        q,
        peaks,
        frames,
        portability,
        d5_rule,
        truth_basis,
        reference_cell,
    )
    published = evaluate_d45_feed(
        q,
        radius,
        peaks,
        frames,
        aggregate_d4,
        portability,
        d5_rule,
        truth_basis,
        reference_cell,
    )
    sign_aware = evaluate_d45_feed(
        q,
        radius,
        peaks,
        frames,
        aggregate_d4_scan_direction_aware,
        portability,
        d5_rule,
        truth_basis,
        reference_cell,
    )

    scan_ids = sorted(
        {frame["scan_configuration_id"] for frame in frames.values()}
    )
    groups = {
        "SCAN_CONFIGURATION_HASH_PARITY_A": {
            scan_id
            for scan_id in scan_ids
            if hashlib.sha256(scan_id.encode("utf-8")).digest()[0] % 2 == 0
        },
        "SCAN_CONFIGURATION_HASH_PARITY_B": {
            scan_id
            for scan_id in scan_ids
            if hashlib.sha256(scan_id.encode("utf-8")).digest()[0] % 2 == 1
        },
    }
    split_results: dict[str, Any] = {}
    for label, allowed_scans in groups.items():
        allowed_frames = {
            frame_id
            for frame_id, frame in frames.items()
            if frame["scan_configuration_id"] in allowed_scans
        }
        local = subset_interface(
            q, radius, peaks, frames, allowed_frames
        )
        split_results[label] = {
            "scan_configuration_count": len(allowed_scans),
            "frame_count": len(allowed_frames),
            "diagnostic": evaluate_d45_feed(
                *local,
                aggregate_d4_scan_direction_aware,
                portability,
                d5_rule,
                truth_basis,
                reference_cell,
            ),
        }

    baseline_basis_equal = (
        baseline.get("basis")
        == historical_result["lattice_candidate"]["metric"][
            "reciprocal_basis_columns_cycles_per_angstrom"
        ]
    )
    return {
        "corpus": "9Z6F_SPECIMEN_A_DEVELOPMENT",
        "input_primary_peak_count": len(q),
        "frame_count": len(frames),
        "scan_configuration_count": len(scan_ids),
        "angle_increment_values_degrees": sorted(
            {
                float(frame["geometry"]["angle_increment_degrees"])
                for frame in frames.values()
            }
        ),
        "unchanged_frozen_d5_per_frame_baseline": baseline,
        "baseline_basis_byte_level_json_equal_to_historical": (
            baseline_basis_equal
        ),
        "published_d45_0_1_0": published,
        "scan_local_increment_sign_aware_diagnostic": sign_aware,
        "diagnostic_split_results": split_results,
        "historical_reference_binding": {
            "historical_result_semantic_sha256": historical_result[
                "canonical_semantic_sha256"
            ],
            "historical_basis": truth_basis.tolist(),
            "reference_cell": reference_cell,
        },
    }


def evaluate_9jq9(
    q: np.ndarray,
    radius: np.ndarray,
    peaks: list[dict[str, Any]],
    frames: dict[str, dict[str, Any]],
    corpus_record: dict[str, Any],
    historical_diagnostic: dict[str, Any],
    deposited_characterization: dict[str, Any],
    portability: Any,
    d5_rule: dict[str, Any],
) -> dict[str, Any]:
    metadata = deposited_characterization["crystallographic_metadata"]
    reference_cell = {
        "lengths_angstrom": metadata["cell_lengths_angstrom"],
        "angles_degrees": metadata["cell_angles_degrees"],
        "volume_angstrom_cubed": metadata["cell_volume_angstrom_cubed"],
        "role": (
            "DEPOSITED_MERGED_REFERENCE_POST_CONSTRUCTION_DIAGNOSTIC_ONLY"
        ),
    }
    baseline = run_frozen_d5(
        q,
        peaks,
        frames,
        portability,
        d5_rule,
        reference_cell=reference_cell,
    )
    published = evaluate_d45_feed(
        q,
        radius,
        peaks,
        frames,
        aggregate_d4,
        portability,
        d5_rule,
        reference_cell=reference_cell,
    )
    sign_aware = evaluate_d45_feed(
        q,
        radius,
        peaks,
        frames,
        aggregate_d4_scan_direction_aware,
        portability,
        d5_rule,
        reference_cell=reference_cell,
    )
    frame_records = sorted(
        corpus_record["frame_records"], key=lambda item: int(item["ordinal"])
    )
    midpoint = len(frame_records) // 2
    split_frame_ids = {
        "FRAMES_1_THROUGH_450": {
            str(item["frame_id"]) for item in frame_records[:midpoint]
        },
        "FRAMES_451_THROUGH_900": {
            str(item["frame_id"]) for item in frame_records[midpoint:]
        },
    }
    split_results: dict[str, Any] = {}
    for label, allowed_frames in split_frame_ids.items():
        local = subset_interface(q, radius, peaks, frames, allowed_frames)
        split_results[label] = evaluate_d45_feed(
            *local,
            aggregate_d4,
            portability,
            d5_rule,
            reference_cell=reference_cell,
        )

    baseline_metric_equal = (
        baseline.get("metric") == historical_diagnostic["d5"]["metric"]
    )
    return {
        "corpus": "9JQ9_PERMANENTLY_EXPLORATORY_ADVERSE_DIAGNOSTIC",
        "input_primary_peak_count": len(q),
        "frame_count": len(frames),
        "scan_configuration_count": len(
            {frame["scan_configuration_id"] for frame in frames.values()}
        ),
        "angle_increment_values_degrees": sorted(
            {
                float(frame["geometry"]["angle_increment_degrees"])
                for frame in frames.values()
            }
        ),
        "unchanged_frozen_d5_per_frame_baseline": baseline,
        "baseline_metric_equal_to_historical": baseline_metric_equal,
        "published_d45_0_1_0": published,
        "scan_local_increment_sign_aware_diagnostic": sign_aware,
        "published_and_sign_aware_identity_comparison": (
            compare_feed_identities(published, sign_aware)
        ),
        "published_d45_independent_half_results": split_results,
        "deposited_reference_binding": {
            "characterization_content_semantic_sha256": semantic_sha256(
                deposited_characterization
            ),
            "reference_cell": reference_cell,
            "provenance_status": "UNRESOLVED_FROM_PUBLIC_EVIDENCE",
        },
    }


def hashed_unit_vector(identifier: str) -> np.ndarray:
    digest = hashlib.sha256(identifier.encode("utf-8")).digest()
    values = np.asarray(
        [
            int.from_bytes(digest[offset : offset + 8], "big")
            / (2**64 - 1)
            * 2.0
            - 1.0
            for offset in (0, 8, 16)
        ],
        dtype=float,
    )
    norm = float(np.linalg.norm(values))
    if norm == 0.0:
        raise RuntimeError("hash-derived direction is zero")
    return values / norm


def hash_persistent(identifier: str, fraction: float) -> bool:
    value = int.from_bytes(
        hashlib.sha256(identifier.encode("utf-8")).digest()[:8], "big"
    )
    return value / (2**64) < fraction


def frame_local_synthetic_case(
    case: dict[str, Any],
    persistence_fraction: float,
) -> tuple[
    np.ndarray,
    np.ndarray,
    list[dict[str, Any]],
    dict[str, dict[str, Any]],
    dict[str, np.ndarray],
    dict[str, int],
]:
    q_rows: list[np.ndarray] = []
    radii: list[float] = []
    peaks: list[dict[str, Any]] = []
    frames: dict[str, dict[str, Any]] = {}
    latent_q: dict[str, np.ndarray] = {}
    scan_by_observation: dict[str, int] = {}

    observations_by_scan: dict[int, list[dict[str, Any]]] = {
        index: [] for index in range(4)
    }
    for observation in case["observations"]:
        observations_by_scan[int(observation["scan_index"])].append(
            observation
        )

    for scan_index in range(4):
        observations = observations_by_scan[scan_index]
        scan_id = f"SYNTHETIC_SCAN_{scan_index}"
        scan_base = float(scan_index * 10000)
        for block_index in range(
            math.ceil(
                len(observations) / SYNTHETIC_OBJECTS_PER_FRAME_BLOCK
            )
        ):
            for offset in range(3):
                frame_id = (
                    f"FRAME_S{scan_index}_B{block_index:04d}_O{offset}"
                )
                frames[frame_id] = {
                    "frame_id": frame_id,
                    "scan_configuration_id": scan_id,
                    "geometry": {
                        "start_angle_degrees": (
                            scan_base + block_index * 4.0 + offset
                        ),
                        "angle_increment_degrees": 1.0,
                    },
                }
            local = observations[
                block_index * SYNTHETIC_OBJECTS_PER_FRAME_BLOCK :
                (block_index + 1) * SYNTHETIC_OBJECTS_PER_FRAME_BLOCK
            ]
            for observation in local:
                observation_id = str(observation["observation_id"])
                center = np.asarray(observation["q"], dtype=float)
                latent_q[observation_id] = center
                scan_by_observation[observation_id] = scan_index
                persistent = hash_persistent(
                    observation_id, persistence_fraction
                )
                offsets = (-1, 0, 1) if persistent else (0,)
                direction = hashed_unit_vector(observation_id)
                signal = float(observation["synthetic_signal_to_noise"])
                for offset in offsets:
                    frame_id = (
                        f"FRAME_S{scan_index}_B{block_index:04d}_"
                        f"O{offset + 1}"
                    )
                    member_token = f"{observation_id}:{offset}"
                    peak_id = (
                        "SYN45_"
                        + hashlib.sha256(
                            member_token.encode("utf-8")
                        ).hexdigest()[:28]
                    )
                    vector = (
                        center
                        + offset
                        * SYNTHETIC_DETECTION_OFFSET
                        * direction
                    )
                    q_rows.append(vector)
                    radii.append(SYNTHETIC_FORMAL_RADIUS)
                    peaks.append(
                        {
                            "peak_id": peak_id,
                            "f": frame_id,
                            "x": float(
                                int.from_bytes(
                                    hashlib.sha256(
                                        (peak_id + ":x").encode("utf-8")
                                    ).digest()[:4],
                                    "big",
                                )
                                % 100000
                            )
                            / 100.0,
                            "y": float(
                                int.from_bytes(
                                    hashlib.sha256(
                                        (peak_id + ":y").encode("utf-8")
                                    ).digest()[:4],
                                    "big",
                                )
                                % 100000
                            )
                            / 100.0,
                            "Shat": signal,
                            "sigma_S": 1.0,
                            "integrated_signal_to_formal_noise": signal,
                            "q": {"validity": VALIDITY},
                            "synthetic_latent_observation_id": observation_id,
                        }
                    )
    return (
        np.asarray(q_rows, dtype=float),
        np.asarray(radii, dtype=float),
        peaks,
        frames,
        latent_q,
        scan_by_observation,
    )


def evaluate_synthetic_subset(
    q: np.ndarray,
    radius: np.ndarray,
    peaks: list[dict[str, Any]],
    frames: dict[str, dict[str, Any]],
    truth_basis: np.ndarray,
    portability: Any,
    d5_rule: dict[str, Any],
) -> tuple[dict[str, Any], np.ndarray, list[dict[str, Any]]]:
    aggregate_q, aggregate_peaks, construction = aggregate_d4(
        q, radius, peaks, frames
    )
    certified_q, certified_peaks = repeat_certified(
        aggregate_q, aggregate_peaks
    )
    result = {
        "construction": construction,
        "all_d45_object_count": len(aggregate_q),
        "repeat_certified_object_count": len(certified_q),
        "frozen_d5": run_frozen_d5(
            certified_q,
            certified_peaks,
            frames,
            portability,
            d5_rule,
            truth_basis,
        ),
    }
    return result, certified_q, certified_peaks


def evaluate_synthetic_case(
    case: dict[str, Any],
    persistence_fraction: float,
    truth_basis: np.ndarray,
    portability: Any,
    d5_rule: dict[str, Any],
) -> dict[str, Any]:
    (
        q,
        radius,
        peaks,
        frames,
        latent_q,
        scan_by_observation,
    ) = frame_local_synthetic_case(case, persistence_fraction)
    full, certified_q, certified_peaks = evaluate_synthetic_subset(
        q, radius, peaks, frames, truth_basis, portability, d5_rule
    )
    member_to_latent = {
        peak["peak_id"]: peak["synthetic_latent_observation_id"]
        for peak in peaks
    }
    centroid_errors: list[float] = []
    mixed_latent_aggregates = 0
    for vector, peak in zip(certified_q, certified_peaks, strict=True):
        latent_ids = {
            member_to_latent[member_id]
            for member_id in peak["d45"]["member_peak_ids"]
        }
        if len(latent_ids) != 1:
            mixed_latent_aggregates += 1
            continue
        latent_id = next(iter(latent_ids))
        centroid_errors.append(
            float(np.linalg.norm(vector - latent_q[latent_id]))
        )

    split_results: dict[str, Any] = {}
    for label, allowed_scans in {
        "SCANS_0_1": {0, 1},
        "SCANS_2_3": {2, 3},
    }.items():
        allowed_latent = {
            observation_id
            for observation_id, scan_index in scan_by_observation.items()
            if scan_index in allowed_scans
        }
        mask = np.asarray(
            [
                peak["synthetic_latent_observation_id"] in allowed_latent
                for peak in peaks
            ],
            dtype=bool,
        )
        local_peaks = [
            peak
            for peak, selected in zip(peaks, mask, strict=True)
            if selected
        ]
        local_frame_ids = {peak["f"] for peak in local_peaks}
        local_frames = {
            frame_id: frames[frame_id] for frame_id in local_frame_ids
        }
        split_results[label], _, _ = evaluate_synthetic_subset(
            q[mask],
            radius[mask],
            local_peaks,
            local_frames,
            truth_basis,
            portability,
            d5_rule,
        )

    split_a = split_results["SCANS_0_1"]["frozen_d5"]
    split_b = split_results["SCANS_2_3"]["frozen_d5"]
    split_equivalence = (
        lattice_equivalence(
            np.asarray(split_a["basis"], dtype=float),
            np.asarray(split_b["basis"], dtype=float),
        )
        if split_a.get("candidate_returned")
        and split_b.get("candidate_returned")
        else {
            "primitive_lattice_equivalent": False,
            "reason": "ONE_OR_BOTH_SPLITS_LACK_A_CANDIDATE",
        }
    )
    required = str(case["required_outcome"])
    if required == "RECOVER_PRIMITIVE_LATTICE":
        success = bool(
            full["frozen_d5"]
            .get("truth_equivalence", {})
            .get("primitive_lattice_equivalent", False)
            and split_a.get("truth_equivalence", {}).get(
                "primitive_lattice_equivalent", False
            )
            and split_b.get("truth_equivalence", {}).get(
                "primitive_lattice_equivalent", False
            )
            and split_equivalence.get(
                "primitive_lattice_equivalent", False
            )
        )
    else:
        success = bool(
            full["frozen_d5"]["execution_status"]
            == "D45_NO_REPEAT_CERTIFIED_OBJECTS"
            or full["frozen_d5"].get(
                "structured_frozen_d5_abstention", False
            )
        )
    return {
        "case_id": case["case_id"],
        "required_outcome": required,
        "latent_observation_count": len(case["observations"]),
        "frame_local_detection_count": len(q),
        "persistence_fraction_rule": persistence_fraction,
        "persistent_latent_observation_count": sum(
            hash_persistent(str(item["observation_id"]), persistence_fraction)
            for item in case["observations"]
        ),
        "origin_labels_used_by_projection": False,
        "full": full,
        "splits": split_results,
        "split_lattice_reproducibility": split_equivalence,
        "certified_centroid_latent_identity": {
            "single_latent_aggregate_count": len(centroid_errors),
            "mixed_latent_aggregate_count": mixed_latent_aggregates,
            "maximum_centroid_error_cycles_per_angstrom": max(
                centroid_errors, default=None
            ),
            "median_centroid_error_cycles_per_angstrom": (
                float(np.median(centroid_errors))
                if centroid_errors
                else None
            ),
        },
        "required_outcome_satisfied": success,
    }


def evaluate_synthetic_suite(
    corpus: dict[str, Any],
    portability: Any,
    d5_rule: dict[str, Any],
) -> dict[str, Any]:
    truth_basis = np.asarray(
        corpus["truth"]["reciprocal_basis_columns_cycles_per_angstrom"],
        dtype=float,
    )
    models: dict[str, Any] = {}
    for label, fraction in {
        "ALL_LATENT_OBJECTS_PERSIST_FOR_THREE_FRAMES": 1.0,
        "ORIGIN_BLIND_HASHED_FIFTY_PERCENT_PERSISTENCE": 0.5,
    }.items():
        cases = []
        for case in corpus["cases"]:
            progress(
                "SYNTHETIC_CASE_START",
                model=label,
                case_id=case["case_id"],
            )
            evaluated = evaluate_synthetic_case(
                case, fraction, truth_basis, portability, d5_rule
            )
            cases.append(evaluated)
            progress(
                "SYNTHETIC_CASE_COMPLETE",
                model=label,
                case_id=case["case_id"],
                required_outcome_satisfied=evaluated[
                    "required_outcome_satisfied"
                ],
                full_execution_status=evaluated["full"]["frozen_d5"][
                    "execution_status"
                ],
            )
        positive = [
            case
            for case in cases
            if case["required_outcome"] == "RECOVER_PRIMITIVE_LATTICE"
        ]
        negative = [
            case
            for case in cases
            if case["required_outcome"] != "RECOVER_PRIMITIVE_LATTICE"
        ]
        models[label] = {
            "persistence_fraction": fraction,
            "cases": cases,
            "positive_recovery_count": sum(
                case["required_outcome_satisfied"] for case in positive
            ),
            "positive_case_count": len(positive),
            "required_abstention_count": sum(
                case["required_outcome_satisfied"] for case in negative
            ),
            "required_abstention_case_count": len(negative),
        }
    return {
        "source_corpus_artifact_id": corpus["artifact_id"],
        "source_corpus_semantic_sha256": corpus[
            "canonical_semantic_sha256"
        ],
        "projection_contract": {
            "input_role": (
                "EACH_EXISTING_D5_BOUNDARY_OBSERVATION_IS_A_LATENT_"
                "THREE_DIMENSIONAL_OBJECT"
            ),
            "frames_per_persistent_object": 3,
            "objects_per_frame_block": SYNTHETIC_OBJECTS_PER_FRAME_BLOCK,
            "formal_radius_cycles_per_angstrom": (
                SYNTHETIC_FORMAL_RADIUS
            ),
            "symmetric_detection_offset_cycles_per_angstrom": (
                SYNTHETIC_DETECTION_OFFSET
            ),
            "offset_direction": "SHA256_DERIVED_FROM_OBSERVATION_ID",
            "origin_or_truth_used_to_choose_persistence": False,
            "centroid_target": (
                "SYMMETRIC_EQUAL_SIGNAL_DETECTIONS_CENTER_ON_THE_ORIGINAL_Q"
            ),
            "scientific_boundary": (
                "D45_NATIVE_ASSOCIATION_CONTROL_NOT_A_DETECTOR_IMAGE_"
                "OR_GEOMETRY_SIMULATION"
            ),
        },
        "models": models,
    }


def public_control_sign_aware_equivalence(
    corpus_paths: list[Path],
    portability: Any,
) -> list[dict[str, Any]]:
    results: list[dict[str, Any]] = []
    for path in corpus_paths:
        corpus = load_gzip_json(path)
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
        old_q, old_peaks, old_record = aggregate_d4(
            q[mask], radius[mask], peaks, frames
        )
        new_q, new_peaks, new_record = (
            aggregate_d4_scan_direction_aware(
                q[mask], radius[mask], peaks, frames
            )
        )
        results.append(
            {
                "input_path_name": path.name,
                "input_sha256": sha256_file(path),
                "frame_count": len(frames),
                "scan_configuration_count": len(
                    {
                        frame["scan_configuration_id"]
                        for frame in frames.values()
                    }
                ),
                "angle_increment_values_degrees": sorted(
                    {
                        float(
                            frame["geometry"][
                                "angle_increment_degrees"
                            ]
                        )
                        for frame in frames.values()
                    }
                ),
                "aggregate_q_array_equal": bool(
                    np.array_equal(old_q, new_q)
                ),
                "aggregate_peak_records_equal": old_peaks == new_peaks,
                "published_membership_sha256": old_record[
                    "output_membership_sha256"
                ],
                "diagnostic_membership_sha256": new_record[
                    "output_membership_sha256"
                ],
                "published_q_sha256": old_record["aggregate_q_sha256"],
                "diagnostic_q_sha256": new_record["aggregate_q_sha256"],
            }
        )
    return results


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--portability-source", type=Path, required=True)
    parser.add_argument("--published-d45-source", type=Path, required=True)
    parser.add_argument("--d4-source", type=Path, required=True)
    parser.add_argument("--d4-rule", type=Path, required=True)
    parser.add_argument("--d5-source", type=Path, required=True)
    parser.add_argument("--d5-rule", type=Path, required=True)
    parser.add_argument("--synthetic-corpus", type=Path, required=True)
    parser.add_argument("--9z6f-d4-corpus", type=Path, required=True)
    parser.add_argument("--9z6f-q-corpus", type=Path, required=True)
    parser.add_argument("--9z6f-historical-result", type=Path, required=True)
    parser.add_argument("--9jq9-array", type=Path, required=True)
    parser.add_argument("--9jq9-record", type=Path, required=True)
    parser.add_argument(
        "--9jq9-historical-diagnostic", type=Path, required=True
    )
    parser.add_argument(
        "--9jq9-deposited-characterization", type=Path, required=True
    )
    parser.add_argument(
        "--public-control-corpus",
        type=Path,
        action="append",
        default=[],
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    frozen_core = args.portability_source.parent / "frozen_core"
    runtime_bindings = verify_runtime(
        args.portability_source,
        args.published_d45_source,
        args.d4_source,
        args.d4_rule,
        args.d5_source,
        args.d5_rule,
    )
    if args.d4_source.resolve() != (
        frozen_core / "d4_raw_peak_pipeline_0_1_0.py"
    ).resolve():
        raise RuntimeError("D4 source is not beside frozen portability")
    if args.d5_source.resolve() != (
        frozen_core / "d5_reciprocal_lattice_pipeline_0_1_0.py"
    ).resolve():
        raise RuntimeError("D5 source is not beside frozen portability")

    portability = load_portability(args.portability_source)
    d5_rule = json.loads(args.d5_rule.read_text(encoding="utf-8"))
    synthetic_corpus = load_gzip_json(args.synthetic_corpus)

    progress("INPUT_LOAD_9Z6F_START")
    q_9z6f, radius_9z6f, peaks_9z6f, frames_9z6f, tiers_9z6f = (
        primary_9z6f(
            args.__dict__["9z6f_d4_corpus"],
            args.__dict__["9z6f_q_corpus"],
        )
    )
    historical_9z6f = json.loads(
        args.__dict__["9z6f_historical_result"].read_text(
            encoding="utf-8"
        )
    )
    progress(
        "INPUT_LOAD_9Z6F_COMPLETE",
        primary_peak_count=len(q_9z6f),
        frame_count=len(frames_9z6f),
    )
    progress("INPUT_LOAD_9JQ9_START")
    q_9jq9, radius_9jq9, peaks_9jq9, frames_9jq9, record_9jq9 = (
        primary_9jq9(
            args.__dict__["9jq9_array"],
            args.__dict__["9jq9_record"],
        )
    )
    historical_9jq9 = json.loads(
        args.__dict__["9jq9_historical_diagnostic"].read_text(
            encoding="utf-8"
        )
    )
    deposited_9jq9 = json.loads(
        args.__dict__["9jq9_deposited_characterization"].read_text(
            encoding="utf-8"
        )
    )
    progress(
        "INPUT_LOAD_9JQ9_COMPLETE",
        primary_peak_count=len(q_9jq9),
        frame_count=len(frames_9jq9),
    )

    progress("SYNTHETIC_SUITE_START")
    synthetic_result = evaluate_synthetic_suite(
        synthetic_corpus, portability, d5_rule
    )
    progress("SYNTHETIC_SUITE_COMPLETE")
    progress("PUBLIC_CONTROL_EQUIVALENCE_START")
    public_equivalence = public_control_sign_aware_equivalence(
        args.public_control_corpus, portability
    )
    progress("PUBLIC_CONTROL_EQUIVALENCE_COMPLETE")
    progress("9Z6F_EVALUATION_START")
    result_9z6f = evaluate_9z6f(
        q_9z6f,
        radius_9z6f,
        peaks_9z6f,
        frames_9z6f,
        historical_9z6f,
        portability,
        d5_rule,
    )
    progress("9Z6F_EVALUATION_COMPLETE")
    progress("9JQ9_EVALUATION_START")
    result_9jq9 = evaluate_9jq9(
        q_9jq9,
        radius_9jq9,
        peaks_9jq9,
        frames_9jq9,
        record_9jq9,
        historical_9jq9,
        deposited_9jq9,
        portability,
        d5_rule,
    )
    progress("9JQ9_EVALUATION_COMPLETE")

    input_paths = {
        "synthetic_corpus": args.synthetic_corpus,
        "9z6f_d4_corpus": args.__dict__["9z6f_d4_corpus"],
        "9z6f_q_corpus": args.__dict__["9z6f_q_corpus"],
        "9z6f_historical_result": args.__dict__[
            "9z6f_historical_result"
        ],
        "9jq9_array": args.__dict__["9jq9_array"],
        "9jq9_record": args.__dict__["9jq9_record"],
        "9jq9_historical_diagnostic": args.__dict__[
            "9jq9_historical_diagnostic"
        ],
        "9jq9_deposited_characterization": args.__dict__[
            "9jq9_deposited_characterization"
        ],
    }
    result = seal(
        {
            "artifact_id": (
                "NFC_CRYST_D45_SYNTHETIC_9Z6F_9JQ9_GENERALIZATION_"
                "AND_SCAN_DIRECTION_DIAGNOSTIC_RESULT_0_1_0"
            ),
            "scientific_scope": (
                "OPEN_EXPLORATORY_D45_AND_UNCHANGED_FROZEN_D5_EVALUATION"
            ),
            "principal_outcome": (
                "D45_SCAN_DIRECTION_DEFECT_ISOLATED_AND_D5_"
                "ABSTENTION_DEFECT_PERSISTS"
            ),
            "evaluation_execution_policy": {
                "per_frozen_d5_fit_budget_seconds": (
                    DEFAULT_D5_FIT_BUDGET_SECONDS
                ),
                "budget_scope": (
                    "EVALUATION_WRAPPER_ONLY_NO_FROZEN_D5_RULE_OR_SOURCE_CHANGE"
                ),
                "budget_exceeded_interpretation": (
                    "OPERATIONAL_NONCOMPLETION_NOT_STRUCTURED_ABSTENTION_"
                    "AND_NOT_A_SCIENTIFIC_PASS_OR_FAIL"
                ),
                "adoption_record": (
                    "ADOPTED_BEFORE_ANY_FROZEN_D5_FIT_IN_THE_CORRECTED_"
                    "STAGED_RUN_TO_ISOLATE_OPERATIONAL_NONCOMPLETION"
                ),
            },
            "evaluation_harness_correction_record": {
                "initial_attempt": (
                    "SINGLE_PROCESS_REMAINED_CPU_ACTIVE_WITHOUT_RESULT_FOR_"
                    "45_MINUTES_AND_WAS_OPERATOR_INTERRUPTED_EXIT_130"
                ),
                "localized_stage": (
                    "9JQ9_NPZ_TO_ROW_INTERFACE_LOADING_BEFORE_SYNTHETIC_"
                    "OR_REAL_FROZEN_D5_EVALUATION"
                ),
                "cause": (
                    "REPEATED_NPZ_MEMBER_DECOMPRESSION_INSIDE_ROW_LOOP"
                ),
                "correction": (
                    "MATERIALIZE_EACH_BOUND_NPZ_ARRAY_EXACTLY_ONCE_BEFORE_"
                    "ROW_CONSTRUCTION"
                ),
                "scientific_values_or_order_changed": False,
                "scientific_method_change": False,
            },
            "runtime_bindings": runtime_bindings,
            "new_diagnostic_source": {
                "path": "d45_scan_direction_diagnostic.py",
                "sha256": sha256_file(
                    Path(__file__).with_name(
                        "d45_scan_direction_diagnostic.py"
                    )
                ),
                "status": (
                    "SEPARATELY_LABELLED_EXPLORATORY_DIAGNOSTIC_NOT_"
                    "PUBLISHED_D45_0_1_0"
                ),
            },
            "input_bindings": {
                key: {
                    "path": path.name,
                    "byte_count": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
                for key, path in input_paths.items()
            },
            "synthetic_d45_native_suite": synthetic_result,
            "public_control_scan_direction_equivalence": public_equivalence,
            "9z6f": {
                **result_9z6f,
                "source_tier_counts": tiers_9z6f,
            },
            "9jq9": result_9jq9,
            "environment": {
                "python": platform.python_version(),
                "implementation": platform.python_implementation(),
                "numpy": np.__version__,
                "scipy": scipy.__version__,
                "platform": platform.platform(),
            },
            "scientific_changes": {
                "frozen_d5_source_or_rule_changed": False,
                "published_d45_source_changed": False,
                "scan_direction_diagnostic_is_a_new_exploratory_variant": True,
                "d6c1_run": False,
                "d6c2_run": False,
                "d7_run": False,
            },
            "claim_boundaries": [
                "NOT_CONFIRMATORY",
                "NOT_INDEPENDENT_VALIDATION",
                "9Z6F_REMAINS_A_DEVELOPMENT_SPECIMEN",
                "9JQ9_REMAINS_PERMANENTLY_EXPLORATORY",
                "SYNTHETIC_FRAME_LOCAL_PROJECTION_IS_NOT_A_DETECTOR_IMAGE_SIMULATION",
                "NO_D6_OR_D7_TRANSFER_RESULT",
                "NO_NFC_TOE_CLAIM",
            ],
        }
    )
    write_json(args.output, result)
    print(
        json.dumps(
            {
                "output": str(args.output),
                "semantic_sha256": result["semantic_sha256"],
                "synthetic_models": {
                    key: {
                        "positive_recovery": (
                            value["positive_recovery_count"],
                            value["positive_case_count"],
                        ),
                        "required_abstention": (
                            value["required_abstention_count"],
                            value["required_abstention_case_count"],
                        ),
                    }
                    for key, value in synthetic_result["models"].items()
                },
                "9z6f_published_d45_repeat_count": result_9z6f[
                    "published_d45_0_1_0"
                ]["repeat_certified_object_count"],
                "9z6f_sign_aware_repeat_count": result_9z6f[
                    "scan_local_increment_sign_aware_diagnostic"
                ]["repeat_certified_object_count"],
                "9jq9_published_d45_repeat_count": result_9jq9[
                    "published_d45_0_1_0"
                ]["repeat_certified_object_count"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
