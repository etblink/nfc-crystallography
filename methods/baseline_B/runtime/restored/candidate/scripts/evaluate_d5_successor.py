#!/usr/bin/env python3
"""Evaluate a truth-free D5 recovery-or-abstention successor candidate.

The candidate generator and robust integer refinement are the unchanged
frozen D5 0.3.0 operations.  The successor contribution is an evidence gate:

* independently construct and fit full, first-half, and second-half feeds;
* require primitive-lattice consensus among all three candidates;
* score each fitted basis on its own feed and across the held-out half;
* subtract an equal-density fractional-origin null;
* return RECOVERED, AMBIGUOUS, or INSUFFICIENT without consulting truth.

Conventional or synthetic truth is loaded only after the decision and used
solely to score whether a recovery decision was physically correct.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import itertools
import json
import math
import platform
import sys
from pathlib import Path
from typing import Any

import numpy as np
import scipy

from d45_aggregation import canonical_bytes, semantic_sha256


VALIDITY = "CERTIFIED_POSITIVE_SEPARATION_INTERIOR"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def load_gzip_json(path: Path) -> Any:
    with gzip.open(path, "rb") as handle:
        return json.loads(handle.read())


def load_portability(path: Path) -> Any:
    source = path.resolve()
    sys.path.insert(0, str(source.parent))
    spec = importlib.util.spec_from_file_location(
        "nfc_portability_d5_successor", source
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load frozen portability source")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def unpack_feed(
    feed: dict[str, Any],
) -> tuple[
    np.ndarray,
    list[dict[str, Any]],
    dict[str, dict[str, Any]],
]:
    q = np.asarray(feed["q_cycles_per_angstrom"], dtype=float)
    objects = feed["objects"]
    if q.shape != (len(objects), 3):
        raise ValueError("fixed-feed q/object mismatch")
    peaks = [
        {
            "peak_id": item["object_id"],
            "f": item["representative_frame_id"],
            "integrated_signal_to_formal_noise": item["signal_to_noise"],
            "q": {"validity": VALIDITY},
        }
        for item in objects
    ]
    frames = feed["frames"]
    if any(peak["f"] not in frames for peak in peaks):
        raise ValueError("fixed-feed representative frame missing")
    observed_q_sha256 = hashlib.sha256(
        np.asarray(q, dtype="<f8", order="C").tobytes(order="C")
    ).hexdigest()
    if observed_q_sha256 != feed["q_sha256"]:
        raise ValueError("fixed-feed q identity mismatch")
    return q, peaks, frames


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
    }


def fit_candidate(
    feed: dict[str, Any],
    portability: Any,
    frozen_rule: dict[str, Any],
    minimum_objects: int,
) -> dict[str, Any]:
    q, peaks, frames = unpack_feed(feed)
    if len(q) < minimum_objects:
        return {
            "status": "INSUFFICIENT_OBJECTS",
            "object_count": len(q),
            "minimum_required": minimum_objects,
            "candidate_returned": False,
        }
    try:
        seed, seed_record = portability.FROZEN_D5.seed_lattice(
            q, peaks, frames, frozen_rule
        )
        basis, refinement = portability.FROZEN_D5.refine_lattice(
            q, seed, frozen_rule
        )
    except Exception as error:
        return {
            "status": "CANDIDATE_KERNEL_FAILURE",
            "object_count": len(q),
            "candidate_returned": False,
            "exception_type": type(error).__name__,
            "exception_message": str(error),
        }
    if not np.all(np.isfinite(basis)) or abs(np.linalg.det(basis)) <= 0:
        return {
            "status": "NONFINITE_OR_SINGULAR_CANDIDATE",
            "object_count": len(q),
            "candidate_returned": False,
        }
    return {
        "status": "CANDIDATE_RETURNED",
        "object_count": len(q),
        "candidate_returned": True,
        "basis": basis.tolist(),
        "absolute_reciprocal_determinant": float(abs(np.linalg.det(basis))),
        "seed": compact_seed(seed_record),
        "refinement": refinement,
    }


def primitive_equivalence(
    left: np.ndarray,
    right: np.ndarray,
    maximum_deviation: float,
    required_abs_determinant: int,
) -> dict[str, Any]:
    transform = np.linalg.inv(left) @ right
    nearest = np.rint(transform).astype(np.int64)
    deviation = float(np.max(np.abs(transform - nearest)))
    determinant = int(round(float(np.linalg.det(nearest))))
    return {
        "nearest_integer_basis_transform": nearest.tolist(),
        "maximum_transform_integer_deviation": deviation,
        "rounded_integer_transform_determinant": determinant,
        "primitive_lattice_equivalent": bool(
            deviation <= maximum_deviation
            and abs(determinant) == required_abs_determinant
        ),
    }


def phase_shifts(rule: dict[str, Any]) -> np.ndarray:
    values = [
        float(value)
        for value in rule["same_density_phase_specificity"][
            "fractional_origin_grid"
        ]
    ]
    shifts = np.asarray(
        [
            item
            for item in itertools.product(values, repeat=3)
            if item != (0.0, 0.0, 0.0)
        ],
        dtype=float,
    )
    expected = int(
        rule["same_density_phase_specificity"]["null_origin_count"]
    )
    if shifts.shape != (expected, 3):
        raise ValueError("phase-null grid cardinality mismatch")
    return shifts


def phase_specificity(
    q: np.ndarray,
    basis: np.ndarray,
    rule: dict[str, Any],
) -> dict[str, Any]:
    cutoff = float(
        rule["same_density_phase_specificity"][
            "residual_cutoff_cycles_per_angstrom"
        ]
    )
    shifts = phase_shifts(rule)
    fractional = q @ np.linalg.inv(basis).T

    def support(shift: np.ndarray) -> float:
        delta = fractional - shift[None, :]
        delta -= np.rint(delta)
        residual = np.linalg.norm(delta @ basis.T, axis=1)
        return float(np.mean(residual <= cutoff))

    zero = support(np.zeros(3, dtype=float))
    null = np.asarray([support(shift) for shift in shifts], dtype=float)
    median = float(np.median(null))
    excess = zero - median
    empirical_p = float(
        (1 + int(np.sum(null >= zero))) / (1 + len(null))
    )
    minimum_excess = float(
        rule["same_density_phase_specificity"][
            "minimum_zero_minus_null_median_support_fraction"
        ]
    )
    maximum_p = float(
        rule["same_density_phase_specificity"][
            "maximum_empirical_rank_p"
        ]
    )
    return {
        "object_count": len(q),
        "residual_cutoff_cycles_per_angstrom": cutoff,
        "zero_origin_support_fraction": zero,
        "null_origin_count": len(null),
        "null_support_fraction": {
            "minimum": float(np.min(null)),
            "median": median,
            "mean": float(np.mean(null)),
            "maximum": float(np.max(null)),
        },
        "zero_minus_null_median_support_fraction": excess,
        "empirical_rank_p": empirical_p,
        "minimum_required_excess": minimum_excess,
        "maximum_allowed_empirical_rank_p": maximum_p,
        "passes": bool(
            excess >= minimum_excess and empirical_p <= maximum_p
        ),
        "density_control": (
            "IDENTICAL_BASIS_AND_DETERMINANT_FOR_ZERO_AND_ALL_NULL_ORIGINS"
        ),
    }


def van_der_corput(index: int, base: int) -> float:
    value = 0.0
    denominator = 1.0
    while index:
        index, remainder = divmod(index, base)
        denominator *= base
        value += remainder / denominator
    return value


def halton_phase_sensitivity(
    q: np.ndarray,
    basis: np.ndarray,
    cutoff: float = 0.003,
    count: int = 127,
) -> dict[str, Any]:
    """Independent continuous-shift sensitivity; never a decision input."""

    shifts = np.asarray(
        [
            [
                van_der_corput(index, 2),
                van_der_corput(index, 3),
                van_der_corput(index, 5),
            ]
            for index in range(1, count + 1)
        ],
        dtype=float,
    )
    fractional = q @ np.linalg.inv(basis).T

    def support(shift: np.ndarray) -> float:
        delta = fractional - shift[None, :]
        delta -= np.rint(delta)
        residual = np.linalg.norm(delta @ basis.T, axis=1)
        return float(np.mean(residual <= cutoff))

    zero = support(np.zeros(3, dtype=float))
    null = np.asarray([support(shift) for shift in shifts], dtype=float)
    return {
        "status": "POST_DECISION_SENSITIVITY_NOT_A_DECISION_INPUT",
        "sequence": "HALTON_BASES_2_3_5_INDICES_1_THROUGH_127",
        "null_origin_count": len(null),
        "zero_origin_support_fraction": zero,
        "null_median_support_fraction": float(np.median(null)),
        "null_maximum_support_fraction": float(np.max(null)),
        "zero_minus_null_median_support_fraction": float(
            zero - np.median(null)
        ),
        "empirical_rank_p": float(
            (1 + int(np.sum(null >= zero))) / (1 + len(null))
        ),
    }


def upper_hnf_coarsenings(maximum_index: int) -> list[np.ndarray]:
    """Enumerate canonical finite-index sublattices through upper HNF."""

    matrices: list[np.ndarray] = []
    for index in range(2, maximum_index + 1):
        for first in range(1, index + 1):
            if index % first:
                continue
            remainder = index // first
            for second in range(1, remainder + 1):
                if remainder % second:
                    continue
                third = remainder // second
                for upper_01 in range(second):
                    for upper_02 in range(third):
                        for upper_12 in range(third):
                            matrices.append(
                                np.asarray(
                                    [
                                        [first, upper_01, upper_02],
                                        [0, second, upper_12],
                                        [0, 0, third],
                                    ],
                                    dtype=np.int64,
                                )
                            )
    return matrices


def zero_origin_support(
    q: np.ndarray,
    basis: np.ndarray,
    cutoff: float,
) -> float:
    indices = np.rint(q @ np.linalg.inv(basis).T)
    residual = np.linalg.norm(q - indices @ basis.T, axis=1)
    return float(np.mean(residual <= cutoff))


def primitive_complexity_audit(
    q: np.ndarray,
    basis: np.ndarray,
    selected_phase: dict[str, Any],
    rule: dict[str, Any],
) -> dict[str, Any]:
    contract = rule["primitive_complexity_audit"]
    maximum_index = int(contract["maximum_countermodel_index"])
    minimum_relative = float(
        contract[
            "minimum_countermodel_support_relative_to_selected_basis"
        ]
    )
    cutoff = float(
        rule["same_density_phase_specificity"][
            "residual_cutoff_cycles_per_angstrom"
        ]
    )
    selected_support = float(
        selected_phase["zero_origin_support_fraction"]
    )
    supported: list[dict[str, Any]] = []
    best_relative = 0.0
    tested = upper_hnf_coarsenings(maximum_index)
    for transform in tested:
        countermodel = basis @ transform
        support = zero_origin_support(q, countermodel, cutoff)
        relative = (
            support / selected_support if selected_support > 0.0 else 0.0
        )
        best_relative = max(best_relative, relative)
        if relative < minimum_relative:
            continue
        phase = phase_specificity(q, countermodel, rule)
        if not phase["passes"]:
            continue
        supported.append(
            {
                "integer_index": int(round(float(np.linalg.det(transform)))),
                "upper_hnf_transform": transform.tolist(),
                "absolute_reciprocal_determinant": float(
                    abs(np.linalg.det(countermodel))
                ),
                "zero_origin_support_fraction": support,
                "support_relative_to_selected_basis": relative,
                "phase_specificity": phase,
            }
        )
    supported.sort(
        key=lambda item: (
            -item["support_relative_to_selected_basis"],
            item["integer_index"],
            item["upper_hnf_transform"],
        )
    )
    allowed = int(
        contract["allowed_supported_coarser_countermodel_count"]
    )
    return {
        "tested_upper_hnf_countermodel_count": len(tested),
        "maximum_integer_index_tested": maximum_index,
        "minimum_required_relative_support": minimum_relative,
        "best_countermodel_relative_support": best_relative,
        "supported_coarser_countermodel_count": len(supported),
        "supported_coarser_countermodels": supported,
        "passes": len(supported) <= allowed,
        "selected_basis_absolute_reciprocal_determinant": float(
            abs(np.linalg.det(basis))
        ),
        "scientific_role": contract["scientific_role"],
    }


def basis_or_none(fit: dict[str, Any]) -> np.ndarray | None:
    if not fit.get("candidate_returned", False):
        return None
    return np.asarray(fit["basis"], dtype=float)


def evaluate_case(
    case: dict[str, Any],
    portability: Any,
    frozen_rule: dict[str, Any],
    successor_rule: dict[str, Any],
) -> dict[str, Any]:
    minimum_objects = int(
        successor_rule["input_representation"][
            "minimum_repeat_certified_objects_per_full_or_split_feed"
        ]
    )
    fits = {
        label: fit_candidate(
            case["feeds"][label],
            portability,
            frozen_rule,
            minimum_objects,
        )
        for label in ("FULL", "HALF_A", "HALF_B")
    }
    q = {
        label: unpack_feed(case["feeds"][label])[0]
        for label in ("FULL", "HALF_A", "HALF_B")
    }
    bases = {label: basis_or_none(fit) for label, fit in fits.items()}
    candidate_failure = any(basis is None for basis in bases.values())
    equivalence: dict[str, Any] = {}
    consensus = False
    own_phase: dict[str, Any] = {}
    cross_phase: dict[str, Any] = {}
    complexity_audit: dict[str, Any] = {
        "status": "NOT_RUN_CANDIDATE_UNAVAILABLE"
    }
    halton_sensitivity: dict[str, Any] = {
        "status": "NOT_RUN_CANDIDATE_UNAVAILABLE"
    }

    if not candidate_failure:
        assert all(basis is not None for basis in bases.values())
        max_deviation = float(
            successor_rule["primitive_lattice_consensus"][
                "maximum_integer_transform_deviation"
            ]
        )
        abs_determinant = int(
            successor_rule["primitive_lattice_consensus"][
                "required_absolute_rounded_transform_determinant"
            ]
        )
        equivalence = {
            "FULL_TO_HALF_A": primitive_equivalence(
                bases["FULL"],
                bases["HALF_A"],
                max_deviation,
                abs_determinant,
            ),
            "FULL_TO_HALF_B": primitive_equivalence(
                bases["FULL"],
                bases["HALF_B"],
                max_deviation,
                abs_determinant,
            ),
            "HALF_A_TO_HALF_B": primitive_equivalence(
                bases["HALF_A"],
                bases["HALF_B"],
                max_deviation,
                abs_determinant,
            ),
        }
        consensus = all(
            item["primitive_lattice_equivalent"]
            for item in equivalence.values()
        )
        own_phase = {
            label: phase_specificity(
                q[label], bases[label], successor_rule
            )
            for label in ("FULL", "HALF_A", "HALF_B")
        }
        cross_phase = {
            "HALF_A_BASIS_ON_HALF_B_HELD_OUT": phase_specificity(
                q["HALF_B"], bases["HALF_A"], successor_rule
            ),
            "HALF_B_BASIS_ON_HALF_A_HELD_OUT": phase_specificity(
                q["HALF_A"], bases["HALF_B"], successor_rule
            ),
        }
        complexity_audit = primitive_complexity_audit(
            q["FULL"],
            bases["FULL"],
            own_phase["FULL"],
            successor_rule,
        )
        halton_sensitivity = {
            "own_feed": {
                label: halton_phase_sensitivity(
                    q[label], bases[label]
                )
                for label in ("FULL", "HALF_A", "HALF_B")
            },
            "held_out": {
                "HALF_A_BASIS_ON_HALF_B_HELD_OUT": (
                    halton_phase_sensitivity(
                        q["HALF_B"], bases["HALF_A"]
                    )
                ),
                "HALF_B_BASIS_ON_HALF_A_HELD_OUT": (
                    halton_phase_sensitivity(
                        q["HALF_A"], bases["HALF_B"]
                    )
                ),
            },
            "status": "POST_DECISION_SENSITIVITY_NOT_A_DECISION_INPUT",
        }

    required_phase_passes = (
        [
            own_phase.get(label, {}).get("passes", False)
            for label in ("FULL", "HALF_A", "HALF_B")
        ]
        + [
            cross_phase.get(label, {}).get("passes", False)
            for label in (
                "HALF_A_BASIS_ON_HALF_B_HELD_OUT",
                "HALF_B_BASIS_ON_HALF_A_HELD_OUT",
            )
        ]
    )
    if candidate_failure:
        decision = "INSUFFICIENT_SIGNAL"
        reason = "ONE_OR_MORE_FULL_OR_SPLIT_CANDIDATES_UNAVAILABLE"
    elif not consensus:
        if own_phase["FULL"]["passes"]:
            decision = "AMBIGUOUS_LATTICE"
            reason = (
                "PHASE_SPECIFIC_POOLED_CANDIDATE_CONFLICTS_WITH_ONE_OR_"
                "MORE_INDEPENDENT_SPLIT_CANDIDATES"
            )
        else:
            decision = "INSUFFICIENT_SIGNAL"
            reason = (
                "SPLIT_CANDIDATES_CONFLICT_AND_POOLED_PHASE_SPECIFICITY_FAILS"
            )
    elif not all(required_phase_passes):
        decision = "INSUFFICIENT_SIGNAL"
        reason = "ONE_OR_MORE_OWN_OR_HELD_OUT_PHASE_SPECIFICITY_GATES_FAIL"
    elif not complexity_audit["passes"]:
        decision = "AMBIGUOUS_LATTICE"
        reason = (
            "ONE_OR_MORE_LOWER_COMPLEXITY_COARSER_RECIPROCAL_LATTICES_"
            "RETAIN_COMPARABLE_PHASE_SPECIFIC_SUPPORT"
        )
    else:
        decision = "LATTICE_RECOVERED"
        reason = (
            "FULL_AND_SPLIT_PRIMITIVE_CONSENSUS_WITH_OWN_AND_HELD_OUT_"
            "SAME_DENSITY_PHASE_SPECIFICITY"
        )

    ablations = {
        "FULL_PHASE_ONLY": (
            "LATTICE_RECOVERED"
            if own_phase.get("FULL", {}).get("passes", False)
            else "INSUFFICIENT_SIGNAL"
        ),
        "SPLIT_CONSENSUS_ONLY": (
            "LATTICE_RECOVERED"
            if not candidate_failure and consensus
            else "ABSTAIN"
        ),
        "COMBINED_SUCCESSOR": decision,
    }

    truth_basis = case["truth"].get("basis")
    truth_scoring: dict[str, Any] = {
        "truth_available": truth_basis is not None,
        "performed_after_decision": True,
    }
    if truth_basis is not None:
        truth = np.asarray(truth_basis, dtype=float)
        scoring = {
            label: (
                primitive_equivalence(
                    truth,
                    basis,
                    float(
                        successor_rule["primitive_lattice_consensus"][
                            "maximum_integer_transform_deviation"
                        ]
                    ),
                    int(
                        successor_rule["primitive_lattice_consensus"][
                            "required_absolute_rounded_transform_determinant"
                        ]
                    ),
                )
                if basis is not None
                else {
                    "primitive_lattice_equivalent": False,
                    "reason": "CANDIDATE_UNAVAILABLE",
                }
            )
            for label, basis in bases.items()
        }
        truth_scoring["candidate_truth_equivalence"] = scoring
        truth_scoring["recovery_physically_correct"] = bool(
            decision == "LATTICE_RECOVERED"
            and all(
                item["primitive_lattice_equivalent"]
                for item in scoring.values()
            )
        )

    required = str(case["required_outcome"])
    if required == "RECOVER_PRIMITIVE_LATTICE":
        required_satisfied = bool(
            truth_scoring.get("recovery_physically_correct", False)
        )
    else:
        required_satisfied = decision in {
            "AMBIGUOUS_LATTICE",
            "INSUFFICIENT_SIGNAL",
        }
    return {
        "case_id": case["case_id"],
        "case_kind": case["case_kind"],
        "role": case.get("role"),
        "persistence_fraction": case.get("persistence_fraction"),
        "required_outcome": required,
        "decision": decision,
        "decision_reason": reason,
        "candidate_fits": fits,
        "primitive_lattice_consensus": {
            "passes": consensus,
            "pairwise": equivalence,
        },
        "own_feed_phase_specificity": own_phase,
        "held_out_phase_specificity": cross_phase,
        "primitive_complexity_audit": complexity_audit,
        "secondary_continuous_phase_null_sensitivity": (
            halton_sensitivity
        ),
        "decision_ablations": ablations,
        "truth_scoring": truth_scoring,
        "required_outcome_satisfied": required_satisfied,
    }


def aggregate_census(cases: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "case_count": len(cases),
        "decision_counts": {
            decision: sum(
                case["decision"] == decision for case in cases
            )
            for decision in (
                "LATTICE_RECOVERED",
                "AMBIGUOUS_LATTICE",
                "INSUFFICIENT_SIGNAL",
            )
        },
        "required_outcome_satisfied_count": sum(
            case["required_outcome_satisfied"] for case in cases
        ),
        "incorrect_recovery_count": sum(
            case["decision"] == "LATTICE_RECOVERED"
            and not case["required_outcome_satisfied"]
            for case in cases
        ),
    }


def sensitivity_decision(
    case: dict[str, Any],
    minimum_excess: float,
    maximum_p: float,
    maximum_deviation: float,
) -> str:
    if any(
        not fit.get("candidate_returned", False)
        for fit in case["candidate_fits"].values()
    ):
        return "INSUFFICIENT_SIGNAL"
    pairwise = case["primitive_lattice_consensus"]["pairwise"]
    consensus = all(
        abs(int(item["rounded_integer_transform_determinant"])) == 1
        and float(item["maximum_transform_integer_deviation"])
        <= maximum_deviation
        for item in pairwise.values()
    )

    def phase_passes(item: dict[str, Any]) -> bool:
        return bool(
            float(
                item["zero_minus_null_median_support_fraction"]
            )
            >= minimum_excess
            and float(item["empirical_rank_p"]) <= maximum_p
        )

    own = case["own_feed_phase_specificity"]
    held_out = case["held_out_phase_specificity"]
    if not consensus:
        return (
            "AMBIGUOUS_LATTICE"
            if phase_passes(own["FULL"])
            else "INSUFFICIENT_SIGNAL"
        )
    if all(
        phase_passes(item)
        for item in list(own.values()) + list(held_out.values())
    ):
        return (
            "LATTICE_RECOVERED"
            if case["primitive_complexity_audit"].get("passes", False)
            else "AMBIGUOUS_LATTICE"
        )
    return "INSUFFICIENT_SIGNAL"


def sensitivity_grid(
    real_cases: list[dict[str, Any]],
    synthetic_results: dict[str, list[dict[str, Any]]],
) -> dict[str, Any]:
    records = []
    for minimum_excess, maximum_p, maximum_deviation in itertools.product(
        (0.05, 0.10, 0.20, 0.30),
        (0.015625, 0.03125, 0.0625),
        (0.02, 0.03, 0.04),
    ):
        decisions: dict[str, str] = {}
        cases_by_key: dict[str, dict[str, Any]] = {}
        for case in real_cases:
            key = f"REAL::{case['case_id']}"
            decisions[key] = sensitivity_decision(
                case, minimum_excess, maximum_p, maximum_deviation
            )
            cases_by_key[key] = case
        for model, cases in synthetic_results.items():
            for case in cases:
                key = f"SYNTHETIC::{model}::{case['case_id']}"
                decisions[key] = sensitivity_decision(
                    case,
                    minimum_excess,
                    maximum_p,
                    maximum_deviation,
                )
                cases_by_key[key] = case
        binary_correct = {
            key: (
                decision == "LATTICE_RECOVERED"
                if case["required_outcome"]
                == "RECOVER_PRIMITIVE_LATTICE"
                else decision != "LATTICE_RECOVERED"
            )
            for key, decision in decisions.items()
            for case in [cases_by_key[key]]
        }
        records.append(
            {
                "minimum_phase_excess": minimum_excess,
                "maximum_empirical_rank_p": maximum_p,
                "maximum_consensus_integer_transform_deviation": (
                    maximum_deviation
                ),
                "recovery_or_abstention_correct_count": sum(
                    binary_correct.values()
                ),
                "case_count": len(binary_correct),
                "all_recovery_or_abstention_classes_correct": all(
                    binary_correct.values()
                ),
                "incorrect_case_keys": sorted(
                    key
                    for key, correct in binary_correct.items()
                    if not correct
                ),
            }
        )
    return {
        "grid": records,
        "configuration_count": len(records),
        "all_case_classes_correct_configuration_count": sum(
            record["all_recovery_or_abstention_classes_correct"]
            for record in records
        ),
        "interpretation": (
            "POST_HOC_ROBUSTNESS_GRID_NOT_ADDITIONAL_MODEL_SELECTION"
        ),
    }


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(canonical_bytes(value))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--feeds", type=Path, required=True)
    parser.add_argument("--portability-source", type=Path, required=True)
    parser.add_argument("--frozen-d5-rule", type=Path, required=True)
    parser.add_argument("--successor-rule", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    feeds = load_gzip_json(args.feeds)
    frozen_rule = json.loads(
        args.frozen_d5_rule.read_text(encoding="utf-8")
    )
    successor_rule = json.loads(
        args.successor_rule.read_text(encoding="utf-8")
    )
    portability = load_portability(args.portability_source)

    real_results = []
    for case in feeds["real_cases"]:
        print(f"EVALUATE REAL {case['case_id']}", flush=True)
        real_results.append(
            evaluate_case(case, portability, frozen_rule, successor_rule)
        )
    synthetic_results: dict[str, list[dict[str, Any]]] = {}
    for model, cases in feeds["synthetic_models"].items():
        synthetic_results[model] = []
        for case in cases:
            print(f"EVALUATE SYNTHETIC {model} {case['case_id']}", flush=True)
            synthetic_results[model].append(
                evaluate_case(
                    case, portability, frozen_rule, successor_rule
                )
            )

    result: dict[str, Any] = {
        "artifact_id": (
            "NFC_CRYST_D45_FIXED_INPUT_D5_CONSENSUS_PHASE_SUCCESSOR_"
            "EVALUATION_0_1_0"
        ),
        "scientific_scope": (
            "OPEN_EXPLORATORY_D45_AND_D5_SUCCESSOR_DEVELOPMENT"
        ),
        "principal_outcome": "PENDING_CENSUS",
        "source_bindings": {
            "fixed_feeds": {
                "sha256": sha256_file(args.feeds),
                "semantic_sha256": feeds["semantic_sha256"],
            },
            "portability_source_sha256": sha256_file(
                args.portability_source
            ),
            "frozen_d5_source_sha256": sha256_file(
                args.portability_source.parent
                / "frozen_core"
                / "d5_reciprocal_lattice_pipeline_0_1_0.py"
            ),
            "frozen_d5_rule_sha256": sha256_file(args.frozen_d5_rule),
            "successor_rule_sha256": sha256_file(args.successor_rule),
        },
        "real_cases": real_results,
        "synthetic_models": synthetic_results,
        "census": {
            "real": aggregate_census(real_results),
            "synthetic": {
                model: aggregate_census(cases)
                for model, cases in synthetic_results.items()
            },
        },
        "decision_threshold_sensitivity": sensitivity_grid(
            real_results, synthetic_results
        ),
        "truth_firewall": successor_rule["truth_firewall"],
        "environment": {
            "python": platform.python_version(),
            "implementation": platform.python_implementation(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "platform": platform.platform(),
        },
        "claim_boundaries": [
            "NOT_CONFIRMATORY",
            "NOT_INDEPENDENT_VALIDATION",
            "ALL_REAL_CORPORA_ARE_DEVELOPMENT_OR_DIAGNOSTIC_CORPORA",
            "FROZEN_METHOD_0_3_0_REMAINS_IMMUTABLE",
            "NO_D6C1_D6C2_OR_D7_TRANSFER_RESULT",
            "NO_NFC_TOE_CLAIM",
        ],
    }
    positive_real = [
        case
        for case in real_results
        if case["required_outcome"] == "RECOVER_PRIMITIVE_LATTICE"
    ]
    adverse_real = [
        case
        for case in real_results
        if case["required_outcome"] != "RECOVER_PRIMITIVE_LATTICE"
    ]
    synthetic_all = [
        case
        for cases in synthetic_results.values()
        for case in cases
    ]
    positive_synthetic = [
        case
        for case in synthetic_all
        if case["required_outcome"] == "RECOVER_PRIMITIVE_LATTICE"
    ]
    adverse_synthetic = [
        case
        for case in synthetic_all
        if case["required_outcome"] != "RECOVER_PRIMITIVE_LATTICE"
    ]
    if (
        all(case["required_outcome_satisfied"] for case in positive_real)
        and all(
            case["required_outcome_satisfied"]
            for case in positive_synthetic
        )
        and all(
            case["required_outcome_satisfied"]
            for case in adverse_synthetic
        )
        and all(
            case["required_outcome_satisfied"] for case in adverse_real
        )
    ):
        result["principal_outcome"] = (
            "D45_FIXED_INPUT_D5_SUCCESSOR_PASSES_CURRENT_EXPLORATORY_"
            "RECOVERY_AND_STRUCTURED_ABSTENTION_CONTROLS"
        )
    else:
        result["principal_outcome"] = (
            "D45_FIXED_INPUT_D5_SUCCESSOR_DOES_NOT_PASS_CURRENT_"
            "EXPLORATORY_CONTROL_SET"
        )
    result["summary"] = {
        "real_positive_recovery": [
            case["case_id"]
            for case in positive_real
            if case["required_outcome_satisfied"]
        ],
        "real_adverse_structured_abstention": [
            case["case_id"]
            for case in adverse_real
            if case["required_outcome_satisfied"]
        ],
        "synthetic_positive_recovery_count": sum(
            case["required_outcome_satisfied"]
            for case in positive_synthetic
        ),
        "synthetic_positive_case_count": len(positive_synthetic),
        "synthetic_adverse_abstention_count": sum(
            case["required_outcome_satisfied"]
            for case in adverse_synthetic
        ),
        "synthetic_adverse_case_count": len(adverse_synthetic),
    }
    body = dict(result)
    result["semantic_sha256"] = semantic_sha256(body)
    write_json(args.output, result)
    print(
        json.dumps(
            {
                "output": str(args.output),
                "semantic_sha256": result["semantic_sha256"],
                "principal_outcome": result["principal_outcome"],
                "summary": result["summary"],
                "real_decisions": {
                    case["case_id"]: case["decision"]
                    for case in real_results
                },
                "synthetic_decisions": {
                    model: {
                        case["case_id"]: case["decision"]
                        for case in cases
                    }
                    for model, cases in synthetic_results.items()
                },
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
