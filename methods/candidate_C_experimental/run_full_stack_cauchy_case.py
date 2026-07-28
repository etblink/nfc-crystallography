#!/usr/bin/env python3
"""Run one fixed full-stack Cauchy-integration case."""

from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.util
import itertools
import json
import platform
import time
from pathlib import Path
from typing import Any

import numpy as np
import scipy


ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "full_stack_cauchy_work"
PROTOCOL_PATH = WORK / "FULL_STACK_CAUCHY_INTEGRATION_PROTOCOL_0_1_0.json"
PROTOCOL_SHA256 = "539d0c190434493b59e4b427213422a02922823e090eea9468a599fceb355da1"
FLOORS = (0.012, 0.009, 0.006)


def canonical_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")


def semantic_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(8 << 20):
            digest.update(block)
    return digest.hexdigest()


def load(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def load_evaluator():
    path = ROOT / "rebuild" / "evaluate.py"
    specification = importlib.util.spec_from_file_location(
        "nfc_full_stack_archived_evaluator", path
    )
    module = importlib.util.module_from_spec(specification)
    assert specification.loader is not None
    specification.loader.exec_module(module)
    return module


def validate_sources(protocol: dict[str, Any]) -> None:
    if sha256_file(PROTOCOL_PATH) != PROTOCOL_SHA256:
        raise RuntimeError("protocol hash mismatch")
    for binding in protocol["source_bindings"].values():
        path = ROOT / binding["path"]
        if sha256_file(path) != binding["sha256"]:
            raise RuntimeError(f"source hash mismatch: {binding['path']}")


def stable_float(value: float) -> float:
    return float(format(float(value), ".15g"))


def stable(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: stable(item) for key, item in value.items()}
    if isinstance(value, list):
        return [stable(item) for item in value]
    if isinstance(value, tuple):
        return [stable(item) for item in value]
    if isinstance(value, (np.floating, float)):
        return stable_float(float(value))
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def origin_rows() -> np.ndarray:
    values = (0.0, 0.25, 0.5, 0.75)
    return np.asarray(
        [(0.0, 0.0, 0.0)]
        + [
            origin
            for origin in itertools.product(values, repeat=3)
            if origin != (0.0, 0.0, 0.0)
        ],
        dtype=float,
    )


ORIGINS = origin_rows()


def cauchy_phase_specificity(
    q: np.ndarray,
    basis: np.ndarray,
    rule: dict[str, Any],
) -> dict[str, Any]:
    del rule
    scale = 0.003
    fractional = q @ np.linalg.inv(basis).T
    scores = np.empty(len(ORIGINS), dtype=float)
    for index, origin in enumerate(ORIGINS):
        delta = fractional - origin[None, :]
        delta -= np.rint(delta)
        residual = np.linalg.norm(delta @ basis.T, axis=1)
        ratio = residual / scale
        scores[index] = float(np.mean(1.0 / (1.0 + ratio * ratio)))

    zero = float(scores[0])
    null = scores[1:]
    best_null_relative = max(
        range(len(null)),
        key=lambda index: (
            float(null[index]),
            tuple(-float(value) for value in ORIGINS[index + 1]),
        ),
    )
    best_null_index = best_null_relative + 1
    best_null = float(scores[best_null_index])
    tolerance = max(1e-15, abs(zero) * 1e-12)
    rank = 1 + int(np.sum(null >= zero - tolerance))
    margin = zero - best_null
    return stable(
        {
            "object_count": len(q),
            "score_name": "CAUCHY_AFFINITY",
            "kernel_scale_cycles_per_angstrom": scale,
            "hard_residual_cutoff": None,
            "zero_origin_score": zero,
            "best_null_origin": ORIGINS[best_null_index].tolist(),
            "best_null_score": best_null,
            "null_origin_count": len(null),
            "null_origin_score": {
                "minimum": float(np.min(null)),
                "median": float(np.median(null)),
                "mean": float(np.mean(null)),
                "maximum": float(np.max(null)),
            },
            "zero_minus_best_null_margin": margin,
            "conservative_zero_origin_rank_among_64": rank,
            "empirical_origin_rank_p": rank / len(ORIGINS),
            "numerical_tie_tolerance": tolerance,
            "passes": bool(rank == 1 and margin > tolerance),
            "probability_or_likelihood_claim": False,
            "calibrated_uncertainty_claim": False,
            "density_control": (
                "IDENTICAL_BASIS_AND_DETERMINANT_FOR_ZERO_AND_ALL_NULL_ORIGINS"
            ),
        }
    )


def build_cauchy_complexity_audit(old):
    def cauchy_complexity_audit(
        q: np.ndarray,
        basis: np.ndarray,
        selected_phase: dict[str, Any],
        rule: dict[str, Any],
    ) -> dict[str, Any]:
        del selected_phase
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
        selected_support = old.zero_origin_support(q, basis, cutoff)
        supported: list[dict[str, Any]] = []
        best_relative = 0.0
        tested = old.upper_hnf_coarsenings(maximum_index)
        for transform in tested:
            countermodel = basis @ transform
            support = old.zero_origin_support(q, countermodel, cutoff)
            relative = (
                support / selected_support if selected_support > 0.0 else 0.0
            )
            best_relative = max(best_relative, relative)
            if relative < minimum_relative:
                continue
            phase = cauchy_phase_specificity(q, countermodel, rule)
            if not phase["passes"]:
                continue
            supported.append(
                {
                    "integer_index": int(
                        round(float(np.linalg.det(transform)))
                    ),
                    "upper_hnf_transform": transform.tolist(),
                    "absolute_reciprocal_determinant": float(
                        abs(np.linalg.det(countermodel))
                    ),
                    "zero_origin_hard_support_fraction": support,
                    "hard_support_relative_to_selected_basis": relative,
                    "phase_specificity": phase,
                }
            )
        supported.sort(
            key=lambda item: (
                -item["hard_support_relative_to_selected_basis"],
                item["integer_index"],
                item["upper_hnf_transform"],
            )
        )
        allowed = int(
            contract["allowed_supported_coarser_countermodel_count"]
        )
        return stable(
            {
                "tested_upper_hnf_countermodel_count": len(tested),
                "maximum_integer_index_tested": maximum_index,
                "minimum_required_hard_support_relative_to_selected_basis": (
                    minimum_relative
                ),
                "selected_basis_zero_origin_hard_support_fraction": (
                    selected_support
                ),
                "best_countermodel_hard_support_relative_to_selected_basis": (
                    best_relative
                ),
                "supported_coarser_countermodel_count": len(supported),
                "supported_coarser_countermodels": supported,
                "passes": len(supported) <= allowed,
                "selected_basis_absolute_reciprocal_determinant": float(
                    abs(np.linalg.det(basis))
                ),
                "phase_gate": "FIXED_CAUCHY_0P003_UNIQUE_FIRST",
                "relative_support_gate": (
                    "UNCHANGED_HARD_0P003_ZERO_ORIGIN_SUPPORT_RATIO"
                ),
                "scientific_role": contract["scientific_role"],
            }
        )

    return cauchy_complexity_audit


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-key", required=True)
    arguments = parser.parse_args()

    protocol = load(PROTOCOL_PATH)
    validate_sources(protocol)
    evaluator = load_evaluator()
    evaluator.old.phase_specificity = cauchy_phase_specificity
    evaluator.old.primitive_complexity_audit = build_cauchy_complexity_audit(
        evaluator.old
    )

    roster = evaluator.load(
        ROOT / protocol["source_bindings"]["case_roster"]["path"]
    )
    descriptors = [
        row
        for row in roster["cases"]
        if row["case_key"] == arguments.case_key
    ]
    if len(descriptors) != 1:
        raise RuntimeError("case key does not resolve uniquely")
    descriptor = descriptors[0]
    case_path = ROOT / descriptor["case_path"]
    if sha256_file(case_path) != descriptor["case_sha256"]:
        raise RuntimeError("case identity mismatch")
    case = evaluator.load(case_path)

    rule = evaluator.load(
        ROOT
        / protocol["source_bindings"]["multiscale_bidirectional_rule"][
            "path"
        ]
    )
    if tuple(rule["seed_floors_cycles_per_angstrom"]) != FLOORS:
        raise RuntimeError("A0 floor ladder mismatch")
    frozen = evaluator.load(
        ROOT / protocol["source_bindings"]["frozen_d5_kernel_rule"]["path"]
    )
    successor = evaluator.load(
        ROOT / protocol["source_bindings"]["fixed_d5_successor_rule"]["path"]
    )
    portability = evaluator.old.load_portability(
        ROOT / protocol["source_bindings"]["portability_source"]["path"]
    )

    started = time.monotonic()
    postdecision_truth_basis = descriptor.get("truth_basis")
    if postdecision_truth_basis is None:
        postdecision_truth_basis = case.get("truth", {}).get("basis")
    result = evaluator.evaluate_case(
        case,
        postdecision_truth_basis,
        portability,
        frozen,
        successor,
        rule,
    )
    result["case_key"] = descriptor["case_key"]
    result["group"] = descriptor["group"]
    result["model"] = descriptor["model"]
    result["protocol_sha256"] = PROTOCOL_SHA256
    result["case_sha256"] = descriptor["case_sha256"]
    result["intervention"] = "FIXED_CAUCHY_0P003_FULL_STACK"
    result["wall_elapsed_seconds"] = time.monotonic() - started
    result["runtime"] = {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "platform": platform.platform(),
    }

    scientific_projection = copy.deepcopy(result)
    scientific_projection.pop("elapsed_seconds", None)
    scientific_projection.pop("wall_elapsed_seconds", None)
    scientific_projection.pop("runtime", None)
    for scale_row in scientific_projection["scales"]:
        scale_row.pop("elapsed_seconds", None)
    result["scientific_semantic_sha256"] = semantic_sha256(
        scientific_projection
    )

    output = WORK / "case_results" / f"{arguments.case_key}.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(canonical_bytes(result))
    print(
        json.dumps(
            {
                "case_key": arguments.case_key,
                "decision": result["decision"]["decision"],
                "reason": result["decision"]["reason"],
                "correct": result["truth_scoring"].get("correct"),
                "required_outcome_satisfied": result[
                    "required_outcome_satisfied"
                ],
                "incorrect_recovery": result["incorrect_recovery"],
                "scientific_semantic_sha256": result[
                    "scientific_semantic_sha256"
                ],
                "output": output.as_posix(),
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
