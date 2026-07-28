from __future__ import annotations

import itertools
import json
import math
from collections.abc import Iterable
from dataclasses import dataclass
from functools import cache
from pathlib import Path
from typing import Any

import numpy as np

from nfc_cryst.canonical import semantic_sha256, sha256_file


@dataclass(frozen=True)
class GateThresholds:
    minimum_assigned_fraction: float = 0.5
    maximum_relative_metric_difference: float = 0.05
    maximum_orientation_difference_degrees: float = 2.0
    basis_search_bound: int = 1


def _determinant_integer(matrix: np.ndarray) -> int:
    return int(round(float(np.linalg.det(matrix))))


@cache
def proper_unimodular_operators(bound: int = 1) -> tuple[np.ndarray, ...]:
    """Enumerate proper unimodular 3x3 matrices inside a finite box.

    The declared finite search includes signed axis permutations, common
    crystallographic indexing ambiguities, point-group operators represented
    in the input basis, and small shear/change-of-basis transforms. The search
    is deliberately bounded and therefore auditable. A result cannot claim
    invariance outside this enumerated set.
    """

    if bound < 1:
        raise ValueError("basis-search bound must be at least one")
    values = range(-bound, bound + 1)
    operators: list[np.ndarray] = []
    for entries in itertools.product(values, repeat=9):
        matrix = np.asarray(entries, dtype=int).reshape(3, 3)
        if _determinant_integer(matrix) == 1:
            operators.append(matrix)
    operators.sort(key=lambda item: tuple(int(x) for x in item.ravel()))
    return tuple(operators)


def reciprocal_basis_from_direct_vectors(
    a: Iterable[float],
    b: Iterable[float],
    c: Iterable[float],
) -> np.ndarray:
    direct = np.column_stack(
        [
            np.asarray(tuple(a), float),
            np.asarray(tuple(b), float),
            np.asarray(tuple(c), float),
        ]
    )
    if direct.shape != (3, 3) or abs(float(np.linalg.det(direct))) < 1e-15:
        raise ValueError("direct-space basis is singular or malformed")
    return np.linalg.inv(direct).T


def reciprocal_basis_from_dials_expt(path: Path) -> tuple[np.ndarray, dict[str, Any]]:
    """Read a single refined crystal model from a DIALS experiment JSON.

    This reader intentionally avoids a DIALS runtime. If an experiment file
    contains multiple nonidentical crystal models, the gate fails closed.
    """

    document = json.loads(path.read_text(encoding="utf-8"))
    crystals = document.get("crystal", [])
    if not crystals:
        raise ValueError(f"{path}: no crystal model")
    bases: list[np.ndarray] = []
    for crystal in crystals:
        try:
            bases.append(
                reciprocal_basis_from_direct_vectors(
                    crystal["real_space_a"],
                    crystal["real_space_b"],
                    crystal["real_space_c"],
                )
            )
        except KeyError as exc:
            raise ValueError(f"{path}: incomplete crystal model") from exc
    reference = bases[0]
    for basis in bases[1:]:
        difference = np.linalg.norm(reference - basis) / max(
            np.linalg.norm(reference), 1e-15
        )
        if difference > 1e-10:
            raise ValueError(f"{path}: multiple nonidentical crystal models")
    metadata = {
        "file_sha256": sha256_file(path),
        "crystal_model_count": len(crystals),
        "space_group_hall_symbols": sorted(
            {str(item.get("space_group_hall_symbol", "UNKNOWN")) for item in crystals}
        ),
    }
    return reference, metadata


def relative_metric_difference(first: np.ndarray, second: np.ndarray) -> float:
    first_metric = first.T @ first
    second_metric = second.T @ second
    return float(
        np.linalg.norm(first_metric - second_metric)
        / max(np.linalg.norm(second_metric), 1e-15)
    )


def _proper_polar_rotation(linear_map: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    left, singular_values, right_t = np.linalg.svd(linear_map)
    rotation = left @ right_t
    if float(np.linalg.det(rotation)) < 0:
        left[:, -1] *= -1
        rotation = left @ right_t
    return rotation, singular_values


def rotation_angle_degrees(rotation: np.ndarray) -> float:
    cosine = float(np.clip((np.trace(rotation) - 1.0) / 2.0, -1.0, 1.0))
    return math.degrees(math.acos(cosine))


def compare_reciprocal_bases(
    reference: np.ndarray,
    other: np.ndarray,
    *,
    metric_threshold: float = 0.05,
    search_bound: int = 1,
) -> dict[str, Any]:
    """Compare two oriented reciprocal bases modulo a finite GL(3,Z) set.

    Every candidate transformation is applied before the physical orientation
    residual is computed. Among metric-admissible transforms, the function
    selects the minimum orientation angle, with metric difference and
    lexicographic matrix order as deterministic tie breakers.
    """

    reference = np.asarray(reference, float)
    other = np.asarray(other, float)
    if reference.shape != (3, 3) or other.shape != (3, 3):
        raise ValueError("reciprocal bases must be 3x3 matrices")

    admissible: list[tuple[tuple[Any, ...], dict[str, Any]]] = []
    nearest_metric: tuple[float, np.ndarray] | None = None
    operators = proper_unimodular_operators(search_bound)
    for operator in operators:
        transformed = reference @ operator
        metric = relative_metric_difference(transformed, other)
        if nearest_metric is None or metric < nearest_metric[0]:
            nearest_metric = (metric, operator)
        if metric > metric_threshold:
            continue
        relative = other @ np.linalg.inv(transformed)
        rotation, singular_values = _proper_polar_rotation(relative)
        angle = rotation_angle_degrees(rotation)
        stretch = float(np.linalg.norm(singular_values - 1.0))
        matrix_key = tuple(int(x) for x in operator.ravel())
        record = {
            "basis_operator": operator.tolist(),
            "operator_determinant": _determinant_integer(operator),
            "relative_metric_difference": metric,
            "orientation_difference_degrees": angle,
            "polar_stretch_frobenius": stretch,
        }
        admissible.append(((angle, metric, stretch, matrix_key), record))

    if not admissible:
        assert nearest_metric is not None
        return {
            "comparison_status": "NO_METRIC_ADMISSIBLE_BASIS_OPERATOR",
            "basis_search_bound": search_bound,
            "operator_count": len(operators),
            "metric_threshold": metric_threshold,
            "nearest_relative_metric_difference": nearest_metric[0],
            "nearest_basis_operator": nearest_metric[1].tolist(),
            "selected": None,
        }

    admissible.sort(key=lambda item: item[0])
    return {
        "comparison_status": "BASIS_INVARIANT_COMPARISON_AVAILABLE",
        "basis_search_bound": search_bound,
        "operator_count": len(operators),
        "metric_threshold": metric_threshold,
        "metric_admissible_operator_count": len(admissible),
        "selected": admissible[0][1],
    }


def qualify_experiments(
    paths: dict[str, Path],
    assigned_fractions: dict[str, float],
    thresholds: GateThresholds | None = None,
) -> dict[str, Any]:
    if thresholds is None:
        thresholds = GateThresholds()
    labels = ("full", "half_a", "half_b")
    if set(paths) != set(labels) or set(assigned_fractions) != set(labels):
        raise ValueError("full, half_a, and half_b are required")
    if any(not 0.0 <= assigned_fractions[label] <= 1.0 for label in labels):
        raise ValueError("assigned fractions must lie in [0, 1]")

    bases: dict[str, np.ndarray] = {}
    records: dict[str, Any] = {}
    for label in labels:
        bases[label], metadata = reciprocal_basis_from_dials_expt(paths[label])
        records[label] = {
            **metadata,
            "assigned_fraction": assigned_fractions[label],
        }

    pairs = (("full", "half_a"), ("full", "half_b"), ("half_a", "half_b"))
    comparisons = {
        f"{first}__{second}": compare_reciprocal_bases(
            bases[first],
            bases[second],
            metric_threshold=thresholds.maximum_relative_metric_difference,
            search_bound=thresholds.basis_search_bound,
        )
        for first, second in pairs
    }
    comparison_available = all(
        value["selected"] is not None for value in comparisons.values()
    )
    maximum_metric = (
        max(
            value["selected"]["relative_metric_difference"]
            for value in comparisons.values()
        )
        if comparison_available
        else None
    )
    maximum_orientation = (
        max(
            value["selected"]["orientation_difference_degrees"]
            for value in comparisons.values()
        )
        if comparison_available
        else None
    )
    minimum_assigned = min(assigned_fractions.values())
    passed = bool(
        comparison_available
        and minimum_assigned >= thresholds.minimum_assigned_fraction
        and maximum_metric is not None
        and maximum_metric <= thresholds.maximum_relative_metric_difference
        and maximum_orientation is not None
        and maximum_orientation <= thresholds.maximum_orientation_difference_degrees
    )
    result = {
        "artifact_type": "BASIS_INVARIANT_CONVENTIONAL_QUALIFICATION",
        "algorithm": {
            "reciprocal_basis_action": "REFERENCE_BASIS_RIGHT_MULTIPLIED_BY_OPERATOR",
            "operator_family": "PROPER_UNIMODULAR_INTEGER_MATRICES_IN_FINITE_BOX",
            "basis_search_bound": thresholds.basis_search_bound,
            "orientation_residual": "PROPER_POLAR_ROTATION_AFTER_BASIS_TRANSFORM",
            "fail_closed": True,
            "selector_or_nfc_method_change": False,
        },
        "records": records,
        "pairwise_comparisons": comparisons,
        "minimum_assigned_fraction": minimum_assigned,
        "maximum_pairwise_relative_metric_difference": maximum_metric,
        "maximum_pairwise_invariant_orientation_difference_degrees": (
            maximum_orientation
        ),
        "thresholds": {
            "minimum_assigned_fraction": thresholds.minimum_assigned_fraction,
            "maximum_relative_metric_difference": (
                thresholds.maximum_relative_metric_difference
            ),
            "maximum_orientation_difference_degrees": (
                thresholds.maximum_orientation_difference_degrees
            ),
        },
        "principal_outcome": (
            "CONVENTIONAL_FULL_AND_SPLIT_INDEXABILITY_STABLE"
            if passed
            else "CONVENTIONAL_FULL_AND_SPLIT_INDEXABILITY_NOT_STABLE"
        ),
    }
    result["semantic_sha256"] = semantic_sha256(result)
    return result
