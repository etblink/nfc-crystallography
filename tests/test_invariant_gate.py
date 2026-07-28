from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from nfc_cryst.conventional_gate.invariant import (
    GateThresholds,
    compare_reciprocal_bases,
    proper_unimodular_operators,
    qualify_experiments,
    reciprocal_basis_from_dials_expt,
)


def rotation_z(degrees: float) -> np.ndarray:
    angle = math.radians(degrees)
    return np.asarray(
        [
            [math.cos(angle), -math.sin(angle), 0.0],
            [math.sin(angle), math.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )


def write_expt(path: Path, reciprocal_basis: np.ndarray) -> None:
    direct = np.linalg.inv(reciprocal_basis).T
    document = {
        "crystal": [
            {
                "real_space_a": direct[:, 0].tolist(),
                "real_space_b": direct[:, 1].tolist(),
                "real_space_c": direct[:, 2].tolist(),
                "space_group_hall_symbol": " P 1",
            }
        ]
    }
    path.write_text(json.dumps(document), encoding="utf-8")


def test_declared_operator_census_is_stable() -> None:
    assert len(proper_unimodular_operators(1)) == 3480


def test_180_degree_indexing_equivalent_is_invariant() -> None:
    basis = np.diag([0.01, 0.015, 0.02])
    operator = np.diag([1, -1, -1])
    result = compare_reciprocal_bases(basis, basis @ operator)
    assert result["selected"] is not None
    assert result["selected"]["orientation_difference_degrees"] < 1e-8


def test_axis_permutation_is_invariant() -> None:
    basis = np.diag([0.01, 0.015, 0.02])
    operator = np.asarray([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
    result = compare_reciprocal_bases(basis, basis @ operator)
    assert result["selected"] is not None
    assert result["selected"]["orientation_difference_degrees"] < 1e-8


def test_genuine_physical_rotation_remains_visible() -> None:
    basis = np.diag([0.01, 0.015, 0.02])
    result = compare_reciprocal_bases(basis, rotation_z(5.0) @ basis)
    assert result["selected"] is not None
    assert result["selected"]["orientation_difference_degrees"] > 4.9


def test_metric_mismatch_fails_closed() -> None:
    first = np.diag([0.01, 0.015, 0.02])
    second = np.diag([0.02, 0.03, 0.04])
    result = compare_reciprocal_bases(first, second, metric_threshold=0.01)
    assert result["selected"] is None
    assert result["comparison_status"] == "NO_METRIC_ADMISSIBLE_BASIS_OPERATOR"


def test_expt_reader_and_full_gate(tmp_path: Path) -> None:
    basis = np.asarray([[0.011, 0.001, 0.0], [0.0, 0.016, 0.002], [0.001, 0.0, 0.021]])
    transforms = {
        "full": np.eye(3, dtype=int),
        "half_a": np.diag([1, -1, -1]),
        "half_b": np.asarray([[0, 1, 0], [0, 0, 1], [1, 0, 0]]),
    }
    paths = {}
    for label, transform in transforms.items():
        path = tmp_path / f"{label}.expt"
        write_expt(path, basis @ transform)
        paths[label] = path
        observed, metadata = reciprocal_basis_from_dials_expt(path)
        assert observed.shape == (3, 3)
        assert metadata["crystal_model_count"] == 1
    result = qualify_experiments(
        paths,
        {"full": 0.9, "half_a": 0.8, "half_b": 0.85},
        GateThresholds(),
    )
    assert (
        result["principal_outcome"] == "CONVENTIONAL_FULL_AND_SPLIT_INDEXABILITY_STABLE"
    )
    assert result["maximum_pairwise_invariant_orientation_difference_degrees"] < 1e-8
