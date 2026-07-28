#!/usr/bin/env python3
"""Exploratory raw-only D4.5 aggregation for public D5 controls.

This module does not modify the frozen D4 or D5 implementation.  It consumes
only D4 peak records, frame geometry, and the reciprocal coordinates/formal
regions produced from them.  Conventional indexing and DIALS reflections are
not construction inputs.

The construction is intentionally conservative:

1. Candidate associations exist only between consecutive observed frames.
2. Two detections are candidates only when their formal reciprocal-space
   regions overlap.
3. Any detection with multiple candidates on either side is explicitly
   unresolved and is not merged.
4. The remaining one-to-one links form chronological paths.
5. A proposed path segment is admissible only when its signal-weighted
   reciprocal centroid lies inside every member's formal region.
6. Each path is partitioned into the minimum number of admissible contiguous
   segments.  If that minimum partition is not unique, the path is explicitly
   unresolved and all of its detections remain separate.

The rule prevents pairwise-overlap chaining from silently becoming a single
three-dimensional object while avoiding any cell, orientation, or expected
answer input.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Any, Iterable

import numpy as np
from scipy.spatial import cKDTree


JSON_KWARGS = {
    "ensure_ascii": False,
    "sort_keys": True,
    "separators": (",", ":"),
}


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, **JSON_KWARGS) + "\n").encode("utf-8")


def semantic_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def _frame_key(frame: dict[str, Any]) -> tuple[float, str]:
    return (
        float(frame["geometry"]["start_angle_degrees"]),
        str(frame["frame_id"]),
    )


def ordered_frame_ids(
    frames: dict[str, dict[str, Any]],
) -> list[str]:
    return [
        frame["frame_id"]
        for frame in sorted(frames.values(), key=_frame_key)
    ]


def _consecutive(
    left: dict[str, Any],
    right: dict[str, Any],
    absolute_tolerance_degrees: float = 1e-10,
) -> bool:
    left_start = float(left["geometry"]["start_angle_degrees"])
    left_increment = float(left["geometry"]["angle_increment_degrees"])
    right_start = float(right["geometry"]["start_angle_degrees"])
    return math.isclose(
        left_start + left_increment,
        right_start,
        rel_tol=0.0,
        abs_tol=absolute_tolerance_degrees,
    )


def _positive_signal(peak: dict[str, Any]) -> float:
    value = float(peak["Shat"])
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(
            f"nonpositive or nonfinite D4 Shat for {peak['peak_id']}"
        )
    return value


def weighted_centroid(
    members: list[int],
    q: np.ndarray,
    peaks: list[dict[str, Any]],
) -> np.ndarray:
    weights = np.asarray(
        [_positive_signal(peaks[index]) for index in members],
        dtype=float,
    )
    return np.sum(q[members] * weights[:, None], axis=0) / np.sum(weights)


def centroid_supported_by_all(
    members: list[int],
    q: np.ndarray,
    radius: np.ndarray,
    peaks: list[dict[str, Any]],
) -> bool:
    center = weighted_centroid(members, q, peaks)
    distances = np.linalg.norm(q[members] - center[None, :], axis=1)
    return bool(np.all(distances <= radius[members]))


@dataclass(frozen=True)
class Partition:
    segments: tuple[tuple[int, ...], ...]
    unique: bool
    minimum_segment_count: int


def minimum_unique_partition(
    path: list[int],
    q: np.ndarray,
    radius: np.ndarray,
    peaks: list[dict[str, Any]],
) -> Partition:
    """Return the unique minimum admissible partition, if one exists.

    Solution counts are capped at two because only uniqueness matters.
    """

    count = len(path)
    best = [count + 1] * (count + 1)
    ways = [0] * (count + 1)
    predecessor: list[int | None] = [None] * (count + 1)
    best[0] = 0
    ways[0] = 1
    valid: dict[tuple[int, int], bool] = {}

    for end in range(1, count + 1):
        for start in range(end):
            key = (start, end)
            if key not in valid:
                valid[key] = centroid_supported_by_all(
                    path[start:end], q, radius, peaks
                )
            if not valid[key]:
                continue
            candidate = best[start] + 1
            if candidate < best[end]:
                best[end] = candidate
                ways[end] = ways[start]
                predecessor[end] = start
            elif candidate == best[end]:
                ways[end] = min(2, ways[end] + ways[start])

    if ways[count] != 1:
        return Partition(
            segments=tuple((index,) for index in path),
            unique=False,
            minimum_segment_count=best[count],
        )

    segments: list[tuple[int, ...]] = []
    end = count
    while end:
        start = predecessor[end]
        if start is None:
            raise RuntimeError("unique partition lacks a predecessor")
        segments.append(tuple(path[start:end]))
        end = start
    segments.reverse()
    return Partition(
        segments=tuple(segments),
        unique=True,
        minimum_segment_count=best[count],
    )


def _overlap_candidates(
    left: list[int],
    right: list[int],
    q: np.ndarray,
    radius: np.ndarray,
) -> tuple[dict[int, list[int]], dict[int, list[int]]]:
    left_to_right = {index: [] for index in left}
    right_to_left = {index: [] for index in right}
    if not left or not right:
        return left_to_right, right_to_left

    right_array = np.asarray(right, dtype=np.int64)
    tree = cKDTree(q[right_array])
    largest_right_radius = float(np.max(radius[right_array]))
    for left_index in left:
        candidates = tree.query_ball_point(
            q[left_index],
            float(radius[left_index]) + largest_right_radius,
        )
        for local_right in sorted(candidates):
            right_index = int(right_array[local_right])
            distance = float(
                np.linalg.norm(q[left_index] - q[right_index])
            )
            if distance <= float(
                radius[left_index] + radius[right_index]
            ):
                left_to_right[left_index].append(right_index)
                right_to_left[right_index].append(left_index)
    for values in left_to_right.values():
        values.sort()
    for values in right_to_left.values():
        values.sort()
    return left_to_right, right_to_left


def _paths_from_edges(
    count: int,
    edges: list[tuple[int, int]],
    frame_ordinal: dict[str, int],
    peaks: list[dict[str, Any]],
) -> list[list[int]]:
    adjacency: list[list[int]] = [[] for _ in range(count)]
    for left, right in edges:
        adjacency[left].append(right)
        adjacency[right].append(left)
    if any(len(neighbors) > 2 for neighbors in adjacency):
        raise RuntimeError("one-to-one consecutive-frame graph branched")

    visited: set[int] = set()
    paths: list[list[int]] = []
    for start in range(count):
        if start in visited:
            continue
        component: list[int] = []
        stack = [start]
        visited.add(start)
        while stack:
            current = stack.pop()
            component.append(current)
            for neighbor in adjacency[current]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    stack.append(neighbor)
        component.sort(
            key=lambda index: (
                frame_ordinal[peaks[index]["f"]],
                peaks[index]["peak_id"],
            )
        )
        if len(
            {frame_ordinal[peaks[index]["f"]] for index in component}
        ) != len(component):
            raise RuntimeError("path contains multiple detections in one frame")
        paths.append(component)
    paths.sort(
        key=lambda path: tuple(peaks[index]["peak_id"] for index in path)
    )
    return paths


def _representative_member(
    members: list[int],
    peaks: list[dict[str, Any]],
    frames: dict[str, dict[str, Any]],
) -> int:
    weights = np.asarray(
        [_positive_signal(peaks[index]) for index in members],
        dtype=float,
    )
    angles = np.asarray(
        [
            float(frames[peaks[index]["f"]]["geometry"][
                "start_angle_degrees"
            ])
            + 0.5
            * float(frames[peaks[index]["f"]]["geometry"][
                "angle_increment_degrees"
            ])
            for index in members
        ],
        dtype=float,
    )
    center = float(np.sum(weights * angles) / np.sum(weights))
    return min(
        members,
        key=lambda index: (
            abs(
                float(
                    frames[peaks[index]["f"]]["geometry"][
                        "start_angle_degrees"
                    ]
                )
                + 0.5
                * float(
                    frames[peaks[index]["f"]]["geometry"][
                        "angle_increment_degrees"
                    ]
                )
                - center
            ),
            peaks[index]["peak_id"],
        ),
    )


def _aggregate_peak(
    members: list[int],
    q: np.ndarray,
    radius: np.ndarray,
    peaks: list[dict[str, Any]],
    frames: dict[str, dict[str, Any]],
    disposition: str,
) -> tuple[np.ndarray, dict[str, Any], dict[str, Any]]:
    member_ids = sorted(peaks[index]["peak_id"] for index in members)
    aggregate_id = (
        "D45_"
        + hashlib.sha256(canonical_bytes(member_ids)).hexdigest()[:28]
    )
    center = weighted_centroid(members, q, peaks)
    representative = _representative_member(members, peaks, frames)
    signal = float(sum(float(peaks[index]["Shat"]) for index in members))
    variance = float(
        sum(float(peaks[index]["sigma_S"]) ** 2 for index in members)
    )
    if variance <= 0.0 or not math.isfinite(variance):
        raise ValueError(f"invalid aggregate variance for {aggregate_id}")
    weights = np.asarray(
        [_positive_signal(peaks[index]) for index in members],
        dtype=float,
    )
    x = float(
        np.sum(
            weights
            * np.asarray([float(peaks[index]["x"]) for index in members])
        )
        / np.sum(weights)
    )
    y = float(
        np.sum(
            weights
            * np.asarray([float(peaks[index]["y"]) for index in members])
        )
        / np.sum(weights)
    )
    residuals = np.linalg.norm(q[members] - center[None, :], axis=1)
    normalized = residuals / radius[members]
    peak = {
        "peak_id": aggregate_id,
        "f": peaks[representative]["f"],
        "x": x,
        "y": y,
        "Shat": signal,
        "sigma_S": math.sqrt(variance),
        "integrated_signal_to_formal_noise": signal / math.sqrt(variance),
        "q": {
            "validity": "CERTIFIED_POSITIVE_SEPARATION_INTERIOR",
        },
        "d45": {
            "disposition": disposition,
            "member_count": len(members),
            "member_peak_ids": member_ids,
            "representative_member_peak_id": peaks[representative][
                "peak_id"
            ],
        },
    }
    record = {
        "aggregate_peak_id": aggregate_id,
        "disposition": disposition,
        "member_count": len(members),
        "member_peak_ids": member_ids,
        "representative_frame_id": peak["f"],
        "aggregate_q": center.tolist(),
        "maximum_member_residual_cycles_per_angstrom": float(
            np.max(residuals, initial=0.0)
        ),
        "maximum_normalized_formal_residual": float(
            np.max(normalized, initial=0.0)
        ),
    }
    return center, peak, record


def aggregate_d4(
    q: np.ndarray,
    radius: np.ndarray,
    peaks: list[dict[str, Any]],
    frames: dict[str, dict[str, Any]],
) -> tuple[np.ndarray, list[dict[str, Any]], dict[str, Any]]:
    """Construct deterministic exploratory D4.5 objects."""

    q = np.asarray(q, dtype=float)
    radius = np.asarray(radius, dtype=float)
    if q.shape != (len(peaks), 3):
        raise ValueError("q/peak cardinality or shape mismatch")
    if radius.shape != (len(peaks),):
        raise ValueError("radius/peak cardinality mismatch")
    if not np.all(np.isfinite(q)) or not np.all(np.isfinite(radius)):
        raise ValueError("nonfinite q or radius")
    if np.any(radius <= 0.0):
        raise ValueError("formal reciprocal radius must be positive")
    if len({peak["peak_id"] for peak in peaks}) != len(peaks):
        raise ValueError("duplicate D4 peak identity")
    if any(peak["f"] not in frames for peak in peaks):
        raise ValueError("D4 peak refers to an absent frame")

    frame_ids = ordered_frame_ids(frames)
    frame_ordinal = {
        frame_id: index for index, frame_id in enumerate(frame_ids)
    }
    by_frame: dict[str, list[int]] = defaultdict(list)
    for index, peak in enumerate(peaks):
        by_frame[peak["f"]].append(index)
    for indices in by_frame.values():
        indices.sort(key=lambda index: peaks[index]["peak_id"])

    candidate_pairs: list[tuple[int, int]] = []
    ambiguous_nodes: set[int] = set()
    candidate_degree_histogram: Counter[int] = Counter()
    nonconsecutive_frame_boundaries = 0
    for left_id, right_id in zip(frame_ids[:-1], frame_ids[1:], strict=True):
        if not _consecutive(frames[left_id], frames[right_id]):
            nonconsecutive_frame_boundaries += 1
            continue
        left_to_right, right_to_left = _overlap_candidates(
            by_frame.get(left_id, []),
            by_frame.get(right_id, []),
            q,
            radius,
        )
        candidate_degree_histogram.update(
            len(values) for values in left_to_right.values()
        )
        candidate_degree_histogram.update(
            len(values) for values in right_to_left.values()
        )
        ambiguous_nodes.update(
            index
            for index, values in left_to_right.items()
            if len(values) > 1
        )
        ambiguous_nodes.update(
            index
            for index, values in right_to_left.items()
            if len(values) > 1
        )
        for left, rights in left_to_right.items():
            if len(rights) != 1:
                continue
            right = rights[0]
            if right_to_left[right] == [left]:
                candidate_pairs.append((left, right))

    edges = [
        pair
        for pair in candidate_pairs
        if pair[0] not in ambiguous_nodes and pair[1] not in ambiguous_nodes
    ]
    paths = _paths_from_edges(
        len(peaks), edges, frame_ordinal, peaks
    )

    segments: list[tuple[list[int], str]] = []
    ambiguous_partition_paths = 0
    ambiguous_partition_members = 0
    for path in paths:
        if any(index in ambiguous_nodes for index in path):
            for index in path:
                disposition = (
                    "UNRESOLVED_MULTIPLE_FORMAL_OVERLAPS"
                    if index in ambiguous_nodes
                    else "UNMERGED_NEIGHBOR_OF_UNRESOLVED_ASSOCIATION"
                )
                segments.append(([index], disposition))
            continue
        partition = minimum_unique_partition(
            path, q, radius, peaks
        )
        if not partition.unique:
            ambiguous_partition_paths += 1
            ambiguous_partition_members += len(path)
            for index in path:
                segments.append(
                    ([index], "UNRESOLVED_NONUNIQUE_MINIMUM_PARTITION")
                )
        else:
            for segment in partition.segments:
                disposition = (
                    "UNIQUE_FORMAL_CONSENSUS_AGGREGATE"
                    if len(segment) > 1
                    else "UNASSOCIATED_SINGLETON"
                )
                segments.append((list(segment), disposition))

    aggregates = [
        _aggregate_peak(
            members, q, radius, peaks, frames, disposition
        )
        for members, disposition in segments
    ]
    aggregates.sort(key=lambda item: item[1]["peak_id"])
    aggregate_q = np.asarray([item[0] for item in aggregates], dtype=float)
    aggregate_peaks = [item[1] for item in aggregates]
    records = [item[2] for item in aggregates]

    covered = [
        peak_id
        for record in records
        for peak_id in record["member_peak_ids"]
    ]
    expected = sorted(peak["peak_id"] for peak in peaks)
    if sorted(covered) != expected or len(covered) != len(expected):
        raise RuntimeError("D4.5 membership is not a total disjoint cover")

    disposition_counts = Counter(
        record["disposition"] for record in records
    )
    member_count_histogram = Counter(
        record["member_count"] for record in records
    )
    construction_record = {
        "construction": (
            "RAW_ONLY_CELL_FREE_UNIQUE_FORMAL_CONSENSUS_AGGREGATION"
        ),
        "input_peak_count": len(peaks),
        "output_object_count": len(aggregate_peaks),
        "consecutive_frame_candidate_pair_count": len(candidate_pairs),
        "accepted_one_to_one_edge_count": len(edges),
        "ambiguous_input_peak_count": len(ambiguous_nodes),
        "path_count": len(paths),
        "maximum_path_member_count": max(map(len, paths), default=0),
        "ambiguous_partition_path_count": ambiguous_partition_paths,
        "ambiguous_partition_member_count": ambiguous_partition_members,
        "nonconsecutive_frame_boundary_count": (
            nonconsecutive_frame_boundaries
        ),
        "candidate_degree_histogram": {
            str(key): value
            for key, value in sorted(candidate_degree_histogram.items())
        },
        "output_disposition_counts": dict(sorted(disposition_counts.items())),
        "output_member_count_histogram": {
            str(key): value
            for key, value in sorted(member_count_histogram.items())
        },
        "input_peak_id_sha256": semantic_sha256(expected),
        "output_membership_sha256": semantic_sha256(
            [
                {
                    "aggregate_peak_id": record["aggregate_peak_id"],
                    "member_peak_ids": record["member_peak_ids"],
                }
                for record in records
            ]
        ),
        "aggregate_q_sha256": hashlib.sha256(
            np.asarray(aggregate_q, dtype="<f8", order="C").tobytes(order="C")
        ).hexdigest(),
        "cell_or_orientation_input": False,
        "dials_input": False,
        "processed_reflection_input": False,
    }
    return aggregate_q, aggregate_peaks, construction_record


def subset_by_frames(
    q: np.ndarray,
    radius: np.ndarray,
    peaks: list[dict[str, Any]],
    frame_ids: Iterable[str],
) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]]]:
    allowed = set(frame_ids)
    keep = np.asarray([peak["f"] in allowed for peak in peaks], dtype=bool)
    return (
        q[keep],
        radius[keep],
        [
            peak
            for peak, selected in zip(peaks, keep, strict=True)
            if selected
        ],
    )
