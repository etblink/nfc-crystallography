#!/usr/bin/env python3
"""Scan-local, increment-sign-aware diagnostic for the exploratory D4.5 rule.

This is not the already evaluated D4.5 implementation.  It changes only the
enumeration of candidate consecutive-frame pairs:

* frames are partitioned by their frozen D4 scan-configuration identity;
* each partition is ordered in the direction declared by its angle increment;
* cross-scan associations are never considered.

Formal-overlap association, ambiguity handling, minimum unique partitioning,
centroid construction, and every downstream frozen-D5 operation remain those
of ``d45_aggregation.py``.
"""

from __future__ import annotations

import hashlib
from collections import Counter, defaultdict
from typing import Any

import numpy as np

from d45_aggregation import (
    _aggregate_peak,
    _consecutive,
    _overlap_candidates,
    _paths_from_edges,
    minimum_unique_partition,
    semantic_sha256,
)


def scan_direction_chains(
    frames: dict[str, dict[str, Any]],
) -> tuple[list[list[str]], dict[str, int], dict[str, Any]]:
    by_scan: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for frame in frames.values():
        by_scan[str(frame["scan_configuration_id"])].append(frame)

    chains: list[list[str]] = []
    ordinal: dict[str, int] = {}
    scan_records: list[dict[str, Any]] = []
    next_ordinal = 0
    for scan_id in sorted(by_scan):
        local = by_scan[scan_id]
        increments = {
            float(frame["geometry"]["angle_increment_degrees"])
            for frame in local
        }
        if not increments or 0.0 in increments:
            raise ValueError(f"zero or absent increment in scan {scan_id}")
        signs = {1 if value > 0.0 else -1 for value in increments}
        if len(signs) != 1:
            raise ValueError(f"mixed increment signs in scan {scan_id}")
        sign = signs.pop()
        ordered = sorted(
            local,
            key=lambda frame: (
                sign * float(frame["geometry"]["start_angle_degrees"]),
                str(frame["frame_id"]),
            ),
        )
        frame_ids = [str(frame["frame_id"]) for frame in ordered]
        if len(frame_ids) != len(set(frame_ids)):
            raise ValueError(f"duplicate frame identity in scan {scan_id}")
        for frame_id in frame_ids:
            ordinal[frame_id] = next_ordinal
            next_ordinal += 1
        chains.append(frame_ids)
        scan_records.append(
            {
                "scan_configuration_id": scan_id,
                "frame_count": len(frame_ids),
                "increment_values_degrees": sorted(increments),
                "increment_sign": "POSITIVE" if sign > 0 else "NEGATIVE",
                "first_start_angle_degrees": float(
                    ordered[0]["geometry"]["start_angle_degrees"]
                ),
                "last_start_angle_degrees": float(
                    ordered[-1]["geometry"]["start_angle_degrees"]
                ),
            }
        )
    if len(ordinal) != len(frames):
        raise RuntimeError("scan chains are not a total frame cover")
    return (
        chains,
        ordinal,
        {
            "scan_chain_count": len(chains),
            "scan_chains": scan_records,
            "cross_scan_candidate_pairs_considered": False,
        },
    )


def aggregate_d4_scan_direction_aware(
    q: np.ndarray,
    radius: np.ndarray,
    peaks: list[dict[str, Any]],
    frames: dict[str, dict[str, Any]],
) -> tuple[np.ndarray, list[dict[str, Any]], dict[str, Any]]:
    """Run the separately labelled scan-direction diagnostic."""

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

    chains, frame_ordinal, chain_record = scan_direction_chains(frames)
    by_frame: dict[str, list[int]] = defaultdict(list)
    for index, peak in enumerate(peaks):
        by_frame[peak["f"]].append(index)
    for indices in by_frame.values():
        indices.sort(key=lambda index: peaks[index]["peak_id"])

    candidate_pairs: list[tuple[int, int]] = []
    ambiguous_nodes: set[int] = set()
    candidate_degree_histogram: Counter[int] = Counter()
    nonconsecutive_frame_boundaries = 0
    declared_frame_boundaries = 0
    for chain in chains:
        for left_id, right_id in zip(chain[:-1], chain[1:], strict=True):
            declared_frame_boundaries += 1
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
    paths = _paths_from_edges(len(peaks), edges, frame_ordinal, peaks)

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
        partition = minimum_unique_partition(path, q, radius, peaks)
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
        _aggregate_peak(members, q, radius, peaks, frames, disposition)
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
            "RAW_ONLY_CELL_FREE_SCAN_LOCAL_INCREMENT_SIGN_AWARE_"
            "UNIQUE_FORMAL_CONSENSUS_AGGREGATION_DIAGNOSTIC"
        ),
        "diagnostic_status": (
            "SEPARATELY_VERSIONED_EXPLORATORY_DIAGNOSTIC_NOT_THE_"
            "PUBLISHED_D45_0_1_0_CONSTRUCTION"
        ),
        "input_peak_count": len(peaks),
        "output_object_count": len(aggregate_peaks),
        "declared_within_scan_frame_boundary_count": (
            declared_frame_boundaries
        ),
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
        **chain_record,
    }
    return aggregate_q, aggregate_peaks, construction_record
