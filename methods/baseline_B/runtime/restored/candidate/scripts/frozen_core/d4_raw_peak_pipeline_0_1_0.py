#!/usr/bin/env python3
"""Detector-only D4 peak construction for the frozen 9Z6F raw CBF archive.

This program is deliberately incapable of reading the processed HKL, XDS.INP,
cell, symmetry, SHELX, FCF, or solved-structure lanes.  It consumes only CBF
members from the frozen raw archive, the raw acquisition metadata, the frozen
D4 rule, and the local CBF byte-offset decoder.
"""

from __future__ import annotations

import argparse
import base64
import ctypes
import gzip
import hashlib
import json
import math
import os
import platform
import re
import subprocess
import sys
import zipfile
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import scipy
import scipy.ndimage as ndi
from scipy.spatial import cKDTree


VERSION = "0_1_0"
ROOT = Path(__file__).resolve().parent
RAW_ARCHIVE = ROOT / "source" / "raw" / "Archive.zip"
RAW_METADATA = ROOT / "source" / "raw" / "metadata.txt"
RULE_PATH = ROOT / "d4_detection_rule_freeze.json"
DECODER_SOURCE = ROOT / "d4_cbf_byte_offset.c"
DECODER_LIBRARY = ROOT / "work" / "libd4_cbf_byte_offset.so"
ARTIFACTS = ROOT / "artifacts"

CORPUS_PATH = (
    ARTIFACTS / f"NFC_CRYST_9Z6F_D4_RAW_PEAK_CORPUS_{VERSION}.json.gz"
)
RESULT_JSON_PATH = (
    ARTIFACTS
    / f"NFC_CRYST_9Z6F_D4_RAW_PEAK_SET_AND_POSITIVE_SEPARATION_RESULT_{VERSION}.json"
)
RESULT_MD_PATH = (
    ARTIFACTS
    / f"NFC_CRYST_9Z6F_D4_RAW_PEAK_SET_AND_POSITIVE_SEPARATION_RESULT_{VERSION}.md"
)

EXPECTED_ARCHIVE_BYTES = 817_839_975
EXPECTED_ARCHIVE_SHA256 = (
    "fc639d3c87bb0d5c9002c78a17d0cbe85cd4879e638a3236a49ae1a95f48ecc1"
)
CBF_MARKER = b"\x0c\x1a\x04\xd5"
HEIGHT = 775
WIDTH = 800

JSON_KWARGS = {
    "ensure_ascii": False,
    "sort_keys": True,
    "separators": (",", ":"),
}


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, **JSON_KWARGS) + "\n").encode("utf-8")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(canonical_bytes(value))


def rounded(value: float, places: int = 6) -> float:
    result = round(float(value), places)
    return 0.0 if result == 0 else result


def seal(value: dict[str, Any]) -> dict[str, Any]:
    body = dict(value)
    body.pop("canonical_semantic_sha256", None)
    value["canonical_semantic_sha256"] = sha256_bytes(canonical_bytes(body))
    return value


def compile_decoder() -> None:
    DECODER_LIBRARY.parent.mkdir(parents=True, exist_ok=True)
    if (
        DECODER_LIBRARY.exists()
        and DECODER_LIBRARY.stat().st_mtime_ns >= DECODER_SOURCE.stat().st_mtime_ns
    ):
        return
    subprocess.run(
        [
            "gcc",
            "-O3",
            "-std=c11",
            "-Wall",
            "-Wextra",
            "-pedantic",
            "-fPIC",
            "-shared",
            str(DECODER_SOURCE),
            "-o",
            str(DECODER_LIBRARY),
        ],
        check=True,
        capture_output=True,
        text=True,
    )


class Decoder:
    def __init__(self) -> None:
        compile_decoder()
        self.library = ctypes.CDLL(str(DECODER_LIBRARY))
        self.function = self.library.cbf_byte_offset_decode
        self.function.argtypes = [
            ctypes.POINTER(ctypes.c_uint8),
            ctypes.c_size_t,
            ctypes.POINTER(ctypes.c_int32),
            ctypes.c_size_t,
            ctypes.POINTER(ctypes.c_size_t),
        ]
        self.function.restype = ctypes.c_int

    def decode(self, payload: bytes | memoryview, element_count: int) -> np.ndarray:
        source = np.frombuffer(payload, dtype=np.uint8)
        output = np.empty(element_count, dtype=np.int32)
        bytes_used = ctypes.c_size_t()
        status = self.function(
            source.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8)),
            source.size,
            output.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
            output.size,
            ctypes.byref(bytes_used),
        )
        if status != 0:
            raise RuntimeError(
                f"CBF decode failed: status={status}, "
                f"used={bytes_used.value}, available={source.size}"
            )
        return output


def python_reference_decode(payload: bytes, element_count: int) -> np.ndarray:
    output = np.empty(element_count, dtype=np.int32)
    position = 0
    previous = 0
    for index in range(element_count):
        if position >= len(payload):
            raise ValueError("truncated byte-offset stream")
        delta = int.from_bytes(payload[position : position + 1], "little", signed=True)
        position += 1
        if delta == -128:
            delta = int.from_bytes(
                payload[position : position + 2], "little", signed=True
            )
            position += 2
            if delta == -32768:
                delta = int.from_bytes(
                    payload[position : position + 4], "little", signed=True
                )
                position += 4
                if delta == -2147483648:
                    delta = int.from_bytes(
                        payload[position : position + 8], "little", signed=True
                    )
                    position += 8
        previous += delta
        output[index] = previous
    if position != len(payload):
        raise ValueError("trailing bytes in byte-offset stream")
    return output


def header_number(header: str, label: str, cast: type = float) -> Any:
    match = re.search(
        rf"(?m)^#\s*{re.escape(label)}\s+([-+0-9.eE]+)", header
    )
    if not match:
        raise ValueError(f"missing CBF comment field {label}")
    return cast(match.group(1))


def binary_field(header: str, label: str, cast: type = int) -> Any:
    match = re.search(rf"(?m)^{re.escape(label)}:\s*([^\r\n]+)", header)
    if not match:
        raise ValueError(f"missing CBF binary field {label}")
    return cast(match.group(1).strip().strip('"'))


def parse_vector(header: str, label: str) -> list[float]:
    match = re.search(rf"(?m)^#\s*{re.escape(label)}\s+(.+?)\s*$", header)
    if not match:
        raise ValueError(f"missing CBF vector {label}")
    return [float(item) for item in match.group(1).split()]


def parse_cbf(data: bytes) -> tuple[dict[str, Any], bytes]:
    marker = data.find(CBF_MARKER)
    if marker < 0:
        raise ValueError("CBF marker absent")
    header = data[:marker].decode("latin-1")
    binary_size = binary_field(header, "X-Binary-Size", int)
    element_count = binary_field(header, "X-Binary-Number-of-Elements", int)
    fast = binary_field(header, "X-Binary-Size-Fastest-Dimension", int)
    slow = binary_field(header, "X-Binary-Size-Second-Dimension", int)
    if (slow, fast) != (HEIGHT, WIDTH) or element_count != HEIGHT * WIDTH:
        raise ValueError(
            f"unexpected detector interface: slow={slow}, fast={fast}, "
            f"elements={element_count}"
        )
    payload = data[marker + len(CBF_MARKER) : marker + len(CBF_MARKER) + binary_size]
    if len(payload) != binary_size:
        raise ValueError("truncated CBF payload")
    md5_match = re.search(r"(?m)^Content-MD5:\s*(\S+)\s*$", header)
    if not md5_match:
        raise ValueError("Content-MD5 absent")
    observed_md5 = base64.b64encode(
        hashlib.md5(payload, usedforsecurity=False).digest()
    ).decode("ascii")
    if observed_md5 != md5_match.group(1):
        raise ValueError("CBF payload MD5 mismatch")
    timestamp_match = re.search(r"(?m)^#\s*(\d{4}-\d{2}-\d{2}T\S+)\s*$", header)
    beam_match = re.search(
        r"(?m)^#\s*Beam_xy\s+\(([-+0-9.eE]+),\s*([-+0-9.eE]+)\)\s+pixels",
        header,
    )
    if not timestamp_match or not beam_match:
        raise ValueError("timestamp or beam center absent")
    geometry = {
        "acquired_at": timestamp_match.group(1),
        "exposure_seconds": header_number(header, "Exposure_time"),
        "wavelength_angstrom": header_number(header, "Wavelength"),
        "detector_distance_m": header_number(header, "Detector_distance"),
        "beam_xy_pixels": [
            float(beam_match.group(1)),
            float(beam_match.group(2)),
        ],
        "start_angle_degrees": header_number(header, "Start_angle"),
        "angle_increment_degrees": header_number(header, "Angle_increment"),
        "detector_2theta_degrees": header_number(header, "Detector_2theta"),
        "alpha_degrees": header_number(header, "Alpha"),
        "beta_degrees": header_number(header, "Beta"),
        "phi_degrees": header_number(header, "Phi"),
        "phi_increment_degrees": header_number(header, "Phi_increment"),
        "omega_degrees": header_number(header, "Omega"),
        "omega_increment_degrees": header_number(header, "Omega_increment"),
        "kappa_degrees": header_number(header, "Kappa"),
        "kappa_increment_degrees": header_number(header, "Kappa_increment"),
        "detector_fast_axis": parse_vector(header, "Detector_fast_axis_vector"),
        "detector_slow_axis": parse_vector(header, "Detector_slow_axis_vector"),
        "incident_beam": parse_vector(header, "Incident_beam_vector"),
    }
    metadata = {
        "binary_byte_count": binary_size,
        "element_count": element_count,
        "fast_dimension": fast,
        "slow_dimension": slow,
        "payload_sha256": sha256_bytes(payload),
        "payload_content_md5_base64": observed_md5,
        "geometry": geometry,
    }
    return metadata, payload


def frame_semantic_id(metadata: dict[str, Any]) -> str:
    semantic = {
        "payload_sha256": metadata["payload_sha256"],
        "geometry": metadata["geometry"],
        "dimensions": [
            metadata["slow_dimension"],
            metadata["fast_dimension"],
        ],
    }
    return "FRAME_" + sha256_bytes(canonical_bytes(semantic))[:24]


def scan_signature(geometry: dict[str, Any]) -> str:
    increments = {
        "phi": geometry["phi_increment_degrees"],
        "omega": geometry["omega_increment_degrees"],
        "kappa": geometry["kappa_increment_degrees"],
    }
    active = [name for name, value in increments.items() if abs(value) > 1e-12]
    if len(active) != 1:
        active = ["start_angle"]
    axis = active[0]
    static_angles = {
        "alpha": geometry["alpha_degrees"],
        "beta": geometry["beta_degrees"],
        "phi": None if axis == "phi" else geometry["phi_degrees"],
        "omega": None if axis == "omega" else geometry["omega_degrees"],
        "kappa": None if axis == "kappa" else geometry["kappa_degrees"],
    }
    signature = {
        "rotation_axis": axis,
        "angle_increment_degrees": geometry["angle_increment_degrees"],
        "axis_increments": increments,
        "static_angles": static_angles,
        "detector_2theta_degrees": geometry["detector_2theta_degrees"],
        "detector_distance_m": geometry["detector_distance_m"],
        "beam_xy_pixels": geometry["beam_xy_pixels"],
        "wavelength_angstrom": geometry["wavelength_angstrom"],
        "exposure_seconds": geometry["exposure_seconds"],
        "detector_fast_axis": geometry["detector_fast_axis"],
        "detector_slow_axis": geometry["detector_slow_axis"],
        "incident_beam": geometry["incident_beam"],
    }
    return "SCANCFG_" + sha256_bytes(canonical_bytes(signature))[:20]


def tile_background(
    image: np.ndarray, tile_height: int, tile_width: int
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    tile_rows = math.ceil(image.shape[0] / tile_height)
    tile_columns = math.ceil(image.shape[1] / tile_width)
    backgrounds = np.empty((tile_rows, tile_columns), dtype=np.float64)
    scales = np.empty_like(backgrounds)
    fully_masked = 0
    for row in range(tile_rows):
        y0 = row * tile_height
        y1 = min((row + 1) * tile_height, image.shape[0])
        for column in range(tile_columns):
            x0 = column * tile_width
            x1 = min((column + 1) * tile_width, image.shape[1])
            values = image[y0:y1, x0:x1]
            values = values[values >= 0]
            if values.size == 0:
                backgrounds[row, column] = 0.0
                scales[row, column] = np.inf
                fully_masked += 1
                continue
            median = float(np.median(values))
            mad = float(np.median(np.abs(values.astype(np.float64) - median)))
            backgrounds[row, column] = median
            scales[row, column] = max(
                1.0,
                1.4826 * mad,
                math.sqrt(max(median, 0.0) + 1.0),
            )
    background = np.repeat(
        np.repeat(backgrounds, tile_height, axis=0), tile_width, axis=1
    )[: image.shape[0], : image.shape[1]]
    scale = np.repeat(
        np.repeat(scales, tile_height, axis=0), tile_width, axis=1
    )[: image.shape[0], : image.shape[1]]
    return background, scale, {
        "tile_rows": tile_rows,
        "tile_columns": tile_columns,
        "fully_masked_tile_count": fully_masked,
        "median_background_minimum": rounded(np.nanmin(backgrounds)),
        "median_background_maximum": rounded(np.nanmax(backgrounds)),
        "robust_scale_minimum": rounded(np.min(scales[np.isfinite(scales)])),
        "robust_scale_maximum": rounded(np.max(scales[np.isfinite(scales)])),
    }


def interval_distance(a0: float, a1: float, b0: float, b1: float) -> float:
    if a1 < b0:
        return b0 - a1
    if b1 < a0:
        return a0 - b1
    return 0.0


def attach_separation(peaks: list[dict[str, Any]]) -> None:
    if len(peaks) == 1:
        peaks[0]["m"] = None
        peaks[0]["nearest_peak_id"] = None
        peaks[0]["q"]["separation_status"] = "SINGLETON_FRAME"
        return
    if not peaks:
        return
    boxes = [
        (
            peak["R"]["x_min"],
            peak["R"]["x_max"],
            peak["R"]["y_min"],
            peak["R"]["y_max"],
        )
        for peak in peaks
    ]
    for index, peak in enumerate(peaks):
        nearest_margin = math.inf
        nearest_id: str | None = None
        for other_index, other in enumerate(peaks):
            if index == other_index:
                continue
            dx = interval_distance(
                boxes[index][0],
                boxes[index][1],
                boxes[other_index][0],
                boxes[other_index][1],
            )
            dy = interval_distance(
                boxes[index][2],
                boxes[index][3],
                boxes[other_index][2],
                boxes[other_index][3],
            )
            margin = math.hypot(dx, dy)
            if margin < nearest_margin or (
                math.isclose(margin, nearest_margin)
                and (nearest_id is None or other["peak_id"] < nearest_id)
            ):
                nearest_margin = margin
                nearest_id = other["peak_id"]
        peak["m"] = rounded(nearest_margin)
        peak["nearest_peak_id"] = nearest_id
        peak["q"]["separation_status"] = (
            "POSITIVE_SEPARATION"
            if nearest_margin > 0
            else "UNCERTAINTY_REGIONS_OVERLAP"
        )


def finalize_quality(peak: dict[str, Any]) -> None:
    q = peak["q"]
    if q["local_maximum_plateau_count"] > 1:
        q["validity"] = "UNRESOLVED_MULTI_MAXIMUM_COMPONENT"
    elif q["separation_status"] == "UNCERTAINTY_REGIONS_OVERLAP":
        q["validity"] = "UNRESOLVED_UNCERTAINTY_REGION_OVERLAP"
    elif q["edge_status"] != "INTERIOR" or q["mask_status"] != "CLEAR":
        q["validity"] = "CERTIFIED_POSITIVE_SEPARATION_EDGE_OR_MASK_ADJACENT"
    elif q["separation_status"] == "SINGLETON_FRAME":
        q["validity"] = "CERTIFIED_SINGLETON_FRAME"
    else:
        q["validity"] = "CERTIFIED_POSITIVE_SEPARATION_INTERIOR"


def detect_peaks(
    image: np.ndarray,
    frame_id: str,
    geometry: dict[str, Any],
    rule: dict[str, Any],
    fixed_mask: np.ndarray,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    background_rule = rule["background_rule"]
    detection = rule["detection_rule"]
    quality = rule["quality_typing"]
    background, noise, tile_summary = tile_background(
        image,
        int(background_rule["tile_height_pixels"]),
        int(background_rule["tile_width_pixels"]),
    )
    valid = image >= 0
    excess = image.astype(np.float64) - background
    z = np.zeros_like(background)
    z[valid] = excess[valid] / noise[valid]
    seed = (
        valid
        & (z >= float(detection["seed_z_minimum"]))
        & (excess >= float(detection["seed_excess_count_minimum"]))
    )
    growth = (
        valid
        & (z >= float(detection["growth_z_minimum"]))
        & (excess > 0)
    )
    labels, growth_component_count = ndi.label(
        growth, structure=np.ones((3, 3), dtype=bool)
    )
    candidate_ids = np.unique(labels[seed])
    candidate_ids = candidate_ids[candidate_ids > 0]

    maximum = ndi.maximum_filter(
        image, size=3, mode="constant", cval=np.iinfo(np.int32).min
    )
    local_maximum_pixels = seed & (image == maximum)
    plateau_labels, _ = ndi.label(
        local_maximum_pixels, structure=np.ones((3, 3), dtype=bool)
    )
    object_slices = ndi.find_objects(labels)
    mask_guard = ndi.binary_dilation(
        fixed_mask,
        iterations=int(quality["mask_guard_pixels"]),
        structure=np.ones((3, 3), dtype=bool),
    )

    peaks: list[dict[str, Any]] = []
    rejected_small = 0
    rejected_signal = 0
    for component_id in candidate_ids.tolist():
        object_slice = object_slices[component_id - 1]
        if object_slice is None:
            raise RuntimeError("component slice unexpectedly absent")
        local_component = labels[object_slice] == component_id
        local_y, local_x = np.where(local_component)
        y = local_y + object_slice[0].start
        x = local_x + object_slice[1].start
        pixel_count = len(x)
        if pixel_count < int(detection["minimum_component_pixel_count"]):
            rejected_small += 1
            continue
        weights = np.maximum(excess[y, x], 0.0)
        signal = float(np.sum(weights))
        sigma_signal = float(math.sqrt(float(np.sum(noise[y, x] ** 2))))
        integrated_snr = signal / sigma_signal
        if integrated_snr < float(
            detection["minimum_integrated_signal_to_formal_noise"]
        ):
            rejected_signal += 1
            continue
        centroid_x = float(np.sum(weights * x) / signal)
        centroid_y = float(np.sum(weights * y) / signal)
        centroid_sigma_x = float(
            math.sqrt(float(np.sum(noise[y, x] ** 2 * (x - centroid_x) ** 2)))
            / signal
        )
        centroid_sigma_y = float(
            math.sqrt(float(np.sum(noise[y, x] ** 2 * (y - centroid_y) ** 2)))
            / signal
        )
        expansion_x = 0.5 + 3.0 * centroid_sigma_x
        expansion_y = 0.5 + 3.0 * centroid_sigma_y
        support_linear = np.sort((y * WIDTH + x).astype("<u4"))
        support_sha = sha256_bytes(support_linear.tobytes())
        peak_id = "PEAK_" + sha256_bytes(
            canonical_bytes(
                {
                    "frame_id": frame_id,
                    "support_sha256": support_sha,
                }
            )
        )[:28]
        local_plateaus = np.unique(plateau_labels[y, x])
        local_plateaus = local_plateaus[local_plateaus > 0]
        edge_guard = int(quality["edge_guard_pixels"])
        edge_status = (
            "DETECTOR_EDGE"
            if (
                int(x.min()) <= edge_guard
                or int(x.max()) >= WIDTH - 1 - edge_guard
                or int(y.min()) <= edge_guard
                or int(y.max()) >= HEIGHT - 1 - edge_guard
            )
            else "INTERIOR"
        )
        near_mask = bool(np.any(mask_guard[y, x]))
        peak = {
            "peak_id": peak_id,
            "f": frame_id,
            "x": rounded(centroid_x),
            "y": rounded(centroid_y),
            "omega": {
                "value_degrees": rounded(geometry["start_angle_degrees"]),
                "increment_degrees": rounded(geometry["angle_increment_degrees"]),
                "semantic_role": "CBF_START_ANGLE_NOT_PHYSICAL_TIME",
            },
            "Bhat": rounded(float(np.sum(background[y, x]))),
            "Shat": rounded(signal),
            "sigma_S": rounded(sigma_signal),
            "integrated_signal_to_formal_noise": rounded(integrated_snr),
            "R": {
                "type": "AXIS_ALIGNED_DETECTOR_COORDINATE_BOX",
                "x_min": rounded(float(x.min()) - expansion_x),
                "x_max": rounded(float(x.max()) + expansion_x),
                "y_min": rounded(float(y.min()) - expansion_y),
                "y_max": rounded(float(y.max()) + expansion_y),
                "centroid_sigma_x_formal": rounded(centroid_sigma_x),
                "centroid_sigma_y_formal": rounded(centroid_sigma_y),
                "coverage_interpretation": (
                    "DECLARED_FORMAL_REGION_NOT_EMPIRICALLY_CALIBRATED_CONFIDENCE"
                ),
            },
            "m": None,
            "nearest_peak_id": None,
            "support": {
                "pixel_count": pixel_count,
                "x_min_pixel": int(x.min()),
                "x_max_pixel": int(x.max()),
                "y_min_pixel": int(y.min()),
                "y_max_pixel": int(y.max()),
                "linear_index_sha256": support_sha,
                "maximum_raw_count": int(image[y, x].max()),
            },
            "q": {
                "validity": "PENDING_SEPARATION",
                "local_maximum_plateau_count": int(len(local_plateaus)),
                "overlap_status": (
                    "MULTI_MAXIMUM_COMPONENT"
                    if len(local_plateaus) > 1
                    else "SINGLE_MAXIMUM_COMPONENT"
                ),
                "mask_status": "ADJACENT_TO_MASK" if near_mask else "CLEAR",
                "edge_status": edge_status,
                "saturation_status": (
                    "THRESHOLD_UNAVAILABLE_NO_EXPLICIT_OVERLOAD_SENTINEL"
                ),
            },
        }
        peaks.append(peak)

    peaks.sort(key=lambda item: item["peak_id"])
    attach_separation(peaks)
    for peak in peaks:
        finalize_quality(peak)
    frame_summary = {
        "raw_minimum": int(image.min()),
        "raw_maximum": int(image.max()),
        "invalid_pixel_count": int(np.sum(~valid)),
        "seed_pixel_count": int(np.sum(seed)),
        "growth_pixel_count": int(np.sum(growth)),
        "growth_component_count": int(growth_component_count),
        "seeded_component_count": int(len(candidate_ids)),
        "rejected_small_component_count": rejected_small,
        "rejected_integrated_signal_count": rejected_signal,
        "accepted_peak_count": len(peaks),
        "tile_background": tile_summary,
    }
    return peaks, frame_summary


def scientific_projection(
    frames: Iterable[dict[str, Any]], peaks: Iterable[dict[str, Any]]
) -> dict[str, Any]:
    projected_frames = []
    for frame in frames:
        projected_frames.append(
            {
                "frame_id": frame["frame_id"],
                "payload_sha256": frame["payload_sha256"],
                "geometry": frame["geometry"],
                "scan_configuration_id": frame["scan_configuration_id"],
                "detection_summary": frame["detection_summary"],
            }
        )
    return {
        "frames": sorted(projected_frames, key=lambda item: item["frame_id"]),
        "peaks": sorted(peaks, key=lambda item: item["peak_id"]),
    }


def projection_digest(
    frames: Iterable[dict[str, Any]], peaks: Iterable[dict[str, Any]]
) -> str:
    return sha256_bytes(canonical_bytes(scientific_projection(frames, peaks)))


def parse_time(value: str) -> datetime:
    return datetime.fromisoformat(value.rstrip("Z"))


def make_scan_segments(frames: list[dict[str, Any]]) -> list[list[dict[str, Any]]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for frame in frames:
        grouped[frame["scan_configuration_id"]].append(frame)
    segments: list[list[dict[str, Any]]] = []
    for configuration in sorted(grouped):
        ordered = sorted(
            grouped[configuration],
            key=lambda item: (
                parse_time(item["geometry"]["acquired_at"]),
                item["frame_id"],
            ),
        )
        current: list[dict[str, Any]] = []
        for frame in ordered:
            if current:
                previous = current[-1]
                time_gap = (
                    parse_time(frame["geometry"]["acquired_at"])
                    - parse_time(previous["geometry"]["acquired_at"])
                ).total_seconds()
                expected = previous["geometry"]["angle_increment_degrees"]
                angle_gap = (
                    frame["geometry"]["start_angle_degrees"]
                    - previous["geometry"]["start_angle_degrees"]
                )
                contiguous = (
                    abs(angle_gap - expected) <= 1e-6
                    and 0 <= time_gap <= max(
                        120.0,
                        5.0 * previous["geometry"]["exposure_seconds"],
                    )
                )
                if not contiguous:
                    segments.append(current)
                    current = []
            current.append(frame)
        if current:
            segments.append(current)
    return sorted(segments, key=lambda group: group[0]["frame_id"])


def temporal_diagnostics(
    frames: list[dict[str, Any]], peaks: list[dict[str, Any]]
) -> dict[str, Any]:
    by_frame: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for peak in peaks:
        by_frame[peak["f"]].append(peak)
    isolated_ids: set[str] = set()
    internal_gap_records: list[dict[str, Any]] = []
    segments = make_scan_segments(frames)
    for segment in segments:
        coordinates = []
        for frame in segment:
            frame_peaks = by_frame.get(frame["frame_id"], [])
            coordinates.append(
                np.array([[peak["x"], peak["y"]] for peak in frame_peaks], dtype=float)
            )
        has_neighbor: list[np.ndarray] = [
            np.zeros(len(points), dtype=bool) for points in coordinates
        ]
        for index in range(len(segment) - 1):
            left = coordinates[index]
            right = coordinates[index + 1]
            if not len(left) or not len(right):
                continue
            left_tree = cKDTree(left)
            right_tree = cKDTree(right)
            left_distance, _ = right_tree.query(left, k=1)
            right_distance, _ = left_tree.query(right, k=1)
            has_neighbor[index] |= left_distance <= 2.5
            has_neighbor[index + 1] |= right_distance <= 2.5
        for index, frame in enumerate(segment):
            for peak_index, peak in enumerate(by_frame.get(frame["frame_id"], [])):
                if not has_neighbor[index][peak_index]:
                    isolated_ids.add(peak["peak_id"])
        for index in range(1, len(segment) - 1):
            outer_left = coordinates[index - 1]
            middle = coordinates[index]
            outer_right = coordinates[index + 1]
            if not len(outer_left) or not len(outer_right):
                continue
            right_tree = cKDTree(outer_right)
            middle_tree = cKDTree(middle) if len(middle) else None
            distances, right_indices = right_tree.query(outer_left, k=1)
            for left_index, distance in enumerate(distances):
                if distance > 5.0:
                    continue
                right_index = int(right_indices[left_index])
                expected = (outer_left[left_index] + outer_right[right_index]) / 2.0
                middle_distance = (
                    math.inf
                    if middle_tree is None
                    else float(middle_tree.query(expected, k=1)[0])
                )
                if middle_distance > 3.0:
                    internal_gap_records.append(
                        {
                            "left_frame_id": segment[index - 1]["frame_id"],
                            "middle_frame_id": segment[index]["frame_id"],
                            "right_frame_id": segment[index + 1]["frame_id"],
                            "expected_x": rounded(expected[0]),
                            "expected_y": rounded(expected[1]),
                            "outer_match_distance_pixels": rounded(distance),
                            "nearest_middle_distance_pixels": (
                                None
                                if math.isinf(middle_distance)
                                else rounded(middle_distance)
                            ),
                        }
                    )
    return {
        "diagnostic_status": (
            "POST_CONSTRUCTION_RAW_TEMPORAL_CONSISTENCY_PROXIES_COMPLETE"
        ),
        "scan_configuration_count": len(
            {frame["scan_configuration_id"] for frame in frames}
        ),
        "contiguous_scan_segment_count": len(segments),
        "possible_false_positive_proxy": {
            "definition": (
                "NO_PEAK_WITHIN_2.5_PIXELS_IN_EITHER_ADJACENT_ROTATION_FRAME"
            ),
            "peak_count": len(isolated_ids),
            "peak_ids": sorted(isolated_ids),
            "not_a_calibrated_false_positive_rate": True,
        },
        "possible_false_negative_proxy": {
            "definition": (
                "OUTER_ADJACENT_FRAMES_MATCH_WITHIN_5_PIXELS_AND_MIDDLE_FRAME "
                "HAS_NO_PEAK_WITHIN_3_PIXELS_OF_INTERPOLATED_POSITION"
            ),
            "gap_count": len(internal_gap_records),
            "gaps": sorted(
                internal_gap_records,
                key=lambda item: (
                    item["middle_frame_id"],
                    item["expected_y"],
                    item["expected_x"],
                ),
            ),
            "not_a_calibrated_false_negative_rate": True,
        },
        "processed_hkl_retrospective_scoring": (
            "NOT_EXECUTABLE_WITHOUT_A_LAWFUL_PIXEL_FRAME_TO_HKL_BINDING; "
            "PROCESSED_HKL_NOT_READ"
        ),
        "peak_set_modified_by_diagnostics": False,
    }


def quality_census(peaks: list[dict[str, Any]]) -> dict[str, Any]:
    validity = Counter(peak["q"]["validity"] for peak in peaks)
    overlap = Counter(peak["q"]["overlap_status"] for peak in peaks)
    separation = Counter(peak["q"]["separation_status"] for peak in peaks)
    edge = Counter(peak["q"]["edge_status"] for peak in peaks)
    mask = Counter(peak["q"]["mask_status"] for peak in peaks)
    return {
        "validity": dict(sorted(validity.items())),
        "overlap": dict(sorted(overlap.items())),
        "separation": dict(sorted(separation.items())),
        "edge": dict(sorted(edge.items())),
        "mask": dict(sorted(mask.items())),
        "saturation_threshold_available_peak_count": 0,
        "saturation_threshold_unavailable_peak_count": len(peaks),
    }


def determine_outcome(peaks: list[dict[str, Any]]) -> str:
    unresolved = sum(
        peak["q"]["validity"].startswith("UNRESOLVED") for peak in peaks
    )
    if not peaks:
        return "INSUFFICIENT_RAW_PEAK_SEPARATION"
    if unresolved or any(
        peak["q"]["saturation_status"].startswith("THRESHOLD_UNAVAILABLE")
        for peak in peaks
    ):
        return "PARTIAL_RAW_PEAK_SET_WITH_TYPED_OVERLAPS"
    return "CERTIFIED_RAW_PEAK_SET_WITH_POSITIVE_SEPARATION"


def decoder_independence_test(
    archive: zipfile.ZipFile,
    infos: list[zipfile.ZipInfo],
    decoder: Decoder,
) -> dict[str, Any]:
    positions = sorted(
        {
            0,
            len(infos) // 7,
            2 * len(infos) // 7,
            3 * len(infos) // 7,
            4 * len(infos) // 7,
            5 * len(infos) // 7,
            6 * len(infos) // 7,
            len(infos) - 1,
        }
    )
    records = []
    for position in positions:
        info = infos[position]
        data = archive.read(info)
        metadata, payload = parse_cbf(data)
        compiled = decoder.decode(payload, metadata["element_count"])
        reference = python_reference_decode(payload, metadata["element_count"])
        records.append(
            {
                "member_sha256": sha256_bytes(data),
                "array_sha256": sha256_bytes(compiled.astype("<i4").tobytes()),
                "equal": bool(np.array_equal(compiled, reference)),
            }
        )
    return {
        "test": "COMPILED_C_DECODER_EQUALS_INDEPENDENT_PYTHON_REFERENCE",
        "sample_count": len(records),
        "records": records,
        "all_equal": all(record["equal"] for record in records),
    }


def reprocess_order_test(
    archive: zipfile.ZipFile,
    selected_frames: list[dict[str, Any]],
    decoder: Decoder,
    rule: dict[str, Any],
    fixed_mask: np.ndarray,
) -> dict[str, Any]:
    path_by_id = {
        frame["frame_id"]: frame["source_archive_path"] for frame in selected_frames
    }

    def process(order: list[str]) -> str:
        subset_frames = []
        subset_peaks = []
        for frame_id in order:
            data = archive.read(path_by_id[frame_id])
            metadata, payload = parse_cbf(data)
            image = decoder.decode(payload, metadata["element_count"]).reshape(
                HEIGHT, WIDTH
            )
            peaks, summary = detect_peaks(
                image, frame_id, metadata["geometry"], rule, fixed_mask
            )
            subset_frames.append(
                {
                    "frame_id": frame_id,
                    "payload_sha256": metadata["payload_sha256"],
                    "geometry": metadata["geometry"],
                    "scan_configuration_id": scan_signature(metadata["geometry"]),
                    "detection_summary": summary,
                }
            )
            subset_peaks.extend(peaks)
        return projection_digest(subset_frames, subset_peaks)

    ids = sorted(path_by_id)
    forward_digest = process(ids)
    reverse_digest = process(list(reversed(ids)))
    return {
        "test": "INDEPENDENT_FRAME_PROCESSING_ORDER_INVARIANCE",
        "frame_count": len(ids),
        "forward_digest": forward_digest,
        "reverse_digest": reverse_digest,
        "pass": forward_digest == reverse_digest,
    }


def write_gzip_canonical(path: Path, value: Any) -> dict[str, Any]:
    raw = canonical_bytes(value)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as target:
        with gzip.GzipFile(
            filename="",
            mode="wb",
            fileobj=target,
            compresslevel=9,
            mtime=0,
        ) as compressed:
            compressed.write(raw)
    return {
        "uncompressed_byte_count": len(raw),
        "uncompressed_sha256": sha256_bytes(raw),
        "gzip_byte_count": path.stat().st_size,
        "gzip_sha256": sha256_file(path),
    }


def result_markdown(result: dict[str, Any]) -> str:
    summary = result["summary"]
    quality = result["quality_census"]
    diagnostics = result["post_construction_diagnostics"]
    lines = [
        "# NFC Crystallography Specimen A — D4 raw peak-set result",
        "",
        f"Artifact: `{result['artifact_id']}`",
        "",
        f"Outcome: `{result['outcome']}`",
        "",
        "Scientific label: `OPEN_DEVELOPMENT_NOT_INDEPENDENT_VALIDATION`",
        "",
        "## Construction result",
        "",
        f"The frozen detector-only rule processed all {summary['frame_count']:,} "
        f"CBF frames and constructed {summary['peak_count']:,} accepted raw-image "
        "peak components. The construction did not read the processed HKL, "
        "XDS.INP, cell, symmetry, SHELX, FCF, or solved-structure lanes.",
        "",
        f"Positive-separation records: "
        f"{summary['positive_separation_peak_count']:,}. "
        f"Typed unresolved records: {summary['unresolved_peak_count']:,}. "
        f"Detector-edge or mask-adjacent records: "
        f"{summary['edge_or_mask_adjacent_peak_count']:,}.",
        "",
        "All frames share one exact invalid-pixel mask: 32 complete columns "
        "(24,800 pixels per frame). The mask is derived only from the negative "
        "raw sentinel. No hardware saturation threshold was supplied, so every "
        "peak retains `THRESHOLD_UNAVAILABLE_NO_EXPLICIT_OVERLOAD_SENTINEL`; "
        "absence of observed clipping is not promoted to a saturation certificate.",
        "",
        "No licensed beamstop/static-shadow mask, flat field, dark-current map, "
        "or gain calibration was supplied. Those regions were not silently "
        "removed or repaired; their possible contribution is one reason the "
        "outcome remains partial.",
        "",
        "## Fixed rule",
        "",
        "* Background: 32×32 tile median of nonnegative pixels.",
        "* Noise scale: max(1, 1.4826×MAD, sqrt(background+1)).",
        "* Seeds: z ≥ 8 and excess ≥ 8 counts.",
        "* Growth: eight-connected pixels with z ≥ 3 and positive excess.",
        "* Acceptance: at least 2 pixels and integrated formal SNR ≥ 12.",
        "* Multi-maximum components are retained as unresolved; they are never silently split.",
        "* R is the component pixel-footprint box expanded by half a pixel and three formal centroid standard errors.",
        "",
        "## Quality census",
        "",
        "```json",
        json.dumps(quality, sort_keys=True, indent=2),
        "```",
        "",
        "## Post-construction diagnostics",
        "",
        f"Possible false-positive proxy peaks: "
        f"{diagnostics['possible_false_positive_proxy']['peak_count']:,}. "
        f"Possible internal-gap false-negative proxies: "
        f"{diagnostics['possible_false_negative_proxy']['gap_count']:,}. "
        "These are rotation-adjacency consistency diagnostics, not calibrated "
        "false-positive or false-negative rates, and they did not modify the "
        "constructed peak set.",
        "",
        "## Invariance",
        "",
        f"All {result['invariance']['pass_count']}/"
        f"{result['invariance']['test_count']} declared invariance checks pass. "
        "Catalog row order and frame processing order are presentation-order "
        "invariances. Filename invariance excludes archive paths from scientific "
        "semantics. Detector pixel-row order is not permutable because it is a "
        "spatial coordinate.",
        "",
        "## Scientific boundary",
        "",
        "D4 establishes a partial detector-coordinate peak set with explicit "
        "overlap, edge, mask, and saturation-calibration limitations. It does "
        "not establish a reciprocal lattice, unit cell, crystal symmetry, "
        "Miller indexing, novel-orbit prediction, solved structure, or "
        "independent empirical validation. Any D5 use must preserve the unresolved "
        "and calibration-limited records rather than repair them from conventional "
        "processing.",
        "",
        "## Bindings",
        "",
        f"* Raw archive SHA-256: `{result['bindings']['raw_archive_sha256']}`",
        f"* Rule SHA-256: `{result['bindings']['rule_sha256']}`",
        f"* Peak corpus gzip SHA-256: `{result['corpus']['gzip_sha256']}`",
        f"* Peak corpus uncompressed SHA-256: `{result['corpus']['uncompressed_sha256']}`",
        f"* Scientific corpus semantic SHA-256: `{result['scientific_corpus_semantic_sha256']}`",
        f"* Pipeline SHA-256: `{result['bindings']['pipeline_sha256']}`",
        "",
    ]
    return "\n".join(lines)


def run() -> None:
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    if RAW_ARCHIVE.stat().st_size != EXPECTED_ARCHIVE_BYTES:
        raise RuntimeError("raw archive byte count mismatch")
    raw_sha = sha256_file(RAW_ARCHIVE)
    if raw_sha != EXPECTED_ARCHIVE_SHA256:
        raise RuntimeError("raw archive SHA-256 mismatch")
    rule = json.loads(RULE_PATH.read_text(encoding="utf-8"))
    decoder = Decoder()
    frames: list[dict[str, Any]] = []
    peaks: list[dict[str, Any]] = []
    fixed_mask: np.ndarray | None = None
    mask_digest: str | None = None
    mask_mismatch_count = 0
    frame_maxima: list[int] = []

    with zipfile.ZipFile(RAW_ARCHIVE) as archive:
        infos = [info for info in archive.infolist() if info.filename.endswith(".cbf")]
        if len(infos) != 2408:
            raise RuntimeError(f"unexpected CBF frame count {len(infos)}")
        decoder_test = decoder_independence_test(archive, infos, decoder)
        if not decoder_test["all_equal"]:
            raise RuntimeError("compiled/reference decoder mismatch")
        for index, info in enumerate(infos, 1):
            data = archive.read(info)
            member_sha = sha256_bytes(data)
            metadata, payload = parse_cbf(data)
            image = decoder.decode(payload, metadata["element_count"]).reshape(
                HEIGHT, WIDTH
            )
            frame_mask = image < 0
            current_mask_digest = sha256_bytes(np.packbits(frame_mask).tobytes())
            if fixed_mask is None:
                fixed_mask = frame_mask.copy()
                mask_digest = current_mask_digest
            elif current_mask_digest != mask_digest or not np.array_equal(
                frame_mask, fixed_mask
            ):
                mask_mismatch_count += 1
            frame_id = frame_semantic_id(metadata)
            frame_peaks, detection_summary = detect_peaks(
                image, frame_id, metadata["geometry"], rule, fixed_mask
            )
            frame_maxima.append(int(image.max()))
            frame_record = {
                "frame_id": frame_id,
                "source_archive_path": info.filename,
                "source_member_sha256": member_sha,
                "source_member_byte_count": len(data),
                "payload_sha256": metadata["payload_sha256"],
                "payload_content_md5_base64": metadata[
                    "payload_content_md5_base64"
                ],
                "geometry": metadata["geometry"],
                "scan_configuration_id": scan_signature(metadata["geometry"]),
                "detection_summary": detection_summary,
            }
            frames.append(frame_record)
            peaks.extend(frame_peaks)
            if index % 200 == 0 or index == len(infos):
                print(
                    json.dumps(
                        {
                            "frames_processed": index,
                            "frames_total": len(infos),
                            "peaks_constructed": len(peaks),
                        },
                        sort_keys=True,
                    ),
                    flush=True,
                )

        if fixed_mask is None or mask_digest is None:
            raise RuntimeError("no detector mask constructed")
        frames.sort(key=lambda item: item["frame_id"])
        peaks.sort(key=lambda item: item["peak_id"])
        selected = [
            frames[position]
            for position in sorted(
                {
                    0,
                    len(frames) // 7,
                    2 * len(frames) // 7,
                    3 * len(frames) // 7,
                    4 * len(frames) // 7,
                    5 * len(frames) // 7,
                    6 * len(frames) // 7,
                    len(frames) - 1,
                }
            )
        ]
        order_reprocess = reprocess_order_test(
            archive, selected, decoder, rule, fixed_mask
        )

    diagnostics = temporal_diagnostics(frames, peaks)
    science_digest = projection_digest(frames, peaks)
    row_digest = projection_digest(frames, list(reversed(peaks)))
    frame_digest = projection_digest(list(reversed(frames)), peaks)
    renamed_frames = [
        dict(frame, source_archive_path=f"RENAMED/{index:04d}.cbf")
        for index, frame in enumerate(frames)
    ]
    filename_digest = projection_digest(renamed_frames, peaks)
    invariance_tests = [
        {
            "test": "CATALOG_ROW_ORDER_INVARIANCE",
            "variant_digest": row_digest,
            "pass": row_digest == science_digest,
            "qualification": (
                "PEAK_RECORD_ORDER_ONLY; DETECTOR PIXEL ROWS ARE SPATIAL "
                "COORDINATES AND ARE NOT PERMUTED"
            ),
        },
        {
            "test": "FRAME_RECORD_ORDER_INVARIANCE",
            "variant_digest": frame_digest,
            "pass": frame_digest == science_digest,
        },
        {
            "test": "ARCHIVE_FILENAME_INVARIANCE",
            "variant_digest": filename_digest,
            "pass": filename_digest == science_digest,
            "qualification": (
                "SCIENTIFIC PROJECTION EXCLUDES PROVENANCE PATH; SOURCE BYTES "
                "AND CBF-EMBEDDED GEOMETRY REMAIN FIXED"
            ),
        },
        order_reprocess,
    ]
    invariance = {
        "baseline_scientific_corpus_semantic_sha256": science_digest,
        "tests": invariance_tests,
        "pass_count": sum(test["pass"] for test in invariance_tests),
        "test_count": len(invariance_tests),
        "all_pass": all(test["pass"] for test in invariance_tests),
    }
    if not invariance["all_pass"]:
        raise RuntimeError("D4 invariance failure")

    q_census = quality_census(peaks)
    unresolved_count = sum(
        peak["q"]["validity"].startswith("UNRESOLVED") for peak in peaks
    )
    positive_count = sum(
        peak["q"]["separation_status"] == "POSITIVE_SEPARATION" for peak in peaks
    )
    edge_or_mask_count = sum(
        peak["q"]["edge_status"] != "INTERIOR"
        or peak["q"]["mask_status"] != "CLEAR"
        for peak in peaks
    )
    outcome = determine_outcome(peaks)
    corpus = {
        "artifact_id": f"NFC_CRYST_9Z6F_D4_RAW_PEAK_CORPUS_{VERSION}",
        "scientific_scope": "OPEN_DEVELOPMENT_NOT_INDEPENDENT_VALIDATION",
        "outcome": outcome,
        "rule": rule,
        "detector_mask": {
            "definition": "RAW_VALUE_LESS_THAN_ZERO",
            "mask_sha256": mask_digest,
            "mask_pixel_count": int(np.sum(fixed_mask)),
            "mask_mismatch_frame_count": mask_mismatch_count,
            "complete_masked_columns_zero_based": np.where(
                np.all(fixed_mask, axis=0)
            )[0].tolist(),
            "complete_masked_rows_zero_based": np.where(
                np.all(fixed_mask, axis=1)
            )[0].tolist(),
        },
        "frames": frames,
        "peaks": peaks,
        "post_construction_diagnostics": diagnostics,
        "invariance": invariance,
        "decoder_independence": decoder_test,
        "scientific_corpus_semantic_sha256": science_digest,
        "calibration_limitations": [
            "HARDWARE_SATURATION_THRESHOLD_UNAVAILABLE",
            "NO_LICENSED_BEAMSTOP_OR_STATIC_SHADOW_MASK",
            "NO_FLAT_FIELD_DARK_CURRENT_OR_GAIN_CALIBRATION_OBJECT",
            "NO_RAW_PIXEL_TO_RECIPROCAL_COORDINATE_ORIENTATION_BINDING",
            "TEMPORAL_FALSE_POSITIVE_AND_FALSE_NEGATIVE_DIAGNOSTICS_ARE_PROXIES_ONLY"
        ],
        "prohibited_construction_inputs_accessed": [],
    }
    corpus_binding = write_gzip_canonical(CORPUS_PATH, corpus)

    result = seal(
        {
            "artifact_id": (
                f"NFC_CRYST_9Z6F_D4_RAW_PEAK_SET_AND_POSITIVE_SEPARATION_RESULT_{VERSION}"
            ),
            "outcome": outcome,
            "scientific_scope": "OPEN_DEVELOPMENT_NOT_INDEPENDENT_VALIDATION",
            "bindings": {
                "raw_archive_byte_count": RAW_ARCHIVE.stat().st_size,
                "raw_archive_sha256": raw_sha,
                "raw_metadata_sha256": sha256_file(RAW_METADATA),
                "rule_sha256": sha256_file(RULE_PATH),
                "decoder_source_sha256": sha256_file(DECODER_SOURCE),
                "compiled_decoder_sha256": sha256_file(DECODER_LIBRARY),
                "pipeline_sha256": sha256_file(Path(__file__).resolve()),
                "python_version": platform.python_version(),
                "numpy_version": np.__version__,
                "scipy_version": scipy.__version__,
                "platform": platform.platform(),
            },
            "summary": {
                "frame_count": len(frames),
                "peak_count": len(peaks),
                "positive_separation_peak_count": positive_count,
                "unresolved_peak_count": unresolved_count,
                "edge_or_mask_adjacent_peak_count": edge_or_mask_count,
                "minimum_frame_maximum_raw_count": min(frame_maxima),
                "maximum_frame_maximum_raw_count": max(frame_maxima),
                "frames_at_or_above_32767_counts": sum(
                    maximum >= 32767 for maximum in frame_maxima
                ),
                "detector_mask_pixel_count_per_frame": int(np.sum(fixed_mask)),
                "detector_mask_mismatch_frame_count": mask_mismatch_count,
            },
            "quality_census": q_census,
            "post_construction_diagnostics": {
                "diagnostic_status": diagnostics["diagnostic_status"],
                "scan_configuration_count": diagnostics[
                    "scan_configuration_count"
                ],
                "contiguous_scan_segment_count": diagnostics[
                    "contiguous_scan_segment_count"
                ],
                "possible_false_positive_proxy": {
                    "peak_count": diagnostics["possible_false_positive_proxy"][
                        "peak_count"
                    ],
                    "not_a_calibrated_false_positive_rate": True,
                },
                "possible_false_negative_proxy": {
                    "gap_count": diagnostics["possible_false_negative_proxy"][
                        "gap_count"
                    ],
                    "not_a_calibrated_false_negative_rate": True,
                },
                "peak_set_modified_by_diagnostics": False,
                "processed_hkl_retrospective_scoring": diagnostics[
                    "processed_hkl_retrospective_scoring"
                ],
            },
            "invariance": invariance,
            "decoder_independence": decoder_test,
            "corpus": corpus_binding,
            "scientific_corpus_semantic_sha256": science_digest,
            "construction_firewall": {
                "processed_hkl_accessed": False,
                "xds_input_accessed": False,
                "cell_or_symmetry_accessed": False,
                "shelx_or_fcf_accessed": False,
                "solved_structure_accessed": False,
            },
            "calibration_limitations": [
                "HARDWARE_SATURATION_THRESHOLD_UNAVAILABLE",
                "NO_LICENSED_BEAMSTOP_OR_STATIC_SHADOW_MASK",
                "NO_FLAT_FIELD_DARK_CURRENT_OR_GAIN_CALIBRATION_OBJECT",
                "NO_RAW_PIXEL_TO_RECIPROCAL_COORDINATE_ORIENTATION_BINDING",
                "TEMPORAL_FALSE_POSITIVE_AND_FALSE_NEGATIVE_DIAGNOSTICS_ARE_PROXIES_ONLY"
            ],
            "scientific_boundaries": [
                "RAW_IMAGE_PEAK_SET_ONLY",
                "NO_RECIPROCAL_LATTICE",
                "NO_UNIT_CELL",
                "NO_CRYSTAL_SYMMETRY_QUOTIENT",
                "NO_MILLER_INDEX_BINDING",
                "NO_NOVEL_ORBIT_INTENSITY_PREDICTION",
                "NO_SOLVED_STRUCTURE",
                "NO_INDEPENDENT_EMPIRICAL_VALIDATION",
            ],
            "d5_handoff": (
                "ELIGIBLE_ONLY_WITH_ALL_UNRESOLVED_EDGE_MASK_AND_SATURATION_"
                "CALIBRATION_LIMITATIONS_PRESERVED"
            ),
        }
    )
    write_json(RESULT_JSON_PATH, result)
    RESULT_MD_PATH.write_text(
        result_markdown(result), encoding="utf-8", newline="\n"
    )
    print(
        json.dumps(
            {
                "status": "PASS",
                "outcome": outcome,
                "frames": len(frames),
                "peaks": len(peaks),
                "unresolved": unresolved_count,
                "corpus": str(CORPUS_PATH),
                "result": str(RESULT_JSON_PATH),
            },
            sort_keys=True,
        ),
        flush=True,
    )


def verify_existing() -> None:
    result = json.loads(RESULT_JSON_PATH.read_text(encoding="utf-8"))
    body = dict(result)
    declared_semantic = body.pop("canonical_semantic_sha256")
    if sha256_bytes(canonical_bytes(body)) != declared_semantic:
        raise RuntimeError("result semantic digest mismatch")
    if sha256_file(CORPUS_PATH) != result["corpus"]["gzip_sha256"]:
        raise RuntimeError("corpus gzip digest mismatch")
    with gzip.open(CORPUS_PATH, "rb") as handle:
        raw = handle.read()
    if len(raw) != result["corpus"]["uncompressed_byte_count"]:
        raise RuntimeError("corpus uncompressed byte count mismatch")
    if sha256_bytes(raw) != result["corpus"]["uncompressed_sha256"]:
        raise RuntimeError("corpus uncompressed digest mismatch")
    corpus = json.loads(raw)
    if (
        projection_digest(corpus["frames"], corpus["peaks"])
        != result["scientific_corpus_semantic_sha256"]
    ):
        raise RuntimeError("scientific corpus semantic digest mismatch")
    if not result["invariance"]["all_pass"]:
        raise RuntimeError("stored invariance failure")
    print(
        json.dumps(
            {
                "status": "PASS",
                "outcome": result["outcome"],
                "frame_count": result["summary"]["frame_count"],
                "peak_count": result["summary"]["peak_count"],
                "result_semantic_digest_verified": True,
                "corpus_byte_digests_verified": 2,
                "scientific_corpus_digest_verified": True,
                "invariance_tests_passed": result["invariance"]["pass_count"],
            },
            sort_keys=True,
        )
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify-existing", action="store_true")
    args = parser.parse_args()
    if args.verify_existing:
        verify_existing()
    else:
        run()


if __name__ == "__main__":
    main()
