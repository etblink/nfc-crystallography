#!/usr/bin/env python3
"""Prospective PILATUS 1.2 portability layer for NFC crystallography 0.3.x.

This module does not modify the frozen 9Z6F 0.1.0 implementation.  It imports
byte-identical copies of the frozen D4 and D5 sources and reuses their
scientific constructors where their semantics remain applicable.  The new
branch is restricted to native single-panel DECTRIS PILATUS 1.2 MiniCBF
frames, literal pixel-domain D4 rules, and a prospectively fixed dxtbx-style
single-axis geometry model.
"""

from __future__ import annotations

import base64
import ctypes
import hashlib
import importlib.util
import json
import math
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

import numpy as np


VERSION = "0_3_0"
ROOT = Path(__file__).resolve().parent
FROZEN_CORE = ROOT / "frozen_core"
WORK = ROOT.parent / "work"
CBF_MARKER = b"\x0c\x1a\x04\xd5"

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


def rounded(value: float, places: int = 9) -> float:
    result = round(float(value), places)
    return 0.0 if result == 0 else result


def rounded_vector(values: Iterable[float], places: int = 9) -> list[float]:
    return [rounded(value, places) for value in values]


def _load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load frozen module {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


FROZEN_D4 = _load_module(
    "nfc_frozen_d4_0_1_0",
    FROZEN_CORE / "d4_raw_peak_pipeline_0_1_0.py",
)
FROZEN_D5 = _load_module(
    "nfc_frozen_d5_0_1_0",
    FROZEN_CORE / "d5_reciprocal_lattice_pipeline_0_1_0.py",
)


def _one_match(pattern: str, text: str, label: str) -> re.Match[str]:
    matches = list(re.finditer(pattern, text))
    if len(matches) != 1:
        raise ValueError(f"{label} must occur exactly once")
    return matches[0]


def _comment_number(header: str, label: str, cast: type = float) -> Any:
    match = _one_match(
        rf"(?m)^#\s*{re.escape(label)}(?:\s*=\s*|\s+)([-+0-9.eE]+)",
        header,
        f"PILATUS comment field {label}",
    )
    return cast(match.group(1))


def _comment_text(header: str, label: str) -> str:
    match = _one_match(
        rf"(?m)^#\s*{re.escape(label)}\s+(.+?)\s*$",
        header,
        f"PILATUS comment field {label}",
    )
    return match.group(1)


def _binary_field(header: str, label: str, cast: type = int) -> Any:
    match = _one_match(
        rf"(?m)^{re.escape(label)}:\s*([^\r\n]+)",
        header,
        f"CBF binary field {label}",
    )
    return cast(match.group(1).strip().strip('"'))


def parse_pilatus_1_2_cbf(data: bytes) -> tuple[dict[str, Any], bytes]:
    """Parse the exact MiniCBF representation admitted by the 0.3.x branch."""

    marker = data.find(CBF_MARKER)
    if marker < 0:
        raise ValueError("CBF binary marker absent")
    header_bytes = data[:marker]
    header = header_bytes.decode("latin-1", errors="strict")

    convention_match = _one_match(
        r'(?m)^_array_data\.header_convention\s+"([^"]+)"',
        header,
        "CBF header convention",
    )
    if convention_match.group(1) != "PILATUS_1.2":
        raise ValueError("header convention is not exactly PILATUS_1.2")
    conversion_match = _one_match(
        r'conversions="([^"]+)"', header, "CBF conversion declaration"
    )
    if conversion_match.group(1) != "x-CBF_BYTE_OFFSET":
        raise ValueError("encoding is not exactly x-CBF_BYTE_OFFSET")

    element_type = _binary_field(header, "X-Binary-Element-Type", str)
    byte_order = _binary_field(header, "X-Binary-Element-Byte-Order", str)
    if element_type != "signed 32-bit integer":
        raise ValueError("element type is not signed 32-bit integer")
    if byte_order != "LITTLE_ENDIAN":
        raise ValueError("element byte order is not LITTLE_ENDIAN")

    binary_size = _binary_field(header, "X-Binary-Size", int)
    binary_padding = _binary_field(header, "X-Binary-Size-Padding", int)
    element_count = _binary_field(header, "X-Binary-Number-of-Elements", int)
    fast = _binary_field(header, "X-Binary-Size-Fastest-Dimension", int)
    slow = _binary_field(header, "X-Binary-Size-Second-Dimension", int)
    if fast <= 0 or slow <= 0 or element_count != fast * slow:
        raise ValueError("invalid CBF dimensions or element count")
    payload = data[marker + len(CBF_MARKER) : marker + len(CBF_MARKER) + binary_size]
    if len(payload) != binary_size:
        raise ValueError("truncated CBF byte-offset payload")
    if binary_padding < 0:
        raise ValueError("negative CBF binary padding")
    expected_trailer = (
        b"\x00" * binary_padding
        + b"\r\n--CIF-BINARY-FORMAT-SECTION----\r\n;\r\n\r\n"
    )
    trailer = data[marker + len(CBF_MARKER) + binary_size :]
    if trailer != expected_trailer:
        raise ValueError("CBF padding or terminal boundary mismatch")

    md5_match = _one_match(
        r"(?m)^Content-MD5:\s*(\S+)\s*$", header, "CBF Content-MD5"
    )
    observed_md5 = base64.b64encode(
        hashlib.md5(payload, usedforsecurity=False).digest()
    ).decode("ascii")
    if observed_md5 != md5_match.group(1):
        raise ValueError("CBF payload Content-MD5 mismatch")

    pixel_match = _one_match(
        r"(?m)^#\s*Pixel_size\s+([-+0-9.eE]+)\s*m\s*x\s*"
        r"([-+0-9.eE]+)\s*m",
        header,
        "PILATUS Pixel_size",
    )
    beam_match = _one_match(
        r"(?m)^#\s*Beam_xy\s+\(([-+0-9.eE]+),\s*([-+0-9.eE]+)\)"
        r"\s+pixels",
        header,
        "PILATUS Beam_xy",
    )
    timestamp_match = _one_match(
        r"(?m)^#\s*(\d{4}-\d{2}-\d{2}T\S+)\s*$",
        header,
        "PILATUS acquisition timestamp",
    )
    detector_match = _one_match(
        r"(?m)^#\s*Detector:\s*(.+?)\s*$", header, "PILATUS detector identity"
    )

    oscillation_axis = _comment_text(header, "Oscillation_axis").strip().upper()
    if oscillation_axis != "PHI":
        raise ValueError("only the prospectively bound PHI rotation is admitted")

    pixel_size = [float(pixel_match.group(1)), float(pixel_match.group(2))]
    count_cutoff = _comment_number(header, "Count_cutoff", int)
    excluded_pixels = _comment_number(header, "N_excluded_pixels", int)
    geometry = {
        "acquired_at": timestamp_match.group(1),
        "exposure_seconds": _comment_number(header, "Exposure_time"),
        "exposure_period_seconds": _comment_number(header, "Exposure_period"),
        "wavelength_angstrom": _comment_number(header, "Wavelength"),
        "detector_distance_m": _comment_number(header, "Detector_distance"),
        "beam_xy_pixels": [
            float(beam_match.group(1)),
            float(beam_match.group(2)),
        ],
        "start_angle_degrees": _comment_number(header, "Start_angle"),
        "angle_increment_degrees": _comment_number(header, "Angle_increment"),
        "detector_2theta_degrees": _comment_number(header, "Detector_2theta"),
        "alpha_degrees": _comment_number(header, "Alpha"),
        "beta_degrees": 0.0,
        "phi_degrees": _comment_number(header, "Phi"),
        "phi_increment_degrees": _comment_number(header, "Phi_increment"),
        "omega_degrees": _comment_number(header, "Omega"),
        "omega_increment_degrees": _comment_number(header, "Omega_increment"),
        "kappa_degrees": _comment_number(header, "Kappa"),
        "kappa_increment_degrees": 0.0,
        "pixel_size_m": pixel_size,
        # Frozen prospectively from the public dxtbx MiniCBF model:
        # detector factory directions "+x", "-y"; simple beam and goniometer.
        "detector_fast_axis": [1.0, 0.0, 0.0],
        "detector_slow_axis": [0.0, -1.0, 0.0],
        "incident_beam": [0.0, 0.0, -1.0],
        "rotation_axis": [1.0, 0.0, 0.0],
        "rotation_axis_name": "PHI",
        "geometry_derivation": "DXTBX_MINICBF_SIMPLE_MODEL_PROSPECTIVELY_BOUND",
    }
    metadata = {
        "header_convention": convention_match.group(1),
        "detector": detector_match.group(1),
        "binary_byte_count": binary_size,
        "binary_padding_byte_count": binary_padding,
        "element_count": element_count,
        "fast_dimension": fast,
        "slow_dimension": slow,
        "element_type": element_type,
        "byte_order": byte_order,
        "conversion": conversion_match.group(1),
        "payload_sha256": sha256_bytes(payload),
        "payload_content_md5_base64": observed_md5,
        "header_sha256": sha256_bytes(header_bytes),
        "count_cutoff": count_cutoff,
        "declared_excluded_pixel_count": excluded_pixels,
        "geometry": geometry,
    }
    return metadata, payload


class ByteOffsetDecoder:
    """ctypes binding to the byte-identical frozen 0.1.0 C decoder."""

    def __init__(self) -> None:
        source = FROZEN_CORE / "d4_cbf_byte_offset_0_1_0.c"
        library_path = WORK / "libnfc_cbf_byte_offset_0_3_0.so"
        WORK.mkdir(parents=True, exist_ok=True)
        if (
            not library_path.exists()
            or library_path.stat().st_mtime_ns < source.stat().st_mtime_ns
        ):
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
                    str(source),
                    "-o",
                    str(library_path),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
        self.library_path = library_path
        self.library = ctypes.CDLL(str(library_path))
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
        used = ctypes.c_size_t()
        status = self.function(
            source.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8)),
            source.size,
            output.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
            output.size,
            ctypes.byref(used),
        )
        if status != 0:
            raise ValueError(
                f"CBF byte-offset decode failed: status={status}, "
                f"used={used.value}, available={source.size}"
            )
        if used.value != source.size:
            raise ValueError(
                "CBF byte-offset payload has trailing undecoded bytes: "
                f"used={used.value}, available={source.size}"
            )
        return output


def python_reference_decode(payload: bytes, element_count: int) -> np.ndarray:
    """Independent, deliberately simple reference decoder for sampled frames."""

    output = np.empty(element_count, dtype=np.int32)
    position = 0
    previous = 0
    for index in range(element_count):
        if position >= len(payload):
            raise ValueError("truncated byte-offset stream")
        delta = int.from_bytes(payload[position : position + 1], "little", signed=True)
        position += 1
        if delta == -128:
            if position + 2 > len(payload):
                raise ValueError("truncated int16 byte-offset escape")
            delta = int.from_bytes(
                payload[position : position + 2], "little", signed=True
            )
            position += 2
            if delta == -32768:
                if position + 4 > len(payload):
                    raise ValueError("truncated int32 byte-offset escape")
                delta = int.from_bytes(
                    payload[position : position + 4], "little", signed=True
                )
                position += 4
                if delta == -2147483648:
                    if position + 8 > len(payload):
                        raise ValueError("truncated int64 byte-offset escape")
                    delta = int.from_bytes(
                        payload[position : position + 8], "little", signed=True
                    )
                    position += 8
        previous += delta
        if previous < -(2**31) or previous > 2**31 - 1:
            raise ValueError("decoded pixel leaves signed-int32 range")
        output[index] = previous
    if position != len(payload):
        raise ValueError("trailing bytes in byte-offset stream")
    return output


def decode_pilatus_file(
    path: Path, decoder: ByteOffsetDecoder
) -> tuple[dict[str, Any], np.ndarray, bytes]:
    data = path.read_bytes()
    metadata, payload = parse_pilatus_1_2_cbf(data)
    flat = decoder.decode(payload, metadata["element_count"])
    image = flat.reshape(
        (metadata["slow_dimension"], metadata["fast_dimension"]), order="C"
    )
    if image.dtype != np.dtype("int32"):
        raise RuntimeError("decoded array is not native signed int32")
    return metadata, image, payload


def pilatus_module_gap_mask(shape: tuple[int, int]) -> np.ndarray:
    """Return the fixed PILATUS module-gap mask for an admitted panel shape.

    PILATUS modules contain 487 fast by 195 slow pixels.  Adjacent modules are
    separated by 7 fast-axis or 17 slow-axis array positions.  The authorized
    1M development panel is therefore 2 x 5 modules with shape 981 x 1043.
    """

    slow, fast = shape
    module_fast, gap_fast = 487, 7
    module_slow, gap_slow = 195, 17
    n_fast, fast_remainder = divmod(fast, module_fast)
    n_slow, slow_remainder = divmod(slow, module_slow)
    if (
        n_fast < 1
        or n_slow < 1
        or fast_remainder != (n_fast - 1) * gap_fast
        or slow_remainder != (n_slow - 1) * gap_slow
    ):
        raise ValueError("array shape is not a supported PILATUS module tiling")
    mask = np.zeros(shape, dtype=bool)
    for module in range(1, n_fast):
        start = module * module_fast + (module - 1) * gap_fast
        mask[:, start : start + gap_fast] = True
    for module in range(1, n_slow):
        start = module * module_slow + (module - 1) * gap_slow
        mask[start : start + gap_slow, :] = True
    return mask


def pilatus_invalid_mask(image: np.ndarray) -> tuple[np.ndarray, dict[str, int]]:
    gap_mask = pilatus_module_gap_mask(tuple(image.shape))
    negative_mask = image < 0
    combined = gap_mask | negative_mask
    return combined, {
        "module_gap_pixel_count": int(np.sum(gap_mask)),
        "negative_sentinel_pixel_count": int(np.sum(negative_mask)),
        "combined_invalid_pixel_count": int(np.sum(combined)),
        "gap_and_negative_overlap_count": int(np.sum(gap_mask & negative_mask)),
    }


def frame_semantic_id(metadata: dict[str, Any]) -> str:
    semantic = {
        "payload_sha256": metadata["payload_sha256"],
        "geometry": metadata["geometry"],
        "dimensions": [
            metadata["slow_dimension"],
            metadata["fast_dimension"],
        ],
        "count_cutoff": metadata["count_cutoff"],
    }
    return "FRAME_" + sha256_bytes(canonical_bytes(semantic))[:24]


def scan_signature(geometry: dict[str, Any]) -> str:
    semantic = {
        "rotation_axis_name": geometry["rotation_axis_name"],
        "rotation_axis": geometry["rotation_axis"],
        "angle_increment_degrees": geometry["angle_increment_degrees"],
        "detector_2theta_degrees": geometry["detector_2theta_degrees"],
        "detector_distance_m": geometry["detector_distance_m"],
        "beam_xy_pixels": geometry["beam_xy_pixels"],
        "pixel_size_m": geometry["pixel_size_m"],
        "wavelength_angstrom": geometry["wavelength_angstrom"],
        "exposure_seconds": geometry["exposure_seconds"],
        "detector_fast_axis": geometry["detector_fast_axis"],
        "detector_slow_axis": geometry["detector_slow_axis"],
        "incident_beam": geometry["incident_beam"],
    }
    return "SCANCFG_" + sha256_bytes(canonical_bytes(semantic))[:20]


def detect_peaks_literal(
    image: np.ndarray,
    frame_id: str,
    geometry: dict[str, Any],
    rule: dict[str, Any],
    fixed_mask: np.ndarray,
    count_cutoff: int,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Apply frozen D4 operators with native dimensions and literal pixel rules."""

    if image.ndim != 2 or image.dtype != np.dtype("int32"):
        raise ValueError("D4 requires one signed-int32 two-dimensional panel")
    if fixed_mask.shape != image.shape:
        raise ValueError("invalid mask shape differs from pixel-array shape")
    detection_image = image.copy()
    detection_image[fixed_mask] = -1
    old_height = FROZEN_D4.HEIGHT
    old_width = FROZEN_D4.WIDTH
    FROZEN_D4.HEIGHT, FROZEN_D4.WIDTH = detection_image.shape
    try:
        peaks, summary = FROZEN_D4.detect_peaks(
            detection_image, frame_id, geometry, rule, fixed_mask
        )
    finally:
        FROZEN_D4.HEIGHT = old_height
        FROZEN_D4.WIDTH = old_width
    overloaded_peak_count = 0
    for peak in peaks:
        if peak["support"]["maximum_raw_count"] >= count_cutoff:
            peak["q"]["saturation_status"] = (
                "EXPLICIT_HEADER_COUNT_CUTOFF_REACHED"
            )
            overloaded_peak_count += 1
        else:
            peak["q"]["saturation_status"] = (
                "EXPLICIT_HEADER_COUNT_CUTOFF_NOT_REACHED"
            )
    summary["count_cutoff"] = count_cutoff
    summary["overloaded_peak_count"] = overloaded_peak_count
    summary["raw_minimum"] = int(image.min())
    summary["raw_maximum"] = int(image.max())
    summary["invalid_pixel_count"] = int(np.sum(fixed_mask))
    summary["invalid_mask_application"] = (
        "RAW_NEGATIVE_SENTINELS_UNION_FIXED_PILATUS_MODULE_GAPS"
    )
    summary["literal_pixel_rules"] = True
    return peaks, summary


def detector_q_lab(
    geometry: dict[str, Any],
    x: np.ndarray,
    y: np.ndarray,
    pixel_offset: float,
) -> np.ndarray:
    fast = np.asarray(geometry["detector_fast_axis"], dtype=float)
    slow = np.asarray(geometry["detector_slow_axis"], dtype=float)
    normal = np.cross(fast, slow)
    normal /= np.linalg.norm(normal)
    pixel_x, pixel_y = map(float, geometry["pixel_size_m"])
    point = (
        geometry["detector_distance_m"] * normal[None, :]
        + pixel_x
        * (np.asarray(x) + pixel_offset - geometry["beam_xy_pixels"][0])[:, None]
        * fast[None, :]
        + pixel_y
        * (np.asarray(y) + pixel_offset - geometry["beam_xy_pixels"][1])[:, None]
        * slow[None, :]
    )
    outgoing = point / np.linalg.norm(point, axis=1)[:, None]
    incident = np.asarray(geometry["incident_beam"], dtype=float)
    incident /= np.linalg.norm(incident)
    return (outgoing - incident[None, :]) / geometry["wavelength_angstrom"]


def q_for_frame_points(
    frame: dict[str, Any],
    x: np.ndarray,
    y: np.ndarray,
    rotation_degrees: float,
    pixel_offset: float = 0.5,
) -> np.ndarray:
    q_lab = detector_q_lab(frame["geometry"], x, y, pixel_offset)
    orientation = FROZEN_D5.axis_rotation(
        np.asarray(frame["geometry"]["rotation_axis"], dtype=float),
        rotation_degrees,
    )
    return q_lab @ orientation


def compute_all_q(
    frames_list: list[dict[str, Any]],
    peaks_input: list[dict[str, Any]],
) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]], dict[str, int]]:
    frames = {frame["frame_id"]: frame for frame in frames_list}
    peaks = sorted(peaks_input, key=lambda item: item["peak_id"])
    positions: dict[str, list[int]] = defaultdict(list)
    for index, peak in enumerate(peaks):
        positions[peak["f"]].append(index)
    q = np.empty((len(peaks), 3), dtype=float)
    radius = np.empty(len(peaks), dtype=float)
    tier_counts: Counter[str] = Counter()
    for frame_id, indices in positions.items():
        frame = frames[frame_id]
        selected = [peaks[index] for index in indices]
        x = np.array([peak["x"] for peak in selected], dtype=float)
        y = np.array([peak["y"] for peak in selected], dtype=float)
        start = frame["geometry"]["start_angle_degrees"]
        increment = frame["geometry"]["angle_increment_degrees"]
        center = q_for_frame_points(frame, x, y, start + 0.5 * increment)
        q[indices] = center
        x_options = [
            np.array([peak["R"]["x_min"] for peak in selected]),
            np.array([peak["R"]["x_max"] for peak in selected]),
        ]
        y_options = [
            np.array([peak["R"]["y_min"] for peak in selected]),
            np.array([peak["R"]["y_max"] for peak in selected]),
        ]
        local_radius = np.zeros(len(selected), dtype=float)
        for corner_x in x_options:
            for corner_y in y_options:
                for rotation in (start, start + increment):
                    corner = q_for_frame_points(
                        frame, corner_x, corner_y, rotation
                    )
                    local_radius = np.maximum(
                        local_radius, np.linalg.norm(corner - center, axis=1)
                    )
        radius[indices] = local_radius
        tier_counts.update(FROZEN_D5.tier_for_peak(peak) for peak in selected)
    return q, radius, peaks, dict(sorted(tier_counts.items()))


def lattice_diagnostic(
    frames_list: list[dict[str, Any]],
    peaks: list[dict[str, Any]],
    d5_rule: dict[str, Any],
) -> dict[str, Any]:
    """Run the unchanged D5 lattice search after the new geometry projection."""

    q, radius, ordered_peaks, tier_counts = compute_all_q(frames_list, peaks)
    frames = {frame["frame_id"]: frame for frame in frames_list}
    primary_indices = np.asarray(
        [
            index
            for index, peak in enumerate(ordered_peaks)
            if FROZEN_D5.tier_for_peak(peak) == "PRIMARY_LATTICE_TIER"
        ],
        dtype=np.int64,
    )
    primary_q = q[primary_indices]
    primary_peaks = [ordered_peaks[index] for index in primary_indices]
    if len(primary_peaks) < 1000:
        raise RuntimeError("insufficient primary peaks for the frozen D5 search")
    seed, seed_record = FROZEN_D5.seed_lattice(
        primary_q, primary_peaks, frames, d5_rule
    )
    basis, refinement = FROZEN_D5.refine_lattice(primary_q, seed, d5_rule)
    direct = np.linalg.inv(basis).T
    metric = FROZEN_D5.cell_metric(direct, basis)
    h, residual, compatible, support = FROZEN_D5.support_census(
        q,
        radius,
        ordered_peaks,
        basis,
        d5_rule["lattice_search"][
            "center_residual_thresholds_cycles_per_angstrom"
        ],
    )
    projection = [
        {
            "peak_id": peak["peak_id"],
            "q": rounded_vector(q[index]),
            "h": h[index].tolist(),
            "residual": rounded(residual[index]),
            "radius": rounded(radius[index]),
            "compatible": bool(compatible[index]),
        }
        for index, peak in enumerate(ordered_peaks)
    ]
    return {
        "status": "NONCONFIRMATORY_DEVELOPMENT_LATTICE_DIAGNOSTIC_COMPLETE",
        "peak_count": len(ordered_peaks),
        "primary_peak_count": len(primary_peaks),
        "tier_counts": tier_counts,
        "single_axis_geometry": {
            "rotation_axis": [1.0, 0.0, 0.0],
            "axis_name": "PHI",
            "pixel_center_offset": 0.5,
            "frame_fraction": 0.5,
            "selection_basis": "PROSPECTIVE_DXTBX_MINICBF_CONVENTION",
            "data_driven_convention_selection": False,
        },
        "seed_search": seed_record,
        "refinement": refinement,
        "metric": metric,
        "support": support,
        "q_projection_sha256": sha256_bytes(canonical_bytes(projection)),
        "claims_not_established": [
            "PHYSICAL_PRIMITIVE_UNIT_CELL",
            "CRYSTAL_SYSTEM",
            "SPACE_GROUP",
            "PROCESSED_MILLER_INDEX_BINDING",
            "D6C1_BRIDGE",
            "D6C2_INTENSITY_SCREEN",
            "D7_PREDICTION",
            "CONFIRMATION",
        ],
    }
