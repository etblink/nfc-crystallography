#!/usr/bin/env python3
"""Build a truth-free fixed D4.5 feed from public PILATUS MiniCBF frames.

This is a narrow interface adapter for the already archived D4.5/D5
successor.  It changes no peak, aggregation, candidate-generation, or
decision rule.  dxtbx supplies the single-panel detector, beam, goniometer,
scan, and pixel array.  A byte-identical copy of the frozen NFC byte-offset
decoder independently verifies every decoded array.

The reciprocal projection deliberately uses dxtbx's SimplePxMmStrategy.
That is the literal-pixel model already fixed by the NFC 0.3.x line; the
optional detector-material parallax correction is not silently introduced.
"""

from __future__ import annotations

import argparse
import base64
import copy
import ctypes
import gzip
import hashlib
import json
import math
import re
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np


CBF_MARKER = bytes.fromhex("0c1a04d5")


def canonical_bytes(value: Any) -> bytes:
    return (
        json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n"
    ).encode("utf-8")


def semantic_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(canonical_bytes(value))


def write_gzip_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(
            filename="", mode="wb", fileobj=raw, mtime=0
        ) as compressed:
            compressed.write(canonical_bytes(value))


def one_match(pattern: str, text: str, label: str) -> re.Match[str]:
    matches = list(re.finditer(pattern, text))
    if len(matches) != 1:
        raise ValueError(f"{label}: expected one match, observed {len(matches)}")
    return matches[0]


def binary_field(header: str, label: str, cast: type = int) -> Any:
    match = one_match(
        rf"(?m)^{re.escape(label)}:\s*([^\r\n]+)",
        header,
        f"CBF binary field {label}",
    )
    return cast(match.group(1).strip().strip('"'))


def comment_number(
    header: str, label: str, cast: type = float
) -> Any:
    match = one_match(
        rf"(?m)^#\s*{re.escape(label)}\s+([-+0-9.eE]+)",
        header,
        f"PILATUS comment field {label}",
    )
    return cast(match.group(1))


def optional_comment_number(
    header: str, label: str, cast: type = float
) -> Any | None:
    matches = list(
        re.finditer(
            rf"(?m)^#\s*{re.escape(label)}\s+([-+0-9.eE]+)",
            header,
        )
    )
    if len(matches) > 1:
        raise ValueError(f"PILATUS comment field {label}: multiple values")
    return cast(matches[0].group(1)) if matches else None


def parse_byte_offset_payload(data: bytes) -> dict[str, Any]:
    """Bind and expose the exact byte-offset payload without geometry parsing."""

    marker = data.find(CBF_MARKER)
    if marker < 0:
        raise ValueError("CBF binary marker absent")
    header_bytes = data[:marker]
    header = header_bytes.decode("latin-1", errors="strict")
    convention = one_match(
        r'(?m)^_array_data\.header_convention\s+"([^"]+)"',
        header,
        "CBF header convention",
    ).group(1)
    if convention != "PILATUS_1.2":
        raise ValueError("header convention is not exactly PILATUS_1.2")
    conversion = one_match(
        r'conversions="([^"]+)"',
        header,
        "CBF conversion declaration",
    ).group(1)
    if conversion != "x-CBF_BYTE_OFFSET":
        raise ValueError("encoding is not x-CBF_BYTE_OFFSET")
    element_type = binary_field(header, "X-Binary-Element-Type", str)
    byte_order = binary_field(header, "X-Binary-Element-Byte-Order", str)
    if element_type != "signed 32-bit integer":
        raise ValueError("element type is not signed 32-bit integer")
    if byte_order != "LITTLE_ENDIAN":
        raise ValueError("element byte order is not LITTLE_ENDIAN")
    binary_size = binary_field(header, "X-Binary-Size", int)
    binary_padding = binary_field(header, "X-Binary-Size-Padding", int)
    element_count = binary_field(
        header, "X-Binary-Number-of-Elements", int
    )
    fast = binary_field(header, "X-Binary-Size-Fastest-Dimension", int)
    slow = binary_field(header, "X-Binary-Size-Second-Dimension", int)
    if fast <= 0 or slow <= 0 or element_count != fast * slow:
        raise ValueError("invalid CBF dimensions or element count")
    start = marker + len(CBF_MARKER)
    payload = data[start : start + binary_size]
    if len(payload) != binary_size:
        raise ValueError("truncated CBF byte-offset payload")
    trailer = data[start + binary_size :]
    expected = (
        b"\x00" * binary_padding
        + b"\r\n--CIF-BINARY-FORMAT-SECTION----\r\n;\r\n\r\n"
    )
    trailer_exact = trailer == expected
    stated_md5 = one_match(
        r"(?m)^Content-MD5:\s*(\S+)\s*$",
        header,
        "CBF Content-MD5",
    ).group(1)
    observed_md5 = base64.b64encode(
        hashlib.md5(payload, usedforsecurity=False).digest()
    ).decode("ascii")
    content_md5_matches = observed_md5 == stated_md5
    timestamp_matches = list(
        re.finditer(r"(?m)^#\s*(\d{4}-\d{2}-\d{2}T\S+)\s*$", header)
    )
    if len(timestamp_matches) != 1:
        raise ValueError("expected one PILATUS acquisition timestamp")
    detector_matches = list(
        re.finditer(r"(?m)^#\s*Detector:\s*(.+?)\s*$", header)
    )
    if len(detector_matches) != 1:
        raise ValueError("expected one PILATUS detector identity")
    return {
        "header": header,
        "header_sha256": sha256_bytes(header_bytes),
        "payload": payload,
        "payload_sha256": sha256_bytes(payload),
        "payload_content_md5_base64": observed_md5,
        "stated_payload_content_md5_base64": stated_md5,
        "payload_content_md5_matches": content_md5_matches,
        "trailer_exact_canonical": trailer_exact,
        "trailer_sha256": sha256_bytes(trailer),
        "header_convention": convention,
        "conversion": conversion,
        "element_type": element_type,
        "byte_order": byte_order,
        "element_count": element_count,
        "fast_dimension": fast,
        "slow_dimension": slow,
        "binary_byte_count": binary_size,
        "binary_padding_byte_count": binary_padding,
        "count_cutoff": comment_number(header, "Count_cutoff", int),
        "declared_excluded_pixel_count": optional_comment_number(
            header, "N_excluded_pixels", int
        ),
        "acquired_at": timestamp_matches[0].group(1),
        "detector": detector_matches[0].group(1),
    }


class BoundByteOffsetDecoder:
    """Compile and bind the explicitly supplied frozen C decoder source."""

    def __init__(self, source: Path, work: Path) -> None:
        self.source = source.resolve()
        work.mkdir(parents=True, exist_ok=True)
        source_sha = sha256_file(self.source)
        self.library_path = (
            work.resolve() / f"libnfc_cbf_byte_offset_{source_sha[:16]}.so"
        )
        if not self.library_path.exists():
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
                    str(self.source),
                    "-o",
                    str(self.library_path),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
        self.library = ctypes.CDLL(str(self.library_path))
        self.function = self.library.cbf_byte_offset_decode
        self.function.argtypes = [
            ctypes.POINTER(ctypes.c_uint8),
            ctypes.c_size_t,
            ctypes.POINTER(ctypes.c_int32),
            ctypes.c_size_t,
            ctypes.POINTER(ctypes.c_size_t),
        ]
        self.function.restype = ctypes.c_int

    def decode(self, payload: bytes, element_count: int) -> np.ndarray:
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
        if status != 0 or used.value != source.size:
            raise ValueError(
                "frozen byte-offset decoder failed or did not consume payload: "
                f"status={status}, used={used.value}, available={source.size}"
            )
        return output


def vector(value: Any) -> list[float]:
    return [float(item) for item in value]


def normalized(value: Any) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    norm = np.linalg.norm(array)
    if not math.isfinite(norm) or norm <= 0:
        raise ValueError("invalid zero or nonfinite vector")
    return array / norm


def scan_series_stem(path: Path) -> str:
    """Return the acquisition-series stem without its trailing image number."""

    match = re.fullmatch(r"(.+?)(\d+)(\.[^.]+)", path.name)
    if match is None:
        raise ValueError(
            f"frame name lacks a terminal numeric image ordinal: {path.name}"
        )
    return match.group(1)


def scan_series_id(
    path: Path, geometry: dict[str, Any], portability: Any
) -> str:
    """Separate same-geometry sweeps without using lattice information."""

    semantic = {
        "geometry_signature": portability.scan_signature(geometry),
        "filename_series_stem": scan_series_stem(path),
    }
    return "SCANSERIES_" + semantic_sha256(semantic)[:20]


def dxtbx_geometry(image: Any, payload: dict[str, Any]) -> dict[str, Any]:
    detector = image.get_detector()
    if len(detector) != 1:
        raise ValueError("adapter admits exactly one detector panel")
    panel = detector[0]
    beam = image.get_beam()
    goniometer = image.get_goniometer()
    scan = image.get_scan()
    if beam is None or goniometer is None or scan is None:
        raise ValueError("beam, single-axis goniometer, and scan are required")
    if tuple(panel.get_image_size()) != (
        payload["fast_dimension"],
        payload["slow_dimension"],
    ):
        raise ValueError("dxtbx image size differs from CBF binary dimensions")
    fixed = np.asarray(goniometer.get_fixed_rotation(), dtype=float).reshape(3, 3)
    setting = np.asarray(
        goniometer.get_setting_rotation(), dtype=float
    ).reshape(3, 3)
    if max(
        float(np.max(np.abs(fixed - np.eye(3)))),
        float(np.max(np.abs(setting - np.eye(3)))),
    ) > 1e-12:
        raise ValueError("nonidentity goniometer fixed/setting rotation unsupported")
    fast = normalized(panel.get_fast_axis())
    slow = normalized(panel.get_slow_axis())
    normal = normalized(np.cross(fast, slow))
    incident = normalized(beam.get_s0())
    axis = normalized(goniometer.get_rotation_axis())
    pixel = np.asarray(panel.get_pixel_size(), dtype=float) / 1000.0
    distance = float(panel.get_distance()) / 1000.0
    beam_xy = np.asarray(
        panel.get_beam_centre_px(beam.get_s0()), dtype=float
    )
    reconstructed_origin = (
        distance * normal
        - pixel[0] * beam_xy[0] * fast
        - pixel[1] * beam_xy[1] * slow
    )
    origin = np.asarray(panel.get_origin(), dtype=float) / 1000.0
    origin_error = float(np.max(np.abs(reconstructed_origin - origin)))
    if origin_error > 1e-12:
        raise ValueError(
            f"panel is not representable by literal flat-panel model: {origin_error}"
        )
    start, increment = map(float, scan.get_oscillation(deg=True))
    if increment == 0.0:
        raise ValueError("zero scan increment")
    exposure_values = list(scan.get_exposure_times())
    exposure = float(exposure_values[0]) if exposure_values else 0.0
    if any(abs(float(item) - exposure) > 1e-12 for item in exposure_values):
        raise ValueError("within-image exposure values are not constant")
    two_theta = math.degrees(
        math.acos(float(np.clip(np.dot(normal, incident), -1.0, 1.0)))
    )
    geometry = {
        "acquired_at": payload["acquired_at"],
        "exposure_seconds": exposure,
        "exposure_period_seconds": exposure,
        "wavelength_angstrom": float(beam.get_wavelength()),
        "detector_distance_m": distance,
        "beam_xy_pixels": vector(beam_xy),
        "start_angle_degrees": start,
        "angle_increment_degrees": increment,
        "detector_2theta_degrees": two_theta,
        "alpha_degrees": 0.0,
        "beta_degrees": 0.0,
        "phi_degrees": start,
        "phi_increment_degrees": increment,
        "omega_degrees": start,
        "omega_increment_degrees": increment,
        "kappa_degrees": 0.0,
        "kappa_increment_degrees": 0.0,
        "pixel_size_m": vector(pixel),
        "detector_fast_axis": vector(fast),
        "detector_slow_axis": vector(slow),
        "incident_beam": vector(incident),
        "rotation_axis": vector(axis),
        "rotation_axis_name": "DXTBX_SINGLE_AXIS",
        "geometry_derivation": (
            "DXTBX_SINGLE_PANEL_SINGLE_AXIS_WITH_LITERAL_SIMPLE_PIXEL_MODEL"
        ),
    }
    return {
        "geometry": geometry,
        "origin_reconstruction_max_abs_m": origin_error,
        "panel_material": str(panel.get_material()),
        "panel_thickness_mm": float(panel.get_thickness()),
        "dxtbx_native_pixel_strategy": type(
            panel.get_px_mm_strategy()
        ).__name__,
        "nfc_pixel_strategy": "SimplePxMmStrategy_LITERAL_PIXEL_CENTERS",
    }


def import_sources(args: argparse.Namespace) -> tuple[Any, Any, Any, Any]:
    candidate_scripts = args.candidate_scripts.resolve()
    sys.path.insert(0, str(candidate_scripts))
    import build_fixed_control_feeds as feed_builder
    import d45_successor
    import pilatus_portability as portability

    if Path(portability.__file__).resolve() != (
        candidate_scripts / "pilatus_portability.py"
    ):
        raise RuntimeError("wrong portability source imported")
    if Path(d45_successor.__file__).resolve() != (
        candidate_scripts / "d45_successor.py"
    ):
        raise RuntimeError("wrong D4.5 source imported")
    return portability, d45_successor, feed_builder, candidate_scripts


def reciprocal_reference_check(
    experiment_path: Path,
    reflection_path: Path,
    portability: Any,
    maximum_rows: int,
) -> dict[str, Any]:
    """Compare NFC q with dxtbx after fixing the same literal pixel strategy."""

    from dials.array_family import flex
    from dxtbx.model import SimplePxMmStrategy
    from dxtbx.model.experiment_list import ExperimentList

    experiments = ExperimentList.from_file(str(experiment_path))
    for experiment in experiments:
        if len(experiment.detector) != 1:
            raise ValueError("reciprocal reference requires one panel")
        experiment.detector[0].set_px_mm_strategy(SimplePxMmStrategy())
    reflections = flex.reflection_table.from_file(str(reflection_path))
    if len(reflections) > maximum_rows:
        step = max(1, len(reflections) // maximum_rows)
        selection = flex.size_t(range(0, len(reflections), step))
        reflections = reflections.select(selection)
    for column in ("rlp", "s1", "xyzobs.mm.value", "xyzobs.mm.variance"):
        if column in reflections:
            del reflections[column]
    reflections.centroid_px_to_mm(experiments)
    reflections.map_centroids_to_reciprocal_space(experiments)
    if len(experiments) != 1:
        raise ValueError("reference check currently requires one experiment")
    experiment = experiments[0]
    payload_stub = {
        "fast_dimension": experiment.detector[0].get_image_size()[0],
        "slow_dimension": experiment.detector[0].get_image_size()[1],
        "acquired_at": "REFERENCE_ONLY",
    }

    class ImageModels:
        def get_detector(self):
            return experiment.detector

        def get_beam(self):
            return experiment.beam

        def get_goniometer(self):
            return experiment.goniometer

        def get_scan(self):
            return experiment.scan

    geometry = dxtbx_geometry(ImageModels(), payload_stub)["geometry"]
    direct = np.empty((len(reflections), 3), dtype="<f8")
    reference = np.empty((len(reflections), 3), dtype="<f8")
    for ordinal in range(len(reflections)):
        x, y, z = reflections["xyzobs.px.value"][ordinal]
        q_lab = portability.detector_q_lab(
            geometry,
            np.asarray([x], dtype=float),
            np.asarray([y], dtype=float),
            0.0,
        )
        angle = (
            geometry["start_angle_degrees"]
            + float(z) * geometry["angle_increment_degrees"]
        )
        direct[ordinal] = (
            q_lab
            @ portability.FROZEN_D5.axis_rotation(
                np.asarray(geometry["rotation_axis"]), angle
            )
        )[0]
        reference[ordinal] = reflections["rlp"][ordinal]
    component = np.abs(direct - reference)
    vector_error = np.linalg.norm(direct - reference, axis=1)
    maximum = float(component.max(initial=0.0))
    if maximum > 2e-12:
        raise RuntimeError(
            f"literal reciprocal projection mismatch against dxtbx: {maximum}"
        )
    return {
        "comparison_count": len(reflections),
        "dxtbx_pixel_strategy_for_reference": "SimplePxMmStrategy",
        "crystal_model_used": False,
        "maximum_absolute_component_difference": maximum,
        "vector_difference_quantiles": {
            "p50": float(np.quantile(vector_error, 0.50)),
            "p90": float(np.quantile(vector_error, 0.90)),
            "p99": float(np.quantile(vector_error, 0.99)),
            "maximum": float(vector_error.max(initial=0.0)),
        },
        "direct_projection_float64_le_sha256": sha256_bytes(
            direct.tobytes(order="C")
        ),
        "dxtbx_projection_float64_le_sha256": sha256_bytes(
            reference.tobytes(order="C")
        ),
        "tolerance": 2e-12,
        "status": "EXACT_LITERAL_PIXEL_RECIPROCAL_PROJECTION_QUALIFIED",
    }


def run(args: argparse.Namespace) -> None:
    from dxtbx.format.Registry import get_format_class_for_file

    started = time.monotonic()
    portability, _d45, feed_builder, candidate_scripts = import_sources(args)
    paths = sorted(args.data_dir.resolve().glob(args.glob))
    if args.limit is not None:
        paths = paths[: args.limit]
    if not paths:
        raise ValueError("no input CBF frames found")
    d4_rule = json.loads(args.d4_rule.read_text(encoding="utf-8"))
    decoder = BoundByteOffsetDecoder(
        args.decoder_c_source, args.decoder_work_directory
    )

    frames: list[dict[str, Any]] = []
    peaks: list[dict[str, Any]] = []
    frame_records: list[dict[str, Any]] = []
    array_digest = hashlib.sha256()
    source_digest = hashlib.sha256()
    mask_hash_counts: Counter[str] = Counter()
    invalid_pixel_counts: list[int] = []
    native_pixel_strategies: Counter[str] = Counter()
    scan_series_counts: Counter[str] = Counter()
    maximum_origin_error = 0.0
    maximum_array_mismatch = 0
    stale_content_md5_count = 0
    noncanonical_trailer_count = 0

    for ordinal, path in enumerate(paths, 1):
        source_bytes = path.read_bytes()
        source_sha = sha256_bytes(source_bytes)
        source_digest.update(bytes.fromhex(source_sha))
        payload = parse_byte_offset_payload(source_bytes)
        stale_content_md5_count += int(
            not payload["payload_content_md5_matches"]
        )
        noncanonical_trailer_count += int(
            not payload["trailer_exact_canonical"]
        )
        format_class = get_format_class_for_file(str(path))
        if format_class is None:
            raise ValueError(f"dxtbx has no format class for {path.name}")
        image_object = format_class(str(path))
        raw = image_object.get_raw_data()
        if isinstance(raw, tuple):
            if len(raw) != 1:
                raise ValueError("multiple raw panels are unsupported")
            raw = raw[0]
        dxtbx_array = np.asarray(raw.as_numpy_array(), dtype=np.int32)
        independent_array = decoder.decode(
            payload["payload"], payload["element_count"]
        ).reshape(payload["slow_dimension"], payload["fast_dimension"])
        mismatch = int(np.sum(dxtbx_array != independent_array))
        maximum_array_mismatch = max(maximum_array_mismatch, mismatch)
        if mismatch:
            raise RuntimeError(
                f"dxtbx/frozen decoder array mismatch in {path.name}: {mismatch}"
            )
        array_sha = sha256_bytes(
            dxtbx_array.astype("<i4", copy=False).tobytes(order="C")
        )
        array_digest.update(bytes.fromhex(array_sha))
        adapter = dxtbx_geometry(image_object, payload)
        geometry = adapter["geometry"]
        native_pixel_strategies[adapter["dxtbx_native_pixel_strategy"]] += 1
        maximum_origin_error = max(
            maximum_origin_error,
            adapter["origin_reconstruction_max_abs_m"],
        )
        metadata = {
            key: value
            for key, value in payload.items()
            if key not in {"header", "payload"}
        }
        metadata["geometry"] = geometry
        frame_id = portability.frame_semantic_id(metadata)
        series_stem = scan_series_stem(path)
        series_id = scan_series_id(path, geometry, portability)
        scan_series_counts[series_stem] += 1
        invalid_mask, local_mask_census = portability.pilatus_invalid_mask(
            dxtbx_array
        )
        local_mask_sha = sha256_bytes(np.packbits(invalid_mask).tobytes())
        mask_hash_counts[local_mask_sha] += 1
        invalid_pixel_counts.append(
            int(local_mask_census["combined_invalid_pixel_count"])
        )
        local_peaks, detection = portability.detect_peaks_literal(
            dxtbx_array,
            frame_id,
            geometry,
            d4_rule,
            invalid_mask,
            payload["count_cutoff"],
        )
        frame = {
            "frame_id": frame_id,
            "source_archive_path": path.name,
            "source_file_sha256": source_sha,
            "payload_sha256": payload["payload_sha256"],
            "array_int32_le_c_order_sha256": array_sha,
            "geometry": geometry,
            "scan_configuration_id": series_id,
            "detection_summary": detection,
        }
        frames.append(frame)
        peaks.extend(local_peaks)
        frame_records.append(
            {
                "ordinal": ordinal,
                "path": path.name,
                "source_sha256": source_sha,
                "payload_sha256": payload["payload_sha256"],
                "array_int32_le_c_order_sha256": array_sha,
                "array_mismatch_count": mismatch,
                "payload_content_md5_matches": payload[
                    "payload_content_md5_matches"
                ],
                "trailer_exact_canonical": payload[
                    "trailer_exact_canonical"
                ],
                "format_class": (
                    f"{format_class.__module__}.{format_class.__name__}"
                ),
                "shape_slow_fast": list(dxtbx_array.shape),
                "d4_peak_count": len(local_peaks),
                "start_angle_degrees": geometry["start_angle_degrees"],
                "angle_increment_degrees": geometry[
                    "angle_increment_degrees"
                ],
                "scan_series_stem": series_stem,
                "scan_configuration_id": series_id,
                "origin_reconstruction_max_abs_m": adapter[
                    "origin_reconstruction_max_abs_m"
                ],
                "invalid_mask_sha256": local_mask_sha,
                "invalid_mask_census": local_mask_census,
            }
        )
        if ordinal % 25 == 0 or ordinal == len(paths):
            print(
                json.dumps(
                    {
                        "stage": "D4",
                        "frames_complete": ordinal,
                        "frames_total": len(paths),
                        "peaks_so_far": len(peaks),
                        "elapsed_seconds": round(
                            time.monotonic() - started, 1
                        ),
                    },
                    sort_keys=True,
                ),
                flush=True,
            )

    if not mask_hash_counts or not invalid_pixel_counts:
        raise RuntimeError("empty mask state")
    q, radius, ordered, tiers = portability.compute_all_q(frames, peaks)
    primary = np.asarray(
        [
            portability.FROZEN_D5.tier_for_peak(item)
            == "PRIMARY_LATTICE_TIER"
            for item in ordered
        ],
        dtype=bool,
    )
    primary_q = q[primary]
    primary_radius = radius[primary]
    primary_peaks = [
        item
        for item, keep in zip(ordered, primary, strict=True)
        if keep
    ]
    frame_map = {frame["frame_id"]: frame for frame in frames}
    groups = feed_builder.chronological_halves(frame_map)
    feeds, constructions = feed_builder.independently_constructed_splits(
        primary_q, primary_radius, primary_peaks, frame_map, groups
    )
    case = {
        "case_id": args.dataset_id,
        "case_kind": "PROSPECTIVE_PUBLIC_CORPUS",
        "role": "OUT_OF_SAMPLE_D45_D5_SUCCESSOR_TRANSFER",
        "required_outcome": "RECOVER_PRIMITIVE_LATTICE",
        "truth": {
            "basis": None,
            "role": "SEALED_UNTIL_AFTER_NFC_DECISION",
            "used_by_construction": False,
            "used_by_candidate_generation_or_decision": False,
        },
        "feeds": feeds,
        "constructions": constructions,
    }
    case["semantic_sha256"] = semantic_sha256(case)

    reference = None
    if args.reference_expt is not None or args.reference_refl is not None:
        if args.reference_expt is None or args.reference_refl is None:
            raise ValueError("both reciprocal reference inputs are required")
        reference = reciprocal_reference_check(
            args.reference_expt,
            args.reference_refl,
            portability,
            args.reference_maximum_rows,
        )

    qualification = {
        "artifact_id": (
            "NFC_CRYST_PILATUS_DXTBX_LITERAL_PIXEL_D45_FEED_"
            f"{args.dataset_id}_0_1_0"
        ),
        "dataset_id": args.dataset_id,
        "scientific_scope": (
            "INTERFACE_QUALIFICATION_AND_TRUTH_FREE_D45_FEED_CONSTRUCTION"
        ),
        "adapter_contract": {
            "single_panel": True,
            "single_axis": True,
            "goniometer_fixed_and_setting_rotations": "IDENTITY_REQUIRED",
            "scan_partition": (
                "GEOMETRY_SIGNATURE_PLUS_TRAILING_ORDINAL_FILENAME_SERIES_STEM"
            ),
            "format": "PILATUS_1.2_x-CBF_BYTE_OFFSET_SIGNED_INT32_LE",
            "pixel_model": "LITERAL_SIMPLE_PIXEL_MODEL_NO_PARALLAX_CORRECTION",
            "array_transform": "NONE",
            "resampling": False,
            "cropping": False,
            "pixel_value_scaling": False,
            "conventional_cell_or_orientation_consumed": False,
            "d4_rule_changed": False,
            "d45_rule_changed": False,
            "d5_rule_changed": False,
        },
        "source_bindings": {
            "adapter_source_sha256": sha256_file(Path(__file__)),
            "candidate_scripts_directory": str(candidate_scripts),
            "pilatus_portability_sha256": sha256_file(
                candidate_scripts / "pilatus_portability.py"
            ),
            "d45_successor_sha256": sha256_file(
                candidate_scripts / "d45_successor.py"
            ),
            "d4_rule_sha256": sha256_file(args.d4_rule),
            "frozen_byte_offset_decoder_c_sha256": sha256_file(
                args.decoder_c_source
            ),
        },
        "frame_count": len(paths),
        "source_hash_chain_sha256": source_digest.hexdigest(),
        "array_hash_chain_sha256": array_digest.hexdigest(),
        "maximum_array_mismatch_count_per_frame": maximum_array_mismatch,
        "stale_content_md5_frame_count": stale_content_md5_count,
        "noncanonical_trailer_frame_count": noncanonical_trailer_count,
        "maximum_origin_reconstruction_abs_m": maximum_origin_error,
        "dxtbx_native_pixel_strategies": dict(
            sorted(native_pixel_strategies.items())
        ),
        "scan_series_counts": dict(sorted(scan_series_counts.items())),
        "nfc_reciprocal_pixel_strategy": (
            "SimplePxMmStrategy_LITERAL_PIXEL_CENTERS"
        ),
        "invalid_mask_stable_across_frames": len(mask_hash_counts) == 1,
        "distinct_invalid_mask_count": len(mask_hash_counts),
        "invalid_mask_hash_counts": dict(sorted(mask_hash_counts.items())),
        "invalid_pixel_count_range": [
            min(invalid_pixel_counts),
            max(invalid_pixel_counts),
        ],
        "d4_peak_count": len(peaks),
        "tier_counts": tiers,
        "primary_peak_count": len(primary_peaks),
        "fixed_d45_counts": {
            label: len(feed["q_cycles_per_angstrom"])
            for label, feed in feeds.items()
        },
        "reciprocal_reference": reference,
        "frame_records": frame_records,
        "elapsed_seconds": time.monotonic() - started,
        "principal_outcome": (
            "PILATUS_DXTBX_LITERAL_PIXEL_D45_FEED_CONSTRUCTION_COMPLETE"
        ),
    }
    qualification["semantic_sha256"] = semantic_sha256(qualification)
    write_gzip_json(args.output_case, case)
    write_json(args.output_qualification, qualification)
    print(
        json.dumps(
            {
                "outcome": qualification["principal_outcome"],
                "case_sha256": sha256_file(args.output_case),
                "case_semantic_sha256": case["semantic_sha256"],
                "qualification_sha256": sha256_file(
                    args.output_qualification
                ),
                "qualification_semantic_sha256": qualification[
                    "semantic_sha256"
                ],
            },
            sort_keys=True,
        ),
        flush=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--glob", default="*.cbf")
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--candidate-scripts", type=Path, required=True)
    parser.add_argument("--d4-rule", type=Path, required=True)
    parser.add_argument("--decoder-c-source", type=Path, required=True)
    parser.add_argument(
        "--decoder-work-directory",
        type=Path,
        default=Path("adapter_dev/work"),
    )
    parser.add_argument("--output-case", type=Path, required=True)
    parser.add_argument("--output-qualification", type=Path, required=True)
    parser.add_argument("--reference-expt", type=Path)
    parser.add_argument("--reference-refl", type=Path)
    parser.add_argument("--reference-maximum-rows", type=int, default=10000)
    parser.add_argument("--limit", type=int)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
