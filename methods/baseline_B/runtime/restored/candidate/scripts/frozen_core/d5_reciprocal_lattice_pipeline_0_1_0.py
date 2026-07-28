#!/usr/bin/env python3
"""D5 detector-to-reciprocal binding and raw-lattice construction for 9Z6F.

The construction lane is deliberately incapable of reading processed HKL,
XDS input/output, known unit-cell or orientation data, SHELX/FCF/CIF records,
or a solved structure.  It consumes only the frozen D4 peak corpus, CBF
geometry headers from the frozen raw archive, and the frozen D5 rule.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import itertools
import json
import math
import re
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import scipy
from scipy.spatial import cKDTree


VERSION = "0_1_0"
ROOT = Path(__file__).resolve().parent
ARTIFACTS = ROOT / "artifacts"
RAW_ARCHIVE = ROOT / "source" / "raw" / "Archive.zip"
D4_CORPUS = ARTIFACTS / f"NFC_CRYST_9Z6F_D4_RAW_PEAK_CORPUS_{VERSION}.json.gz"
D4_RESULT = (
    ARTIFACTS
    / f"NFC_CRYST_9Z6F_D4_RAW_PEAK_SET_AND_POSITIVE_SEPARATION_RESULT_{VERSION}.json"
)
RULE_PATH = ROOT / "d5_reciprocal_rule_freeze.json"

D5A_PATH = (
    ARTIFACTS
    / f"NFC_CRYST_9Z6F_D5A_DETECTOR_TO_RECIPROCAL_BINDING_{VERSION}.json"
)
D5B_PATH = (
    ARTIFACTS
    / f"NFC_CRYST_9Z6F_D5B_RAW_RECIPROCAL_LATTICE_CANDIDATE_{VERSION}.json"
)
Q_CORPUS_PATH = (
    ARTIFACTS / f"NFC_CRYST_9Z6F_D5_RECIPROCAL_COORDINATE_CORPUS_{VERSION}.json.gz"
)
RESULT_JSON_PATH = (
    ARTIFACTS / f"NFC_CRYST_9Z6F_D5_RAW_RECIPROCAL_LATTICE_RESULT_{VERSION}.json"
)
RESULT_MD_PATH = (
    ARTIFACTS / f"NFC_CRYST_9Z6F_D5_RAW_RECIPROCAL_LATTICE_RESULT_{VERSION}.md"
)

EXPECTED_ARCHIVE_BYTES = 817_839_975
EXPECTED_ARCHIVE_SHA256 = (
    "fc639d3c87bb0d5c9002c78a17d0cbe85cd4879e638a3236a49ae1a95f48ecc1"
)
EXPECTED_D4_GZIP_SHA256 = (
    "94e38a7236d5e769ca4a2e3a3952cda9e06c1cff6be349053532dc6be28936e7"
)
EXPECTED_D4_UNCOMPRESSED_SHA256 = (
    "86c0ba5d0269249185788a9e6dc563bca9feb6ae83a9386b376a51fd05287480"
)
EXPECTED_D4_SCIENCE_SHA256 = (
    "191e9f59ac7c1eea9776c2b1d0f3581d40be176751ba42d94a69c4d3c5876641"
)

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


def rounded_matrix(values: np.ndarray, places: int = 9) -> list[list[float]]:
    return [rounded_vector(row, places) for row in np.asarray(values)]


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(canonical_bytes(value))


def seal(value: dict[str, Any]) -> dict[str, Any]:
    body = dict(value)
    body.pop("canonical_semantic_sha256", None)
    value["canonical_semantic_sha256"] = sha256_bytes(canonical_bytes(body))
    return value


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


def load_d4() -> tuple[dict[str, Any], bytes]:
    if sha256_file(D4_CORPUS) != EXPECTED_D4_GZIP_SHA256:
        raise RuntimeError("D4 corpus gzip digest mismatch")
    with gzip.open(D4_CORPUS, "rb") as handle:
        raw = handle.read()
    if sha256_bytes(raw) != EXPECTED_D4_UNCOMPRESSED_SHA256:
        raise RuntimeError("D4 corpus uncompressed digest mismatch")
    d4 = json.loads(raw)
    if d4["scientific_corpus_semantic_sha256"] != EXPECTED_D4_SCIENCE_SHA256:
        raise RuntimeError("D4 scientific projection digest mismatch")
    return d4, raw


def header_float(header: str, label: str) -> float:
    match = re.search(rf"(?m)^#\s*{re.escape(label)}\s+([-+0-9.eE]+)", header)
    if not match:
        raise ValueError(f"missing CBF header value {label}")
    return float(match.group(1))


def header_vector(header: str, label: str) -> np.ndarray:
    match = re.search(rf"(?m)^#\s*{re.escape(label)}\s+(.+?)\s*$", header)
    if not match:
        raise ValueError(f"missing CBF header vector {label}")
    values = np.array([float(item) for item in match.group(1).split()], dtype=float)
    if values.shape != (3,):
        raise ValueError(f"invalid CBF vector {label}")
    return values


def parse_cbf_geometry_header(header_bytes: bytes) -> dict[str, Any]:
    header = header_bytes.decode("latin-1", errors="strict")
    convention_match = re.search(
        r'(?m)^_array_data\.header_convention\s+"([^"]+)"', header
    )
    pixel_match = re.search(
        r"(?m)^#\s*Pixel_size\s+([-+0-9.eE]+)\s*m\s*x\s*"
        r"([-+0-9.eE]+)\s*m",
        header,
    )
    beam_match = re.search(
        r"(?m)^#\s*Beam_xy\s+\(([-+0-9.eE]+),\s*([-+0-9.eE]+)\)"
        r"\s+pixels",
        header,
    )
    if not convention_match or not pixel_match or not beam_match:
        raise ValueError("CBF geometry convention, pixel size, or Beam_xy absent")
    disclaimer = (
        "the miniCBF header values are based on right-handed rotations"
        in header
        and "CAP ESPERANTO SECTION" in header
    )
    return {
        "header_convention": convention_match.group(1),
        "pixel_size_m": [float(pixel_match.group(1)), float(pixel_match.group(2))],
        "beam_xy_pixels": [
            float(beam_match.group(1)),
            float(beam_match.group(2)),
        ],
        "detector_distance_m": header_float(header, "Detector_distance"),
        "wavelength_angstrom": header_float(header, "Wavelength"),
        "start_angle_degrees": header_float(header, "Start_angle"),
        "angle_increment_degrees": header_float(header, "Angle_increment"),
        "detector_2theta_degrees": header_float(header, "Detector_2theta"),
        "omega_degrees": header_float(header, "Omega"),
        "kappa_degrees": header_float(header, "Kappa"),
        "phi_degrees": header_float(header, "Phi"),
        "detector_fast_axis": header_vector(
            header, "Detector_fast_axis_vector"
        ).tolist(),
        "detector_slow_axis": header_vector(
            header, "Detector_slow_axis_vector"
        ).tolist(),
        "incident_beam": header_vector(header, "Incident_beam_vector").tolist(),
        "omega_axis": header_vector(header, "Omega_axis_vector").tolist(),
        "kappa_axis": header_vector(header, "Kappa_axis_vector").tolist(),
        "phi_axis": header_vector(header, "Phi_axis_vector").tolist(),
        "two_theta_axis": header_vector(header, "2Theta_axis_vector").tolist(),
        "rotation_axis": header_vector(header, "Rotation_axis_vector").tolist(),
        "right_handed_surrogate_disclaimer_present": disclaimer,
        "cap_esperanto_payload_present": (
            header.count("CAP ESPERANTO SECTION") > 1
        ),
    }


def close_list(left: Iterable[float], right: Iterable[float], atol: float = 1e-8) -> bool:
    return bool(np.allclose(list(left), list(right), rtol=0.0, atol=atol))


def acquire_header_geometry(
    d4_frames: list[dict[str, Any]],
) -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
    by_path = {frame["source_archive_path"]: frame for frame in d4_frames}
    if len(by_path) != 2408:
        raise RuntimeError("D4 frame-path universe is not exactly 2,408")
    headers: dict[str, dict[str, Any]] = {}
    with zipfile.ZipFile(RAW_ARCHIVE) as archive:
        cbf_infos = [
            info for info in archive.infolist() if info.filename.endswith(".cbf")
        ]
        if len(cbf_infos) != 2408:
            raise RuntimeError("raw CBF frame count mismatch")
        for info in cbf_infos:
            with archive.open(info) as handle:
                prefix = handle.read(8192)
            geometry = parse_cbf_geometry_header(prefix)
            d4_geometry = by_path[info.filename]["geometry"]
            scalar_pairs = [
                ("detector_distance_m", "detector_distance_m"),
                ("wavelength_angstrom", "wavelength_angstrom"),
                ("start_angle_degrees", "start_angle_degrees"),
                ("angle_increment_degrees", "angle_increment_degrees"),
                ("detector_2theta_degrees", "detector_2theta_degrees"),
                ("omega_degrees", "omega_degrees"),
                ("kappa_degrees", "kappa_degrees"),
                ("phi_degrees", "phi_degrees"),
            ]
            for header_key, d4_key in scalar_pairs:
                if not math.isclose(
                    geometry[header_key],
                    d4_geometry[d4_key],
                    rel_tol=0.0,
                    abs_tol=1e-8,
                ):
                    raise RuntimeError(
                        f"CBF/D4 geometry mismatch: {info.filename} {header_key}"
                    )
            for key in (
                "beam_xy_pixels",
                "detector_fast_axis",
                "detector_slow_axis",
                "incident_beam",
            ):
                if not close_list(geometry[key], d4_geometry[key], 1e-8):
                    raise RuntimeError(
                        f"CBF/D4 geometry mismatch: {info.filename} {key}"
                    )
            headers[info.filename] = geometry
    axes = {
        key: sorted(
            {
                tuple(rounded_vector(geometry[key], 8))
                for geometry in headers.values()
            }
        )
        for key in (
            "omega_axis",
            "kappa_axis",
            "phi_axis",
            "two_theta_axis",
            "rotation_axis",
        )
    }
    census = {
        "frame_count": len(headers),
        "header_conventions": sorted(
            {geometry["header_convention"] for geometry in headers.values()}
        ),
        "pixel_sizes_m": sorted(
            {tuple(geometry["pixel_size_m"]) for geometry in headers.values()}
        ),
        "right_handed_surrogate_disclaimer_frame_count": sum(
            geometry["right_handed_surrogate_disclaimer_present"]
            for geometry in headers.values()
        ),
        "cap_esperanto_payload_frame_count": sum(
            geometry["cap_esperanto_payload_present"]
            for geometry in headers.values()
        ),
        "axis_universes": axes,
    }
    if census["header_conventions"] != ["RIGAKU_1.3"]:
        raise RuntimeError("unexpected CBF header convention")
    if census["pixel_sizes_m"] != [(0.0001, 0.0001)]:
        raise RuntimeError("unexpected detector pixel size")
    if census["right_handed_surrogate_disclaimer_frame_count"] != 2408:
        raise RuntimeError("right-handed miniCBF disclaimer is not universal")
    if any(len(values) != 1 for values in axes.values()):
        raise RuntimeError("goniometer-axis vectors vary across the raw headers")
    return headers, census


def axis_rotation(axis: np.ndarray, degrees: float) -> np.ndarray:
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    theta = math.radians(float(degrees))
    cross = np.array(
        [
            [0.0, -axis[2], axis[1]],
            [axis[2], 0.0, -axis[0]],
            [-axis[1], axis[0], 0.0],
        ]
    )
    return (
        np.eye(3)
        + math.sin(theta) * cross
        + (1.0 - math.cos(theta)) * (cross @ cross)
    )


def batch_axis_rotation(axis: np.ndarray, degrees: np.ndarray) -> np.ndarray:
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    theta = np.deg2rad(np.asarray(degrees, dtype=float))
    cosine = np.cos(theta)
    sine = np.sin(theta)
    complement = 1.0 - cosine
    x, y, z = axis
    result = np.empty((len(theta), 3, 3), dtype=float)
    result[:, 0, 0] = cosine + x * x * complement
    result[:, 0, 1] = x * y * complement - z * sine
    result[:, 0, 2] = x * z * complement + y * sine
    result[:, 1, 0] = y * x * complement + z * sine
    result[:, 1, 1] = cosine + y * y * complement
    result[:, 1, 2] = y * z * complement - x * sine
    result[:, 2, 0] = z * x * complement - y * sine
    result[:, 2, 1] = z * y * complement + x * sine
    result[:, 2, 2] = cosine + z * z * complement
    return result


def detector_q_lab(
    geometry: dict[str, Any],
    x: np.ndarray,
    y: np.ndarray,
    pixel_offset: float,
) -> np.ndarray:
    fast = np.asarray(geometry["detector_fast_axis"], dtype=float)
    slow = np.asarray(geometry["detector_slow_axis"], dtype=float)
    normal = np.cross(fast, slow)
    normal = normal / np.linalg.norm(normal)
    point = (
        geometry["detector_distance_m"] * normal[None, :]
        + 0.0001
        * (np.asarray(x) + pixel_offset - geometry["beam_xy_pixels"][0])[:, None]
        * fast[None, :]
        + 0.0001
        * (np.asarray(y) + pixel_offset - geometry["beam_xy_pixels"][1])[:, None]
        * slow[None, :]
    )
    outgoing = point / np.linalg.norm(point, axis=1)[:, None]
    incident = np.asarray(geometry["incident_beam"], dtype=float)
    incident = incident / np.linalg.norm(incident)
    return (outgoing - incident[None, :]) / geometry["wavelength_angstrom"]


def sample_orientation(
    geometry: dict[str, Any],
    header: dict[str, Any],
    omega_degrees: float,
) -> np.ndarray:
    return (
        axis_rotation(np.asarray(header["omega_axis"]), omega_degrees)
        @ axis_rotation(np.asarray(header["kappa_axis"]), geometry["kappa_degrees"])
        @ axis_rotation(np.asarray(header["phi_axis"]), geometry["phi_degrees"])
    )


def q_for_frame_points(
    frame: dict[str, Any],
    header: dict[str, Any],
    x: np.ndarray,
    y: np.ndarray,
    omega_degrees: float,
    pixel_offset: float = 0.5,
) -> np.ndarray:
    q_lab = detector_q_lab(frame["geometry"], x, y, pixel_offset)
    orientation = sample_orientation(
        frame["geometry"], header, omega_degrees
    )
    return q_lab @ orientation


def tier_for_peak(peak: dict[str, Any]) -> str:
    validity = peak["q"]["validity"]
    if validity == "CERTIFIED_POSITIVE_SEPARATION_INTERIOR":
        return "PRIMARY_LATTICE_TIER"
    if validity == "CERTIFIED_POSITIVE_SEPARATION_EDGE_OR_MASK_ADJACENT":
        return "SENSITIVITY_TIER"
    return "UNRESOLVED_CHALLENGE_TIER"


def compute_all_q(
    d4: dict[str, Any],
    headers: dict[str, dict[str, Any]],
) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]], dict[str, int]]:
    frames = {frame["frame_id"]: frame for frame in d4["frames"]}
    peaks = sorted(d4["peaks"], key=lambda item: item["peak_id"])
    positions = defaultdict(list)
    for index, peak in enumerate(peaks):
        positions[peak["f"]].append(index)
    q = np.empty((len(peaks), 3), dtype=float)
    radius = np.empty(len(peaks), dtype=float)
    tier_counts: Counter[str] = Counter()
    for frame_id, indices in positions.items():
        frame = frames[frame_id]
        header = headers[frame["source_archive_path"]]
        selected = [peaks[index] for index in indices]
        x = np.array([peak["x"] for peak in selected], dtype=float)
        y = np.array([peak["y"] for peak in selected], dtype=float)
        start = frame["geometry"]["omega_degrees"]
        increment = frame["geometry"]["omega_increment_degrees"]
        middle = start + 0.5 * increment
        center = q_for_frame_points(frame, header, x, y, middle)
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
        for corner_x, corner_y, omega in itertools.product(
            x_options, y_options, (start, start + increment)
        ):
            corner = q_for_frame_points(
                frame, header, corner_x, corner_y, omega
            )
            local_radius = np.maximum(
                local_radius, np.linalg.norm(corner - center, axis=1)
            )
        radius[indices] = local_radius
        tier_counts.update(tier_for_peak(peak) for peak in selected)
    return q, radius, peaks, dict(sorted(tier_counts.items()))


def balanced_peaks(
    peaks: list[dict[str, Any]],
    frames: dict[str, dict[str, Any]],
    validity: str,
    cap_per_scan: int,
) -> list[dict[str, Any]]:
    by_scan: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for peak in peaks:
        if peak["q"]["validity"] == validity:
            scan = frames[peak["f"]]["scan_configuration_id"]
            by_scan[scan].append(peak)
    selected: list[dict[str, Any]] = []
    for scan in sorted(by_scan):
        ordered = sorted(
            by_scan[scan],
            key=lambda peak: (
                -peak["integrated_signal_to_formal_noise"],
                peak["peak_id"],
            ),
        )
        selected.extend(ordered[:cap_per_scan])
    return selected


def q_lab_angle_arrays(
    selected: list[dict[str, Any]],
    frames: dict[str, dict[str, Any]],
    pixel_offset: float,
    frame_fraction: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    q_lab = []
    angles = []
    scan_labels = []
    scan_ids = {
        scan: index
        for index, scan in enumerate(
            sorted(
                {
                    frames[peak["f"]]["scan_configuration_id"]
                    for peak in selected
                }
            )
        )
    }
    for peak in selected:
        frame = frames[peak["f"]]
        geometry = frame["geometry"]
        q_lab.append(
            detector_q_lab(
                geometry,
                np.array([peak["x"]]),
                np.array([peak["y"]]),
                pixel_offset,
            )[0]
        )
        angles.append(
            [
                geometry["omega_degrees"]
                + frame_fraction * geometry["omega_increment_degrees"],
                geometry["kappa_degrees"],
                geometry["phi_degrees"],
            ]
        )
        scan_labels.append(scan_ids[frame["scan_configuration_id"]])
    return np.asarray(q_lab), np.asarray(angles), np.asarray(scan_labels)


def transformed_q(
    q_lab: np.ndarray,
    angles: np.ndarray,
    axes: dict[str, np.ndarray],
    order: tuple[str, str, str],
    signs: tuple[int, int, int],
) -> np.ndarray:
    angle_column = {"omega": 0, "kappa": 1, "phi": 2}
    matrices = [
        batch_axis_rotation(
            axes[key], signs[index] * angles[:, angle_column[key]]
        )
        for index, key in enumerate(order)
    ]
    orientation = matrices[0] @ matrices[1] @ matrices[2]
    return np.einsum(
        "nij,nj->ni", np.transpose(orientation, (0, 2, 1)), q_lab
    )


def cross_scan_coincidence(
    q: np.ndarray, scan_labels: np.ndarray, nearest_k: int = 16
) -> dict[str, Any]:
    tree = cKDTree(q)
    distances, indices = tree.query(q, k=nearest_k, workers=-1)
    best = np.full(len(q), np.inf)
    for column in range(1, nearest_k):
        other_scan = scan_labels[indices[:, column]] != scan_labels
        best = np.minimum(
            best, np.where(other_scan, distances[:, column], np.inf)
        )
    finite = best[np.isfinite(best)]
    return {
        "count_below_0_002": int(np.sum(finite < 0.002)),
        "count_below_0_003": int(np.sum(finite < 0.003)),
        "count_below_0_005": int(np.sum(finite < 0.005)),
        "distance_quantiles": {
            str(probability): rounded(np.quantile(finite, probability))
            for probability in (0.01, 0.05, 0.1, 0.25, 0.5)
        },
    }


def convention_census(
    d4: dict[str, Any],
    headers: dict[str, dict[str, Any]],
    rule: dict[str, Any],
) -> dict[str, Any]:
    frames = {frame["frame_id"]: frame for frame in d4["frames"]}
    selected = balanced_peaks(
        d4["peaks"],
        frames,
        "CERTIFIED_POSITIVE_SEPARATION_INTERIOR",
        rule["convention_census"]["balanced_peak_cap_per_scan"],
    )
    first_header = headers[frames[selected[0]["f"]]["source_archive_path"]]
    axes = {
        "omega": np.asarray(first_header["omega_axis"]),
        "kappa": np.asarray(first_header["kappa_axis"]),
        "phi": np.asarray(first_header["phi_axis"]),
    }
    q_lab, angles, scans = q_lab_angle_arrays(selected, frames, 0.5, 0.5)
    rotation_records = []
    for order in itertools.permutations(("omega", "kappa", "phi")):
        for signs in itertools.product((1, -1), repeat=3):
            candidate_q = transformed_q(q_lab, angles, axes, order, signs)
            score = cross_scan_coincidence(candidate_q, scans)
            rotation_records.append(
                {
                    "order": list(order),
                    "signs": list(signs),
                    **score,
                }
            )
    rotation_records.sort(
        key=lambda item: (
            -item["count_below_0_003"],
            -item["count_below_0_002"],
            item["order"],
            item["signs"],
        )
    )
    expected = {
        "order": ["omega", "kappa", "phi"],
        "signs": [1, 1, 1],
    }
    if {
        "order": rotation_records[0]["order"],
        "signs": rotation_records[0]["signs"],
    } != expected:
        raise RuntimeError("raw cross-scan census did not select Ωκφ right-handed order")
    pixel_frame_records = []
    for pixel_offset, frame_fraction in itertools.product(
        rule["convention_census"]["pixel_center_offsets_tested"],
        rule["convention_census"]["frame_fractions_tested"],
    ):
        local_lab, local_angles, local_scans = q_lab_angle_arrays(
            selected, frames, pixel_offset, frame_fraction
        )
        candidate_q = transformed_q(
            local_lab,
            local_angles,
            axes,
            ("omega", "kappa", "phi"),
            (1, 1, 1),
        )
        pixel_frame_records.append(
            {
                "pixel_center_offset": pixel_offset,
                "frame_fraction": frame_fraction,
                **cross_scan_coincidence(candidate_q, local_scans),
            }
        )
    pixel_frame_records.sort(
        key=lambda item: (
            -item["count_below_0_003"],
            -item["count_below_0_002"],
            item["pixel_center_offset"],
            item["frame_fraction"],
        )
    )
    if (
        pixel_frame_records[0]["pixel_center_offset"],
        pixel_frame_records[0]["frame_fraction"],
    ) != (0.5, 0.5):
        raise RuntimeError("raw cross-scan census did not select pixel/frame midpoints")
    return {
        "balanced_peak_count": len(selected),
        "selection_uses_processed_or_reference_data": False,
        "rotation_candidate_count": len(rotation_records),
        "rotation_candidates_ranked": rotation_records,
        "selected_rotation": rotation_records[0],
        "runner_up_rotation": rotation_records[1],
        "pixel_frame_candidate_count": len(pixel_frame_records),
        "pixel_frame_candidates_ranked": pixel_frame_records,
        "selected_pixel_frame_convention": pixel_frame_records[0],
    }


def canonical_vector_sign(vectors: np.ndarray) -> np.ndarray:
    vectors = np.asarray(vectors, dtype=float).copy()
    for index, vector in enumerate(vectors):
        first = int(np.argmax(np.abs(vector) > 1e-12))
        if vector[first] < 0:
            vectors[index] *= -1
    return vectors


def seed_lattice(
    primary_q: np.ndarray,
    primary_peaks: list[dict[str, Any]],
    frames: dict[str, dict[str, Any]],
    rule: dict[str, Any],
) -> tuple[np.ndarray, dict[str, Any]]:
    peak_to_q = {
        peak["peak_id"]: primary_q[index]
        for index, peak in enumerate(primary_peaks)
    }
    selected = balanced_peaks(
        primary_peaks,
        frames,
        "CERTIFIED_POSITIVE_SEPARATION_INTERIOR",
        rule["lattice_search"]["balanced_peak_cap_per_scan"],
    )
    q = np.asarray([peak_to_q[peak["peak_id"]] for peak in selected])
    nearest_k = rule["lattice_search"]["nearest_neighbors_per_peak"]
    tree = cKDTree(q)
    _, neighbors = tree.query(q, k=nearest_k, workers=-1)
    differences = []
    norm_min = rule["lattice_search"][
        "difference_norm_min_cycles_per_angstrom"
    ]
    norm_max = rule["lattice_search"][
        "difference_norm_max_cycles_per_angstrom"
    ]
    for column in range(1, nearest_k):
        vectors = q[neighbors[:, column]] - q
        norms = np.linalg.norm(vectors, axis=1)
        vectors = vectors[(norms > norm_min) & (norms < norm_max)]
        differences.append(canonical_vector_sign(vectors))
    difference_array = np.concatenate(differences, axis=0)
    width = rule["lattice_search"][
        "difference_voxel_width_cycles_per_angstrom"
    ]
    keys = np.rint(difference_array / width).astype(np.int64)
    unique_keys, counts = np.unique(keys, axis=0, return_counts=True)
    threshold = int(
        math.ceil(
            counts.max()
            * rule["lattice_search"]["mode_support_fraction_of_maximum"]
        )
    )
    order = np.argsort(counts)[::-1]
    raw_modes = []
    local_radius = width * rule["lattice_search"]["mode_local_radius_voxels"]
    for key_index in order:
        if counts[key_index] < threshold:
            break
        center = unique_keys[key_index] * width
        local = (
            np.linalg.norm(difference_array - center[None, :], axis=1)
            < local_radius
        )
        median = np.median(difference_array[local], axis=0)
        raw_modes.append(
            {
                "bin_count": int(counts[key_index]),
                "local_count": int(np.sum(local)),
                "vector": median,
                "norm": float(np.linalg.norm(median)),
            }
        )
    deduplication = rule["lattice_search"][
        "mode_deduplication_radius_cycles_per_angstrom"
    ]
    modes = []
    for mode in sorted(raw_modes, key=lambda item: (-item["bin_count"], item["norm"])):
        if any(
            np.linalg.norm(mode["vector"] - kept["vector"]) < deduplication
            for kept in modes
        ):
            continue
        modes.append(mode)
    by_norm = sorted(modes, key=lambda item: (item["norm"], -item["bin_count"]))
    first = by_norm[0]
    second = next(
        mode
        for mode in by_norm[1:]
        if np.linalg.norm(np.cross(first["vector"], mode["vector"]))
        / (first["norm"] * mode["norm"])
        > rule["lattice_search"]["minimum_pair_sine"]
    )
    third = next(
        mode
        for mode in by_norm[1:]
        if abs(
            np.linalg.det(
                np.column_stack(
                    [first["vector"], second["vector"], mode["vector"]]
                )
            )
        )
        / (first["norm"] * second["norm"] * mode["norm"])
        > rule["lattice_search"]["minimum_triple_normalized_determinant"]
    )
    seed = np.column_stack(
        [first["vector"], second["vector"], third["vector"]]
    )
    if np.linalg.det(seed) < 0:
        seed[:, 2] *= -1
    record = {
        "balanced_peak_count": len(selected),
        "difference_vector_count": len(difference_array),
        "voxel_width_cycles_per_angstrom": width,
        "maximum_voxel_count": int(counts.max()),
        "mode_minimum_bin_count": threshold,
        "deduplicated_mode_count": len(modes),
        "deduplicated_modes": [
            {
                **{key: value for key, value in mode.items() if key != "vector"},
                "vector": rounded_vector(mode["vector"]),
            }
            for mode in by_norm
        ],
        "selected_seed_vectors": rounded_matrix(seed.T),
        "selected_seed_determinant": rounded(abs(np.linalg.det(seed)), 12),
    }
    return seed, record


def refine_lattice(
    q: np.ndarray, seed: np.ndarray, rule: dict[str, Any]
) -> tuple[np.ndarray, dict[str, Any]]:
    basis = np.asarray(seed, dtype=float).copy()
    schedule_records = []
    for cutoff in rule["lattice_search"][
        "integer_refinement_cut_schedule_cycles_per_angstrom"
    ]:
        h = np.rint(q @ np.linalg.inv(basis).T).astype(np.int64)
        residual = np.linalg.norm(q - h @ basis.T, axis=1)
        accepted = (residual < cutoff) & np.any(h != 0, axis=1)
        updated = np.linalg.lstsq(h[accepted], q[accepted], rcond=None)[0].T
        schedule_records.append(
            {
                "cutoff": cutoff,
                "accepted_count": int(np.sum(accepted)),
                "median_residual": rounded(np.median(residual[accepted])),
                "p90_residual": rounded(np.quantile(residual[accepted], 0.9)),
                "maximum_basis_delta": rounded(np.max(np.abs(updated - basis))),
            }
        )
        basis = updated
    tukey_cutoff = rule["lattice_search"][
        "final_tukey_cutoff_cycles_per_angstrom"
    ]
    tukey_records = []
    for iteration in range(
        rule["lattice_search"]["final_tukey_iterations"]
    ):
        h = np.rint(q @ np.linalg.inv(basis).T).astype(np.int64)
        residual = np.linalg.norm(q - h @ basis.T, axis=1)
        scaled = residual / tukey_cutoff
        weights = np.where(scaled < 1.0, (1.0 - scaled**2) ** 2, 0.0)
        root = np.sqrt(weights)
        updated = np.linalg.lstsq(
            h * root[:, None], q * root[:, None], rcond=None
        )[0].T
        tukey_records.append(
            {
                "iteration": iteration + 1,
                "positive_weight_count": int(np.sum(weights > 0)),
                "maximum_basis_delta": rounded(np.max(np.abs(updated - basis))),
            }
        )
        basis = updated
    direct = np.linalg.inv(basis).T
    permutation = np.argsort(np.linalg.norm(direct, axis=0))
    direct = direct[:, permutation]
    basis = basis[:, permutation]
    if np.linalg.det(direct) < 0:
        direct[:, 2] *= -1
        basis[:, 2] *= -1
    if not np.allclose(basis.T @ direct, np.eye(3), atol=1e-9, rtol=0.0):
        raise RuntimeError("reciprocal/direct duality failure")
    return basis, {
        "integer_refinement_schedule": schedule_records,
        "tukey_refinement": {
            "cutoff": tukey_cutoff,
            "iterations": tukey_records,
        },
        "final_direct_length_sort_permutation": permutation.tolist(),
    }


def cell_metric(direct: np.ndarray, reciprocal: np.ndarray) -> dict[str, Any]:
    lengths = np.linalg.norm(direct, axis=0)
    pairs = ((1, 2), (0, 2), (0, 1))
    angles = [
        math.degrees(
            math.acos(
                float(
                    np.clip(
                        np.dot(direct[:, left], direct[:, right])
                        / (lengths[left] * lengths[right]),
                        -1.0,
                        1.0,
                    )
                )
            )
        )
        for left, right in pairs
    ]
    return {
        "direct_basis_columns_angstrom": rounded_matrix(direct),
        "reciprocal_basis_columns_cycles_per_angstrom": rounded_matrix(reciprocal),
        "direct_lengths_angstrom_length_sorted": rounded_vector(lengths),
        "direct_angles_degrees_alpha_beta_gamma": rounded_vector(angles),
        "direct_metric_tensor_angstrom_squared": rounded_matrix(direct.T @ direct),
        "reciprocal_metric_tensor_cycles_squared_per_angstrom_squared": rounded_matrix(
            reciprocal.T @ reciprocal
        ),
        "direct_cell_volume_angstrom_cubed": rounded(abs(np.linalg.det(direct))),
        "reciprocal_cell_volume_cycles_cubed_per_angstrom_cubed": rounded(
            abs(np.linalg.det(reciprocal)), 12
        ),
        "duality_maximum_absolute_error": rounded(
            np.max(np.abs(reciprocal.T @ direct - np.eye(3))), 12
        ),
        "metric_symmetry_assignment": "NOT_PERFORMED",
        "angle_snapping_performed": False,
    }


def support_census(
    q: np.ndarray,
    radius: np.ndarray,
    peaks: list[dict[str, Any]],
    basis: np.ndarray,
    thresholds: list[float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    h = np.rint(q @ np.linalg.inv(basis).T).astype(np.int64)
    residual = np.linalg.norm(q - h @ basis.T, axis=1)
    compatible = residual <= radius + 1e-12
    tier_indices: dict[str, list[int]] = defaultdict(list)
    for index, peak in enumerate(peaks):
        tier_indices[tier_for_peak(peak)].append(index)
    tiers = {}
    for tier, indices_list in sorted(tier_indices.items()):
        indices = np.asarray(indices_list, dtype=np.int64)
        local_residual = residual[indices]
        local_compatible = compatible[indices]
        local_h = h[indices]
        tiers[tier] = {
            "peak_count": len(indices),
            "formal_region_compatible_count": int(np.sum(local_compatible)),
            "formal_region_compatible_fraction": rounded(
                np.mean(local_compatible)
            ),
            "center_residual_quantiles": {
                str(probability): rounded(
                    np.quantile(local_residual, probability)
                )
                for probability in (0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1.0)
            },
            "formal_region_radius_quantiles": {
                str(probability): rounded(np.quantile(radius[indices], probability))
                for probability in (0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1.0)
            },
            "center_threshold_counts": {
                str(threshold): int(np.sum(local_residual < threshold))
                for threshold in thresholds
            },
            "unique_nearest_integer_triplet_count_formal_compatible": int(
                len(np.unique(local_h[local_compatible], axis=0))
            ),
        }
    return h, residual, compatible, {"tiers": tiers}


def ambiguity_record(
    direct: np.ndarray,
    reciprocal: np.ndarray,
    q: np.ndarray,
    residual: np.ndarray,
    compatible: np.ndarray,
    h: np.ndarray,
) -> dict[str, Any]:
    transforms = [
        np.eye(3, dtype=int),
        np.array([[1, 1, 0], [0, 1, 0], [0, 0, 1]], dtype=int),
        np.array([[1, 0, 0], [0, 1, 1], [0, 0, 1]], dtype=int),
        np.array([[1, 0, 1], [0, 1, 0], [0, 0, 1]], dtype=int),
    ]
    equivalent = []
    for transform in transforms:
        alternative_direct = direct @ transform
        alternative_reciprocal = reciprocal @ np.linalg.inv(transform).T
        alternative_h = np.rint(
            q @ np.linalg.inv(alternative_reciprocal).T
        ).astype(np.int64)
        alternative_residual = np.linalg.norm(
            q - alternative_h @ alternative_reciprocal.T, axis=1
        )
        equivalent.append(
            {
                "unimodular_direct_basis_transform": transform.tolist(),
                "determinant": int(round(np.linalg.det(transform))),
                "direct_basis_columns_angstrom": rounded_matrix(alternative_direct),
                "reciprocal_basis_columns_cycles_per_angstrom": rounded_matrix(
                    alternative_reciprocal
                ),
                "maximum_residual_change": rounded(
                    np.max(np.abs(alternative_residual - residual)), 12
                ),
            }
        )
    supercells = []
    for axis in range(3):
        transform = np.eye(3)
        transform[axis, axis] = 2.0
        alternative_direct = direct @ transform
        alternative_reciprocal = reciprocal @ np.linalg.inv(transform).T
        alternative_h = np.rint(
            q @ np.linalg.inv(alternative_reciprocal).T
        ).astype(np.int64)
        alternative_residual = np.linalg.norm(
            q - alternative_h @ alternative_reciprocal.T, axis=1
        )
        supercells.append(
            {
                "doubled_direct_axis_zero_based": axis,
                "direct_transform": transform.astype(int).tolist(),
                "inequivalent_lattice_relation": (
                    "RECIPROCAL_SUPERLATTICE_CONTAINS_THE_OBSERVED_CANDIDATE; "
                    "ADDITIONAL_HALF_STEP_NODES_ARE_NOT_CERTIFIED_ABSENT"
                ),
                "maximum_residual_change_for_observed_centers": rounded(
                    np.max(np.abs(alternative_residual - residual)), 12
                ),
                "formal_compatible_count_for_observed_centers": int(
                    np.sum(compatible)
                ),
                "not_excluded_without_completeness_or_systematic_absence_evidence": True,
            }
        )
    primary_h = h[compatible]
    parity = {
        f"axis_{axis}": {
            "even": int(np.sum(primary_h[:, axis] % 2 == 0)),
            "odd": int(np.sum(primary_h[:, axis] % 2 != 0)),
        }
        for axis in range(3)
    }
    return {
        "basis_equivalence_relation": "GL_3_Z_UNIMODULAR_CHANGE_OF_BASIS",
        "equivalent_basis_witnesses": equivalent,
        "physical_primitivity_status": "NOT_ESTABLISHED",
        "real_space_axis_doubling_countermodels": supercells,
        "compatible_index_parity_census": parity,
        "decisive_limitation": (
            "D4 HAS NO CALIBRATED FALSE_NEGATIVE OR COMPLETE_ABSENCE MODEL; "
            "OBSERVED PEAK POSITIONS ALONE DO NOT EXCLUDE REAL_SPACE_SUPERCELLS"
        ),
    }


def diagnostic_census(
    d4: dict[str, Any],
    q: np.ndarray,
    radius: np.ndarray,
    peaks: list[dict[str, Any]],
    residual: np.ndarray,
    compatible: np.ndarray,
    basis: np.ndarray,
    headers: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    isolation_ids = set(
        d4["post_construction_diagnostics"]["possible_false_positive_proxy"][
            "peak_ids"
        ]
    )
    primary = np.array(
        [
            tier_for_peak(peak) == "PRIMARY_LATTICE_TIER"
            for peak in peaks
        ]
    )
    isolation = np.array([peak["peak_id"] in isolation_ids for peak in peaks])

    def group(mask: np.ndarray) -> dict[str, Any]:
        count = int(np.sum(mask))
        return {
            "peak_count": count,
            "formal_compatible_count": int(np.sum(compatible[mask])),
            "formal_compatible_fraction": rounded(np.mean(compatible[mask]))
            if count
            else None,
            "center_within_0_010_count": int(np.sum(residual[mask] < 0.01)),
            "center_within_0_010_fraction": rounded(
                np.mean(residual[mask] < 0.01)
            )
            if count
            else None,
        }

    frames = {frame["frame_id"]: frame for frame in d4["frames"]}
    gap_residuals = []
    for gap in d4["post_construction_diagnostics"][
        "possible_false_negative_proxy"
    ]["gaps"]:
        frame = frames[gap["middle_frame_id"]]
        header = headers[frame["source_archive_path"]]
        middle = (
            frame["geometry"]["omega_degrees"]
            + 0.5 * frame["geometry"]["omega_increment_degrees"]
        )
        gap_q = q_for_frame_points(
            frame,
            header,
            np.array([gap["expected_x"]]),
            np.array([gap["expected_y"]]),
            middle,
        )[0]
        gap_h = np.rint(np.linalg.inv(basis) @ gap_q).astype(np.int64)
        gap_residuals.append(float(np.linalg.norm(gap_q - basis @ gap_h)))
    gap_residuals_array = np.asarray(gap_residuals)
    return {
        "peak_set_or_lattice_modified_by_diagnostics": False,
        "possible_false_positive_proxy": {
            "all_primary": group(primary),
            "isolation_proxy_primary": group(primary & isolation),
            "non_isolation_primary": group(primary & ~isolation),
            "not_a_calibrated_false_positive_rate": True,
        },
        "possible_false_negative_gap_proxy": {
            "gap_count": len(gap_residuals),
            "candidate_lattice_center_threshold_counts": {
                str(threshold): int(np.sum(gap_residuals_array < threshold))
                for threshold in (0.003, 0.005, 0.008, 0.01)
            },
            "residual_quantiles": {
                str(probability): rounded(
                    np.quantile(gap_residuals_array, probability)
                )
                for probability in (0.0, 0.25, 0.5, 0.75, 0.9, 1.0)
            },
            "not_a_calibrated_false_negative_rate": True,
        },
    }


def q_scientific_projection(records: list[dict[str, Any]]) -> str:
    projection = [
        {
            "peak_id": record["peak_id"],
            "tier": record["tier"],
            "q": record["q_cycles_per_angstrom"],
            "h": record["nearest_integer_triplet"],
            "residual": record["center_residual_cycles_per_angstrom"],
            "radius": record["formal_q_region_corner_radius"],
            "compatible": record["formal_region_compatible"],
        }
        for record in sorted(records, key=lambda item: item["peak_id"])
    ]
    return sha256_bytes(canonical_bytes(projection))


def result_markdown(result: dict[str, Any]) -> str:
    metric = result["lattice_candidate"]["metric"]
    support = result["lattice_candidate"]["support"]["tiers"]
    primary = support["PRIMARY_LATTICE_TIER"]
    sensitivity = support["SENSITIVITY_TIER"]
    unresolved = support["UNRESOLVED_CHALLENGE_TIER"]
    lines = [
        "# NFC Crystallography Specimen A — D5 raw reciprocal-lattice result",
        "",
        f"Artifact: `{result['artifact_id']}`",
        "",
        f"Outcome: `{result['outcome']}`",
        "",
        "Scientific label: `OPEN_DEVELOPMENT_NOT_INDEPENDENT_VALIDATION`",
        "",
        "## D5a detector-to-reciprocal binding",
        "",
        "The construction reads only the frozen D4 peak corpus and geometry "
        "embedded in the 2,408 raw CBF headers. It reads no processed HKL, "
        "XDS input/output, known cell, orientation matrix, space group, "
        "SHELX/FCF/CIF, or solved structure.",
        "",
        "The fixed mapping uses pixel centers `(x+0.5,y+0.5)`, the reported "
        "fast/slow detector axes, `fast×slow` detector normal, 34 mm distance, "
        "1.541838 Å wavelength, frame midpoint, and the right-handed "
        "`Omega*Kappa*Phi` sample orientation. Reciprocal vectors are reported "
        "in cycles/Å without a 2π factor.",
        "",
        f"The raw-only convention census tested "
        f"{result['d5a']['convention_census']['rotation_candidate_count']} "
        "rotation-order/sign candidates and "
        f"{result['d5a']['convention_census']['pixel_frame_candidate_count']} "
        "pixel/frame conventions. The declared mapping ranks first by the "
        "predeclared 0.003 cycles/Å cross-scan coincidence count.",
        "",
        "This is a binding to the CBF-declared right-handed miniCBF geometry. "
        "The headers contain the CAP-model warning but no CAP Esperanto payload; "
        "therefore no stronger independently calibrated instrument-model claim "
        "is made.",
        "",
        "## D5b lattice candidate",
        "",
        "A deterministic nearest-neighbor difference-mode search and robust "
        "integer refinement recovered one minimal observed reciprocal "
        "translation-lattice candidate without prior cell information.",
        "",
        f"* Length-sorted direct lengths (Å): "
        f"`{metric['direct_lengths_angstrom_length_sorted']}`",
        f"* Direct angles α, β, γ (degrees): "
        f"`{metric['direct_angles_degrees_alpha_beta_gamma']}`",
        f"* Direct volume (Å³): `{metric['direct_cell_volume_angstrom_cubed']}`",
        f"* Primary formal-region compatibility: "
        f"{primary['formal_region_compatible_count']:,}/"
        f"{primary['peak_count']:,} "
        f"({primary['formal_region_compatible_fraction']:.6f})",
        f"* Sensitivity-tier compatibility: "
        f"{sensitivity['formal_region_compatible_count']:,}/"
        f"{sensitivity['peak_count']:,} "
        f"({sensitivity['formal_region_compatible_fraction']:.6f})",
        f"* Unresolved challenge-tier compatibility: "
        f"{unresolved['formal_region_compatible_count']:,}/"
        f"{unresolved['peak_count']:,} "
        f"({unresolved['formal_region_compatible_fraction']:.6f})",
        "",
        "The candidate metric is not snapped to an orthogonal or other crystal "
        "system. Crystal symmetry and a physical primitive unit cell remain "
        "unestablished.",
        "",
        "## Explicit ambiguity",
        "",
        "Unimodular basis changes give exactly the same lattice and are recorded "
        "as equivalent witnesses. In addition, doubling any direct-space axis "
        "creates an inequivalent supercell whose reciprocal superlattice still "
        "contains every observed candidate node. Because D4 has no calibrated "
        "false-negative or complete-absence model, the extra half-step reciprocal "
        "nodes cannot yet be certified absent. D5 therefore establishes a "
        "minimal observed lattice candidate, not physical primitivity.",
        "",
        "## Scientific boundary",
        "",
        "D5 does not establish crystal symmetry, a conventional space group, "
        "a Miller-index binding to processed data, novel-orbit intensity "
        "prediction, a solved structure, or independent empirical validation. "
        "D6 may analyze symmetry only on this raw candidate while preserving "
        "the CAP-model, calibration, contamination, unresolved-peak, and "
        "supercell limitations.",
        "",
        "## Bindings",
        "",
        f"* Raw archive SHA-256: `{result['bindings']['raw_archive_sha256']}`",
        f"* D4 corpus gzip SHA-256: `{result['bindings']['d4_corpus_gzip_sha256']}`",
        f"* D5 rule SHA-256: `{result['bindings']['d5_rule_sha256']}`",
        f"* Reciprocal corpus gzip SHA-256: `{result['reciprocal_corpus']['gzip_sha256']}`",
        f"* Reciprocal scientific projection SHA-256: "
        f"`{result['reciprocal_corpus']['scientific_projection_sha256']}`",
        f"* Pipeline SHA-256: `{result['bindings']['pipeline_sha256']}`",
        "",
    ]
    return "\n".join(lines)


def run() -> None:
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    if RAW_ARCHIVE.stat().st_size != EXPECTED_ARCHIVE_BYTES:
        raise RuntimeError("raw archive byte count mismatch")
    raw_archive_sha256 = sha256_file(RAW_ARCHIVE)
    if raw_archive_sha256 != EXPECTED_ARCHIVE_SHA256:
        raise RuntimeError("raw archive digest mismatch")
    d4, _ = load_d4()
    rule = json.loads(RULE_PATH.read_text(encoding="utf-8"))
    headers, header_census = acquire_header_geometry(d4["frames"])
    print(json.dumps({"stage": "headers", **header_census}, sort_keys=True))

    convention = convention_census(d4, headers, rule)
    print(
        json.dumps(
            {
                "stage": "convention_census",
                "selected_rotation": convention["selected_rotation"],
                "runner_up_rotation": convention["runner_up_rotation"],
                "selected_pixel_frame": convention[
                    "selected_pixel_frame_convention"
                ],
            },
            sort_keys=True,
        )
    )

    q, radius, peaks, tier_counts = compute_all_q(d4, headers)
    frames = {frame["frame_id"]: frame for frame in d4["frames"]}
    primary_indices = np.asarray(
        [
            index
            for index, peak in enumerate(peaks)
            if tier_for_peak(peak) == "PRIMARY_LATTICE_TIER"
        ],
        dtype=np.int64,
    )
    primary_q = q[primary_indices]
    primary_peaks = [peaks[index] for index in primary_indices]
    if len(primary_peaks) != 76_166:
        raise RuntimeError("primary D5 tier count mismatch")
    seed, seed_record = seed_lattice(
        primary_q, primary_peaks, frames, rule
    )
    basis, refinement = refine_lattice(primary_q, seed, rule)
    direct = np.linalg.inv(basis).T
    metric = cell_metric(direct, basis)
    print(json.dumps({"stage": "lattice", "metric": metric}, sort_keys=True))

    h, residual, compatible, support = support_census(
        q,
        radius,
        peaks,
        basis,
        rule["lattice_search"][
            "center_residual_thresholds_cycles_per_angstrom"
        ],
    )
    ambiguity = ambiguity_record(
        direct,
        basis,
        q[primary_indices],
        residual[primary_indices],
        compatible[primary_indices],
        h[primary_indices],
    )
    diagnostics = diagnostic_census(
        d4, q, radius, peaks, residual, compatible, basis, headers
    )

    d5a = seal(
        {
            "artifact_id": (
                f"NFC_CRYST_9Z6F_D5A_DETECTOR_TO_RECIPROCAL_BINDING_{VERSION}"
            ),
            "status": "CBF_DECLARED_DETECTOR_TO_RECIPROCAL_BINDING_ESTABLISHED",
            "scientific_scope": "OPEN_DEVELOPMENT_NOT_INDEPENDENT_VALIDATION",
            "header_geometry_census": header_census,
            "mapping": rule["coordinate_convention"],
            "convention_census": convention,
            "geometric_checks": {
                "detector_fast_slow_orthonormality_maximum_error": rounded(
                    max(
                        abs(
                            np.dot(
                                frame["geometry"]["detector_fast_axis"],
                                frame["geometry"]["detector_slow_axis"],
                            )
                        )
                        for frame in d4["frames"]
                    ),
                    12,
                ),
                "detector_normal_two_theta_maximum_error_degrees": rounded(
                    max(
                        abs(
                            math.degrees(
                                math.acos(
                                    float(
                                        np.clip(
                                            np.dot(
                                                np.cross(
                                                    frame["geometry"][
                                                        "detector_fast_axis"
                                                    ],
                                                    frame["geometry"][
                                                        "detector_slow_axis"
                                                    ],
                                                )
                                                / np.linalg.norm(
                                                    np.cross(
                                                        frame["geometry"][
                                                            "detector_fast_axis"
                                                        ],
                                                        frame["geometry"][
                                                            "detector_slow_axis"
                                                        ],
                                                    )
                                                ),
                                                np.asarray(
                                                    frame["geometry"][
                                                        "incident_beam"
                                                    ]
                                                )
                                                / np.linalg.norm(
                                                    frame["geometry"][
                                                        "incident_beam"
                                                    ]
                                                ),
                                            ),
                                            -1.0,
                                            1.0,
                                        )
                                    )
                                )
                            )
                            - abs(frame["geometry"]["detector_2theta_degrees"])
                        )
                        for frame in d4["frames"]
                    )
                ),
            },
            "qualification": (
                "CBF_DECLARED_RIGHT_HANDED_MINICBF_GEOMETRY_ONLY; "
                "CAP_ESPERANTO_PAYLOAD_ABSENT; NO_INDEPENDENT_INSTRUMENT_CALIBRATION"
            ),
            "prohibited_construction_inputs_accessed": [],
        }
    )
    write_json(D5A_PATH, d5a)

    d5b = seal(
        {
            "artifact_id": (
                f"NFC_CRYST_9Z6F_D5B_RAW_RECIPROCAL_LATTICE_CANDIDATE_{VERSION}"
            ),
            "status": (
                "MINIMAL_OBSERVED_RECIPROCAL_TRANSLATION_LATTICE_CANDIDATE_ESTABLISHED"
            ),
            "scientific_scope": "OPEN_DEVELOPMENT_NOT_INDEPENDENT_VALIDATION",
            "input_tier_counts": tier_counts,
            "seed_search": seed_record,
            "refinement": refinement,
            "metric": metric,
            "support": support,
            "ambiguity": ambiguity,
            "contamination_diagnostics": diagnostics,
            "claims_not_established": [
                "PHYSICAL_PRIMITIVE_UNIT_CELL",
                "CRYSTAL_SYSTEM",
                "CRYSTAL_SYMMETRY_QUOTIENT",
                "SPACE_GROUP",
                "PROCESSED_MILLER_INDEX_BINDING",
                "NOVEL_ORBIT_INTENSITY_PREDICTION",
                "SOLVED_STRUCTURE",
                "INDEPENDENT_EMPIRICAL_VALIDATION"
            ],
            "prohibited_construction_inputs_accessed": [],
        }
    )
    write_json(D5B_PATH, d5b)

    isolation_ids = set(
        d4["post_construction_diagnostics"]["possible_false_positive_proxy"][
            "peak_ids"
        ]
    )
    q_records = []
    for index, peak in enumerate(peaks):
        q_records.append(
            {
                "peak_id": peak["peak_id"],
                "frame_id": peak["f"],
                "tier": tier_for_peak(peak),
                "d4_validity": peak["q"]["validity"],
                "q_cycles_per_angstrom": rounded_vector(q[index]),
                "q_two_pi_per_angstrom": rounded_vector(2.0 * math.pi * q[index]),
                "nearest_integer_triplet": h[index].tolist(),
                "center_residual_cycles_per_angstrom": rounded(residual[index]),
                "formal_q_region_corner_radius": rounded(radius[index]),
                "formal_region_compatible": bool(compatible[index]),
                "isolation_proxy_flag": peak["peak_id"] in isolation_ids,
                "saturation_status": peak["q"]["saturation_status"],
            }
        )
    q_projection = q_scientific_projection(q_records)
    q_corpus = {
        "artifact_id": (
            f"NFC_CRYST_9Z6F_D5_RECIPROCAL_COORDINATE_CORPUS_{VERSION}"
        ),
        "scientific_scope": "OPEN_DEVELOPMENT_NOT_INDEPENDENT_VALIDATION",
        "mapping_artifact_semantic_sha256": d5a["canonical_semantic_sha256"],
        "lattice_artifact_semantic_sha256": d5b["canonical_semantic_sha256"],
        "record_count": len(q_records),
        "tier_counts": tier_counts,
        "scientific_projection_sha256": q_projection,
        "records": q_records,
    }
    q_binding = write_gzip_canonical(Q_CORPUS_PATH, q_corpus)

    outcome = (
        "RAW_RECIPROCAL_LATTICE_CANDIDATE_ESTABLISHED_WITH_EXPLICIT_AMBIGUITY"
    )
    result = seal(
        {
            "artifact_id": (
                f"NFC_CRYST_9Z6F_D5_RAW_RECIPROCAL_LATTICE_RESULT_{VERSION}"
            ),
            "outcome": outcome,
            "scientific_scope": "OPEN_DEVELOPMENT_NOT_INDEPENDENT_VALIDATION",
            "d5a": d5a,
            "lattice_candidate": d5b,
            "reciprocal_corpus": {
                **q_binding,
                "scientific_projection_sha256": q_projection,
                "record_count": len(q_records),
            },
            "bindings": {
                "raw_archive_byte_count": RAW_ARCHIVE.stat().st_size,
                "raw_archive_sha256": raw_archive_sha256,
                "d4_corpus_gzip_sha256": sha256_file(D4_CORPUS),
                "d4_result_sha256": sha256_file(D4_RESULT),
                "d5_rule_sha256": sha256_file(RULE_PATH),
                "pipeline_sha256": sha256_file(Path(__file__).resolve()),
                "numpy_version": np.__version__,
                "scipy_version": scipy.__version__,
            },
            "construction_firewall": {
                "permitted_inputs": rule["construction_inputs"],
                "prohibited_inputs": rule["prohibited_construction_inputs"],
                "prohibited_inputs_accessed": [],
                "retrospective_comparison_in_this_artifact": False,
            },
            "scientific_status": {
                "D5A_DETECTOR_TO_RECIPROCAL_BINDING": (
                    "ESTABLISHED_RELATIVE_TO_CBF_DECLARED_RIGHT_HANDED_GEOMETRY"
                ),
                "D5B_MINIMAL_OBSERVED_RECIPROCAL_LATTICE_CANDIDATE": "ESTABLISHED",
                "D5_PERIODICITY_AND_LATTICE_EQUIVALENCE_CLASS": (
                    "PARTIAL_CANDIDATE_WITH_EXPLICIT_SUPERCELL_AMBIGUITY"
                ),
                "PHYSICAL_PRIMITIVE_UNIT_CELL": "NOT_ESTABLISHED",
                "CRYSTAL_SYMMETRY_QUOTIENT": "NOT_ESTABLISHED",
                "INDEPENDENT_EMPIRICAL_VALIDATION": "NOT_YET_PERFORMED",
            },
            "d6_handoff": (
                "ELIGIBLE_ON_MINIMAL_OBSERVED_LATTICE_CANDIDATE_ONLY_WITH_CAP_"
                "MODEL_CALIBRATION_CONTAMINATION_UNRESOLVED_PEAK_AND_SUPERCELL_"
                "LIMITATIONS_PRESERVED"
            ),
        }
    )
    write_json(RESULT_JSON_PATH, result)
    RESULT_MD_PATH.write_text(result_markdown(result), encoding="utf-8")
    print(
        json.dumps(
            {
                "outcome": outcome,
                "result_semantic_sha256": result["canonical_semantic_sha256"],
                "d5a_semantic_sha256": d5a["canonical_semantic_sha256"],
                "d5b_semantic_sha256": d5b["canonical_semantic_sha256"],
                "reciprocal_corpus": q_binding,
                "metric": metric,
                "primary_support": support["tiers"]["PRIMARY_LATTICE_TIER"],
            },
            sort_keys=True,
        )
    )


def verify_existing() -> None:
    result = json.loads(RESULT_JSON_PATH.read_text(encoding="utf-8"))
    expected = result["canonical_semantic_sha256"]
    body = dict(result)
    body.pop("canonical_semantic_sha256")
    if sha256_bytes(canonical_bytes(body)) != expected:
        raise RuntimeError("result semantic digest mismatch")
    with gzip.open(Q_CORPUS_PATH, "rb") as handle:
        q_corpus = json.loads(handle.read())
    projection = q_scientific_projection(q_corpus["records"])
    if projection != q_corpus["scientific_projection_sha256"]:
        raise RuntimeError("reciprocal corpus projection mismatch")
    if result["reciprocal_corpus"]["scientific_projection_sha256"] != projection:
        raise RuntimeError("result/corpus projection binding mismatch")
    for embedded, path in (
        (result["d5a"], D5A_PATH),
        (result["lattice_candidate"], D5B_PATH),
    ):
        disk = json.loads(path.read_text(encoding="utf-8"))
        if canonical_bytes(disk) != canonical_bytes(embedded):
            raise RuntimeError(f"embedded component mismatch: {path.name}")
    print(
        json.dumps(
            {
                "verification": "PASS",
                "result_semantic_sha256": expected,
                "reciprocal_scientific_projection_sha256": projection,
                "record_count": len(q_corpus["records"]),
            },
            sort_keys=True,
        )
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--verify-existing", action="store_true", help="verify stored D5 artifacts"
    )
    arguments = parser.parse_args()
    if arguments.verify_existing:
        verify_existing()
    else:
        run()


if __name__ == "__main__":
    main()
