"""Resolve release records in a checkout or an installed distribution."""

from __future__ import annotations

from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parent
SOURCE_CHECKOUT_ROOT = Path(__file__).resolve().parents[2]
INSTALLED_RELEASE_ROOT = PACKAGE_ROOT / "_release"


def release_root() -> Path:
    """Return the root containing methods, results, and data manifests."""
    required = ("methods", "results", "data_manifests")
    for candidate in (SOURCE_CHECKOUT_ROOT, INSTALLED_RELEASE_ROOT):
        if all((candidate / name).is_dir() for name in required):
            return candidate
    raise FileNotFoundError(
        "release payload is unavailable: expected methods, results, and "
        "data_manifests in either the source checkout or installed wheel"
    )
