#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import urllib.request
from pathlib import Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--destination", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    download = manifest.get("download", {})
    if download.get("verification") != "sha256":
        raise SystemExit(
            "manifest does not provide a whole-object SHA-256; refusing "
            "unverified convenience download"
        )
    url = download["url"]
    expected = download["sha256"]
    args.destination.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(url, args.destination)
    actual = sha256_file(args.destination)
    if actual != expected:
        raise SystemExit(f"SHA-256 mismatch: expected {expected}, received {actual}")
    print(actual)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
