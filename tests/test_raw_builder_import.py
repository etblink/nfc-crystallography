from __future__ import annotations

import subprocess
import sys
import textwrap
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def test_raw_builder_imports_frozen_support_modules_from_clean_checkout() -> None:
    script = textwrap.dedent(
        """
        import argparse
        import importlib.util
        import sys
        from pathlib import Path

        root = Path.cwd()
        builder_path = (
            root
            / "methods/baseline_B/raw_pipeline/build_pilatus_d45_feed.py"
        )
        candidate_scripts = (
            root
            / "methods/baseline_B/runtime/restored/candidate/scripts"
        ).resolve()

        spec = importlib.util.spec_from_file_location(
            "nfc_raw_builder_import_smoke",
            builder_path,
        )
        if spec is None or spec.loader is None:
            raise RuntimeError("cannot load raw builder")
        builder = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(builder)

        portability, d45, feed_builder, imported_directory = (
            builder.import_sources(
                argparse.Namespace(candidate_scripts=candidate_scripts)
            )
        )

        expected = {
            "build_fixed_control_feeds": "build_fixed_control_feeds.py",
            "d45_scan_direction_diagnostic": (
                "d45_scan_direction_diagnostic.py"
            ),
            "d45_successor": "d45_successor.py",
            "pilatus_portability": "pilatus_portability.py",
            "prior_run_generalization": "prior_run_generalization.py",
        }
        assert imported_directory == candidate_scripts
        assert Path(portability.__file__).resolve() == (
            candidate_scripts / expected["pilatus_portability"]
        )
        assert Path(d45.__file__).resolve() == (
            candidate_scripts / expected["d45_successor"]
        )
        assert Path(feed_builder.__file__).resolve() == (
            candidate_scripts / expected["build_fixed_control_feeds"]
        )
        for module_name, filename in expected.items():
            assert Path(sys.modules[module_name].__file__).resolve() == (
                candidate_scripts / filename
            )
        """
    )
    subprocess.run(
        [sys.executable, "-I", "-c", script],
        cwd=REPOSITORY_ROOT,
        check=True,
    )
