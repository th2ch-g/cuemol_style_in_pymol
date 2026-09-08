"""Isolated real-PyMOL validation in addition to geometric unit tests."""

import json
from pathlib import Path
import subprocess
import sys

import pytest


def test_real_pymol_headless_lifecycle(tmp_path):
    script = Path(__file__).with_name("check_cuemol_style.py")
    result = subprocess.run(
        [sys.executable, str(script), "--output", str(tmp_path)],
        capture_output=True,
        text=True,
        timeout=180,
    )
    if result.returncode == 77:
        pytest.skip("Real PyMOL is unavailable")
    assert result.returncode == 0, result.stdout + result.stderr
    report = json.loads((tmp_path / "report.json").read_text())
    assert "richardson" in report["headless_profiles"]
    assert len(report["headless_profiles"]) >= 25
