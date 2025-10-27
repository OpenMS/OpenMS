import os
import sys
import re
import importlib
import importlib.util
from pathlib import Path

import pytest


cli = importlib.import_module("pyopenms.pytopp.tools.pytopp_openswath_chromatogram_extractor_cli")

TEST_MZML = Path(__file__).parent / 'data' / 'OpenSwathChromatogramExtractor_input.mzML'
TEST_TRAML = Path(__file__).parent / 'data' / 'OpenSwathChromatogramExtractor_input.TraML'


def _run_cli_with_args(monkeypatch, args):
    """
    Invoke CLI like from the shell: sets sys.argv then calls main().
    Handles both styles:
      - main() returning an int
      - main() calling sys.exit(code)
    """
    monkeypatch.setattr(sys, "argv", ["pytopp_openswath_chromatogram_extractor_cli.py", *args], raising=False)
    try:
        rc = cli.main()
    except SystemExit as e:
        # pytest treats SystemExit as an exception; normalize to an int rc
        rc = e.code if isinstance(e.code, int) else 1
    return 0 if rc is None else rc

def test_cli_pyprophet_end_to_end(monkeypatch, tmp_path):
    """
    Test a minimal end-to-end run of the CLI.
    """
    if importlib.util.find_spec("pyopenms") is None:
        pytest.skip("pyopenms not importable")
    if not TEST_MZML.exists():
        pytest.skip(f"Test data missing: {TEST_MZML}")
    if not TEST_TRAML.exists():
        pytest.skip(f"Test data missing: {TEST_TRAML}")

    out_mzml = tmp_path / "OpenSwathChromatogramExtractor_output.mzML"

    rc = _run_cli_with_args(monkeypatch, [
        "-in", str(TEST_MZML),
        "-tr", str(TEST_TRAML),
        "-out", str(out_mzml),
    ])
    assert rc == 0

    assert out_mzml.exists() and out_mzml.stat().st_size > 0

