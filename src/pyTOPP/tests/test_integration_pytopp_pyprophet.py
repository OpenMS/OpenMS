"""
A very light "integration" check against the CLI module.

We verify the internal helpers that
massage CLI strings (levels/contexts) continue to behave as expected and that
the module can be imported in a dev environment (thanks to the stubs in
conftest.py).
"""

import os
import sys
import re
import importlib
import importlib.util
from pathlib import Path

import pytest


cli = importlib.import_module("pyopenms.pytopp.tools.pytopp_pyprophet_cli")

TEST_OSW = Path(__file__).parent / 'data' / 'pyprophet' / 'test_data.osw'


def _run_cli_with_args(monkeypatch, args):
    """
    Invoke CLI like from the shell: sets sys.argv then calls main().
    Handles both styles:
      - main() returning an int
      - main() calling sys.exit(code)
    """
    monkeypatch.setattr(sys, "argv", ["pytopp_pyprophet_cli.py", *args], raising=False)
    try:
        rc = cli.main()
    except SystemExit as e:
        # pytest treats SystemExit as an exception; normalize to an int rc
        rc = e.code if isinstance(e.code, int) else 1
    return 0 if rc is None else rc


def test_canonical_levels_valid_sets():
    c = cli._canonical_levels
    assert c("ms1") == ["ms1"]
    assert c("ms2") == ["ms2"]
    assert c("ms1,ms2") == ["ms1", "ms2"]
    assert c("ms1ms2") == ["ms1ms2"]

    # order-insensitive; de-duplication
    assert c("ms2,ms1,ms2") == ["ms1", "ms2"]


def test_canonical_levels_invalid():
    import pytest
    c = cli._canonical_levels

    with pytest.raises(ValueError):
        c("foo")

    with pytest.raises(ValueError):
        c("ms1ms2,ms2")  # mixed combos not allowed

    # 'transition' can only be used with one of ms1/ms2/ms1ms2
    with pytest.raises(ValueError):
        c("transition")


def test_parse_contexts():
    p = cli._parse_contexts
    assert p("run-specific") == ["run-specific"]
    assert p("global") == ["global"]
    assert p("global,run-specific") == ["global", "run-specific"]

    # Invalid values raise
    import pytest
    with pytest.raises(ValueError):
        p("other")


def test_cli_pyprophet_dry_run_cmds(monkeypatch, tmp_path, capsys, regtest):
    """
    Snapshot the exact pyprophet invocations (score/export) the CLI would run.
    No pyprophet install required thanks to '-dry_run true'.
    """
    if not TEST_OSW.exists():
        pytest.skip("tests/test_data.osw not present")

    out_osw = tmp_path / "scored.osw"

    rc = _run_cli_with_args(monkeypatch, [
        "-in", str(TEST_OSW),
        "-out", str(out_osw),
        "-scoring:level", "ms1",
        "-dry_run", "true",
        # keep exports on to capture export cmds too; to reduce lines add:
        # "-export:matrix", "false"
    ])
    assert rc == 0

    out = capsys.readouterr().out.splitlines()
    cmds = [ln for ln in out if ln.startswith("CMD:")]

    # Stable snapshot of the command lines we would execute.
    regtest.write("\n".join(cmds))

    # More tolerant assertions: look for 'pyprophet' and a subcommand token later
    pat_score  = re.compile(r"\bpyprophet\b.*\bscore\b", re.IGNORECASE)
    pat_export = re.compile(r"\bpyprophet\b.*\bexport\b", re.IGNORECASE)

    assert any(pat_score.search(ln) for ln in cmds), "No score command found.\nGot:\n" + "\n".join(cmds)
    assert any(pat_export.search(ln) for ln in cmds), "No export command found.\nGot:\n" + "\n".join(cmds)



def test_cli_pyprophet_end_to_end(monkeypatch, tmp_path):
    """
    Test a minimal end-to-end run of the CLI with pyprophet scoring.
    """
    if importlib.util.find_spec("pyprophet") is None:
        pytest.skip("pyprophet not importable")
    if not TEST_OSW.exists():
        pytest.skip("tests/test_data.osw missing")

    out_osw = tmp_path / "scored.osw"

    rc = _run_cli_with_args(monkeypatch, [
        "-in", str(TEST_OSW),
        "-out", str(out_osw),
        "-scoring:level", "ms2",
        "-scoring:classifier", "LDA",
        "-scoring:ss_num_iter", "3",
        "-scoring:xeval_num_iter", "3",
        "-scoring:extra", "+--test +--pi0_lambda 0 0 0",
        "-infer:run_protein", "false",
        "-export:run_matrix", "false",
        "-export:score_report", "false",
        "-dry_run", "false",
    ])
    assert rc == 0

    assert out_osw.exists() and out_osw.stat().st_size > 0
    out_tsv = out_osw.with_suffix(".tsv")
    assert out_tsv.exists() and out_tsv.stat().st_size > 0

    header = out_tsv.read_text(encoding="utf-8", errors="ignore").splitlines()[0]
    plausible_cols = {"m_score", "peptide", "run_id", "transition_group_id"}
    assert any(col in header for col in plausible_cols), f"Unexpected TSV header: {header!r}"

