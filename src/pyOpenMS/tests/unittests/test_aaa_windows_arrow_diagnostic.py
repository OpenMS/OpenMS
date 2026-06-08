"""
DIAGNOSTIC (temporary) — pin the Windows heap corruption (exit code 0xC0000374,
STATUS_HEAP_CORRUPTION) seen in the pyOpenMS wheel test for
``ProteinIdentificationArrowIO.exportSearchParamsToParquet``.

Background
----------
In the cibuildwheel Windows job, ``test_arrow_io_classes.py`` crashes the whole
pytest process at ``test_export_search_params_to_parquet``. The identical C++
code path (``ProteinIdentificationArrowIO_test``) passes in ctest on Windows, and
the zero-copy Arrow bridge + the proteins/groups Parquet exports pass in the same
wheel run. The crash only happens when **pyarrow's Arrow runtime is co-loaded**
with OpenMS's own (contrib) Arrow runtime, which strongly suggests an ABI /
process-global-state interaction rather than a logic bug in the export.

This module reproduces the scenario in **fresh subprocesses** so that a crash in
one scenario does not abort the others. We inspect each subprocess return code to
build a truth table that pins:

  * does the crash need pyarrow loaded at all?
  * is it specific to the search_params export (vs proteins / groups)?
  * does pyarrow import order (before vs after pyopenms) matter?
  * does running the full in-process test sequence reproduce it?

The single test prints the table (as a warning, visible on pass) and asserts that
every scenario exited cleanly (visible on fail). This file is named ``test_aaa_*``
so it runs *before* ``test_arrow_io_classes.py`` — otherwise the in-process crash
would kill the session before this diagnostic is collected.

Remove once the root cause is identified and fixed.
"""

import os
import subprocess
import sys
import textwrap
import warnings

import pytest

# Windows STATUS_HEAP_CORRUPTION, as reported by subprocess.returncode
# (unsigned 3221226356 == signed -1073740940).
_HEAP_CORRUPTION = (3221226356, -1073740940)

# Minimal equivalent of the ``protein_identifications`` fixture from
# test_arrow_io_classes.py, emitted as source so it can run in a subprocess.
_FIXTURE_SRC = textwrap.dedent(
    """
    prot_id = oms.ProteinIdentification()
    prot_id.setIdentifier("PI_0")
    prot_id.setSearchEngine("Comet")
    prot_id.setScoreType("q-value")
    sp = oms.SearchParameters()
    sp.db = "uniprot_human.fasta"
    sp.missed_cleavages = 2
    prot_id.setSearchParameters(sp)
    ph1 = oms.ProteinHit(); ph1.setAccession("PROT_A"); ph1.setScore(0.99)
    ph2 = oms.ProteinHit(); ph2.setAccession("PROT_B"); ph2.setScore(0.95)
    prot_id.setHits([ph1, ph2])
    pg = oms.ProteinGroup(); pg.probability = 0.01
    pg.accessions = [b"PROT_A", b"PROT_B"]
    prot_id.insertProteinGroup(pg)
    prot_ids = [prot_id]
    """
)

_EXPORT_SRC = {
    "proteins": 'oms.ProteinIdentificationArrowIO.exportProteinsToParquet('
                'prot_ids, os.path.join(tmp, "proteins.parquet"))',
    "groups": 'oms.ProteinIdentificationArrowIO.exportProteinGroupsToParquet('
              'prot_ids, os.path.join(tmp, "groups.parquet"))',
    "search_params": 'oms.ProteinIdentificationArrowIO.exportSearchParamsToParquet('
                     'prot_ids, os.path.join(tmp, "params.parquet"))',
}


def _build_snippet(export, *, pyarrow_before, pyarrow_after, do_proteins_groups_first):
    lines = ["import os, tempfile"]
    if pyarrow_before:
        lines.append("import pyarrow")
    lines.append("import pyopenms as oms")
    if pyarrow_after:
        lines.append("import pyarrow")
    lines.append("tmp = tempfile.mkdtemp()")
    lines.append(_FIXTURE_SRC)
    if do_proteins_groups_first:
        lines.append("assert " + _EXPORT_SRC["proteins"])
        lines.append("assert " + _EXPORT_SRC["groups"])
    lines.append("assert " + _EXPORT_SRC[export])
    lines.append("print('DIAG_EXPORT_OK')")
    return "\n".join(lines)


def _run(snippet):
    proc = subprocess.run(
        [sys.executable, "-c", snippet],
        capture_output=True,
        text=True,
        timeout=300,
    )
    return proc


def _classify(rc):
    if rc == 0:
        return "OK"
    if rc in _HEAP_CORRUPTION:
        return "HEAP_CORRUPTION(0xC0000374)"
    return f"FAIL(rc={rc})"


# Each scenario: (label, snippet)
def _scenarios():
    yield (
        "search_params | no pyarrow",
        _build_snippet("search_params", pyarrow_before=False, pyarrow_after=False,
                       do_proteins_groups_first=False),
    )
    yield (
        "search_params | pyarrow BEFORE pyopenms",
        _build_snippet("search_params", pyarrow_before=True, pyarrow_after=False,
                       do_proteins_groups_first=False),
    )
    yield (
        "search_params | pyarrow AFTER pyopenms",
        _build_snippet("search_params", pyarrow_before=False, pyarrow_after=True,
                       do_proteins_groups_first=False),
    )
    yield (
        "proteins (control) | pyarrow BEFORE pyopenms",
        _build_snippet("proteins", pyarrow_before=True, pyarrow_after=False,
                       do_proteins_groups_first=False),
    )
    yield (
        "groups (control) | pyarrow BEFORE pyopenms",
        _build_snippet("groups", pyarrow_before=True, pyarrow_after=False,
                       do_proteins_groups_first=False),
    )
    yield (
        "proteins+groups+search_params | pyarrow BEFORE pyopenms",
        _build_snippet("search_params", pyarrow_before=True, pyarrow_after=False,
                       do_proteins_groups_first=True),
    )


def test_windows_arrow_runtime_diagnostic():
    """Subprocess truth table for the search_params heap corruption."""
    pytest.importorskip("pyarrow")
    oms = pytest.importorskip("pyopenms")
    if not hasattr(oms, "ProteinIdentificationArrowIO"):
        pytest.skip("Arrow IO classes not available")

    import pyarrow

    rows = [
        f"platform={sys.platform} python={sys.version.split()[0]} "
        f"pyarrow={pyarrow.__version__} pyopenms={oms.__version__}"
    ]
    failures = []
    for label, snippet in _scenarios():
        proc = _run(snippet)
        verdict = _classify(proc.returncode)
        rows.append(f"  [{verdict:<28}] {label}")
        if verdict != "OK":
            tail = (proc.stderr or "").strip().splitlines()[-3:]
            rows.append(f"        stderr: {' | '.join(tail)}")
            failures.append((label, verdict))

    table = "\n".join(["", "=== Windows Arrow runtime diagnostic ===", *rows, ""])
    # Visible even when the test passes (pytest prints the warnings summary).
    warnings.warn(UserWarning(table))

    assert not failures, table
