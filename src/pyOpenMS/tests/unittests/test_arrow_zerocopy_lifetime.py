"""
Lifetime tests for the Arrow zero-copy exporters (OpenMS issue #9460, Section 2).

``pyopenms._arrow_zerocopy.spectra_to_arrow()`` / ``chromatograms_to_arrow()``
hand the peak buffers to PyArrow over the C Data Interface; the buffers' lifetime
rides on Arrow's ``release`` callback, **not** on the source ``MSExperiment``. The
existing ``test_arrow_zerocopy.py`` only checks shapes/columns; this pins the
harder contract: the exported table stays byte-identical after the source object
is deleted and garbage-collected.

Note: the snapshot is taken *before* deleting the source, so the comparison
cannot pass vacuously (comparing corrupted data against itself).
"""

import gc

import pytest

import pyopenms

np = pytest.importorskip("numpy")
pa = pytest.importorskip("pyarrow")

from pyopenms._arrow_zerocopy import spectra_to_arrow, chromatograms_to_arrow


def _make_spectra_experiment():
    exp = pyopenms.MSExperiment()
    for i in range(3):
        spec = pyopenms.MSSpectrum()
        spec.setRT(10.0 + i * 5.0)
        spec.setMSLevel(1)
        spec.setNativeID(f"scan={i}")
        mz = np.array([100.0 + j for j in range(4)], dtype=np.float64)
        intensity = np.array([1000.0 * (j + 1) for j in range(4)], dtype=np.float32)
        spec.set_peaks(mz, intensity)
        exp.addSpectrum(spec)
    return exp


def _make_chromatogram_experiment():
    exp = pyopenms.MSExperiment()
    for i in range(2):
        chrom = pyopenms.MSChromatogram()
        chrom.setNativeID(f"chrom={i}")
        rt = np.array([1.0 + j for j in range(5)], dtype=np.float64)
        intensity = np.array([10.0 * (j + 1) for j in range(5)], dtype=np.float32)
        chrom.set_peaks(rt, intensity)
        exp.addChromatogram(chrom)
    return exp


def _assert_survives_owner_gc(table, snapshot):
    """The table (whose source is already deleted by the caller) must re-read
    byte-identically to the pre-deletion snapshot, even after a forced GC."""
    gc.collect()
    for col in table.column_names:
        assert table.column(col).to_pylist() == snapshot[col], (
            f"column {col!r} changed after the source experiment was GC'd"
        )
    # the owner-less buffers must also still rebuild into a fresh Arrow table
    rebuilt = pa.table({col: table.column(col) for col in table.column_names})
    assert rebuilt.num_rows == table.num_rows


def test_spectra_to_arrow_survives_source_gc():
    exp = _make_spectra_experiment()
    table = spectra_to_arrow(exp, format="long")
    assert table.num_rows > 0

    # snapshot WHILE the source experiment is still alive
    snapshot = {col: table.column(col).to_pylist() for col in table.column_names}

    del exp  # drop the source; the Arrow buffers must outlive it
    _assert_survives_owner_gc(table, snapshot)


def test_chromatograms_to_arrow_survives_source_gc():
    exp = _make_chromatogram_experiment()
    table = chromatograms_to_arrow(exp, format="long")
    assert table.num_rows > 0

    snapshot = {col: table.column(col).to_pylist() for col in table.column_names}

    del exp
    _assert_survives_owner_gc(table, snapshot)


if __name__ == "__main__":
    test_spectra_to_arrow_survives_source_gc()
    test_chromatograms_to_arrow_survives_source_gc()
    print("All arrow zero-copy lifetime tests passed!")
