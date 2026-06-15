"""
Lifetime tests for the numpy peak-data getters (OpenMS issue #9460, Section 2).

MSSpectrum.get_peaks() and MSExperiment.get2DPeakData()/get2DPeakDataLong() are
nanobind wrappers that return their peak data as numpy arrays whose buffer is
owned by an nb::capsule (delete[]), independent of the C++ source object. They
are safe by construction today; this pins that contract so a future change to a
*view*-returning implementation (sharing the source buffer) is caught.

test_spectrum_zerocopy.py already covers get_peaks_struct(); these cover the
get_peaks()/get2DPeakData() getters named in the issue. The snapshot is taken
*before* deleting the source, so the comparison cannot pass vacuously.
"""

import gc

import pytest

import pyopenms

np = pytest.importorskip("numpy")


def _make_spectrum():
    spec = pyopenms.MSSpectrum()
    spec.setRT(42.0)
    spec.setMSLevel(1)
    spec.set_peaks((np.array([300.0, 301.0, 302.0], dtype=np.float64),
                    np.array([9.0, 8.0, 7.0], dtype=np.float32)))
    return spec


def _make_experiment():
    exp = pyopenms.MSExperiment()
    for i in range(3):
        spec = pyopenms.MSSpectrum()
        spec.setRT(10.0 + i)
        spec.setMSLevel(1)
        spec.set_peaks((np.array([100.0 + j for j in range(4)], dtype=np.float64),
                        np.array([5.0 * (j + 1) for j in range(4)], dtype=np.float32)))
        exp.addSpectrum(spec)
    return exp


def test_msspectrum_get_peaks_survives_owner_gc():
    spec = _make_spectrum()
    mz, intensity = spec.get_peaks()
    assert len(mz) > 0

    # snapshot WHILE the source spectrum is still alive
    mz_snap, in_snap = mz.copy(), intensity.copy()

    del spec
    gc.collect()

    assert np.array_equal(mz, mz_snap), "get_peaks() mz array changed after owner GC"
    assert np.array_equal(intensity, in_snap), "get_peaks() intensity array changed after owner GC"


def test_msexperiment_get2dpeakdata_survives_source_gc():
    exp = _make_experiment()
    rt, mz, intensity = exp.get2DPeakData(0.0, 1.0e9, 0.0, 1.0e9, 1)
    assert len(rt) > 0

    snap = (rt.copy(), mz.copy(), intensity.copy())

    del exp
    gc.collect()

    assert np.array_equal(rt, snap[0]), "get2DPeakData() rt changed after source GC"
    assert np.array_equal(mz, snap[1]), "get2DPeakData() mz changed after source GC"
    assert np.array_equal(intensity, snap[2]), "get2DPeakData() intensity changed after source GC"


def test_msexperiment_get2dpeakdatalong_survives_source_gc():
    exp = _make_experiment()
    rt, mz, intensity = exp.get2DPeakDataLong(0.0, 1.0e9, 0.0, 1.0e9, 1)
    assert len(rt) > 0

    snap = (rt.copy(), mz.copy(), intensity.copy())

    del exp
    gc.collect()

    assert np.array_equal(rt, snap[0]), "get2DPeakDataLong() rt changed after source GC"
    assert np.array_equal(mz, snap[1]), "get2DPeakDataLong() mz changed after source GC"
    assert np.array_equal(intensity, snap[2]), "get2DPeakDataLong() intensity changed after source GC"


if __name__ == "__main__":
    test_msspectrum_get_peaks_survives_owner_gc()
    test_msexperiment_get2dpeakdata_survives_source_gc()
    test_msexperiment_get2dpeakdatalong_survives_source_gc()
    print("All numpy-getter lifetime tests passed!")
