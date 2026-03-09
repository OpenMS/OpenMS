"""Tests for Mobilogram zero-copy structured array access."""
import numpy as np
import pyopenms


def test_mobilogram_get_peaks_struct():
    """get_peaks_struct() should return zero-copy structured array {mobility: float64, intensity: float32}."""
    mob = pyopenms.Mobilogram()
    for i in range(5):
        p = pyopenms.MobilityPeak1D()
        p.setMobility(float(i) * 0.1)
        p.setIntensity(float(i) * 100.0)
        mob.push_back(p)

    arr = mob.get_peaks_struct()
    assert arr.dtype['mobility'] == np.float64
    assert arr.dtype['intensity'] == np.float32
    assert len(arr) == 5
    np.testing.assert_allclose(arr['mobility'], [0.0, 0.1, 0.2, 0.3, 0.4])
    np.testing.assert_allclose(arr['intensity'], [0.0, 100.0, 200.0, 300.0, 400.0], rtol=1e-6)


def test_mobilogram_zero_copy_modification():
    """Modifications through struct view should reflect in original Mobilogram."""
    mob = pyopenms.Mobilogram()
    for i in range(3):
        p = pyopenms.MobilityPeak1D()
        p.setMobility(float(i))
        p.setIntensity(float(i * 10))
        mob.push_back(p)

    arr = mob.get_peaks_struct()
    arr['intensity'][1] = np.float32(999.0)
    assert mob[1].getIntensity() == np.float32(999.0), "Zero-copy write should be reflected in original"


def test_mobilogram_empty():
    """get_peaks_struct() on empty Mobilogram should return empty array."""
    mob = pyopenms.Mobilogram()
    arr = mob.get_peaks_struct()
    assert len(arr) == 0
    assert arr.dtype['mobility'] == np.float64
    assert arr.dtype['intensity'] == np.float32


def test_mobilogram_memory_safety():
    """Struct view should remain valid when Mobilogram is kept alive."""
    def get_view():
        mob = pyopenms.Mobilogram()
        p = pyopenms.MobilityPeak1D()
        p.setMobility(1.5)
        p.setIntensity(42.0)
        mob.push_back(p)
        return mob.get_peaks_struct(), mob  # return mob to keep it alive

    arr, mob = get_view()
    assert arr['mobility'][0] == 1.5
    assert arr['intensity'][0] == np.float32(42.0)
