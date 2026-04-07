"""
Tests for MSSpectrum class bindings.
"""

import pytest
import numpy as np


def test_msspectrum_creation():
    """Test MSSpectrum constructor."""
    from pyopenms import MSSpectrum

    # Default constructor
    spec = MSSpectrum()
    assert spec.size() == 0
    assert len(spec) == 0

    # Copy constructor
    spec.setRT(123.45)
    spec.setMSLevel(2)
    spec2 = MSSpectrum(spec)
    assert spec2.getRT() == pytest.approx(123.45)
    assert spec2.getMSLevel() == 2


def test_msspectrum_basic_properties():
    """Test MSSpectrum basic property getters/setters."""
    from pyopenms import MSSpectrum

    spec = MSSpectrum()

    # RT
    spec.setRT(456.78)
    assert spec.getRT() == pytest.approx(456.78)

    # MS Level
    spec.setMSLevel(3)
    assert spec.getMSLevel() == 3

    # Name
    spec.setName("test_spectrum")
    assert spec.getName() == "test_spectrum"

    # Drift time
    spec.setDriftTime(25.5)
    assert spec.getDriftTime() == pytest.approx(25.5)


def test_msspectrum_peaks():
    """Test MSSpectrum peak operations."""
    from pyopenms import MSSpectrum, Peak1D

    spec = MSSpectrum()

    # Add peaks using push_back
    for i in range(5):
        p = Peak1D()
        p.setMZ(100.0 + i * 10)
        p.setIntensity(1000.0 + i * 100)
        spec.push_back(p)

    assert spec.size() == 5
    assert len(spec) == 5

    # Access individual peaks
    assert spec[0].getMZ() == pytest.approx(100.0)
    assert spec[4].getMZ() == pytest.approx(140.0)

    # Index out of range
    with pytest.raises(IndexError):
        _ = spec[10]


def test_msspectrum_get_set_peaks():
    """Test MSSpectrum get_peaks/set_peaks methods."""
    from pyopenms import MSSpectrum

    spec = MSSpectrum()

    # Set peaks from numpy arrays
    mz = np.array([100.0, 200.0, 300.0, 400.0, 500.0])
    intensity = np.array([1000.0, 2000.0, 3000.0, 2000.0, 1000.0])

    spec.set_peaks(mz, intensity)

    assert spec.size() == 5

    # Get peaks back
    mz_out, int_out = spec.get_peaks()

    np.testing.assert_array_almost_equal(mz_out, mz)
    np.testing.assert_array_almost_equal(int_out, intensity, decimal=3)


def test_msspectrum_get_set_peaks_tuple():
    """Test MSSpectrum set_peaks with tuple argument."""
    from pyopenms import MSSpectrum

    spec = MSSpectrum()

    mz = np.array([100.0, 200.0, 300.0])
    intensity = np.array([1000.0, 2000.0, 3000.0])

    # Set using tuple
    spec.set_peaks((mz, intensity))

    mz_out, int_out = spec.get_peaks()
    np.testing.assert_array_almost_equal(mz_out, mz)


def test_msspectrum_get_peaks_dtypes():
    """Test get_peaks returns correct dtypes."""
    from pyopenms import MSSpectrum

    spec = MSSpectrum()
    mz_exp = [100.0, 200.0, 300.0]
    int_exp = [1000.0, 2000.0, 500.0]
    spec.set_peaks((mz_exp, int_exp))

    mz, intensities = spec.get_peaks()

    assert len(mz) == 3
    assert len(intensities) == 3
    assert mz.dtype == np.float64
    assert intensities.dtype == np.float32


def test_msspectrum_get_peaks_empty():
    """Test get_peaks on empty spectrum."""
    from pyopenms import MSSpectrum

    spec = MSSpectrum()
    mz, intensities = spec.get_peaks()

    assert len(mz) == 0
    assert len(intensities) == 0
    assert mz.dtype == np.float64
    assert intensities.dtype == np.float32


def test_msspectrum_get_mz_array():
    """Test get_mz_array method."""
    from pyopenms import MSSpectrum

    spec = MSSpectrum()
    mz_exp = np.array([100.0, 200.0, 300.0])
    int_exp = np.array([1000.0, 2000.0, 500.0])
    spec.set_peaks((mz_exp, int_exp))

    mz = spec.get_mz_array()
    assert len(mz) == 3
    assert mz.dtype == np.float64
    np.testing.assert_array_almost_equal(mz, mz_exp)


def test_msspectrum_get_intensity_array():
    """Test get_intensity_array method."""
    from pyopenms import MSSpectrum

    spec = MSSpectrum()
    mz_exp = np.array([100.0, 200.0, 300.0])
    int_exp = np.array([1000.0, 2000.0, 500.0])
    spec.set_peaks((mz_exp, int_exp))

    intensities = spec.get_intensity_array()
    assert len(intensities) == 3
    assert intensities.dtype == np.float32
    np.testing.assert_array_almost_equal(intensities, int_exp, decimal=1)


def test_msspectrum_update_ranges():
    """Test MSSpectrum updateRanges and range getters."""
    from pyopenms import MSSpectrum, Peak1D

    spec = MSSpectrum()
    p = Peak1D()
    p.setMZ(500.0)
    p.setIntensity(1e5)
    spec.push_back(p)

    spec.updateRanges()
    assert isinstance(spec.getMinMZ(), float)
    assert isinstance(spec.getMaxMZ(), float)
    assert isinstance(spec.getMinIntensity(), float)
    assert isinstance(spec.getMaxIntensity(), float)
    assert spec.getMinIntensity() == pytest.approx(1e5)
    assert spec.getMaxIntensity() == pytest.approx(1e5)


def test_msspectrum_iteration():
    """Test MSSpectrum iteration."""
    from pyopenms import MSSpectrum

    spec = MSSpectrum()
    mz = np.array([100.0, 200.0, 300.0])
    intensity = np.array([1000.0, 2000.0, 3000.0])
    spec.set_peaks(mz, intensity)

    # Iterate over peaks
    peaks = list(spec)
    assert len(peaks) == 3

    for i, peak in enumerate(spec):
        assert peak.getMZ() == pytest.approx(mz[i])


def test_msspectrum_sorting():
    """Test MSSpectrum sorting methods."""
    from pyopenms import MSSpectrum

    spec = MSSpectrum()
    mz = np.array([300.0, 100.0, 200.0])
    intensity = np.array([3000.0, 1000.0, 2000.0])
    spec.set_peaks(mz, intensity)

    # Should not be sorted initially
    assert not spec.isSorted()

    # Sort by position (m/z)
    spec.sortByPosition()
    assert spec.isSorted()

    mz_sorted, _ = spec.get_peaks()
    np.testing.assert_array_almost_equal(mz_sorted, [100.0, 200.0, 300.0])


def test_msspectrum_find_nearest():
    """Test MSSpectrum findNearest methods."""
    from pyopenms import MSSpectrum

    spec = MSSpectrum()
    mz = np.array([100.0, 200.0, 300.0, 400.0, 500.0])
    intensity = np.array([1000.0, 2000.0, 3000.0, 2000.0, 1000.0])
    spec.set_peaks(mz, intensity)
    spec.sortByPosition()

    # Find nearest
    idx = spec.findNearest(195.0)
    assert idx == 1  # 200.0 is nearest

    # Find nearest with tolerance
    idx = spec.findNearest(195.0, 10.0)
    assert idx == 1

    # No match within tolerance
    idx = spec.findNearest(195.0, 1.0)
    assert idx == -1


def test_msspectrum_tic():
    """Test MSSpectrum TIC calculation."""
    from pyopenms import MSSpectrum

    spec = MSSpectrum()
    mz = np.array([100.0, 200.0, 300.0])
    intensity = np.array([1000.0, 2000.0, 3000.0])
    spec.set_peaks(mz, intensity)

    tic = spec.calculateTIC()
    assert tic == pytest.approx(6000.0)


def test_msspectrum_comparison():
    """Test MSSpectrum comparison operators."""
    from pyopenms import MSSpectrum

    spec1 = MSSpectrum()
    spec1.setRT(100.0)

    spec2 = MSSpectrum()
    spec2.setRT(100.0)

    spec3 = MSSpectrum()
    spec3.setRT(200.0)

    # Same RT should be equal (before adding peaks)
    assert spec1 == spec2
    assert spec1 != spec3


def test_msspectrum_repr():
    """Test MSSpectrum string representation."""
    from pyopenms import MSSpectrum

    spec = MSSpectrum()
    spec.setRT(123.45)
    spec.setMSLevel(2)
    mz = np.array([100.0, 200.0, 300.0])
    intensity = np.array([1000.0, 2000.0, 3000.0])
    spec.set_peaks(mz, intensity)

    repr_str = repr(spec)
    assert "MSSpectrum" in repr_str
    assert "3" in repr_str  # n_peaks
    assert "123" in repr_str  # RT


def test_msspectrum_reserve_resize():
    """Test MSSpectrum reserve and resize."""
    from pyopenms import MSSpectrum

    spec = MSSpectrum()

    # Reserve space (should not change size)
    spec.reserve(100)
    assert spec.size() == 0

    # Resize
    spec.resize(10)
    assert spec.size() == 10


def test_msspectrum_clear():
    """Test MSSpectrum clear method."""
    from pyopenms import MSSpectrum

    spec = MSSpectrum()
    mz = np.array([100.0, 200.0, 300.0])
    intensity = np.array([1000.0, 2000.0, 3000.0])
    spec.set_peaks(mz, intensity)
    spec.setRT(123.0)
    spec.setMSLevel(2)

    assert spec.size() == 3

    # Clear without clearing metadata
    spec.clear(False)
    assert spec.size() == 0
    assert spec.getRT() == pytest.approx(123.0)

    # Restore and clear with metadata
    spec.set_peaks(mz, intensity)
    spec.clear(True)
    assert spec.size() == 0


def test_msspectrum_drift_time_no_im():
    """Test drift time methods when no IM data present."""
    from pyopenms import MSSpectrum

    spec = MSSpectrum()
    spec.set_peaks(([100.0, 200.0], [1000.0, 2000.0]))

    assert not spec.containsIMData()
    assert spec.get_drift_time_array() is None
    assert spec.get_drift_time_array_view() is None
    assert spec.get_drift_time_unit() is None


def test_msspectrum_drift_time_with_im():
    """Test drift time methods with ion mobility data."""
    from pyopenms import MSSpectrum, FloatDataArray

    spec = MSSpectrum()
    spec.set_peaks(([100.0, 200.0, 300.0], [1000.0, 2000.0, 500.0]))

    # Add ion mobility data via FloatDataArray
    fda = FloatDataArray()
    fda.setName("Ion Mobility")
    fda.push_back(1.5)
    fda.push_back(2.5)
    fda.push_back(3.5)
    spec.setFloatDataArrays([fda])

    assert spec.containsIMData()

    drift = spec.get_drift_time_array()
    assert drift is not None
    assert len(drift) == 3
    assert drift.dtype == np.float32
    assert drift[0] == pytest.approx(1.5, abs=0.1)
    assert drift[1] == pytest.approx(2.5, abs=0.1)
    assert drift[2] == pytest.approx(3.5, abs=0.1)

    assert spec.get_drift_time_unit() is not None


def test_float_data_array_get_data():
    """Test FloatDataArray get_data method (returns copy - safe)."""
    from pyopenms import FloatDataArray

    fda = FloatDataArray()
    fda.push_back(1.0)
    fda.push_back(2.0)
    fda.push_back(3.0)

    # Get a copy (safe default)
    data_copy = fda.get_data()
    assert len(data_copy) == 3
    assert data_copy.dtype == np.float32
    assert data_copy[0] == pytest.approx(1.0, abs=0.1)
    assert data_copy[1] == pytest.approx(2.0, abs=0.1)
    assert data_copy[2] == pytest.approx(3.0, abs=0.1)

    # Modify copy - original should be unchanged
    data_copy[0] = 100.0
    assert fda[0] == pytest.approx(1.0, abs=0.1)


def test_float_data_array_get_data_view():
    """Test FloatDataArray get_data_view method (memory view - fast, unsafe)."""
    from pyopenms import FloatDataArray

    fda = FloatDataArray()
    fda.push_back(1.0)
    fda.push_back(2.0)
    fda.push_back(3.0)

    # Get a view (fast but unsafe)
    data_view = fda.get_data_view()
    assert len(data_view) == 3
    assert data_view[0] == pytest.approx(1.0, abs=0.1)
    assert data_view[1] == pytest.approx(2.0, abs=0.1)
    assert data_view[2] == pytest.approx(3.0, abs=0.1)

    # Modify view - original SHOULD change (it's a view)
    data_view[0] = 100.0
    assert fda[0] == pytest.approx(100.0, abs=0.1)


def test_float_data_array_empty():
    """Test FloatDataArray methods on empty array."""
    from pyopenms import FloatDataArray

    fda = FloatDataArray()

    data_copy = fda.get_data()
    assert len(data_copy) == 0

    data_view = fda.get_data_view()
    assert isinstance(data_view, np.ndarray)
    assert len(data_view) == 0


def test_msspectrum_get_type_with_query_data():
    """Test getType(bool) method with query_data parameter."""
    from pyopenms import MSSpectrum, SpectrumSettings

    spec = MSSpectrum()

    # Set type explicitly
    spec.setType(SpectrumSettings.SpectrumType.CENTROID)
    assert spec.getType(False) == SpectrumSettings.SpectrumType.CENTROID

    # Test with query_data=True on spectrum with peaks
    spec2 = MSSpectrum()
    spec2.set_peaks(([100.0, 200.0, 300.0], [1000.0, 2000.0, 500.0]))

    # Without setting type, getType(False) should return UNKNOWN
    assert spec2.getType(False) == SpectrumSettings.SpectrumType.UNKNOWN

    # With query_data=True, it may estimate the type based on peak data
    spec_type_with_query = spec2.getType(True)
    assert spec_type_with_query in [
        SpectrumSettings.SpectrumType.UNKNOWN,
        SpectrumSettings.SpectrumType.CENTROID,
        SpectrumSettings.SpectrumType.PROFILE,
    ]
