"""Verify zero-copy API naming conventions across all pyOpenMS types.

Convention:
  _mv    = typed zero-copy column view, returns ndarray<T> or None if empty
  _struct = structured zero-copy record view, returns structured ndarray or empty array

No public methods should use _view or _as_view suffixes.
"""
import numpy as np
import pyopenms


def test_mv_returns_none_on_empty():
    """All _mv methods return None when container is empty."""
    # FloatDataArray
    fda = pyopenms.FloatDataArray()
    assert fda.get_data_mv() is None

    # IntegerDataArray
    ida = pyopenms.IntegerDataArray()
    assert ida.get_data_mv() is None

    # MSSpectrum drift time (no IM data)
    spec = pyopenms.MSSpectrum()
    assert spec.get_drift_time_array_mv() is None

    # MatrixDouble (empty)
    mat = pyopenms.MatrixDouble()
    assert mat.get_matrix_mv() is None


def test_mv_returns_typed_array():
    """All _mv methods return typed ndarray (not structured, not raw bytes)."""
    fda = pyopenms.FloatDataArray()
    fda.push_back(1.0)
    fda.push_back(2.0)
    arr = fda.get_data_mv()
    assert arr.dtype == np.float32
    assert arr.shape == (2,)
    # Direct arithmetic works on _mv arrays
    result = arr * 2.0
    assert np.isclose(result[0], 2.0)


def test_struct_returns_empty_array_on_empty():
    """All _struct methods return empty structured array (not None) when container is empty."""
    spec = pyopenms.MSSpectrum()
    arr = spec.get_peaks_struct()
    assert isinstance(arr, np.ndarray)
    assert len(arr) == 0
    assert 'mz' in arr.dtype.names

    chrom = pyopenms.MSChromatogram()
    arr = chrom.get_peaks_struct()
    assert isinstance(arr, np.ndarray)
    assert len(arr) == 0
    assert 'rt' in arr.dtype.names

    mob = pyopenms.Mobilogram()
    arr = mob.get_peaks_struct()
    assert isinstance(arr, np.ndarray)
    assert len(arr) == 0
    assert 'mobility' in arr.dtype.names


def test_struct_returns_structured_array():
    """All _struct methods return structured arrays with named fields."""
    spec = pyopenms.MSSpectrum()
    p = pyopenms.Peak1D()
    p.setMZ(100.0)
    p.setIntensity(50.0)
    spec.push_back(p)
    arr = spec.get_peaks_struct()
    assert arr.dtype.names is not None
    assert 'mz' in arr.dtype.names
    assert 'intensity' in arr.dtype.names


def test_no_public_view_or_as_view_methods():
    """No public methods should use _view or _as_view suffix (they should be _ prefixed or removed)."""
    for cls_name in ['MSSpectrum', 'MSChromatogram', 'Mobilogram', 'MatrixDouble']:
        cls = getattr(pyopenms, cls_name)
        public_methods = [m for m in dir(cls) if not m.startswith('_')]
        view_methods = [m for m in public_methods if m.endswith('_view') or '_as_view' in m]
        assert view_methods == [], \
            f"{cls_name} has public _view/_as_view methods: {view_methods}. " \
            f"These should be prefixed with _ (internal) or removed."


def test_mv_writeback():
    """All _mv methods return writable views (modifications affect C++ object)."""
    fda = pyopenms.FloatDataArray()
    fda.push_back(1.0)
    arr = fda.get_data_mv()
    arr[0] = 42.0
    assert fda[0] == 42.0


def test_struct_writeback():
    """All _struct methods return writable views (modifications affect C++ object)."""
    spec = pyopenms.MSSpectrum()
    p = pyopenms.Peak1D()
    p.setMZ(100.0)
    p.setIntensity(50.0)
    spec.push_back(p)
    arr = spec.get_peaks_struct()
    arr['intensity'][0] = np.float32(999.0)
    assert spec[0].getIntensity() == np.float32(999.0)

    chrom = pyopenms.MSChromatogram()
    cp = pyopenms.ChromatogramPeak()
    cp.setRT(1.0)
    cp.setIntensity(50.0)
    chrom.push_back(cp)
    arr = chrom.get_peaks_struct()
    arr['intensity'][0] = np.float32(999.0)
    _, ints = chrom.get_peaks()
    assert np.isclose(ints[0], np.float32(999.0))
