"""Verify zero-copy API naming conventions across all pyOpenMS types.

Convention:
  _view  = typed zero-copy column view, returns empty ndarray<T> if empty
  _struct = structured zero-copy record view, returns empty structured ndarray if empty

No public methods should use _as_view suffix.
"""
import numpy as np
import pyopenms


def test_view_returns_empty_array_on_empty():
    """All _view methods return empty arrays (not None) when container is empty."""
    # FloatDataArray
    fda = pyopenms.FloatDataArray()
    arr = fda.get_data_view()
    assert isinstance(arr, np.ndarray)
    assert len(arr) == 0
    assert arr.dtype == np.float32

    # IntegerDataArray
    ida = pyopenms.IntegerDataArray()
    arr = ida.get_data_view()
    assert isinstance(arr, np.ndarray)
    assert len(arr) == 0
    assert arr.dtype == np.int32

    # MSSpectrum drift time (no IM data — None is correct here, data doesn't exist)
    spec = pyopenms.MSSpectrum()
    assert spec.get_drift_time_array_view() is None

    # MatrixDouble (empty)
    mat = pyopenms.MatrixDouble()
    arr = mat.get_matrix_view()
    assert isinstance(arr, np.ndarray)
    assert arr.size == 0


def test_view_returns_typed_array():
    """All _view methods return typed ndarray (not structured, not raw bytes)."""
    fda = pyopenms.FloatDataArray()
    fda.push_back(1.0)
    fda.push_back(2.0)
    arr = fda.get_data_view()
    assert arr.dtype == np.float32
    assert arr.shape == (2,)
    # Direct arithmetic works on _view arrays
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


def test_no_public_as_view_methods():
    """No public methods should use _as_view suffix (removed convention)."""
    # get_matrix_as_view is an intentional deprecated alias for get_matrix_view
    deprecated_aliases = {
        'MatrixDouble': {'get_matrix_as_view'},
    }
    for cls_name in ['MSSpectrum', 'MSChromatogram', 'Mobilogram', 'MatrixDouble',
                     'OSBinaryDataArray', 'OSSpectrum', 'OSChromatogram']:
        cls = getattr(pyopenms, cls_name)
        allowed = deprecated_aliases.get(cls_name, set())
        public_methods = [m for m in dir(cls) if not m.startswith('_')]
        bad_methods = [m for m in public_methods if '_as_view' in m and m not in allowed]
        assert bad_methods == [], \
            f"{cls_name} has public _as_view methods: {bad_methods}. " \
            f"Use _view suffix instead."


def test_view_writeback():
    """All _view methods return writable views (modifications affect C++ object)."""
    fda = pyopenms.FloatDataArray()
    fda.push_back(1.0)
    arr = fda.get_data_view()
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


def test_deprecated_mv_aliases_emit_warning():
    """Deprecated _mv aliases should still work but emit DeprecationWarning."""
    import warnings

    # FloatDataArray
    fda = pyopenms.FloatDataArray()
    fda.push_back(1.0)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        arr = fda.get_data_mv()
        assert len(w) == 1
        assert issubclass(w[0].category, DeprecationWarning)
        assert "get_data_view" in str(w[0].message)
    assert arr is not None
    assert arr.dtype == np.float32

    # MatrixDouble
    mat = pyopenms.MatrixDouble(2, 2, 1.0)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        view = mat.get_matrix_mv()
        assert len(w) == 1
        assert issubclass(w[0].category, DeprecationWarning)
        assert "get_matrix_view" in str(w[0].message)
    assert view is not None
