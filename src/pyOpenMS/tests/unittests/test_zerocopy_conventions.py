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
    arr = fda.data_view()
    assert isinstance(arr, np.ndarray)
    assert len(arr) == 0
    assert arr.dtype == np.float32

    # IntegerDataArray
    ida = pyopenms.IntegerDataArray()
    arr = ida.data_view()
    assert isinstance(arr, np.ndarray)
    assert len(arr) == 0
    assert arr.dtype == np.int32

    # MSSpectrum drift time (no IM data — None is correct here, data doesn't exist)
    spec = pyopenms.MSSpectrum()
    assert spec.drift_time_array_view() is None

    # MatrixDouble (empty)
    mat = pyopenms.MatrixDouble()
    arr = mat.matrix_view()
    assert isinstance(arr, np.ndarray)
    assert arr.size == 0


def test_view_returns_typed_array():
    """All _view methods return typed ndarray (not structured, not raw bytes)."""
    fda = pyopenms.FloatDataArray()
    fda.push_back(1.0)
    fda.push_back(2.0)
    arr = fda.data_view()
    assert arr.dtype == np.float32
    assert arr.shape == (2,)
    # Direct arithmetic works on _view arrays
    result = arr * 2.0
    assert np.isclose(result[0], 2.0)


def test_struct_returns_empty_array_on_empty():
    """All _struct methods return empty structured array (not None) when container is empty."""
    spec = pyopenms.MSSpectrum()
    arr = spec.peaks_struct()
    assert isinstance(arr, np.ndarray)
    assert len(arr) == 0
    assert 'mz' in arr.dtype.names

    chrom = pyopenms.MSChromatogram()
    arr = chrom.peaks_struct()
    assert isinstance(arr, np.ndarray)
    assert len(arr) == 0
    assert 'rt' in arr.dtype.names

    mob = pyopenms.Mobilogram()
    arr = mob.peaks_struct()
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
    arr = spec.peaks_struct()
    assert arr.dtype.names is not None
    assert 'mz' in arr.dtype.names
    assert 'intensity' in arr.dtype.names


def test_struct_dtype_matches_cpp_peak_layout():
    """_struct dtypes must describe the real AoS layout, not a plausible-looking one.

    The dtype overlays Peak1D/ChromatogramPeak/MobilityPeak1D storage directly, so a
    wrong offset silently reinterprets bytes rather than failing. C++ static_asserts in
    bindings/peak_layout.h pin the offsets at build time; this pins what reaches Python.
    """
    spec = pyopenms.MSSpectrum()
    spec.set_peaks(([100.5, 200.25], [1.0, 2.0]))
    chrom = pyopenms.MSChromatogram()
    chrom.set_peaks(([10.5, 20.25], [3.0, 4.0]))
    mob = pyopenms.Mobilogram()
    mob.set_peaks(([0.5, 0.75], [5.0, 6.0]))

    for obj, pos_field in [(spec, 'mz'), (chrom, 'rt'), (mob, 'mobility')]:
        dt = obj.peaks_struct().dtype
        assert dt.names == (pos_field, 'intensity')
        assert dt.itemsize == 16
        assert dt.fields[pos_field][1] == 0
        assert dt.fields['intensity'][1] == 8
        assert dt.fields[pos_field][0] == np.float64
        assert dt.fields['intensity'][0] == np.float32


def test_struct_view_reads_the_same_values_as_the_scalar_accessors():
    """End-to-end check that the dtype offsets land on the right members.

    Values are chosen so a swapped or shifted field would not coincidentally match.
    """
    spec = pyopenms.MSSpectrum()
    spec.set_peaks(([100.5, 200.25, 300.125], [11.0, 22.0, 33.0]))
    arr = spec.peaks_struct()
    assert [p.getMZ() for p in spec] == list(arr['mz'])
    assert [p.getIntensity() for p in spec] == list(arr['intensity'])

    chrom = pyopenms.MSChromatogram()
    chrom.set_peaks(([10.5, 20.25, 30.125], [44.0, 55.0, 66.0]))
    arr = chrom.peaks_struct()
    assert [p.getRT() for p in chrom] == list(arr['rt'])
    assert [p.getIntensity() for p in chrom] == list(arr['intensity'])

    mob = pyopenms.Mobilogram()
    mob.set_peaks(([0.5, 0.75, 1.25], [77.0, 88.0, 99.0]))
    arr = mob.peaks_struct()
    assert [p.getMobility() for p in mob] == list(arr['mobility'])
    assert [p.getIntensity() for p in mob] == list(arr['intensity'])


def test_no_public_as_view_methods():
    """No public methods should use _as_view suffix (removed convention)."""
    # get_matrix_as_view is an intentional deprecated alias for matrix_view
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
    arr = fda.data_view()
    arr[0] = 42.0
    assert fda[0] == 42.0


def test_struct_writeback():
    """All _struct methods return writable views (modifications affect C++ object)."""
    spec = pyopenms.MSSpectrum()
    p = pyopenms.Peak1D()
    p.setMZ(100.0)
    p.setIntensity(50.0)
    spec.push_back(p)
    arr = spec.peaks_struct()
    arr['intensity'][0] = np.float32(999.0)
    assert spec[0].getIntensity() == np.float32(999.0)

    chrom = pyopenms.MSChromatogram()
    cp = pyopenms.ChromatogramPeak()
    cp.setRT(1.0)
    cp.setIntensity(50.0)
    chrom.push_back(cp)
    arr = chrom.peaks_struct()
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
        assert "data_view" in str(w[0].message)
    assert arr is not None
    assert arr.dtype == np.float32

    # MatrixDouble
    mat = pyopenms.MatrixDouble(2, 2, 1.0)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        view = mat.get_matrix_mv()
        assert len(w) == 1
        assert issubclass(w[0].category, DeprecationWarning)
        assert "matrix_view" in str(w[0].message)
    assert view is not None
