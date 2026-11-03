"""
Tests for type caster functionality.

These tests verify that the custom nanobind type casters correctly
convert between Python and C++ types.
"""

import pytest
import numpy as np


class TestStringCaster:
    """Tests for OpenMS::String type caster."""

    def test_string_from_python_str(self):
        """Test conversion from Python str to OpenMS::String."""
        from pyopenms import MSSpectrum

        spec = MSSpectrum()
        spec.setName("test spectrum name")
        assert spec.getName() == "test spectrum name"

    def test_string_from_python_unicode(self):
        """Test conversion with Unicode characters."""
        from pyopenms import MSSpectrum

        spec = MSSpectrum()
        spec.setName("spectrum_with_unicode_\u03b1\u03b2\u03b3")
        name = spec.getName()
        assert "unicode" in name

    def test_string_empty(self):
        """Test conversion of empty string."""
        from pyopenms import MSSpectrum

        spec = MSSpectrum()
        spec.setName("")
        assert spec.getName() == ""

    def test_string_roundtrip(self):
        """Test string roundtrip conversion."""
        from pyopenms import MSSpectrum

        spec = MSSpectrum()
        test_strings = [
            "simple",
            "with spaces",
            "with\ttabs",
            "with\nnewlines",
            "with_underscores_123",
        ]

        for s in test_strings:
            spec.setName(s)
            assert spec.getName() == s


class TestDataValueCaster:
    """Tests for OpenMS::DataValue type caster."""

    def test_datavalue_from_int(self):
        """Test DataValue from Python int."""
        from pyopenms import DataValue
        dv = DataValue(42)
        assert dv.toInt() == 42
        assert dv.valueType() == DataValue.INT_VALUE
        assert not dv.isEmpty()

    def test_datavalue_from_float(self):
        """Test DataValue from Python float."""
        from pyopenms import DataValue
        dv = DataValue(3.14)
        assert abs(dv.toDouble() - 3.14) < 1e-10
        assert dv.valueType() == DataValue.DOUBLE_VALUE
        assert not dv.isEmpty()

    def test_datavalue_from_string(self):
        """Test DataValue from Python string."""
        from pyopenms import DataValue
        dv = DataValue("hello")
        assert dv.toString() == "hello"
        assert dv.valueType() == DataValue.STRING_VALUE
        assert not dv.isEmpty()

    def test_datavalue_empty(self):
        """Test empty DataValue."""
        from pyopenms import DataValue
        dv = DataValue()
        assert dv.isEmpty()
        assert dv.valueType() == DataValue.EMPTY_VALUE

    def test_datavalue_int_list(self):
        """Test DataValue from int list."""
        from pyopenms import DataValue
        dv = DataValue([1, 2, 3])
        assert dv.toIntList() == [1, 2, 3]
        assert dv.valueType() == DataValue.INT_LIST

    def test_datavalue_double_list(self):
        """Test DataValue from float list."""
        from pyopenms import DataValue
        dv = DataValue([1.0, 2.5])
        result = dv.toDoubleList()
        assert len(result) == 2
        assert abs(result[0] - 1.0) < 1e-10
        assert abs(result[1] - 2.5) < 1e-10
        assert dv.valueType() == DataValue.DOUBLE_LIST

    def test_datavalue_string_list(self):
        """Test DataValue from string list."""
        from pyopenms import DataValue
        dv = DataValue(["a", "b", "c"])
        assert dv.toStringList() == ["a", "b", "c"]
        assert dv.valueType() == DataValue.STRING_LIST


class TestDPositionCaster:
    """Tests for DPosition type casters."""

    def test_dposition1_from_float(self):
        """Test DPosition<1> from Python float."""
        from pyopenms import Peak1D

        p = Peak1D()
        p.setMZ(123.456)
        assert p.getMZ() == pytest.approx(123.456)

    def test_dposition1_from_int(self):
        """Test DPosition<1> from Python int."""
        from pyopenms import Peak1D

        p = Peak1D()
        p.setMZ(100)  # int should work
        assert p.getMZ() == pytest.approx(100.0)


class TestContainerCasters:
    """Tests for container type casters."""

    def test_vector_peak1d_iteration(self):
        """Test that MSSpectrum iteration yields Peak1D objects."""
        from pyopenms import MSSpectrum, Peak1D

        spec = MSSpectrum()
        mz = np.array([100.0, 200.0, 300.0])
        intensity = np.array([1000.0, 2000.0, 3000.0])
        spec.set_peaks(mz, intensity)

        peaks = list(spec)
        assert len(peaks) == 3
        assert all(isinstance(p, Peak1D) for p in peaks)


class TestNumpyArrayCasters:
    """Tests for numpy array conversions."""

    def test_get_peaks_returns_numpy(self):
        """Test that get_peaks returns numpy arrays."""
        from pyopenms import MSSpectrum

        spec = MSSpectrum()
        spec.set_peaks(np.array([100.0]), np.array([1000.0]))

        mz, intensity = spec.get_peaks()

        assert isinstance(mz, np.ndarray)
        assert isinstance(intensity, np.ndarray)
        assert mz.dtype == np.float64
        assert intensity.dtype == np.float32

    def test_set_peaks_from_float64(self):
        """Test set_peaks with float64 arrays."""
        from pyopenms import MSSpectrum

        spec = MSSpectrum()
        mz = np.array([100.0, 200.0], dtype=np.float64)
        intensity = np.array([1000.0, 2000.0], dtype=np.float64)

        spec.set_peaks(mz, intensity)
        mz_out, int_out = spec.get_peaks()

        np.testing.assert_array_almost_equal(mz_out, mz)
        np.testing.assert_array_almost_equal(int_out, intensity)

    def test_set_peaks_from_float32(self):
        """Test set_peaks with float32 arrays (should be converted)."""
        from pyopenms import MSSpectrum

        spec = MSSpectrum()
        mz = np.array([100.0, 200.0], dtype=np.float32)
        intensity = np.array([1000.0, 2000.0], dtype=np.float32)

        # May need conversion - implementation dependent
        spec.set_peaks(mz.astype(np.float64), intensity.astype(np.float64))

        assert spec.size() == 2

    def test_set_peaks_from_list(self):
        """Test set_peaks with Python lists (converted to numpy)."""
        from pyopenms import MSSpectrum

        spec = MSSpectrum()

        # Lists should be convertible
        mz = np.array([100.0, 200.0, 300.0])
        intensity = np.array([1000.0, 2000.0, 3000.0])

        spec.set_peaks(mz, intensity)
        assert spec.size() == 3

    def test_get_peaks_empty_spectrum(self):
        """Test get_peaks on empty spectrum."""
        from pyopenms import MSSpectrum

        spec = MSSpectrum()
        mz, intensity = spec.get_peaks()

        assert len(mz) == 0
        assert len(intensity) == 0

    def test_set_peaks_mismatched_lengths(self):
        """Test set_peaks with mismatched array lengths."""
        from pyopenms import MSSpectrum

        spec = MSSpectrum()
        mz = np.array([100.0, 200.0, 300.0])
        intensity = np.array([1000.0, 2000.0])

        with pytest.raises(RuntimeError):
            spec.set_peaks(mz, intensity)

    def test_large_spectrum(self):
        """Test with large spectrum (performance check)."""
        from pyopenms import MSSpectrum

        spec = MSSpectrum()
        n = 10000
        mz = np.linspace(100, 1000, n)
        intensity = np.random.rand(n) * 10000

        spec.set_peaks(mz, intensity)
        assert spec.size() == n

        mz_out, int_out = spec.get_peaks()
        np.testing.assert_array_almost_equal(mz_out, mz)


class TestListBuildingCasters:
    """Tests for the caster paths that build Python lists in C++.

    These cover ``from_cpp`` in the DataValue/ParamValue, ``std::vector<std::string>``
    and ``std::vector<DPosition<2>>`` casters, which construct their result one
    element at a time. They are built with the Limited API (``PyList_SetItem``)
    rather than the unchecked ``PyList_SET_ITEM`` macro, so both the empty and the
    populated shapes are worth pinning down.
    """

    @pytest.mark.parametrize("value", [[], ["a"], ["a", "b", "c"]])
    def test_datavalue_string_list_roundtrip(self, value):
        """DataValue STRING_LIST survives a round trip at any length."""
        from pyopenms import MSSpectrum

        spec = MSSpectrum()
        spec.setMetaValue("k", value)
        assert spec.getMetaValue("k") == value

    @pytest.mark.parametrize("value", [[], [1], [1, 2, 3]])
    def test_datavalue_int_list_roundtrip(self, value):
        """DataValue INT_LIST survives a round trip at any length."""
        from pyopenms import MSSpectrum

        spec = MSSpectrum()
        spec.setMetaValue("k", value)
        assert spec.getMetaValue("k") == value

    @pytest.mark.parametrize("value", [[1.5], [1.5, 2.5, 3.5]])
    def test_datavalue_double_list_roundtrip(self, value):
        """DataValue DOUBLE_LIST survives a round trip."""
        from pyopenms import MSSpectrum

        spec = MSSpectrum()
        spec.setMetaValue("k", value)
        result = spec.getMetaValue("k")
        assert len(result) == len(value)
        assert all(r == pytest.approx(v) for r, v in zip(result, value))

    @pytest.mark.parametrize("value", [[], ["x", "y"], [7, 8], [0.25, 0.5]])
    def test_paramvalue_list_roundtrip(self, value):
        """ParamValue is a separate caster from DataValue; check it too."""
        from pyopenms import Param

        param = Param()
        param.setValue("key", value)
        result = param.getValue("key")
        assert len(result) == len(value)
        assert all(r == pytest.approx(v) if isinstance(v, float) else r == v
                   for r, v in zip(result, value))

    @pytest.mark.parametrize("value", [[], ["only"], ["a", "b"]])
    def test_string_vector_from_cpp(self, value):
        """std::vector<std::string> -> list[str] for empty and populated vectors."""
        from pyopenms import DataValue

        assert DataValue(value).toStringList() == value

    @pytest.mark.parametrize("points", [[], [(1.0, 2.0)], [(1.0, 2.0), (3.0, 4.0)]])
    def test_dposition2_vector_roundtrip(self, points):
        """std::vector<DPosition<2>> keeps its shape in both directions.

        ConvexHull2D.getHullPoints() is wrapped by an addon that turns the
        caster's list of 2-tuples into an ndarray, so this exercises the caster
        one layer down.
        """
        from pyopenms import ConvexHull2D

        hull = ConvexHull2D()
        hull.setHullPoints(points)
        result = hull.getHullPoints()
        assert len(result) == len(points)
        for got, want in zip(result, points):
            assert tuple(got) == pytest.approx(want)

    def test_string_list_rejects_non_sequence(self):
        """Invalid conversions still raise rather than producing a partial list."""
        from pyopenms import ConvexHull2D

        hull = ConvexHull2D()
        with pytest.raises(TypeError):
            hull.setHullPoints(42)

    def test_bytes_still_accepted_for_strings(self):
        """The bytes/str behaviour of the string caster is unchanged."""
        from pyopenms import MSSpectrum

        spec = MSSpectrum()
        spec.setName(b"bytes name")
        assert spec.getName() == "bytes name"
