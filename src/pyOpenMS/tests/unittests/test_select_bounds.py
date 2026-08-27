## ----------------------------------------------------------------------------
## $Maintainer: $
## $Authors: Alivia Hossain $
## ----------------------------------------------------------------------------
"""Regression tests for MSSpectrum.select() bounds/size validation (issue #9807).

``MSSpectrum::select()`` is exposed to Python.  Before #9807 it performed no index-range
check and validated data-array sizes only *after* permuting the peaks, so from Python::

    spec.select([999999])          # out-of-bounds read -> segfault

crashed the interpreter, and a mismatched data array left the spectrum half-modified.
The checked ``select()`` now validates every index up front and leaves the spectrum
unchanged on rejection.  Duplicate indices are undefined behaviour in C++ (moved-from
entries), so the nanobind ``select()`` binding rejects them with a ``ValueError`` before
the call ever reaches C++.
"""

import numpy as np
import pytest

import pyopenms


def _spectrum(n):
    """A spectrum with `n` peaks and one aligned string data array ("0".."n-1")."""
    s = pyopenms.MSSpectrum()
    mz = np.arange(1.0, n + 1.0, dtype=np.float64)
    inten = (np.arange(1.0, n + 1.0) * 10.0).astype(np.float32)
    s.set_peaks((mz, inten))
    sda = pyopenms.StringDataArray()
    sda.set_data([str(i) for i in range(n)])
    s.setStringDataArrays([sda])
    return s


def _snapshot(s):
    """Full observable state: peak arrays plus the string data-array values."""
    mz, inten = s.get_peaks()
    return mz.tolist(), inten.tolist(), [list(a) for a in s.getStringDataArrays()]


def test_out_of_range_index_raises_and_leaves_spectrum_untouched():
    s = _spectrum(3)
    before = _snapshot(s)
    with pytest.raises(Exception):
        s.select([0, 999999])  # used to segfault
    assert _snapshot(s) == before  # not just sizes: no value was altered before the raise


def test_valid_subset_reorders_and_keeps_arrays_aligned():
    s = _spectrum(3)
    s.select([2, 0])
    assert s.size() == 2
    assert list(s.getStringDataArrays()[0]) == ["2", "0"]


def test_duplicate_indices_rejected():
    s = _spectrum(3)
    before = _snapshot(s)
    with pytest.raises(ValueError):
        s.select([0, 0])
    assert _snapshot(s) == before


def test_mis_sized_data_array_rejected_before_mutation():
    s = _spectrum(3)
    sda = pyopenms.StringDataArray()
    sda.set_data(["a", "b"])  # 2 entries for 3 peaks
    s.setStringDataArrays([sda])
    before = _snapshot(s)
    with pytest.raises(Exception):
        s.select([0, 1, 2])
    assert _snapshot(s) == before  # unchanged, values included
