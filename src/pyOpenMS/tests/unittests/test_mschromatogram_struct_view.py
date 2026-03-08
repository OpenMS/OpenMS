import gc
import numpy as np
import pyopenms as poms


def test_chromatogram_zero_copy_modification():
    chrom = poms.MSChromatogram()

    chrom.set_peaks([[1.0, 2.0], [10.0, 20.0]])

    arr = chrom.get_peaks_struct()

    # Value correctness
    assert np.isclose(arr["rt"][0], 1.0)
    assert np.isclose(arr["intensity"][1], 20.0)
    assert arr["intensity"].dtype == np.float32

    # PROVE ZERO-COPY
    arr["intensity"][0] = 999.0

    _, int_copy = chrom.get_peaks()
    assert np.isclose(int_copy[0], 999.0), "Zero-copy modification failed!"


def test_chromatogram_memory_safety():
    def get_view():
        chrom = poms.MSChromatogram()
        chrom.set_peaks([[5.0], [42.0]])
        return chrom.get_peaks_struct()

    arr = get_view()

    # Force garbage collection
    gc.collect()

    # If lifetime policy is broken, this will segfault
    assert np.isclose(arr["intensity"][0], 42.0), "Memory safety/lifetime sharing failed!"


def test_empty_chromatogram():
    chrom = poms.MSChromatogram()
    arr = chrom.get_peaks_struct()

    assert isinstance(arr, np.ndarray)
    assert arr.size == 0
    assert arr.dtype["rt"] == np.dtype(np.float64)
    assert arr.dtype["intensity"] == np.dtype(np.float32)


def test_consistency_with_get_peaks():
    rts = [0.5, 1.5, 2.5]
    ints = [100.0, 200.0, 300.0]

    chrom = poms.MSChromatogram()
    chrom.set_peaks([rts, ints])

    arr = chrom.get_peaks_struct()
    rt_copy, int_copy = chrom.get_peaks()

    assert int_copy.dtype == np.float32
    np.testing.assert_array_equal(arr["rt"], rt_copy)
    np.testing.assert_array_equal(arr["intensity"], int_copy)
