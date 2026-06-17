import gc
import numpy as np
import pytest
import pyopenms as poms


def test_chromatogram_zero_copy_modification():
    chrom = poms.MSChromatogram()

    chrom.set_peaks([[1.0, 2.0], [10.0, 20.0]])

    arr = chrom.get_peaks_struct()

    # Value correctness
    assert np.isclose(arr['rt'][0], 1.0)   # RT
    assert np.isclose(arr['intensity'][1], np.float32(20.0))  # intensity

    # PROVE ZERO-COPY
    arr['intensity'][0] = np.float32(999.0)

    _, int_copy = chrom.get_peaks()
    assert np.isclose(int_copy[0], np.float32(999.0)), "Zero-copy modification failed!"


def test_chromatogram_memory_safety():
    def get_view():
        chrom = poms.MSChromatogram()
        chrom.set_peaks([[5.0], [42.0]])
        return chrom.get_peaks_struct()

    arr = get_view()

    # Force garbage collection
    gc.collect()

    # If lifetime policy is broken, this will segfault
    assert np.isclose(arr['intensity'][0], np.float32(42.0)), "Memory safety/lifetime sharing failed!"


def test_empty_chromatogram():
    chrom = poms.MSChromatogram()
    arr = chrom.get_peaks_struct()

    assert isinstance(arr, np.ndarray)
    assert arr.shape == (0,)
    assert 'rt' in arr.dtype.names
    assert 'intensity' in arr.dtype.names
    assert arr.dtype['rt'] == np.float64
    assert arr.dtype['intensity'] == np.float32


def test_consistency_with_get_peaks():
    rts = [0.5, 1.5, 2.5]
    ints = [100.0, 200.0, 300.0]

    chrom = poms.MSChromatogram()
    chrom.set_peaks([rts, ints])

    arr = chrom.get_peaks_struct()
    rt_copy, int_copy = chrom.get_peaks()

    np.testing.assert_array_equal(arr['rt'], rt_copy)
    np.testing.assert_array_equal(arr['intensity'], int_copy)


def test_chromatogram_intensity_type_is_float32():
    """ChromatogramPeak intensity should be float32, consistent with Peak1D (issue #3379)."""
    chrom = poms.MSChromatogram()
    p = poms.ChromatogramPeak()
    p.setRT(1.0)
    p.setIntensity(100.0)
    chrom.push_back(p)

    rt, intensity = chrom.get_peaks()
    assert rt.dtype == np.float64, f"RT should be float64, got {rt.dtype}"
    assert intensity.dtype == np.float32, f"Intensity should be float32 (like Peak1D), got {intensity.dtype}"


def test_chromatogram_struct_view_dtype():
    """get_peaks_struct() should return {rt: float64, intensity: float32} matching Peak1D layout."""
    chrom = poms.MSChromatogram()
    for i in range(3):
        p = poms.ChromatogramPeak()
        p.setRT(float(i))
        p.setIntensity(float(i * 100))
        chrom.push_back(p)

    arr = chrom.get_peaks_struct()
    assert arr.dtype['rt'] == np.float64, f"RT dtype should be float64, got {arr.dtype['rt']}"
    assert arr.dtype['intensity'] == np.float32, f"Intensity dtype should be float32, got {arr.dtype['intensity']}"
    assert arr['rt'][0] == 0.0
    assert arr['rt'][2] == 2.0
    assert arr['intensity'][1] == np.float32(100.0)


def test_peak_intensity_type_consistency():
    """All peak types should produce float32 intensity in get_peaks()."""
    # MSSpectrum (Peak1D)
    spec = poms.MSSpectrum()
    p1 = poms.Peak1D()
    p1.setMZ(100.0)
    p1.setIntensity(50.0)
    spec.push_back(p1)
    _, spec_int = spec.get_peaks()

    # MSChromatogram (ChromatogramPeak)
    chrom = poms.MSChromatogram()
    cp = poms.ChromatogramPeak()
    cp.setRT(1.0)
    cp.setIntensity(50.0)
    chrom.push_back(cp)
    _, chrom_int = chrom.get_peaks()

    # Mobilogram (MobilityPeak1D)
    mob = poms.Mobilogram()
    mp = poms.MobilityPeak1D()
    mp.setMobility(0.5)
    mp.setIntensity(50.0)
    mob.push_back(mp)
    _, mob_int = mob.get_peaks()

    assert spec_int.dtype == chrom_int.dtype == mob_int.dtype == np.float32, \
        f"All intensity dtypes should be float32: spec={spec_int.dtype}, chrom={chrom_int.dtype}, mob={mob_int.dtype}"
