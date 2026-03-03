"""
Tests for MSSpectrum.rasterizeIMFrame() method in pyOpenMS.

These tests cover:
- Basic rasterization with SUM and MAX aggregation
- Edge cases (spectrum without IM data, out-of-range peaks)
- Invalid input validation (non-contiguous, invalid aggregation)
- In-place modification and zero-initialization
- Sum vs max aggregation behavior
"""

import pytest
import numpy as np
import pyopenms as poms


def create_im_spectrum():
    """Create a test MSSpectrum with ion mobility data."""
    spec = poms.MSSpectrum()
    spec.setRT(100.0)
    spec.setMSLevel(1)

    # 5 peaks at different m/z values
    mz = np.array([400.0, 450.0, 500.0, 550.0, 600.0], dtype=np.float64)
    intensity = np.array([100.0, 200.0, 300.0, 400.0, 500.0], dtype=np.float32)
    spec.set_peaks(mz, intensity)

    # Add ion mobility float data array
    fda = poms.FloatDataArray()
    fda.setName("Ion Mobility")
    for v in [0.1, 0.3, 0.5, 0.7, 0.9]:
        fda.push_back(v)
    spec.setFloatDataArrays([fda])

    assert spec.containsIMData(), "Spectrum should contain IM data"
    return spec


def test_rasterizeIMFrame_basic_sum():
    """Test basic rasterization with SUM aggregation."""
    spec = create_im_spectrum()

    mz_bins, im_bins = 10, 10
    output = np.zeros((mz_bins, im_bins), dtype=np.float32)

    spec.rasterizeIMFrame(output, 0.0, 1.0, 400.0, 600.0, "sum")

    assert np.any(output > 0), "Output should contain non-zero intensity values"
    assert output.dtype == np.float32


def test_rasterizeIMFrame_basic_max():
    """Test basic rasterization with MAX aggregation."""
    spec = create_im_spectrum()

    mz_bins, im_bins = 10, 10
    output = np.zeros((mz_bins, im_bins), dtype=np.float32)

    spec.rasterizeIMFrame(output, 0.0, 1.0, 400.0, 600.0, "max")

    assert np.any(output > 0), "Output should contain non-zero intensity values"


def test_rasterizeIMFrame_default_aggregation():
    """Test that default aggregation is 'sum'."""
    spec = create_im_spectrum()

    mz_bins, im_bins = 10, 10
    output_explicit = np.zeros((mz_bins, im_bins), dtype=np.float32)
    output_default = np.zeros((mz_bins, im_bins), dtype=np.float32)

    spec.rasterizeIMFrame(output_explicit, 0.0, 1.0, 400.0, 600.0, "sum")
    spec.rasterizeIMFrame(output_default, 0.0, 1.0, 400.0, 600.0)

    np.testing.assert_array_equal(output_explicit, output_default)


def test_rasterizeIMFrame_sum_vs_max():
    """Test that SUM and MAX aggregation produce different results when peaks overlap."""
    spec = poms.MSSpectrum()

    # Multiple peaks that will land in the same bin
    mz = np.array([500.0, 500.1, 500.2], dtype=np.float64)
    intensity = np.array([1000.0, 1000.0, 1000.0], dtype=np.float32)
    spec.set_peaks(mz, intensity)

    fda = poms.FloatDataArray()
    fda.setName("Ion Mobility")
    for v in [0.5, 0.5, 0.5]:
        fda.push_back(v)
    spec.setFloatDataArrays([fda])

    mz_bins, im_bins = 2, 2
    output_sum = np.zeros((mz_bins, im_bins), dtype=np.float32)
    output_max = np.zeros((mz_bins, im_bins), dtype=np.float32)

    spec.rasterizeIMFrame(output_sum, 0.0, 1.0, 490.0, 510.0, "sum")
    spec.rasterizeIMFrame(output_max, 0.0, 1.0, 490.0, 510.0, "max")

    assert np.max(output_sum) > np.max(output_max), \
        "SUM should accumulate higher than MAX for overlapping peaks"
    assert np.isclose(np.max(output_max), 1000.0, rtol=0.01)


def test_rasterizeIMFrame_out_of_range():
    """Test rasterization when data is outside requested range."""
    spec = create_im_spectrum()

    mz_bins, im_bins = 10, 10
    output = np.ones((mz_bins, im_bins), dtype=np.float32) * 999.0

    # Request m/z range outside spectrum data
    spec.rasterizeIMFrame(output, 0.0, 1.0, 1000.0, 2000.0, "sum")

    assert np.all(output == 0.0), "Out-of-range request should produce all-zero output"


def test_rasterizeIMFrame_case_insensitive_aggregation():
    """Test that aggregation mode accepts uppercase."""
    spec = create_im_spectrum()

    mz_bins, im_bins = 5, 5
    output = np.zeros((mz_bins, im_bins), dtype=np.float32)

    spec.rasterizeIMFrame(output, 0.0, 1.0, 400.0, 600.0, "SUM")
    assert np.any(output > 0)

    output[:] = 0
    spec.rasterizeIMFrame(output, 0.0, 1.0, 400.0, 600.0, "MAX")
    assert np.any(output > 0)


def test_rasterizeIMFrame_in_place_modification():
    """Test that the function modifies the array in-place."""
    spec = create_im_spectrum()

    output = np.zeros((10, 10), dtype=np.float32)
    original_id = id(output)

    spec.rasterizeIMFrame(output, 0.0, 1.0, 400.0, 600.0, "sum")

    assert id(output) == original_id, "Array should be modified in-place"
    assert np.any(output > 0), "Array should contain non-zero values"


def test_rasterizeIMFrame_no_im_data_raises():
    """Test that calling on a spectrum without IM data raises RuntimeError."""
    spec = poms.MSSpectrum()
    spec.set_peaks(np.array([500.0]), np.array([100.0], dtype=np.float32))

    output = np.zeros((5, 5), dtype=np.float32)

    with pytest.raises(RuntimeError):
        spec.rasterizeIMFrame(output, 0.0, 1.0, 400.0, 600.0)


def test_rasterizeIMFrame_exception_non_contiguous():
    """Test that non-contiguous arrays raise ValueError."""
    spec = create_im_spectrum()

    output = np.zeros((10, 20), dtype=np.float32)
    non_contiguous = output.T
    assert not non_contiguous.flags['C_CONTIGUOUS']

    with pytest.raises(ValueError, match="C-contiguous"):
        spec.rasterizeIMFrame(non_contiguous, 0.0, 1.0, 400.0, 600.0)


def test_rasterizeIMFrame_exception_invalid_aggregation():
    """Test that invalid aggregation mode raises ValueError."""
    spec = create_im_spectrum()
    output = np.zeros((10, 10), dtype=np.float32)

    with pytest.raises(ValueError, match="Invalid aggregation mode"):
        spec.rasterizeIMFrame(output, 0.0, 1.0, 400.0, 600.0, "invalid")


def test_rasterizeIMFrame_shape_validation():
    """Test various output shapes are correctly interpreted."""
    spec = create_im_spectrum()

    for mz_bins, im_bins in [(5, 10), (20, 30), (100, 200)]:
        output = np.zeros((mz_bins, im_bins), dtype=np.float32)
        spec.rasterizeIMFrame(output, 0.0, 1.0, 400.0, 600.0)
        assert output.shape == (mz_bins, im_bins)


def test_rasterizeIMFrame_total_intensity_conserved():
    """Test that sum aggregation conserves total intensity for non-overlapping bins."""
    spec = poms.MSSpectrum()

    mz = np.array([500.0], dtype=np.float64)
    intensity = np.array([1234.0], dtype=np.float32)
    spec.set_peaks(mz, intensity)

    fda = poms.FloatDataArray()
    fda.setName("Ion Mobility")
    fda.push_back(0.5)
    spec.setFloatDataArrays([fda])

    output = np.zeros((100, 100), dtype=np.float32)
    spec.rasterizeIMFrame(output, 0.0, 1.0, 400.0, 600.0, "sum")

    assert np.isclose(output.sum(), 1234.0, rtol=0.01), \
        f"Total intensity should be conserved, got {output.sum()}"
