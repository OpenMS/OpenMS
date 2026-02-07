"""
Tests for MSExperiment.rasterizeRTMZ() method in pyOpenMS.

These tests cover:
- Basic rasterization with SUM and MAX aggregation
- Edge cases (empty experiment, single spectrum)
- MS level filtering
- Invalid input validation (exceptions)
- Array contiguity requirements
"""

import pytest
import numpy as np
import pyopenms as poms


def create_test_experiment():
    """Create a test MSExperiment with known data."""
    exp = poms.MSExperiment()
    
    # Create 3 spectra with known RT and m/z values
    for i in range(3):
        spec = poms.MSSpectrum()
        spec.setRT(100.0 + i * 50.0)  # RT: 100, 150, 200
        spec.setMSLevel(1)
        
        # Add 3 peaks per spectrum
        spec.set_peaks((
            np.array([500.0, 550.0, 600.0], dtype=np.float64),
            np.array([1000.0, 2000.0, 3000.0], dtype=np.float32)
        ))
        
        exp.addSpectrum(spec)
    
    exp.updateRanges()
    return exp


def test_rasterizeRTMZ_basic_sum():
    """Test basic rasterization with SUM aggregation."""
    exp = create_test_experiment()
    
    rt_bins, mz_bins = 10, 10
    output = np.empty((mz_bins, rt_bins), dtype=np.float32)
    
    exp.rasterizeRTMZ(output, 100.0, 200.0, 500.0, 600.0, 1, "sum")
    
    # Output should contain non-zero values
    assert np.any(output > 0), "Output should contain non-zero intensity values"
    assert output.dtype == np.float32, "Output dtype should be float32"


def test_rasterizeRTMZ_basic_max():
    """Test basic rasterization with MAX aggregation."""
    exp = create_test_experiment()
    
    rt_bins, mz_bins = 10, 10
    output = np.empty((mz_bins, rt_bins), dtype=np.float32)
    
    exp.rasterizeRTMZ(output, 100.0, 200.0, 500.0, 600.0, 1, "max")
    
    # Output should contain non-zero values
    assert np.any(output > 0), "Output should contain non-zero intensity values"


def test_rasterizeRTMZ_empty_experiment():
    """Test rasterization with empty experiment."""
    exp = poms.MSExperiment()
    
    rt_bins, mz_bins = 5, 5
    output = np.ones((mz_bins, rt_bins), dtype=np.float32)  # Initialize with non-zero
    
    exp.rasterizeRTMZ(output, 0.0, 100.0, 0.0, 1000.0, 1, "sum")
    
    # Output should be all zeros
    assert np.all(output == 0.0), "Empty experiment should produce all-zero output"


def test_rasterizeRTMZ_single_spectrum():
    """Test rasterization with single spectrum."""
    exp = poms.MSExperiment()
    spec = poms.MSSpectrum()
    spec.setRT(50.0)
    spec.setMSLevel(1)
    spec.set_peaks((
        np.array([500.0], dtype=np.float64),
        np.array([100.0], dtype=np.float32)
    ))
    exp.addSpectrum(spec)
    
    rt_bins, mz_bins = 10, 10
    output = np.empty((mz_bins, rt_bins), dtype=np.float32)
    
    exp.rasterizeRTMZ(output, 0.0, 100.0, 400.0, 600.0, 1, "sum")
    
    # Should have exactly one bin with intensity
    nonzero_count = np.count_nonzero(output)
    assert nonzero_count == 1, f"Expected 1 non-zero bin, got {nonzero_count}"


def test_rasterizeRTMZ_ms_level_filtering():
    """Test that MS level filtering works correctly."""
    exp = poms.MSExperiment()
    
    # Add alternating MS1 and MS2 spectra
    for i in range(5):
        spec = poms.MSSpectrum()
        spec.setRT(100.0 + i * 10.0)
        spec.setMSLevel(1 if i % 2 == 0 else 2)
        spec.set_peaks((
            np.array([500.0], dtype=np.float64),
            np.array([1000.0], dtype=np.float32)
        ))
        exp.addSpectrum(spec)
    
    rt_bins, mz_bins = 10, 10
    output_ms1 = np.empty((mz_bins, rt_bins), dtype=np.float32)
    output_ms2 = np.empty((mz_bins, rt_bins), dtype=np.float32)
    
    # Rasterize MS1 only
    exp.rasterizeRTMZ(output_ms1, 100.0, 150.0, 400.0, 600.0, 1, "sum")
    
    # Rasterize MS2 only
    exp.rasterizeRTMZ(output_ms2, 100.0, 150.0, 400.0, 600.0, 2, "sum")
    
    # Both should have non-zero values
    assert np.sum(output_ms1) > 0, "MS1 rasterization should produce non-zero output"
    assert np.sum(output_ms2) > 0, "MS2 rasterization should produce non-zero output"


def test_rasterizeRTMZ_out_of_range():
    """Test rasterization when data is outside requested range."""
    exp = create_test_experiment()
    
    rt_bins, mz_bins = 10, 10
    output = np.empty((mz_bins, rt_bins), dtype=np.float32)
    
    # Request range outside experiment data
    exp.rasterizeRTMZ(output, 1000.0, 2000.0, 500.0, 600.0, 1, "sum")
    
    # Output should be all zeros
    assert np.all(output == 0.0), "Out-of-range request should produce all-zero output"


def test_rasterizeRTMZ_sum_vs_max():
    """Test that SUM and MAX aggregation produce different results."""
    exp = poms.MSExperiment()
    
    # Create multiple spectra with overlapping peaks (same RT bin)
    for i in range(5):
        spec = poms.MSSpectrum()
        spec.setRT(100.0 + i * 2.0)  # Close together
        spec.setMSLevel(1)
        spec.set_peaks((
            np.array([500.0], dtype=np.float64),
            np.array([1000.0], dtype=np.float32)
        ))
        exp.addSpectrum(spec)
    
    rt_bins, mz_bins = 2, 2
    output_sum = np.empty((mz_bins, rt_bins), dtype=np.float32)
    output_max = np.empty((mz_bins, rt_bins), dtype=np.float32)
    
    exp.rasterizeRTMZ(output_sum, 100.0, 110.0, 490.0, 510.0, 1, "sum")
    exp.rasterizeRTMZ(output_max, 100.0, 110.0, 490.0, 510.0, 1, "max")
    
    sum_value = np.max(output_sum)
    max_value = np.max(output_max)
    
    # SUM should accumulate (5 * 1000 = 5000), MAX should be 1000
    assert sum_value > max_value, "SUM aggregation should produce higher values than MAX"
    assert np.isclose(max_value, 1000.0, rtol=0.01), f"MAX value should be ~1000, got {max_value}"


def test_rasterizeRTMZ_exception_non_contiguous():
    """Test that non-contiguous arrays raise ValueError."""
    exp = create_test_experiment()
    
    # Create non-contiguous array using transpose
    output = np.empty((10, 20), dtype=np.float32)
    non_contiguous = output.T  # Transpose makes it non-contiguous
    
    assert not non_contiguous.flags['C_CONTIGUOUS'], "Array should be non-contiguous"
    
    # Accept any message that mentions 'C-contiguous' (different builds may phrase it)
    with pytest.raises(ValueError, match="C-contiguous"):
        exp.rasterizeRTMZ(non_contiguous, 100.0, 200.0, 500.0, 600.0, 1, "sum")


def test_rasterizeRTMZ_exception_invalid_aggregation():
    """Test that invalid aggregation mode raises ValueError."""
    exp = create_test_experiment()
    output = np.empty((10, 10), dtype=np.float32)
    
    with pytest.raises(ValueError, match="Invalid aggregation mode"):
        exp.rasterizeRTMZ(output, 100.0, 200.0, 500.0, 600.0, 1, "invalid")


def test_rasterizeRTMZ_exception_inverted_rt_range():
    """Test that inverted RT range raises exception."""
    exp = create_test_experiment()
    output = np.empty((10, 10), dtype=np.float32)

    # Equal RT values
    with pytest.raises(RuntimeError):
        exp.rasterizeRTMZ(output, 100.0, 100.0, 500.0, 600.0, 1, "sum")

    # Inverted RT range
    with pytest.raises(RuntimeError):
        exp.rasterizeRTMZ(output, 200.0, 100.0, 500.0, 600.0, 1, "sum")


def test_rasterizeRTMZ_exception_inverted_mz_range():
    """Test that inverted m/z range raises exception."""
    exp = create_test_experiment()
    output = np.zeros((10, 10), dtype=np.float32)

    # Equal m/z values
    with pytest.raises(RuntimeError):
        exp.rasterizeRTMZ(output, 100.0, 200.0, 600.0, 600.0, 1, "sum")

    # Inverted m/z range
    with pytest.raises(RuntimeError):
        exp.rasterizeRTMZ(output, 100.0, 200.0, 700.0, 500.0, 1, "sum")


def test_rasterizeRTMZ_shape_validation():
    """Test that output array shape is correctly interpreted."""
    exp = create_test_experiment()
    
    # Test various shapes
    for mz_bins, rt_bins in [(5, 10), (20, 30), (100, 200)]:
        output = np.empty((mz_bins, rt_bins), dtype=np.float32)
        
        # Should not raise exception
        exp.rasterizeRTMZ(output, 100.0, 200.0, 500.0, 600.0, 1, "sum")
        
        assert output.shape == (mz_bins, rt_bins), "Output shape should be preserved"


def test_rasterizeRTMZ_in_place_modification():
    """Test that the function modifies the array in-place."""
    exp = create_test_experiment()
    
    output = np.empty((10, 10), dtype=np.float32)
    original_id = id(output)
    
    exp.rasterizeRTMZ(output, 100.0, 200.0, 500.0, 600.0, 1, "sum")
    
    # Should modify in-place, not create new array
    assert id(output) == original_id, "Array should be modified in-place"
    assert np.any(output > 0), "Array should be modified with non-zero values"


def test_rasterizeRTMZ_zero_initialization():
    """Test that output is zero-initialized even with pre-existing values."""
    exp = create_test_experiment()
    
    output = np.ones((10, 10), dtype=np.float32) * 999.0  # Pre-fill with large values
    
    # Request range outside data (should zero-initialize)
    exp.rasterizeRTMZ(output, 1000.0, 2000.0, 500.0, 600.0, 1, "sum")
    
    # All values should be zero (not 999)
    assert np.all(output == 0.0), "Output should be zero-initialized"
