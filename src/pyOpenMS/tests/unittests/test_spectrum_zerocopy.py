import gc
import numpy as np
import pytest
import pyopenms as poms

def test_spectrum_zero_copy_modification():
    spec = poms.MSSpectrum()
    
    # Populate with some dummy data
    spec.push_back(poms.Peak1D(100.0, 50.0))
    spec.push_back(poms.Peak1D(200.0, 75.0))
    
    # Get the zero-copy structured array view
    peaks_array = spec.get_peaks_struct()
    
    # 1. Check if values match
    assert np.isclose(peaks_array['mz'][0], 100.0)
    assert np.isclose(peaks_array['intensity'][1], 75.0)
    
    # 2. PROVE ZERO-COPY: Modify the numpy array and check the C++ object
    peaks_array['intensity'][0] = 999.0
    assert np.isclose(spec[0].getIntensity(), 999.0), "Zero-copy modification failed!"

def test_spectrum_memory_safety():
    # Create the spectrum and get the view inside a local scope
    def get_view():
        spec = poms.MSSpectrum()
        spec.push_back(poms.Peak1D(300.0, 10.0))
        return spec.get_peaks_struct()
    
    # 'spec' goes out of scope here and would normally be destroyed, 
    # but the numpy array view should keep the C++ memory alive.
    peaks_array = get_view()
    
    # Force Python's garbage collection
    gc.collect()
    
    # If nanobind's reference_internal policy failed, this next line will segfault (crash).
    # If it passes, the memory was safely kept alive by the view.
    assert np.isclose(peaks_array['mz'][0], 300.0), "Memory safety/lifetime sharing failed!"

def test_spectrum_empty_struct():
    """get_peaks_struct() on empty MSSpectrum should return empty structured array."""
    spec = poms.MSSpectrum()
    arr = spec.get_peaks_struct()
    assert isinstance(arr, np.ndarray)
    assert arr.shape == (0,)
    assert arr.dtype['mz'] == np.float64
    assert arr.dtype['intensity'] == np.float32
