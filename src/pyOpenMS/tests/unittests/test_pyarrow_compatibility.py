"""
Test pyarrow compatibility with pyopenms.

This module tests that pyopenms and pyarrow can be imported in any order
without causing import conflicts or symbol collisions. This is important because
both libraries may depend on shared C++ libraries (like abseil-cpp) that could
have conflicting symbols if not properly isolated.

Background:
- pyopenms uses Arrow/Parquet for file I/O
- pyarrow provides Python bindings for Apache Arrow
- Both may link against shared libraries that need proper symbol isolation
"""

import pytest


def test_import_pyarrow_then_pyopenms():
    """Test importing pyarrow before pyopenms works correctly."""
    # Skip if pyarrow is not installed
    pytest.importorskip('pyarrow')
    
    # Import in order: pyarrow first, then pyopenms
    import pyarrow
    import pyopenms
    
    # Basic sanity checks
    assert hasattr(pyarrow, '__version__')
    assert hasattr(pyopenms, 'MSSpectrum')


def test_import_pyopenms_then_pyarrow():
    """Test importing pyopenms before pyarrow works correctly."""
    # Skip if pyarrow is not installed
    pytest.importorskip('pyarrow')
    
    # Import in order: pyopenms first, then pyarrow
    import pyopenms
    import pyarrow
    
    # Basic sanity checks
    assert hasattr(pyopenms, 'MSSpectrum')
    assert hasattr(pyarrow, '__version__')


def test_both_libraries_functional():
    """Test that both libraries remain functional after being imported together."""
    pytest.importorskip('pyarrow')
    
    import pyopenms
    import pyarrow as pa
    
    # Test pyopenms functionality
    spec = pyopenms.MSSpectrum()
    spec.setMSLevel(1)
    assert spec.getMSLevel() == 1
    
    # Test pyarrow functionality
    arr = pa.array([1, 2, 3, 4])
    assert len(arr) == 4
    assert arr[0].as_py() == 1
