"""
Test pyarrow compatibility with pyopenms.

This module rigorously tests that pyopenms and pyarrow can be imported in any order
without causing import conflicts or symbol collisions. Each test spawns a separate
Python process to avoid module caching effects (sys.modules), ensuring true
import-order independence.

This is important because both libraries may depend on shared C++ libraries
(like abseil-cpp) that could have conflicting symbols if not properly isolated.

Background:
- pyopenms uses Arrow/Parquet for file I/O
- pyarrow provides Python bindings for Apache Arrow
- Both may link against shared libraries that need proper symbol isolation
"""

import sys
import subprocess
import pytest


def test_import_pyarrow_then_pyopenms():
    """Test importing pyarrow before pyopenms works correctly in a fresh process."""
    pytest.importorskip('pyarrow')
    
    code = """
import sys
import pyarrow
import pyopenms
assert hasattr(pyarrow, '__version__'), 'pyarrow missing __version__'
assert hasattr(pyopenms, 'MSSpectrum'), 'pyopenms missing MSSpectrum'
print('SUCCESS: pyarrow then pyopenms')
"""
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        timeout=30
    )
    assert result.returncode == 0, f"Import failed:\n{result.stderr}\n{result.stdout}"


def test_import_pyopenms_then_pyarrow():
    """Test importing pyopenms before pyarrow works correctly in a fresh process."""
    pytest.importorskip('pyarrow')
    
    code = """
import sys
import pyopenms
import pyarrow
assert hasattr(pyopenms, 'MSSpectrum'), 'pyopenms missing MSSpectrum'
assert hasattr(pyarrow, '__version__'), 'pyarrow missing __version__'
print('SUCCESS: pyopenms then pyarrow')
"""
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        timeout=30
    )
    assert result.returncode == 0, f"Import failed:\n{result.stderr}\n{result.stdout}"


def test_both_libraries_functional():
    """Test that both libraries remain functional after being imported together in a fresh process."""
    pytest.importorskip('pyarrow')
    
    code = """
import sys
import pyopenms
import pyarrow as pa

# Test pyopenms functionality
spec = pyopenms.MSSpectrum()
spec.setMSLevel(1)
assert spec.getMSLevel() == 1, 'pyopenms MSSpectrum not functional'

# Test pyarrow functionality
arr = pa.array([1, 2, 3, 4])
assert len(arr) == 4, 'pyarrow array wrong length'
assert arr[0].as_py() == 1, 'pyarrow array wrong value'

print('SUCCESS: both libraries functional')
"""
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        timeout=30
    )
    assert result.returncode == 0, f"Functionality test failed:\n{result.stderr}\n{result.stdout}"
