"""
pytest configuration for pyOpenMS unit tests.

Path setup is handled by the parent tests/conftest.py (runs first).
This file only provides fixtures specific to unittests.
"""

import os
import sys
import pytest

# test_data_paths sits next to this file; make it importable once for the whole session
_UNITTESTS_DIR = os.path.dirname(os.path.abspath(__file__))
if _UNITTESTS_DIR not in sys.path:
    sys.path.insert(0, _UNITTESTS_DIR)


@pytest.fixture(scope='session')
def openms_test_data_dir():
    """Pytest fixture providing the path to the OpenMS test data directory."""
    # Check environment variable for override
    env_path = os.environ.get('OPENMS_TEST_DATA_PATH')
    if env_path and os.path.isdir(env_path):
        return env_path

    # Use relative path from this conftest.py file
    tests_dir = os.path.dirname(os.path.abspath(__file__))
    relative_path = os.path.join(tests_dir, '..', '..', '..', '..', 'src', 'tests', 'topp')
    relative_path = os.path.normpath(relative_path)
    if os.path.isdir(relative_path):
        return relative_path

    raise RuntimeError(
        "Could not locate OpenMS test data directory. "
        "Please set the OPENMS_TEST_DATA_PATH environment variable."
    )


@pytest.fixture(scope='session')
def imzml_continuous_path() -> str:
    """Path to the continuous imzML example of the C++ class tests; skips when it is unavailable."""
    from test_data_paths import get_class_test_data_dir

    try:
        path = os.path.join(get_class_test_data_dir(), "ImzMLFile_1_Example_Continuous.imzML")
    except RuntimeError as e:
        pytest.skip(str(e))
    if not os.path.isfile(path):
        pytest.skip(f"imzML test data not found: {path}")
    return path
