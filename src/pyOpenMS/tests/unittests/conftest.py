"""
pytest configuration for pyOpenMS unit tests.

Path setup is handled by the parent tests/conftest.py (runs first).
This file only provides fixtures specific to unittests.
"""

import os
import pytest


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
