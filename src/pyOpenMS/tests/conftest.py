"""
Pytest configuration file for pyOpenMS tests.

This file provides fixtures and configuration for test data paths,
replacing the old env.py import pattern.
"""

import os

# Import pytest only if available (needed for fixtures)
try:
    import pytest
    _HAS_PYTEST = True
except ImportError:
    _HAS_PYTEST = False


def _get_test_data_dir():
    """
    Get the path to the OpenMS test data directory.
    
    This function tries multiple strategies to find the test data:
    1. Check the OPENMS_TEST_DATA_PATH environment variable (for custom locations)
    2. Use relative path from this file (works for regular source clones)
    3. Try to import env.py if available (for CMake builds, backwards compatibility)
    
    Returns:
        str: Path to the test data directory (src/tests/topp)
    """
    # Strategy 1: Check environment variable for override
    env_path = os.environ.get('OPENMS_TEST_DATA_PATH')
    if env_path and os.path.isdir(env_path):
        return env_path
    
    # Strategy 2: Use relative path from this conftest.py file
    # This works for regular clones where the structure is:
    # src/pyOpenMS/tests/conftest.py  (this file)
    # src/tests/topp/                 (test data)
    tests_dir = os.path.dirname(os.path.abspath(__file__))
    relative_path = os.path.join(tests_dir, '..', '..', '..', 'src', 'tests', 'topp')
    relative_path = os.path.normpath(relative_path)
    if os.path.isdir(relative_path):
        return relative_path
    
    # Strategy 3: Try to import env.py (backwards compatibility for CMake builds)
    try:
        import env
        if hasattr(env, 'PYOPENMS_SRC_DIR'):
            env_based_path = os.path.join(env.PYOPENMS_SRC_DIR, "..", "..", "src", "tests", "topp")
            env_based_path = os.path.normpath(env_based_path)
            if os.path.isdir(env_based_path):
                return env_based_path
    except ImportError:
        pass
    
    # If all strategies fail, raise an error
    raise RuntimeError(
        "Could not locate OpenMS test data directory. "
        "Please set the OPENMS_TEST_DATA_PATH environment variable to point to "
        "the 'src/tests/topp' directory, or ensure you're running tests from a "
        "complete OpenMS source tree."
    )


# Only define pytest fixtures if pytest is available
if _HAS_PYTEST:
    @pytest.fixture(scope='session')
    def openms_test_data_dir():
        """
        Pytest fixture providing the path to the OpenMS test data directory.
        
        Usage in tests:
            def test_something(openms_test_data_dir):
                data_file = os.path.join(openms_test_data_dir, "test_file.mzML")
        
        Returns:
            str: Path to the test data directory
        """
        return _get_test_data_dir()

