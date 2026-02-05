"""
pytest configuration for legacy pyOpenMS tests running against nanobind build.

Ensures the nanobind build is imported when running legacy tests.
"""

import sys
from pathlib import Path
import importlib
import importlib.util
import os


def _has_compiled_extensions(path: Path) -> bool:
    """Check if pyopenms directory has compiled extensions."""
    pyopenms_dir = path / "pyopenms"
    if not pyopenms_dir.exists():
        return False
    # Look for any _pyopenms_*.so or _pyopenms2_*.so file
    return bool(list(pyopenms_dir.glob("_pyopenms_*.so")) or list(pyopenms_dir.glob("_pyopenms2_*.so")))


def _setup_pyopenms():
    """Ensure the nanobind build is used by pre-loading it into sys.modules."""
    build_path = None

    # Check PYTHONPATH for build directory
    for p in sys.path:
        if ("pyOpenMS2-build" in p or "pyOpenMS" in p) and _has_compiled_extensions(Path(p)):
            build_path = Path(p)
            break

    # Fallback: relative to source tree
    if not build_path:
        # src/pyOpenMS/tests/unittests/conftest.py -> ../../../.. = repo root
        repo_root = Path(__file__).parent.parent.parent.parent.parent
        for candidate in [
            repo_root.parent / "OpenMS-build" / "pyOpenMS",
            repo_root.parent / "OpenMS-build" / "pyOpenMS2-build",
            repo_root / "build" / "pyOpenMS",
            repo_root / "build" / "pyOpenMS2-build",
        ]:
            if candidate.exists() and _has_compiled_extensions(candidate):
                build_path = candidate
                break

    if not build_path:
        # Skip setup if no nanobind build found (will use system/source pyopenms)
        return

    # Clear any existing pyopenms from sys.modules
    to_remove = [k for k in list(sys.modules.keys()) if k.startswith("pyopenms")]
    for k in to_remove:
        del sys.modules[k]

    # Remove source pyopenms directories from sys.path to prevent accidental import
    source_markers = ["src/pyOpenMS/pyopenms", "src/pyOpenMS2/pyopenms"]
    for p in list(sys.path):
        if any(marker in str(Path(p).resolve()) for marker in source_markers):
            if "build" not in str(p).lower():
                sys.path.remove(p)

    # Add build path to front of sys.path
    build_path_str = str(build_path)
    if build_path_str in sys.path:
        sys.path.remove(build_path_str)
    sys.path.insert(0, build_path_str)

    # Force import from the correct location
    pyopenms_init = build_path / "pyopenms" / "__init__.py"
    spec = importlib.util.spec_from_file_location("pyopenms", pyopenms_init)
    if spec and spec.loader:
        pyopenms = importlib.util.module_from_spec(spec)
        sys.modules["pyopenms"] = pyopenms
        spec.loader.exec_module(pyopenms)


# Set up pyopenms immediately when conftest is loaded
_setup_pyopenms()


# Import fixtures from parent conftest
from pathlib import Path
import pytest


def _get_test_data_dir():
    """Get the path to the OpenMS test data directory."""
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
def openms_test_data_dir():
    """Pytest fixture providing the path to the OpenMS test data directory."""
    return _get_test_data_dir()
