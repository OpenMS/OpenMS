"""
pytest configuration for pyOpenMS2 tests.

Ensures the built pyopenms package is imported rather than the source package.
"""

import sys
from pathlib import Path
import importlib
import importlib.util


def _has_compiled_extensions(path: Path) -> bool:
    """Check if pyopenms directory has compiled extensions."""
    pyopenms_dir = path / "pyopenms"
    if not pyopenms_dir.exists():
        return False
    # Look for any _pyopenms2_*.so file (handles different Python versions and domain names)
    return bool(list(pyopenms_dir.glob("_pyopenms2_*.so")))


def _setup_pyopenms():
    """Ensure the built pyopenms is used by pre-loading it into sys.modules."""
    build_path = Path("/home/sachsenb/Development/tmp/OpenMS/OpenMS-build/pyOpenMS2-build")

    if not build_path.exists() or not _has_compiled_extensions(build_path):
        raise ImportError(
            f"Built pyopenms not found at {build_path}.\n"
            "Please build pyopenms2 first with: "
            "cmake --build OpenMS-build --target pyopenms2"
        )

    # Clear any existing pyopenms from sys.modules
    to_remove = [k for k in list(sys.modules.keys()) if k.startswith("pyopenms")]
    for k in to_remove:
        del sys.modules[k]

    # Add build path to front of sys.path
    build_path_str = str(build_path)
    if build_path_str in sys.path:
        sys.path.remove(build_path_str)
    sys.path.insert(0, build_path_str)

    # Remove source paths to prevent accidental import
    source_markers = ["src/pyOpenMS2"]
    for p in list(sys.path):
        if any(marker in str(Path(p).resolve()) for marker in source_markers):
            if "build" not in str(p).lower():
                sys.path.remove(p)

    # Force import from the correct location
    pyopenms_init = build_path / "pyopenms" / "__init__.py"
    spec = importlib.util.spec_from_file_location("pyopenms", pyopenms_init)
    if spec and spec.loader:
        pyopenms = importlib.util.module_from_spec(spec)
        sys.modules["pyopenms"] = pyopenms
        spec.loader.exec_module(pyopenms)


# Set up pyopenms immediately when conftest is loaded
_setup_pyopenms()
