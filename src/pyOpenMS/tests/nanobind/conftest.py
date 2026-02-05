"""
pytest configuration for pyOpenMS nanobind tests.

Ensures the built pyopenms package is imported rather than the source package,
and that the generator package is importable.
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
    # Look for any _pyopenms_*.so or _pyopenms2_*.so file (handles different Python versions and domain names)
    return bool(list(pyopenms_dir.glob("_pyopenms_*.so")) or list(pyopenms_dir.glob("_pyopenms2_*.so")))


def _setup_pyopenms():
    """Ensure the built pyopenms is used by pre-loading it into sys.modules."""
    # Try to find build path from PYTHONPATH or relative to source
    build_path = None

    # Check PYTHONPATH for build directory
    for p in sys.path:
        if ("pyOpenMS2-build" in p or "pyOpenMS" in p) and _has_compiled_extensions(Path(p)):
            build_path = Path(p)
            break

    # Fallback: relative to source tree
    if not build_path:
        # src/pyOpenMS/tests/nanobind/conftest.py -> ../../../.. = repo root
        repo_root = Path(__file__).parent.parent.parent.parent
        for candidate in [
            repo_root.parent / "OpenMS-build" / "pyOpenMS",
            repo_root.parent / "OpenMS-build" / "pyOpenMS2-build",
            repo_root / "build" / "pyOpenMS",
            repo_root / "build" / "pyOpenMS2-build",
        ]:
            if candidate.exists() and _has_compiled_extensions(candidate):
                build_path = candidate
                break

    if not build_path or not build_path.exists() or not _has_compiled_extensions(build_path):
        raise ImportError(
            "Built pyopenms not found.\n"
            "Please build pyopenms first with: "
            "cmake --build OpenMS-build --target pyopenms -j$(nproc)"
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

    # Remove source paths to prevent accidental import of source pyopenms
    source_markers = ["src/pyOpenMS2", "src/pyOpenMS"]
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


def _setup_generator_path():
    """Add src/pyOpenMS/ to sys.path so 'from generator...' imports work."""
    # src/pyOpenMS/tests/nanobind/conftest.py -> src/pyOpenMS/
    pyopenms_src = Path(__file__).parent.parent.parent
    pyopenms_src_str = str(pyopenms_src)
    if pyopenms_src_str not in sys.path:
        sys.path.insert(0, pyopenms_src_str)


# Set up paths immediately when conftest is loaded
_setup_pyopenms()
_setup_generator_path()
