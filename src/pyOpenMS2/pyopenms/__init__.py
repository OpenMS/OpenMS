"""
pyOpenMS2: Python bindings for OpenMS using nanobind

This package provides Python bindings for the OpenMS C++ library,
enabling access to mass spectrometry data structures and algorithms.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

__version__ = "3.4.0"

# Set up OpenMS data path if not already set
if "OPENMS_DATA_PATH" not in os.environ:
    # Try to find the data path relative to this package
    _possible_paths = [
        Path(__file__).parent / "share" / "OpenMS",
        Path(__file__).parent.parent / "share" / "OpenMS",
        Path(sys.prefix) / "share" / "OpenMS",
    ]
    for _path in _possible_paths:
        if _path.exists():
            os.environ["OPENMS_DATA_PATH"] = str(_path)
            break


def _import_submodules():
    """Import all nanobind submodules and merge into this namespace."""
    import importlib
    import importlib.util

    # Dynamically discover all _pyopenms2_N modules
    # This handles any value of PY_NUM_MODULES without hardcoding
    _imported_modules = []
    i = 1
    while True:
        module_name = f"_pyopenms2_{i}"
        # Check if module exists before trying to import
        spec = importlib.util.find_spec(f".{module_name}", package=__name__)
        if spec is None:
            break  # No more modules
        try:
            mod = importlib.import_module(f".{module_name}", package=__name__)
            _imported_modules.append(mod)
            # Re-export all public symbols from each submodule
            for name in dir(mod):
                if not name.startswith("_"):
                    globals()[name] = getattr(mod, name)
        except ImportError:
            # Module exists but failed to load - continue to next
            pass
        i += 1

    if not _imported_modules:
        raise ImportError(
            "Failed to import any pyOpenMS2 bindings.\n"
            "Make sure the package was built correctly with nanobind."
        )


# Import bindings
_import_submodules()

# Import pure Python modules
from ._dataframes import DataFrameMixin

# Apply Python addons to classes
from .addons import apply_addons

apply_addons(globals())

# Clean up namespace
del _import_submodules, apply_addons, DataFrameMixin
