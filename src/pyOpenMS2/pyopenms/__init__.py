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

    # Import the main module which brings in all bindings
    try:
        from . import _pyopenms2

        # Re-export all public symbols from the main module
        for name in dir(_pyopenms2):
            if not name.startswith("_"):
                globals()[name] = getattr(_pyopenms2, name)
    except ImportError as e:
        raise ImportError(
            f"Failed to import pyOpenMS2 bindings: {e}\n"
            "Make sure the package was built correctly with nanobind."
        ) from e


# Import bindings
_import_submodules()

# Import pure Python modules
from ._dataframes import DataFrameMixin

# Apply Python addons to classes
from .addons import apply_addons

apply_addons(globals())

# Clean up namespace
del _import_submodules, apply_addons, DataFrameMixin
