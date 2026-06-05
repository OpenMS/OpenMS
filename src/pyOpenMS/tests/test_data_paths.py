"""Shared helpers for locating OpenMS C++ class test fixtures."""

from __future__ import annotations

import os
from pathlib import Path


def get_class_test_data_dir() -> str:
    """Return path to ``src/tests/class_tests/openms/data``.

    Resolution order:
    1. ``OPENMS_CLASS_TEST_DATA_PATH`` environment variable
    2. Source-tree layout (``src/pyOpenMS/tests/`` → ``src/tests/...``)
    3. Build-tree layout (``OpenMS-build/pyOpenMS/tests/`` → ``../src/tests/...``)
    """
    env_path = os.environ.get("OPENMS_CLASS_TEST_DATA_PATH")
    if env_path and os.path.isdir(env_path):
        return env_path

    here = Path(__file__).resolve().parent
    candidates = (
        here.parent.parent / "tests" / "class_tests" / "openms" / "data",
        here.parents[2] / "src" / "tests" / "class_tests" / "openms" / "data",
    )
    for candidate in candidates:
        resolved = candidate.resolve()
        if resolved.is_dir():
            return str(resolved)

    raise RuntimeError(
        "Could not locate class test data directory "
        "(src/tests/class_tests/openms/data). "
        "Set OPENMS_CLASS_TEST_DATA_PATH to override."
    )
