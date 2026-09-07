"""Shared helpers for locating OpenMS C++ class test fixtures."""

from __future__ import annotations

import os
from pathlib import Path


def get_class_test_data_dir() -> str:
    """Return path to ``src/tests/class_tests/openms/data``.

    Resolution order:
    1. ``OPENMS_CLASS_TEST_DATA_PATH`` environment variable
    2. Walk upward from this file looking for ``src/tests/class_tests/openms/data``.
       Depth-independent, so it works both in the source tree
       (``src/pyOpenMS/tests/unittests/``) and in the build tree
       (``OpenMS-build/pyOpenMS/tests/unittests/``, where ``src`` sits one level
       above the build directory).
    """
    env_path = os.environ.get("OPENMS_CLASS_TEST_DATA_PATH")
    if env_path and os.path.isdir(env_path):
        return env_path

    rel = Path("src") / "tests" / "class_tests" / "openms" / "data"
    for base in Path(__file__).resolve().parents:
        candidate = base / rel
        if candidate.is_dir():
            return str(candidate)

    raise RuntimeError(
        "Could not locate class test data directory "
        "(src/tests/class_tests/openms/data). "
        "Set OPENMS_CLASS_TEST_DATA_PATH to override."
    )
