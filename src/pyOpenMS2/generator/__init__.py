"""
pyOpenMS2 Generator

This package contains the code generator that:
1. Parses existing .pxd files from pyOpenMS
2. Uses libclang for accurate C++ type information
3. Emits nanobind C++ binding code
4. Processes addon files for Python injection

Usage:
    python -m generator --pxd-dir /path/to/pxds --openms-include-dir /path/to/include --output-dir /path/to/output
"""

from .pxd_parser import PxdParser
from .nanobind_emitter_v2 import NanobindEmitterV2
from .type_registry import TypeRegistry
from .addon_processor import AddonProcessor
from .cpp_parser import CppHeaderParser, MergedClass, MergedMethod

__all__ = [
    "PxdParser",
    "NanobindEmitterV2",
    "TypeRegistry",
    "AddonProcessor",
    "CppHeaderParser",
    "MergedClass",
    "MergedMethod",
]
