"""
pyOpenMS2 Generator

This package contains the code generator that:
1. Parses existing .pxd files from pyOpenMS
2. Emits nanobind C++ binding code
3. Processes addon files for Python injection

Usage:
    python -m generator --pxd-dir /path/to/pxds --output-dir /path/to/output
"""

from .pxd_parser import PxdParser
from .nanobind_emitter import NanobindEmitter
from .type_registry import TypeRegistry
from .addon_processor import AddonProcessor

__all__ = ["PxdParser", "NanobindEmitter", "TypeRegistry", "AddonProcessor"]
