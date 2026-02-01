"""Cython declaration for the Parquet filter builder wrapper.

This pxd exposes the cdef class for the standalone `_parquet_filter` module.
It intentionally lives outside the autowrapped `_pyopenms_*` modules to provide
stable cimports for addons and query helpers.
"""

from libcpp.memory cimport shared_ptr
from ParquetFilter cimport ParquetFilterBuilder as _ParquetFilterBuilder

cdef class ParquetFilterBuilder:
    cdef shared_ptr[_ParquetFilterBuilder] inst
