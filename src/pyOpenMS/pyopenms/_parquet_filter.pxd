"""Cython declaration for the Parquet filter builder wrapper.

This pxd exposes the cdef class for the standalone `_parquet_filter` module.
It intentionally lives outside the autowrapped `_pyopenms_*` modules to provide
stable cimports for addons and query helpers.
"""

from libc.stdint cimport int64_t
from libcpp.vector cimport vector as libcpp_vector
from libcpp.memory cimport shared_ptr

# Inline C++ declarations - must match the pyx file
cdef extern from "<OpenMS/DATASTRUCTURES/String.h>" namespace "OpenMS":
    cdef cppclass _String "OpenMS::String":
        _String() except +
        _String(const char*) except +

cdef extern from "<OpenMS/FORMAT/ParquetFilter.h>" namespace "OpenMS":
    cdef cppclass _ParquetFilterBuilder "OpenMS::ParquetFilterBuilder":
        _ParquetFilterBuilder() except +
        _ParquetFilterBuilder& andNext()
        _ParquetFilterBuilder& orNext()
        _ParquetFilterBuilder& eq(const _String& column, int64_t value)
        _ParquetFilterBuilder& ne(const _String& column, int64_t value)
        _ParquetFilterBuilder& lt(const _String& column, int64_t value)
        _ParquetFilterBuilder& le(const _String& column, int64_t value)
        _ParquetFilterBuilder& gt(const _String& column, int64_t value)
        _ParquetFilterBuilder& ge(const _String& column, int64_t value)
        _ParquetFilterBuilder& in_ "in"(const _String& column, const libcpp_vector[int64_t]& values)
        _ParquetFilterBuilder& eq(const _String& column, const _String& value)
        _ParquetFilterBuilder& ne(const _String& column, const _String& value)
        _ParquetFilterBuilder& lt(const _String& column, const _String& value)
        _ParquetFilterBuilder& le(const _String& column, const _String& value)
        _ParquetFilterBuilder& gt(const _String& column, const _String& value)
        _ParquetFilterBuilder& ge(const _String& column, const _String& value)
        _ParquetFilterBuilder& in_ "in"(const _String& column, const libcpp_vector[_String]& values)
        bint empty() const

cdef class ParquetFilterBuilder:
    cdef shared_ptr[_ParquetFilterBuilder] inst
