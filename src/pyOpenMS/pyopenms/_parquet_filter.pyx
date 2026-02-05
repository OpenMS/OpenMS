# distutils: language = c++
# cython: language_level=3
"""Typed Parquet filter builder.

This module provides a small, stable wrapper around the C++
OpenMS::ParquetFilterBuilder for use in pyOpenMS query builders.

Why this exists:
  - The autowrapped pyOpenMS modules are split across `_pyopenms_*` files,
    which makes importing a generated ParquetFilter wrapper by module index
    brittle.
  - This standalone module gives a stable import path
    (`pyopenms._parquet_filter.ParquetFilterBuilder`) for predicate pushdown
    builders.

Usage:
    from pyopenms._parquet_filter import ParquetFilterBuilder
    f = ParquetFilterBuilder().eq("precursor_id", 1318).andNext().eq("annotation", "y3^1")

The resulting builder can be passed to XICParquetFile query helpers or used
internally by higher-level Python query classes.
"""

from cython.operator cimport dereference as deref

# Import declarations from the pxd file
from pyopenms._parquet_filter cimport (
    _String, _ParquetFilterBuilder, ParquetFilterBuilder,
    shared_ptr, libcpp_vector, int64_t
)


cdef inline _String _to_string(object value):
    cdef bytes b
    if isinstance(value, bytes):
        b = value
    elif hasattr(value, "toString"):
        b = (<str>value.toString()).encode("utf-8")
    else:
        b = (<str>value).encode("utf-8")
    return _String(<char *>b)


cdef class ParquetFilterBuilder:
    """Typed filter builder for Parquet-backed readers.

    This class mirrors OpenMS::ParquetFilterBuilder and is intended to be used
    by higher-level query helpers (e.g., ChromatogramQuery). It is not tied to
    any specific `_pyopenms_*` split module.
    """

    def __cinit__(self):
        self.inst = shared_ptr[_ParquetFilterBuilder](new _ParquetFilterBuilder())

    def andNext(self):
        """andNext(self) -> ParquetFilterBuilder"""
        self.inst.get().andNext()
        return self

    def orNext(self):
        """orNext(self) -> ParquetFilterBuilder"""
        self.inst.get().orNext()
        return self

    def eq(self, column, value):
        if isinstance(value, (str, bytes)):
            self.inst.get().eq(_to_string(column), _to_string(value))
        else:
            self.inst.get().eq(_to_string(column), <int64_t>value)
        return self

    def ne(self, column, value):
        if isinstance(value, (str, bytes)):
            self.inst.get().ne(_to_string(column), _to_string(value))
        else:
            self.inst.get().ne(_to_string(column), <int64_t>value)
        return self

    def lt(self, column, value):
        if isinstance(value, (str, bytes)):
            self.inst.get().lt(_to_string(column), _to_string(value))
        else:
            self.inst.get().lt(_to_string(column), <int64_t>value)
        return self

    def le(self, column, value):
        if isinstance(value, (str, bytes)):
            self.inst.get().le(_to_string(column), _to_string(value))
        else:
            self.inst.get().le(_to_string(column), <int64_t>value)
        return self

    def gt(self, column, value):
        if isinstance(value, (str, bytes)):
            self.inst.get().gt(_to_string(column), _to_string(value))
        else:
            self.inst.get().gt(_to_string(column), <int64_t>value)
        return self

    def ge(self, column, value):
        if isinstance(value, (str, bytes)):
            self.inst.get().ge(_to_string(column), _to_string(value))
        else:
            self.inst.get().ge(_to_string(column), <int64_t>value)
        return self

    def in_(self, column, values):
        cdef libcpp_vector[_String] v1
        cdef libcpp_vector[int64_t] v2
        if values is None:
            raise TypeError("values must be a list or tuple")
        if len(values) == 0:
            raise ValueError("values must be non-empty")
        if isinstance(values[0], (str, bytes)):
            for item1 in values:
                v1.push_back(_to_string(item1))
            self.inst.get().in_(_to_string(column), v1)
        else:
            for item2 in values:
                v2.push_back(<int64_t>item2)
            self.inst.get().in_(_to_string(column), v2)
        return self

    def empty(self):
        return self.inst.get().empty()
