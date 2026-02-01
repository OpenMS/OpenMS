from libc.stdint cimport int64_t
from libcpp.vector cimport vector as libcpp_vector
from libcpp.memory cimport shared_ptr

from String cimport String as _String
from ParquetFilter cimport ParquetFilter as _ParquetFilter
from ParquetFilter cimport PyParquetFilter as _PyParquetFilter
from ._pyopenms_1 cimport convString


cdef class PyParquetFilter:
    """
    Typed filter builder for Parquet-backed readers.
    """
    cdef shared_ptr[_ParquetFilter] inst

    def __cinit__(self):
        self.inst = shared_ptr[_ParquetFilter](new _ParquetFilter())

    def andNext(self):
        """
        andNext(self) -> PyParquetFilter
        """
        self.inst.get().andNext()
        return self

    def orNext(self):
        """
        orNext(self) -> PyParquetFilter
        """
        self.inst.get().orNext()
        return self

    def eq(self, column, value):
        if isinstance(value, (str, bytes, String)):
            self.inst.get().eq(deref((convString(column)).get()), deref((convString(value)).get()))
        else:
            self.inst.get().eq(deref((convString(column)).get()), <int64_t>value)
        return self

    def ne(self, column, value):
        if isinstance(value, (str, bytes, String)):
            self.inst.get().ne(deref((convString(column)).get()), deref((convString(value)).get()))
        else:
            self.inst.get().ne(deref((convString(column)).get()), <int64_t>value)
        return self

    def lt(self, column, value):
        if isinstance(value, (str, bytes, String)):
            self.inst.get().lt(deref((convString(column)).get()), deref((convString(value)).get()))
        else:
            self.inst.get().lt(deref((convString(column)).get()), <int64_t>value)
        return self

    def le(self, column, value):
        if isinstance(value, (str, bytes, String)):
            self.inst.get().le(deref((convString(column)).get()), deref((convString(value)).get()))
        else:
            self.inst.get().le(deref((convString(column)).get()), <int64_t>value)
        return self

    def gt(self, column, value):
        if isinstance(value, (str, bytes, String)):
            self.inst.get().gt(deref((convString(column)).get()), deref((convString(value)).get()))
        else:
            self.inst.get().gt(deref((convString(column)).get()), <int64_t>value)
        return self

    def ge(self, column, value):
        if isinstance(value, (str, bytes, String)):
            self.inst.get().ge(deref((convString(column)).get()), deref((convString(value)).get()))
        else:
            self.inst.get().ge(deref((convString(column)).get()), <int64_t>value)
        return self

    def in_(self, column, values):
        cdef libcpp_vector[_String] v1
        cdef libcpp_vector[int64_t] v2
        if values is None:
            raise TypeError("values must be a list or tuple")
        if len(values) == 0:
            raise ValueError("values must be non-empty")
        if isinstance(values[0], (str, bytes, String)):
            for item1 in values:
                v1.push_back(deref((convString(item1)).get()))
            self.inst.get().in_(deref((convString(column)).get()), v1)
        else:
            for item2 in values:
                v2.push_back(<int64_t>item2)
            self.inst.get().in_(deref((convString(column)).get()), v2)
        return self

    def empty(self):
        return self.inst.get().empty()
