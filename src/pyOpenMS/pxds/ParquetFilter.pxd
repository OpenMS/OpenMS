from libcpp.vector cimport vector as libcpp_vector
from libcpp.memory cimport shared_ptr
from libc.stdint cimport int64_t
from String cimport String

cdef extern from "<OpenMS/FORMAT/ParquetFilter.h>" namespace "OpenMS":
    cdef cppclass ParquetFilter:
        ParquetFilter() except +
        ParquetFilter(const ParquetFilter&) except +

        ParquetFilter& andNext()
        ParquetFilter& orNext()

        ParquetFilter& eq(const String& column, int64_t value)
        ParquetFilter& ne(const String& column, int64_t value)
        ParquetFilter& lt(const String& column, int64_t value)
        ParquetFilter& le(const String& column, int64_t value)
        ParquetFilter& gt(const String& column, int64_t value)
        ParquetFilter& ge(const String& column, int64_t value)
        ParquetFilter& in_ "in"(const String& column, const libcpp_vector[int64_t]& values)

        ParquetFilter& eq(const String& column, const String& value)
        ParquetFilter& ne(const String& column, const String& value)
        ParquetFilter& lt(const String& column, const String& value)
        ParquetFilter& le(const String& column, const String& value)
        ParquetFilter& gt(const String& column, const String& value)
        ParquetFilter& ge(const String& column, const String& value)
        ParquetFilter& in_ "in"(const String& column, const libcpp_vector[String]& values)

        bint empty() const

cdef class PyParquetFilter:
    cdef shared_ptr[ParquetFilter] inst
