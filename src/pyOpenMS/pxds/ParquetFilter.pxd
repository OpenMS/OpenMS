from libcpp.vector cimport vector as libcpp_vector
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

    cdef cppclass ParquetFilterBuilder:
        ParquetFilterBuilder() except +
        ParquetFilterBuilder(const ParquetFilterBuilder&) except +

        ParquetFilterBuilder& andNext()
        ParquetFilterBuilder& orNext()

        ParquetFilterBuilder& eq(const String& column, int64_t value)
        ParquetFilterBuilder& ne(const String& column, int64_t value)
        ParquetFilterBuilder& lt(const String& column, int64_t value)
        ParquetFilterBuilder& le(const String& column, int64_t value)
        ParquetFilterBuilder& gt(const String& column, int64_t value)
        ParquetFilterBuilder& ge(const String& column, int64_t value)
        ParquetFilterBuilder& in_ "in"(const String& column, const libcpp_vector[int64_t]& values)

        ParquetFilterBuilder& eq(const String& column, const String& value)
        ParquetFilterBuilder& ne(const String& column, const String& value)
        ParquetFilterBuilder& lt(const String& column, const String& value)
        ParquetFilterBuilder& le(const String& column, const String& value)
        ParquetFilterBuilder& gt(const String& column, const String& value)
        ParquetFilterBuilder& ge(const String& column, const String& value)
        ParquetFilterBuilder& in_ "in"(const String& column, const libcpp_vector[String]& values)

        const ParquetFilter& filter() const
        bint empty() const
