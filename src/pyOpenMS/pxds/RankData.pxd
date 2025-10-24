from libcpp.vector cimport vector

cdef extern from "<OpenMS/MATH/STATISTICS/RankData.h>" namespace "OpenMS::Math":

    # C++ class (no enums inside the class)
    cdef cppclass RankData "OpenMS::Math::RankData":
        RankData() except + nogil
        RankData(const RankData&) except + nogil

    # Treat C++ scoped enums as ints (robust with autowrap)
    ctypedef int RankData_Method "OpenMS::Math::RankData::Method"
    ctypedef int RankData_NaNPolicy "OpenMS::Math::RankData::NaNPolicy"

    # Export constants
    cdef const int RANKDATA_AVERAGE  "OpenMS::Math::RankData::Average"
    cdef const int RANKDATA_MIN      "OpenMS::Math::RankData::Min"
    cdef const int RANKDATA_MAX      "OpenMS::Math::RankData::Max"
    cdef const int RANKDATA_DENSE    "OpenMS::Math::RankData::Dense"
    cdef const int RANKDATA_ORDINAL  "OpenMS::Math::RankData::Ordinal"

    cdef const int RANKDATA_PROPAGATE "OpenMS::Math::RankData::Propagate"
    cdef const int RANKDATA_OMIT      "OpenMS::Math::RankData::Omit"
    cdef const int RANKDATA_RAISE     "OpenMS::Math::RankData::Raise"

    # Free functions (BY VALUE) — match the C++ free wrappers above
    vector[double] rankdata_double(vector[double] a,
                                   RankData_Method method,
                                   RankData_NaNPolicy policy) except + nogil

    vector[double] rankdata_float(vector[float] a,
                                  RankData_Method method,
                                  RankData_NaNPolicy policy) except + nogil

    vector[double] rankdata_int(vector[int] a,
                                RankData_Method method,
                                RankData_NaNPolicy policy) except + nogil
