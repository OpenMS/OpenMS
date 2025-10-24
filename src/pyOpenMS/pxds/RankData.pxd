from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/MATH/STATISTICS/RankData.h>" namespace "OpenMS::Math":

    cdef cppclass RankData "OpenMS::Math::RankData":
        RankData() except + nogil
        RankData(const RankData&) except + nogil

    # Treat C++ scoped enums as ints
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

    # Free functions (BY VALUE) — use libcpp_vector, not vector
    libcpp_vector[double] rankdata_double(libcpp_vector[double] a,
                                          RankData_Method method,
                                          RankData_NaNPolicy policy) except + nogil

    libcpp_vector[double] rankdata_float(libcpp_vector[float] a,
                                         RankData_Method method,
                                         RankData_NaNPolicy policy) except + nogil

    libcpp_vector[double] rankdata_int(libcpp_vector[int] a,
                                       RankData_Method method,
                                       RankData_NaNPolicy policy) except + nogil
