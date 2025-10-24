from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/MATH/STATISTICS/RankData.h>" namespace "OpenMS::Math":

    cdef cppclass RankData "OpenMS::Math::RankData":
        RankData() except + nogil
        RankData(const RankData&) except + nogil

    # You may keep these if you like, but we won't *use* them below:
    ctypedef int RankData_Method "OpenMS::Math::RankData::Method"
    ctypedef int RankData_NaNPolicy "OpenMS::Math::RankData::NaNPolicy"

    # Correct, fully-qualified enum-class constants:
    cdef const int RANKDATA_AVERAGE   "OpenMS::Math::RankData::Method::Average"
    cdef const int RANKDATA_MIN       "OpenMS::Math::RankData::Method::Min"
    cdef const int RANKDATA_MAX       "OpenMS::Math::RankData::Method::Max"
    cdef const int RANKDATA_DENSE     "OpenMS::Math::RankData::Method::Dense"
    cdef const int RANKDATA_ORDINAL   "OpenMS::Math::RankData::Method::Ordinal"

    cdef const int RANKDATA_PROPAGATE "OpenMS::Math::RankData::NaNPolicy::Propagate"
    cdef const int RANKDATA_OMIT      "OpenMS::Math::RankData::NaNPolicy::Omit"
    cdef const int RANKDATA_RAISE     "OpenMS::Math::RankData::NaNPolicy::Raise"
