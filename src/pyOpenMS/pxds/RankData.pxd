from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/MATH/STATISTICS/RankData.h>" namespace "OpenMS::Math":
    
    cdef cppclass RankData:
        # wrap-ignore
        # Dummy class to anchor the RankData.pyx addon module
        # The actual functionality is provided via Python functions in the addon
        RankData() except + nogil  # wrap-ignore
    
    # Enum constants for ranking and NaN policies
    cdef const int RANKDATA_AVERAGE   "OpenMS::Math::RankData::Method::Average"
    cdef const int RANKDATA_MIN       "OpenMS::Math::RankData::Method::Min"
    cdef const int RANKDATA_MAX       "OpenMS::Math::RankData::Method::Max"
    cdef const int RANKDATA_DENSE     "OpenMS::Math::RankData::Method::Dense"
    cdef const int RANKDATA_ORDINAL   "OpenMS::Math::RankData::Method::Ordinal"

    cdef const int RANKDATA_PROPAGATE "OpenMS::Math::RankData::NaNPolicy::Propagate"
    cdef const int RANKDATA_OMIT      "OpenMS::Math::RankData::NaNPolicy::Omit"
    cdef const int RANKDATA_RAISE     "OpenMS::Math::RankData::NaNPolicy::Raise"
