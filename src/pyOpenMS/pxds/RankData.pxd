from libcpp.vector cimport vector as libcpp_vector_as_np


cdef extern from "<OpenMS/MATH/STATISTICS/RankData.h>" namespace "OpenMS::Math":

    cdef cppclass RankData:
        # wrap-doc:
        #  Rank items (1-based) with SciPy-like tie and NaN handling.
        #
        #  This wrapper exposes the concrete forwarders rankdata_double/float/int
        #  and the enum classes RankData.Method and RankData.NaNPolicy.

        RankData() except + nogil  # wrap-ignore
        RankData(RankData &) except + nogil  # wrap-ignore

        @staticmethod
        libcpp_vector_as_np[double] rankdata_double(const libcpp_vector_as_np[double]& a,
                                                    Method method,
                                                    NaNPolicy nan_policy) except + nogil

        @staticmethod
        libcpp_vector_as_np[double] rankdata_float(const libcpp_vector_as_np[float]& a,
                                                   Method method,
                                                   NaNPolicy nan_policy) except + nogil

        @staticmethod
        libcpp_vector_as_np[double] rankdata_int(const libcpp_vector_as_np[int]& a,
                                                 Method method,
                                                 NaNPolicy nan_policy) except + nogil


cdef extern from "<OpenMS/MATH/STATISTICS/RankData.h>" namespace "OpenMS::Math::RankData":

    cdef enum class Method "OpenMS::Math::RankData::Method":
        # wrap-attach:
        #    RankData
        # wrap-doc:
        #  Tie handling method for ranking.
        Average
        Min
        Max
        Dense
        Ordinal

    cdef enum class NaNPolicy "OpenMS::Math::RankData::NaNPolicy":
        # wrap-attach:
        #    RankData
        # wrap-doc:
        #  NaN handling policy for ranking.
        Propagate
        Omit
        Raise
