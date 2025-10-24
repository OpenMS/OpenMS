from libcpp.vector cimport vector

cdef extern from "<OpenMS/MATH/STATISTICS/RankData.h>" namespace "OpenMS::Math":
    cdef cppclass RankData "OpenMS::Math::RankData":
        # default & copy ctors (recommended by pyOpenMS docs even if unused)
        RankData() nogil except +
        RankData(RankData&) nogil except +

        cdef enum Method:
            Average
            Min
            Max
            Dense
            Ordinal

        cdef enum NaNPolicy:
            Propagate
            Omit
            Raise

    # Expose concrete overloads and *attach* them to RankData as static methods
    vector[double] rankdata_double(vector[double] a,
                                   RankData.Method method,
                                   RankData.NaNPolicy policy) nogil except +  # wrap-attach:RankData
    vector[double] rankdata_float(vector[float] a,
                                  RankData.Method method,
                                  RankData.NaNPolicy policy) nogil except +   # wrap-attach:RankData
    vector[int]    _identity_vec_int(vector[int] a) nogil except +  # wrap-ignore

    vector[double] rankdata_int(vector[int] a,
                                RankData.Method method,
                                RankData.NaNPolicy policy) nogil except +     # wrap-attach:RankData
