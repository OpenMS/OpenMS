

from libcpp.vector cimport vector as libcpp_vector
import numpy as np
cimport cython

# Directly declare the C++ API we need. No cimport from pyopenms.pyopenms.
cdef extern from "<OpenMS/MATH/STATISTICS/RankData.h>" namespace "OpenMS::Math":
    ctypedef int RankData_Method    "OpenMS::Math::RankData::Method"
    ctypedef int RankData_NaNPolicy "OpenMS::Math::RankData::NaNPolicy"

    # Scoped-enum constants (as ints)
    cdef const int RANKDATA_AVERAGE   "OpenMS::Math::RankData::Average"
    cdef const int RANKDATA_MIN       "OpenMS::Math::RankData::Min"
    cdef const int RANKDATA_MAX       "OpenMS::Math::RankData::Max"
    cdef const int RANKDATA_DENSE     "OpenMS::Math::RankData::Dense"
    cdef const int RANKDATA_ORDINAL   "OpenMS::Math::RankData::Ordinal"

    cdef const int RANKDATA_PROPAGATE "OpenMS::Math::RankData::Propagate"
    cdef const int RANKDATA_OMIT      "OpenMS::Math::RankData::Omit"
    cdef const int RANKDATA_RAISE     "OpenMS::Math::RankData::Raise"

    # Free functions taking/returning libcpp_vector (BY VALUE)
    libcpp_vector[double] rankdata_double(libcpp_vector[double] a,
                                          RankData_Method method,
                                          RankData_NaNPolicy policy) except +
    libcpp_vector[double] rankdata_float (libcpp_vector[float]  a,
                                          RankData_Method method,
                                          RankData_NaNPolicy policy) except +
    libcpp_vector[double] rankdata_int   (libcpp_vector[int]    a,
                                          RankData_Method method,
                                          RankData_NaNPolicy policy) except +

cdef RankData_Method _parse_method(object m):
    if isinstance(m, (int, np.integer)):
        return <RankData_Method> int(m)
    s = str(m).lower()
    if s == "average": return <RankData_Method> RANKDATA_AVERAGE
    if s == "min":     return <RankData_Method> RANKDATA_MIN
    if s == "max":     return <RankData_Method> RANKDATA_MAX
    if s == "dense":   return <RankData_Method> RANKDATA_DENSE
    if s == "ordinal": return <RankData_Method> RANKDATA_ORDINAL
    raise ValueError(f"Unknown method: {m!r}")

cdef RankData_NaNPolicy _parse_policy(object p):
    if isinstance(p, (int, np.integer)):
        return <RankData_NaNPolicy> int(p)
    s = str(p).lower()
    if s == "propagate": return <RankData_NaNPolicy> RANKDATA_PROPAGATE
    if s == "omit":      return <RankData_NaNPolicy> RANKDATA_OMIT
    if s == "raise":     return <RankData_NaNPolicy> RANKDATA_RAISE
    raise ValueError(f"Unknown nan_policy: {p!r}")

@cython.boundscheck(False)
@cython.wraparound(False)
def _to_vec_double(arr):
    cdef libcpp_vector[double] v
    a = np.asarray(arr, dtype=float)
    cdef double[:] mv = a.ravel()
    v.reserve(mv.shape[0])
    for i in range(mv.shape[0]):
        v.push_back(mv[i])
    return v

@cython.boundscheck(False)
@cython.wraparound(False)
def _to_vec_float(arr):
    cdef libcpp_vector[float] v
    a = np.asarray(arr, dtype=float)
    cdef double[:] mv = a.ravel()
    v.reserve(mv.shape[0])
    for i in range(mv.shape[0]):
        v.push_back(<float> mv[i])
    return v

@cython.boundscheck(False)
@cython.wraparound(False)
def _to_vec_int(arr):
    cdef libcpp_vector[int] v
    a = np.asarray(arr, dtype=int)
    cdef long[:] mv = a.ravel()
    v.reserve(mv.shape[0])
    for i in range(mv.shape[0]):
        v.push_back(<int> mv[i])
    return v

def rankdata(values, method="average", nan_policy="propagate", dtype="float64"):
    cdef RankData_Method m = _parse_method(method)
    cdef RankData_NaNPolicy p = _parse_policy(nan_policy)

    if dtype == "float64":
        out = rankdata_double(_to_vec_double(values), m, p)
    elif dtype == "float32":
        out = rankdata_float(_to_vec_float(values), m, p)
    elif dtype == "int32":
        out = rankdata_int(_to_vec_int(values), m, p)
    else:
        raise ValueError(f"Unsupported dtype: {dtype!r}")

    return np.asarray(out, dtype=float)
