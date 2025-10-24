

from libcpp.vector cimport vector as libcpp_vector
import numpy as np
cimport cython

from pyopenms.pyopenms cimport (
    RankData,
    RankData_Method, RankData_NaNPolicy,
    RANKDATA_AVERAGE, RANKDATA_MIN, RANKDATA_MAX, RANKDATA_DENSE, RANKDATA_ORDINAL,
    RANKDATA_PROPAGATE, RANKDATA_OMIT, RANKDATA_RAISE,
    rankdata_double, rankdata_float, rankdata_int
)

cdef RankData_Method _parse_method(object m):
    if isinstance(m, (int, np.integer)):
        return <RankData_Method> int(m)
    s = str(m).lower()
    if s == "average": return <RankData_Method> RankData.Average
    if s == "min":     return <RankData_Method> RankData.Min
    if s == "max":     return <RankData_Method> RankData.Max
    if s == "dense":   return <RankData_Method> RankData.Dense
    if s == "ordinal": return <RankData_Method> RankData.Ordinal
    raise ValueError(f"Unknown method: {m!r}")

cdef RankData_NaNPolicy _parse_policy(object p):
    if isinstance(p, (int, np.integer)):
        return <RankData_NaNPolicy> int(p)
    s = str(p).lower()
    if s == "propagate": return <RankData_NaNPolicy> RankData.Propagate
    if s == "omit":      return <RankData_NaNPolicy> RankData.Omit
    if s == "raise":     return <RankData_NaNPolicy> RankData.Raise
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
