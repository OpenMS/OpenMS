

from libcpp.vector cimport vector
import numpy as np
cimport cython

# Bring C++ API in
from pyopenms.pyopenms cimport RankData, rankdata_double, rankdata_float, rankdata_int

cdef RankData.Method _parse_method(object m):
    if isinstance(m, RankData.Method):
        return <RankData.Method>m
    s = str(m).lower()
    if s == "average": return RankData.Average
    if s == "min":     return RankData.Min
    if s == "max":     return RankData.Max
    if s == "dense":   return RankData.Dense
    if s == "ordinal": return RankData.Ordinal
    raise ValueError("Unknown method: %r" % m)

cdef RankData.NaNPolicy _parse_policy(object p):
    if isinstance(p, RankData.NaNPolicy):
        return <RankData.NaNPolicy>p
    s = str(p).lower()
    if s == "propagate": return RankData.Propagate
    if s == "omit":      return RankData.Omit
    if s == "raise":     return RankData.Raise
    raise ValueError("Unknown nan_policy: %r" % p)

@cython.boundscheck(False)
@cython.wraparound(False)
def _to_vec_double(arr):
    cdef vector[double] v
    cdef Py_ssize_t i, n
    a = np.asarray(arr, dtype=float)  # accepts list/tuple/np.ndarray
    n = a.size
    v.reserve(n)
    cdef double[:] mv = a.ravel()  # memoryview for speed
    for i in range(n):
        v.push_back(mv[i])
    return v

@cython.boundscheck(False)
@cython.wraparound(False)
def _to_vec_float(arr):
    cdef vector[float] v
    cdef Py_ssize_t i, n
    a = np.asarray(arr, dtype=float)
    n = a.size
    v.reserve(n)
    cdef double[:] mv = a.ravel()
    for i in range(n):
        v.push_back(<float>mv[i])
    return v

@cython.boundscheck(False)
@cython.wraparound(False)
def _to_vec_int(arr):
    cdef vector[int] v
    cdef Py_ssize_t i, n
    a = np.asarray(arr, dtype=int)
    n = a.size
    v.reserve(n)
    cdef long[:] mv = a.ravel()
    for i in range(n):
        v.push_back(<int>mv[i])
    return v

def rankdata(values, method="average", nan_policy="propagate", dtype="float64"):
    """
    Rank values (1-based) with SciPy-like tie handling.

    Parameters
    ----------
    values : array-like
        Input sequence
    method : {'average','min','max','dense','ordinal'} or RankData.Method
    nan_policy : {'propagate','omit','raise'} or RankData.NaNPolicy
    dtype : {'float64','float32','int32'}
        Dispatches to the underlying C++ specialization.

    Returns
    -------
    numpy.ndarray (float64)
        Ranks (1-based)
    """
    m = _parse_method(method)
    p = _parse_policy(nan_policy)

    if dtype == "float64":
        vd = _to_vec_double(values)
        out = rankdata_double(vd, m, p)
    elif dtype == "float32":
        vf = _to_vec_float(values)
        out = rankdata_float(vf, m, p)
    elif dtype == "int32":
        vi = _to_vec_int(values)
        out = rankdata_int(vi, m, p)
    else:
        raise ValueError("Unsupported dtype: %r" % dtype)

    # autowrap converts vector[double] -> Python list; cast to np.array
    return np.asarray(out, dtype=float)
