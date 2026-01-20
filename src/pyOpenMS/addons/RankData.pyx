

# -*- coding: utf-8 -*-
# cython: language_level=3, warn.unreachable=False, warn.undeclared=False

from libcpp.vector cimport vector as libcpp_vector
import numpy as np
cimport cython

# Pull the scoped-enum constants straight from the header
cdef extern from "<OpenMS/MATH/STATISTICS/RankData.h>" namespace "OpenMS::Math":
    cdef const int RANKDATA_AVERAGE   "OpenMS::Math::RankData::Method::Average"
    cdef const int RANKDATA_MIN       "OpenMS::Math::RankData::Method::Min"
    cdef const int RANKDATA_MAX       "OpenMS::Math::RankData::Method::Max"
    cdef const int RANKDATA_DENSE     "OpenMS::Math::RankData::Method::Dense"
    cdef const int RANKDATA_ORDINAL   "OpenMS::Math::RankData::Method::Ordinal"

    cdef const int RANKDATA_PROPAGATE "OpenMS::Math::RankData::NaNPolicy::Propagate"
    cdef const int RANKDATA_OMIT      "OpenMS::Math::RankData::NaNPolicy::Omit"
    cdef const int RANKDATA_RAISE     "OpenMS::Math::RankData::NaNPolicy::Raise"

# Inline shims that take ints (so we don't expose enum classes to Python)
cdef extern from * namespace "OpenMS::Math":
    """
    #include <OpenMS/MATH/STATISTICS/RankData.h>
    namespace OpenMS { namespace Math {

      inline std::vector<double>
      rankdata_double_i(std::vector<double> a, int m, int p) {
        return rankdata_double(
          a,
          static_cast<RankData::Method>(m),
          static_cast<RankData::NaNPolicy>(p)
        );
      }

      inline std::vector<double>
      rankdata_float_i(std::vector<float> a, int m, int p) {
        return rankdata_float(
          a,
          static_cast<RankData::Method>(m),
          static_cast<RankData::NaNPolicy>(p)
        );
      }

      inline std::vector<double>
      rankdata_int_i(std::vector<int> a, int m, int p) {
        return rankdata_int(
          a,
          static_cast<RankData::Method>(m),
          static_cast<RankData::NaNPolicy>(p)
        );
      }

    }} // namespace OpenMS::Math
    """
    libcpp_vector[double] rankdata_double_i(libcpp_vector[double] a, int method, int policy) except +
    libcpp_vector[double] rankdata_float_i (libcpp_vector[float]  a, int method, int policy) except +
    libcpp_vector[double] rankdata_int_i   (libcpp_vector[int]    a, int method, int policy) except +

cdef int _parse_method(object m):
    if isinstance(m, (int, np.integer)):
        return int(m)
    s = str(m).lower()
    if s == "average": return <int> RANKDATA_AVERAGE
    if s == "min":     return <int> RANKDATA_MIN
    if s == "max":     return <int> RANKDATA_MAX
    if s == "dense":   return <int> RANKDATA_DENSE
    if s == "ordinal": return <int> RANKDATA_ORDINAL
    raise ValueError(f"Unknown method: {m!r}")

cdef int _parse_policy(object p):
    if isinstance(p, (int, np.integer)):
        return int(p)
    s = str(p).lower()
    if s == "propagate": return <int> RANKDATA_PROPAGATE
    if s == "omit":      return <int> RANKDATA_OMIT
    if s == "raise":     return <int> RANKDATA_RAISE
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
    a = np.asarray(arr, dtype=np.int32)
    cdef int[:] mv = a.ravel()
    v.reserve(mv.shape[0])
    for i in range(mv.shape[0]):
        v.push_back(mv[i])
    return v

def rankdata(values, method="average", nan_policy="propagate", dtype="float64"):
    """
    Rank items (1-based) with SciPy-like tie handling.

    Mirrors scipy.stats.rankdata semantics on a flattened 1D view of the input.
    Ranks start at 1. Dispatches to OpenMS C++ implementation (OpenMS::Math::RankData).

    Args:
        values (array-like): Input data. Any shape; flattened before ranking.
        method (str): Tie handling: "average" (default), "min", "max", "dense", "ordinal"
        nan_policy (str): NaN handling: "propagate" (default), "omit", or "raise"
        dtype (str): Execution type: "float64" (default), "float32", or "int32"

    Returns:
        ndarray: 1-based ranks with shape (values.size,) and dtype float64

    Example:
        >>> import pyopenms as poms
        >>> poms.rankdata([0, 2, 3, 2], "average")
        array([1. , 2.5, 4. , 2.5])
        >>> poms.rankdata([0, 2, 3, 2], "min")
        array([1., 2., 4., 2.])
        >>> poms.rankdata([0, 2, 3, 2], "dense")
        array([1., 2., 3., 2.])
    """
    cdef int m = _parse_method(method)
    cdef int p = _parse_policy(nan_policy)

    if dtype == "float64":
        out = rankdata_double_i(_to_vec_double(values), m, p)
    elif dtype == "float32":
        out = rankdata_float_i(_to_vec_float(values), m, p)
    elif dtype == "int32":
        out = rankdata_int_i(_to_vec_int(values), m, p)
    else:
        raise ValueError(f"Unsupported dtype: {dtype!r}")

    return np.asarray(out, dtype=float)