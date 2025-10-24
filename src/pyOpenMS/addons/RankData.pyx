

from libcpp.vector cimport vector as libcpp_vector
import numpy as np
cimport cython

# We still want the constant values (correctly qualified)
cdef extern from "<OpenMS/MATH/STATISTICS/RankData.h>" namespace "OpenMS::Math":
    cdef const int RANKDATA_AVERAGE   "OpenMS::Math::RankData::Method::Average"
    cdef const int RANKDATA_MIN       "OpenMS::Math::RankData::Method::Min"
    cdef const int RANKDATA_MAX       "OpenMS::Math::RankData::Method::Max"
    cdef const int RANKDATA_DENSE     "OpenMS::Math::RankData::Method::Dense"
    cdef const int RANKDATA_ORDINAL   "OpenMS::Math::RankData::Method::Ordinal"

    cdef const int RANKDATA_PROPAGATE "OpenMS::Math::RankData::NaNPolicy::Propagate"
    cdef const int RANKDATA_OMIT      "OpenMS::Math::RankData::NaNPolicy::Omit"
    cdef const int RANKDATA_RAISE     "OpenMS::Math::RankData::NaNPolicy::Raise"

# Inline wrappers that take ints and cast to enum class on the C++ side
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
    if s == "average": return RANKDATA_AVERAGE
    if s == "min":     return RANKDATA_MIN
    if s == "max":     return RANKDATA_MAX
    if s == "dense":   return RANKDATA_DENSE
    if s == "ordinal": return RANKDATA_ORDINAL
    raise ValueError(f"Unknown method: {m!r}")

cdef int _parse_policy(object p):
    if isinstance(p, (int, np.integer)):
        return int(p)
    s = str(p).lower()
    if s == "propagate": return RANKDATA_PROPAGATE
    if s == "omit":      return RANKDATA_OMIT
    if s == "raise":     return RANKDATA_RAISE
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
