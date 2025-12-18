# cython: language_level=3

from libcpp.vector cimport vector as libcpp_vector
import numpy as np
cimport cython

cdef extern from * namespace "OpenMS::Math":
    """
    #include <OpenMS/MATH/STATISTICS/MultipleTesting.h>
    namespace OpenMS { namespace Math {

      inline std::vector<double> compute_model_fdr_double_i(std::vector<double> a) {
        return computeModelFdr<double>(a);
      }
      inline std::vector<double> compute_model_fdr_float_i(std::vector<float> a) {
        return computeModelFdr<float>(a);
      }
      inline std::vector<double> compute_model_fdr_int_i(std::vector<int> a) {
        return computeModelFdr<int>(a);
      }

      inline std::vector<double> pemp_double_i(std::vector<double> s, std::vector<double> s0) { return pEmp<double>(s, s0); }
      inline std::vector<double> pemp_float_i(std::vector<float> s, std::vector<float> s0) { return pEmp<float>(s, s0); }
      inline std::vector<double> pemp_int_i(std::vector<int> s, std::vector<int> s0) { return pEmp<int>(s, s0); }

    inline std::vector<double> qvalue_c(std::vector<double> p_values, double pi0, bool pfdr) { return qValue(p_values, pi0, pfdr); }
    inline std::vector<double> pnorm_c(std::vector<double> stat, std::vector<double> stat0) { return pNorm(stat, stat0); }
    inline Pi0Result pi0est_c(std::vector<double> p_values, std::vector<double> lambda_) { return pi0Est(p_values, lambda_); }
    inline std::vector<double> lfdr_c(std::vector<double> p_values, double pi0, bool trunc, bool monotone, std::string transf, double adj, double eps, size_t gridsize, double cut) { 
      return lfdr(p_values, pi0, trunc, monotone, transf, adj, eps, gridsize, cut); 
    }

    }}
    """
    libcpp_vector[double] compute_model_fdr_double_i(libcpp_vector[double] a) except +
    libcpp_vector[double] compute_model_fdr_float_i(libcpp_vector[float] a) except +
    libcpp_vector[double] compute_model_fdr_int_i(libcpp_vector[int] a) except +

    libcpp_vector[double] pemp_double_i(libcpp_vector[double] s, libcpp_vector[double] s0) except +
    libcpp_vector[double] pemp_float_i(libcpp_vector[float] s, libcpp_vector[float] s0) except +
    libcpp_vector[double] pemp_int_i(libcpp_vector[int] s, libcpp_vector[int] s0) except +
    libcpp_vector[double] qvalue_c(libcpp_vector[double] p_values, double pi0, bint pfdr) except +
    libcpp_vector[double] pnorm_c(libcpp_vector[double] stat, libcpp_vector[double] stat0) except +
    # Pi0Result declared later; forward declare return as object using same name in next extern block
    Pi0Result pi0est_c(libcpp_vector[double] p_values, libcpp_vector[double] lambda_) except +
    libcpp_vector[double] lfdr_c(libcpp_vector[double] p_values, double pi0, bint trunc, bint monotone, const char* transf, double adj, double eps, size_t gridsize, double cut) except +


cdef extern from "<OpenMS/MATH/STATISTICS/MultipleTesting.h>" namespace "OpenMS::Math":
    libcpp_vector[double] qValue(libcpp_vector[double] p_values, double pi0, bint pfdr) except +
    libcpp_vector[double] pNorm(libcpp_vector[double] stat, libcpp_vector[double] stat0) except +
    libcpp_vector[double] lfdr(libcpp_vector[double] p_values, double pi0, bint trunc, bint monotone, const char* transf, double adj, double eps, size_t gridsize, double cut) except +

    cdef cppclass Pi0Result:
        double pi0
        libcpp_vector[double] pi0_lambda
        libcpp_vector[double] lambda_
        bint pi0_smooth

    Pi0Result pi0Est(libcpp_vector[double] p_values, libcpp_vector[double] lambda_) except +


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

def _vec_to_numpy(libcpp_vector[double] v):
    cdef Py_ssize_t n = v.size()
    if n == 0:
        return np.empty((0,), dtype=float)
    cdef np.ndarray out = np.empty((n,), dtype=float)
    cdef double[:] outv = out
    for i in range(n):
        outv[i] = v[i]
    return out


def compute_model_fdr(values, dtype="float64"):
    """Compute model-based FDR estimates from posterior error probabilities."""
    if dtype == "float64":
        out = compute_model_fdr_double_i(_to_vec_double(values))
    elif dtype == "float32":
        out = compute_model_fdr_float_i(_to_vec_float(values))
    elif dtype == "int32":
        out = compute_model_fdr_int_i(_to_vec_int(values))
    else:
        raise ValueError(f"Unsupported dtype: {dtype!r}")
    return np.asarray(out, dtype=float)


def pemp(stat, stat0, dtype="float64"):
    """Empirical p-values (qvalue::empPvals translation)."""
    if dtype == "float64":
        out = pemp_double_i(_to_vec_double(stat), _to_vec_double(stat0))
    elif dtype == "float32":
        out = pemp_float_i(_to_vec_float(stat), _to_vec_float(stat0))
    elif dtype == "int32":
        out = pemp_int_i(_to_vec_int(stat), _to_vec_int(stat0))
    else:
        raise ValueError(f"Unsupported dtype: {dtype!r}")
    return np.asarray(out, dtype=float)


def qvalue(p_values, pi0=1.0, pfdr=False):
    """Compute q-values from p-values and pi0."""
    out = qvalue_c(_to_vec_double(p_values), <double> pi0, <bint> pfdr)
    return np.asarray(out, dtype=float)


def pnorm(stat, stat0):
    """Compute tail probabilities under a normal distribution fitted to stat0.
    Returns array of P(X > stat_i) where X ~ N(mu, sigma^2) with mu/sigma estimated from stat0."""
    out = pnorm_c(_to_vec_double(stat), _to_vec_double(stat0))
    return np.asarray(out, dtype=float)


def pi0est(p_values, lambda_=None):
    """Estimate pi0. Returns dict {pi0, pi0_lambda, lambda, pi0_smooth}."""
    if lambda_ is None:
        lamb = libcpp_vector[double]()
    else:
        lamb = _to_vec_double(lambda_)

    res = pi0est_c(_to_vec_double(p_values), lamb)

    return {
        'pi0': res.pi0,
        'pi0_lambda': _vec_to_numpy(res.pi0_lambda),
        'lambda': _vec_to_numpy(res.lambda_),
        'pi0_smooth': res.pi0_smooth != 0
    }


def lfdr(p_values, pi0, trunc=True, monotone=True, transf="probit", adj=1.5, eps=1e-8, gridsize=512, cut=3.0):
    """
    Estimate local false discovery rate (local FDR) from p-values.
    
    Parameters
    ----------
    p_values : array-like
        Vector of p-values from hypothesis tests
    pi0 : float
        Estimated proportion of true null hypotheses (from pi0est())
    trunc : bool, optional
        If True, truncate lfdr values to [0,1] range (default True)
    monotone : bool, optional
        If True, enforce monotonicity constraint (default True)
    transf : str, optional
        Transformation to apply: "probit" (inverse normal) or "logit" (log-odds) (default "probit")
    adj : float, optional
        Bandwidth adjustment factor (multiplied by automatic bandwidth selection) (default 1.5)
    eps : float, optional
        Small constant added to density estimates to avoid division by zero (default 1e-8)
    gridsize : int, optional
        Number of FFT grid points for KDE (default 512)
    cut : float, optional
        Grid extension factor in units of bandwidth (default 3.0)
    
    Returns
    -------
    ndarray
        Vector of local FDR values, one for each input p-value
    """
    transf_bytes = transf.encode('utf-8')
    out = lfdr_c(
        _to_vec_double(p_values), 
        <double> pi0, 
        <bint> trunc, 
        <bint> monotone,
        <const char*> transf_bytes,
        <double> adj,
        <double> eps,
        <size_t> gridsize,
        <double> cut
    )
    return np.asarray(out, dtype=float)
