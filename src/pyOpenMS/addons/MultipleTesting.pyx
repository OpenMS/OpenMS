# cython: language_level=3, warn.unreachable=False, warn.undeclared=False

from libcpp.vector cimport vector as libcpp_vector
import numpy as np
cimport cython

# First define the types we need
cdef extern from "<OpenMS/MATH/STATISTICS/MultipleTesting.h>" namespace "OpenMS::Math":
    cdef cppclass Pi0Result:
        double pi0
        libcpp_vector[double] pi0_lambda
        libcpp_vector[double] lambda_
        bint pi0_smooth

    libcpp_vector[double] qValue(libcpp_vector[double] p_values, double pi0, bint pfdr) except +
    libcpp_vector[double] pNorm(libcpp_vector[double] stat, libcpp_vector[double] stat0) except +
    libcpp_vector[double] lfdr(libcpp_vector[double] p_values, double pi0, bint trunc, bint monotone, const char* transf, double adj, double eps, size_t gridsize, double cut) except +
    Pi0Result pi0Est(libcpp_vector[double] p_values, libcpp_vector[double] lambda_, const char* pi0_method, int smooth_df, bint smooth_log_pi0) except +

# Now inline wrapper functions
cdef extern from * namespace "OpenMS::Math":
    """
    #include <OpenMS/MATH/STATISTICS/MultipleTesting.h>
    namespace OpenMS { namespace Math {

      inline std::vector<double> compute_model_fdr_double_i(std::vector<double> a) {
        return computeModelFDR<double>(a);
      }
      inline std::vector<double> compute_model_fdr_float_i(std::vector<float> a) {
        return computeModelFDR<float>(a);
      }
      inline std::vector<double> compute_model_fdr_int_i(std::vector<int> a) {
        return computeModelFDR<int>(a);
      }

      inline std::vector<double> pemp_double_i(std::vector<double> s, std::vector<double> s0) { return pEmp<double>(s, s0); }
      inline std::vector<double> pemp_float_i(std::vector<float> s, std::vector<float> s0) { return pEmp<float>(s, s0); }
      inline std::vector<double> pemp_int_i(std::vector<int> s, std::vector<int> s0) { return pEmp<int>(s, s0); }

    inline std::vector<double> qvalue_c(std::vector<double> p_values, double pi0, bool pfdr) { return qValue(p_values, pi0, pfdr); }
    inline std::vector<double> pnorm_c(std::vector<double> stat, std::vector<double> stat0) { return pNorm(stat, stat0); }
    inline Pi0Result pi0est_c(std::vector<double> p_values, std::vector<double> lambda_, std::string pi0_method, int smooth_df, bool smooth_log_pi0) { 
      return pi0Est(p_values, lambda_, pi0_method, smooth_df, smooth_log_pi0); 
    }
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
    Pi0Result pi0est_c(libcpp_vector[double] p_values, libcpp_vector[double] lambda_, const char* pi0_method, int smooth_df, bint smooth_log_pi0) except +
    libcpp_vector[double] lfdr_c(libcpp_vector[double] p_values, double pi0, bint trunc, bint monotone, const char* transf, double adj, double eps, size_t gridsize, double cut) except +


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
    a = np.asarray(arr, dtype=np.float32)
    cdef float[:] mv = a.ravel()
    v.reserve(mv.shape[0])
    for i in range(mv.shape[0]):
        v.push_back(mv[i])
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
    """
    Compute model-based FDR estimates from posterior error probabilities.

    Given a vector of posterior error probabilities (PEPs), calculates
    the expected FDR at each threshold using a model-based approach.

    Args:
        values (array-like): Vector of posterior error probabilities (PEPs)
        dtype (str): Data type for computation: "float64" (default), "float32", or "int32"

    Returns:
        ndarray: Vector of FDR estimates corresponding to each input value

    Example:
        >>> peps = [0.01, 0.05, 0.1, 0.2, 0.5]
        >>> fdrs = compute_model_fdr(peps)
    """
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
    """
    Compute empirical p-values from test statistics and null distribution.

    Calculates empirical p-values by comparing test statistics against a null
    distribution. For each test statistic, computes the proportion of null
    statistics that are more extreme (greater or equal). Equivalent to the
    qvalue::empPvals function from the R qvalue package.

    Args:
        stat (array-like): Vector of test statistics from actual data
        stat0 (array-like): Vector of null statistics (e.g., from permutations)
        dtype (str): Data type: "float64" (default), "float32", or "int32"

    Returns:
        ndarray: Vector of empirical p-values, one for each test statistic

    Example:
        >>> test_stats = [2.5, 3.0, 1.5]
        >>> null_stats = [0.5, 1.0, 1.5, 2.0, 2.5]
        >>> pvals = pemp(test_stats, null_stats)
    """
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
    """
    Compute q-values from p-values using the Storey-Tibshirani method.

    The q-value is the expected proportion of false positives among all features
    as significant as or more significant than the observed one. Implements the
    method from Storey and Tibshirani (2003), PNAS 100:9440-9445.

    Args:
        p_values (array-like): Vector of p-values from hypothesis tests
        pi0 (float): Proportion of true nulls (0-1). Use pi0est() or 1.0 (default).
        pfdr (bool): If False (default), compute q-values. If True, compute pFDR.

    Returns:
        ndarray: Vector of q-values corresponding to each input p-value

    Example:
        >>> p = [0.001, 0.01, 0.05, 0.1, 0.5]
        >>> q = qvalue(p, pi0=0.8)
    """
    out = qvalue_c(_to_vec_double(p_values), <double> pi0, <bint> pfdr)
    return np.asarray(out, dtype=float)


def pnorm(stat, stat0):
    """
    Compute parametric p-values under a normal distribution fitted to null statistics.

    Fits a normal distribution to the null statistics (stat0) and computes tail
    probabilities for the test statistics (stat). Estimates mu and sigma from
    stat0, then returns P(X > stat_i) where X ~ N(mu, sigma^2).

    This method assumes the null distribution is approximately normal. For
    non-normal nulls, consider using pemp() for empirical p-values.

    Args:
        stat (array-like): Vector of test statistics from actual data
        stat0 (array-like): Vector of null statistics used to fit the normal model

    Returns:
        ndarray: Vector of tail probabilities (p-values) for each test statistic

    Example:
        >>> test_stats = [2.5, 3.0, 1.5]
        >>> null_stats = np.random.randn(1000)
        >>> pvals = pnorm(test_stats, null_stats)
    """
    out = pnorm_c(_to_vec_double(stat), _to_vec_double(stat0))
    return np.asarray(out, dtype=float)


def pi0est(p_values, lambda_=None, pi0_method="smoother", smooth_df=3, smooth_log_pi0=False):
    """
    Estimate the proportion of true null hypotheses (pi0) using the Storey method.

    Estimates pi0, the proportion of hypotheses that are truly null, which is
    a key parameter for q-value computation. Based on Storey (2002), JRSS-B 64:479-498.

    Args:
        p_values (array-like): Vector of p-values from hypothesis tests
        lambda_ (array-like): Lambda thresholds. If None, uses 0.05 to 0.95.
        pi0_method (str): "smoother" (default) or "bootstrap"
        smooth_df (int): Degrees of freedom for smoothing spline (default: 3)
        smooth_log_pi0 (bool): If True, smooth on log(pi0) scale (default: False)

    Returns:
        dict: Dictionary with 'pi0', 'pi0_lambda', 'lambda', 'pi0_smooth' keys

    Example:
        >>> import numpy as np
        >>> p = np.concatenate([np.random.uniform(0, 0.01, 100),
        ...                     np.random.uniform(0, 1, 900)])
        >>> result = pi0est(p)
        >>> print(f"Estimated pi0: {result['pi0']:.2f}")
    """
    if lambda_ is None:
        lamb = libcpp_vector[double]()
    else:
        lamb = _to_vec_double(lambda_)

    pi0_method_bytes = pi0_method.encode('utf-8')
    res = pi0est_c(_to_vec_double(p_values), lamb, <const char*> pi0_method_bytes, 
                   <int> smooth_df, <bint> smooth_log_pi0)

    return {
        'pi0': res.pi0,
        'pi0_lambda': _vec_to_numpy(res.pi0_lambda),
        'lambda': _vec_to_numpy(res.lambda_),
        'pi0_smooth': res.pi0_smooth != 0
    }


def lfdr(p_values, pi0, trunc=True, monotone=True, transf="probit", adj=1.5, eps=1e-8, gridsize=512, cut=3.0):
    """
    Estimate local false discovery rate (local FDR) from p-values.

    Args:
        p_values (array-like): Vector of p-values from hypothesis tests
        pi0 (float): Estimated proportion of true null hypotheses (from pi0est())
        trunc (bool): If True, truncate lfdr values to [0,1] range (default: True)
        monotone (bool): If True, enforce monotonicity constraint (default: True)
        transf (str): Transformation: "probit" (inverse normal) or "logit" (default: "probit")
        adj (float): Bandwidth adjustment factor (default: 1.5)
        eps (float): Small constant to avoid division by zero (default: 1e-8)
        gridsize (int): Number of FFT grid points for KDE (default: 512)
        cut (float): Grid extension factor in bandwidth units (default: 3.0)

    Returns:
        ndarray: Vector of local FDR values, one for each input p-value
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
