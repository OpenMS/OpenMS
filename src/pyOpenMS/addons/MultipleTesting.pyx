# cython: language_level=3

from libcpp.vector cimport vector as libcpp_vector
import numpy as np
cimport cython

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
    # Pi0Result declared later; forward declare return as object using same name in next extern block
    Pi0Result pi0est_c(libcpp_vector[double] p_values, libcpp_vector[double] lambda_, const char* pi0_method, int smooth_df, bint smooth_log_pi0) except +
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

    Pi0Result pi0Est(libcpp_vector[double] p_values, libcpp_vector[double] lambda_, const char* pi0_method, int smooth_df, bint smooth_log_pi0) except +


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


cpdef compute_model_fdr(values, dtype="float64"):
    """
    Compute model-based FDR estimates from posterior error probabilities.
    
    This function computes false discovery rate (FDR) estimates using a model-based
    approach. Given a vector of posterior error probabilities (PEPs), it calculates
    the expected FDR at each threshold.
    
    Parameters
    ----------
    values : array-like
        Vector of posterior error probabilities (PEPs) or similar scores
    dtype : str, optional
        Data type for computation: "float64" (default), "float32", or "int32"
    
    Returns
    -------
    ndarray
        Vector of FDR estimates corresponding to each input value
    
    Examples
    --------
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


cpdef pemp(stat, stat0, dtype="float64"):
    """
    Compute empirical p-values from test statistics and null distribution.
    
    This function calculates empirical p-values by comparing test statistics
    against a null distribution. For each test statistic, it computes the
    proportion of null statistics that are more extreme (greater or equal).
    This is equivalent to the qvalue::empPvals function from the R qvalue package.
    
    Parameters
    ----------
    stat : array-like
        Vector of test statistics from actual data
    stat0 : array-like
        Vector of test statistics from null distribution (e.g., from permutations)
    dtype : str, optional
        Data type for computation: "float64" (default), "float32", or "int32"
    
    Returns
    -------
    ndarray
        Vector of empirical p-values, one for each test statistic in `stat`
    
    Examples
    --------
    >>> test_stats = [2.5, 3.0, 1.5]
    >>> null_stats = [0.5, 1.0, 1.5, 2.0, 2.5]  # null distribution
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


cpdef qvalue(p_values, pi0=1.0, pfdr=False):
    """
    Compute q-values from p-values using the Storey-Tibshirani method.
    
    The q-value is the expected proportion of false positives among all features
    as significant as or more significant than the observed one. This function
    implements the q-value method from Storey and Tibshirani (2003).
    
    Parameters
    ----------
    p_values : array-like
        Vector of p-values from hypothesis tests
    pi0 : float, optional
        Estimated proportion of true null hypotheses. Should be between 0 and 1.
        Use pi0est() to estimate this value, or set to 1.0 (default) for 
        conservative estimates (equivalent to BH-FDR).
    pfdr : bool, optional
        If False (default), compute q-values (FDR). If True, compute positive
        FDR (pFDR) which conditions on rejecting at least one hypothesis.
    
    Returns
    -------
    ndarray
        Vector of q-values corresponding to each input p-value
    
    References
    ----------
    Storey, J. D., and Tibshirani, R. (2003). Statistical significance for
    genome-wide studies. Proceedings of the National Academy of Sciences,
    100: 9440-9445.
    
    Examples
    --------
    >>> p = [0.001, 0.01, 0.05, 0.1, 0.5]
    >>> q = qvalue(p, pi0=0.8)  # Assume 80% nulls
    """
    out = qvalue_c(_to_vec_double(p_values), <double> pi0, <bint> pfdr)
    return np.asarray(out, dtype=float)


cpdef pnorm(stat, stat0):
    """
    Compute parametric p-values under a normal distribution fitted to null statistics.
    
    This function fits a normal distribution to the null statistics (stat0) and
    computes tail probabilities for the test statistics (stat). Specifically, it
    estimates mu and sigma from stat0, then returns P(X > stat_i) where X ~ N(mu, sigma^2).
    This provides a parametric alternative to empirical p-value computation.
    
    Parameters
    ----------
    stat : array-like
        Vector of test statistics from actual data
    stat0 : array-like
        Vector of test statistics from null distribution used to fit the normal model.
        The mean and standard deviation are estimated from these values.
    
    Returns
    -------
    ndarray
        Vector of tail probabilities (p-values) for each test statistic,
        assuming a normal null distribution
    
    Notes
    -----
    This method assumes the null distribution is approximately normal. For
    non-normal null distributions, consider using pemp() for empirical p-values
    or applying an appropriate transformation first.
    
    Examples
    --------
    >>> test_stats = [2.5, 3.0, 1.5]
    >>> null_stats = np.random.randn(1000)  # standard normal null
    >>> pvals = pnorm(test_stats, null_stats)
    """
    out = pnorm_c(_to_vec_double(stat), _to_vec_double(stat0))
    return np.asarray(out, dtype=float)


cpdef pi0est(p_values, lambda_=None, pi0_method="smoother", smooth_df=3, smooth_log_pi0=False):
    """
    Estimate the proportion of true null hypotheses (pi0) using the Storey method.
    
    This function estimates pi0, the proportion of hypotheses that are truly null,
    which is a key parameter for q-value computation. The method evaluates the
    proportion of p-values above various thresholds (lambda) and may apply smoothing.
    
    Parameters
    ----------
    p_values : array-like
        Vector of p-values from hypothesis tests. Should be uniformly distributed
        under the null hypothesis.
    lambda_ : array-like, optional
        Vector of lambda thresholds for estimating pi0. If None (default), a
        sequence of values from 0.05 to 0.95 will be used automatically.
    pi0_method : str, optional
        Method for pi0 estimation: "smoother" (default) applies smoothing spline,
        "bootstrap" uses bootstrap resampling (not yet implemented).
    smooth_df : int, optional
        Degrees of freedom for the smoothing spline (default: 3).
        Only used when pi0_method="smoother".
    smooth_log_pi0 : bool, optional
        If True, apply smoothing on log(pi0) scale (default: False).
        Only used when pi0_method="smoother".
    
    Returns
    -------
    dict
        Dictionary containing:
        - 'pi0' : float
            The final estimated proportion of true nulls (may be smoothed)
        - 'pi0_lambda' : ndarray
            Raw pi0 estimates at each lambda threshold
        - 'lambda' : ndarray
            The lambda thresholds used for estimation
        - 'pi0_smooth' : bool
            Whether smoothing was applied to obtain the final pi0 estimate
    
    Notes
    -----
    The pi0 estimate is crucial for q-value calculation. Conservative analysis
    uses pi0=1.0, while estimated pi0 < 1 can increase power if most hypotheses
    are truly null.
    
    References
    ----------
    Storey, J. D. (2002). A direct approach to false discovery rates.
    Journal of the Royal Statistical Society, Series B, 64: 479-498.
    
    Examples
    --------
    >>> p = np.concatenate([np.random.uniform(0, 0.01, 100),  # 100 true positives
    ...                      np.random.uniform(0, 1, 900)])    # 900 true nulls
    >>> result = pi0est(p)
    >>> print(f"Estimated pi0: {result['pi0']:.2f}")  # Should be ~0.9
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


cpdef lfdr(p_values, pi0, trunc=True, monotone=True, transf="probit", adj=1.5, eps=1e-8, gridsize=512, cut=3.0):
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
