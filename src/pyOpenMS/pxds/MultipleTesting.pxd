from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
from libcpp.string cimport string as libcpp_utf8_string
from libcpp.string cimport string as libcpp_utf8_output_string
from Types cimport *

cdef extern from "<OpenMS/MATH/STATISTICS/MultipleTesting.h>" namespace "OpenMS::Math":

    cdef cppclass Pi0Result:
        # wrap-doc:
        #  Result of pi0 estimation for multiple testing correction.
        #
        #  Contains the estimated proportion of true null hypotheses (pi0),
        #  which is a key parameter in FDR estimation.
        #
        #  Attributes:
        #    pi0: Estimated proportion of true null hypotheses (0-1)
        #    pi0_lambda: Vector of pi0 estimates at each lambda threshold
        #    lambda_: Vector of lambda threshold values used
        #    pi0_smooth: Whether smoothing was successfully applied

        Pi0Result() except + nogil
        Pi0Result(Pi0Result &) except + nogil

        double pi0
        libcpp_vector[double] pi0_lambda
        libcpp_vector[double] lambda_
        bool pi0_smooth

    cdef cppclass MultipleTesting:
        # wrap-doc:
        #  Statistical functions for multiple testing correction.
        #
        #  Provides FDR estimation, q-value calculation, and local FDR computation
        #  based on the Storey-Tibshirani method.
        #
        #  Example:
        #    >>> import pyopenms
        #    >>> p_values = [0.001, 0.01, 0.05, 0.1, 0.5]
        #    >>> # Estimate pi0
        #    >>> pi0_result = pyopenms.MultipleTesting.pi0Est(p_values, [])
        #    >>> # Calculate q-values
        #    >>> q_values = pyopenms.MultipleTesting.qValue(p_values, pi0_result.pi0, False)
        #    >>> # Calculate local FDR
        #    >>> lfdr_values = pyopenms.MultipleTesting.lfdr(
        #    ...     p_values, pi0_result.pi0, True, True,
        #    ...     pyopenms.MultipleTesting.LfdrTransform.Probit)

        MultipleTesting() except + nogil  # wrap-ignore
        MultipleTesting(MultipleTesting &) except + nogil  # wrap-ignore

        # Enum conversion static methods
        @staticmethod
        libcpp_utf8_output_string pi0MethodToString(Pi0Method m) except + nogil
        # wrap-doc:Convert Pi0Method enum to string representation

        @staticmethod
        Pi0Method toPi0Method(const libcpp_utf8_string& s) except + nogil
        # wrap-doc:Convert string to Pi0Method enum (case-insensitive)

        @staticmethod
        libcpp_utf8_output_string lfdrTransformToString(LfdrTransform t) except + nogil
        # wrap-doc:Convert LfdrTransform enum to string representation

        @staticmethod
        LfdrTransform toLfdrTransform(const libcpp_utf8_string& s) except + nogil
        # wrap-doc:Convert string to LfdrTransform enum (case-insensitive)

        # Statistical methods
        @staticmethod
        libcpp_vector[double] qValue(libcpp_vector[double] p_values,
                                     double pi0,
                                     bool pfdr) except + nogil
        # wrap-doc:Calculate q-values (FDR-adjusted p-values) for multiple testing correction

        @staticmethod
        Pi0Result pi0Est(libcpp_vector[double] p_values,
                         libcpp_vector[double] lambda_,
                         Pi0Method method,
                         int smooth_df,
                         bool smooth_log_pi0) except + nogil
        # wrap-doc:Estimate the proportion of true null hypotheses (pi0)

        @staticmethod
        libcpp_vector[double] lfdr(libcpp_vector[double] p_values,
                                   double pi0,
                                   bool trunc,
                                   bool monotone,
                                   LfdrTransform transf,
                                   double adj,
                                   double eps,
                                   size_t gridsize,
                                   double cut) except + nogil
        # wrap-doc:Estimate local false discovery rate (local FDR) from p-values

        @staticmethod
        libcpp_vector[double] pNorm(libcpp_vector[double] stat,
                                    libcpp_vector[double] stat0) except + nogil
        # wrap-doc:Compute tail probabilities under a fitted normal distribution


cdef extern from "<OpenMS/MATH/STATISTICS/MultipleTesting.h>" namespace "OpenMS::Math::MultipleTesting":

    cdef enum class Pi0Method "OpenMS::Math::MultipleTesting::Pi0Method":
        # wrap-attach:
        #    MultipleTesting
        # wrap-doc:
        #  Method for estimating proportion of true null hypotheses (pi0).
        #
        #  - Smoother: Spline smoothing method (default)
        #  - Bootstrap: Bootstrap resampling method
        Smoother     # Spline smoothing method (default)
        Bootstrap    # Bootstrap resampling method

    cdef enum class LfdrTransform "OpenMS::Math::MultipleTesting::LfdrTransform":
        # wrap-attach:
        #    MultipleTesting
        # wrap-doc:
        #  Transformation for local FDR estimation.
        #
        #  - Probit: Inverse normal (Gaussian) transformation (default)
        #  - Logit: Log-odds transformation
        Probit       # Inverse normal (Gaussian) transformation (default)
        Logit        # Log-odds transformation
