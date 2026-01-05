from Types cimport *
from libcpp.string cimport string as libcpp_utf8_string

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MRMRTNormalizer.h>" namespace "OpenMS":

    cdef cppclass MRMRTNormalizer:

        MRMRTNormalizer() except + nogil
        MRMRTNormalizer(MRMRTNormalizer &) except + nogil  # compiler

        @staticmethod
        libcpp_vector[libcpp_pair[double,double]] removeOutliersIterative(
                libcpp_vector[libcpp_pair[double,double]] & pairs,
                double rsq_limit,
                double coverage_limit,
                bool use_chauvenet,
                libcpp_utf8_string outlier_detection_method
                ) except + nogil

        @staticmethod
        libcpp_vector[libcpp_pair[double,double]] removeOutliersRANSAC(
                libcpp_vector[libcpp_pair[double,double]] & pairs,
                double rsq_limit,
                double coverage_limit,
                size_t max_iterations,
                double max_rt_threshold,
                size_t sampling_size
                ) except + nogil

        @staticmethod
        double chauvenet_probability(libcpp_vector[ double ] residuals, int pos) except + nogil

        @staticmethod
        bool chauvenet(libcpp_vector[ double ] residuals, int pos) except + nogil

        @staticmethod
        bool computeBinnedCoverage(libcpp_pair[double,double] rtRange,
                                   libcpp_vector[libcpp_pair[double,double]] & pairs,
                                   int nrBins,
                                   int minPeptidesPerBin,
                                   int minBinsFilled) except + nogil
