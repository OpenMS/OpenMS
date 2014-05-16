from Types cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MRMRTNormalizer.h>" namespace "OpenMS":

    cdef cppclass MRMRTNormalizer:
        pass

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MRMRTNormalizer.h>" namespace "OpenMS::MRMRTNormalizer":

    libcpp_vector[libcpp_pair[double,double]] removeOutliersIterative(
            libcpp_vector[libcpp_pair[double,double]] & pairs,
            double rsq_limit,
            double coverage_limit, 
            bool use_chauvenet,
            libcpp_string outlier_detection_method 
            ) nogil except + # wrap-attach:MRMRTNormalizer

    libcpp_vector[libcpp_pair[double,double]] removeOutliersRANSAC(
            libcpp_vector[libcpp_pair[double,double]] & pairs,
            double rsq_limit,
            double coverage_limit, 
            size_t max_iterations,
            double max_rt_threshold
            ) nogil except + # wrap-attach:MRMRTNormalizer

    libcpp_vector[libcpp_pair[double,double]] ransac(
            libcpp_vector[libcpp_pair[double,double]] & pairs,
            size_t n, size_t k, double t, size_t d, bool test
            ) nogil except + # wrap-attach:MRMRTNormalizer

