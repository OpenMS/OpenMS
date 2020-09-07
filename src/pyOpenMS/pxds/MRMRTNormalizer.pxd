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
            double max_rt_threshold,
            size_t sampling_size
            ) nogil except + # wrap-attach:MRMRTNormalizer

    double chauvenet_probability(libcpp_vector[ double ] residuals, int pos) nogil except + # wrap-attach:MRMRTNormalizer
    bool chauvenet(libcpp_vector[ double ] residuals, int pos) nogil except + # wrap-attach:MRMRTNormalizer

    bool computeBinnedCoverage(libcpp_pair[double,double] rtRange, 
                               libcpp_vector[libcpp_pair[double,double]] & pairs,
                               int nrBins, 
                               int minPeptidesPerBin,
                               int minBinsFilled) nogil except + # wrap-attach:MRMRTNormalizer

