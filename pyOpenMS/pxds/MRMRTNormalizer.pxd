from libcpp.vector cimport vector as libcpp_vector
from libcpp.pair cimport pair as libcpp_pair
from libcpp cimport bool

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MRMRTNormalizer.h>" namespace "OpenMS":

    cdef cppclass MRMRTNormalizer:
        pass

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MRMRTNormalizer.h>" namespace "OpenMS::MRMRTNormalizer":


    # @brief This function computes a candidate outlier peptide by iteratively
    #  leaving one peptide out to find the one which results in the maximum R^2
    #  of a first order linear regression of the remaining ones. The datapoints
    #  are submitted as two vectors of doubles (x- and y-coordinates).
    #
    # @return The position of the candidate outlier peptide as supplied by the
    #  vector is returned.
    #
    # @exception Exception::UnableToFit is thrown if fitting cannot be performed
    int outlier_candidate(libcpp_vector[double] x, libcpp_vector[double] y
                                ) nogil except + # wrap-attach:MRMRTNormalizer

    # @brief This function removes potential outliers from a set of paired points.
    #  Two thresholds need to be defined, first a lower R^2 limit to accept the
    #  regression for the RT normalization and second, the lower limit of peptide
    #  coverage. The algorithms then selects candidate outlier peptides and applies
    #  the Chauvenet's criterion on the assumption that the residuals are normal
    #  distributed to determine whether the peptides can be removed. This is done
    #  iteratively until both limits are reached.
    #
    # @return A vector of pairs is returned if the R^2 limit was reached without
    #  reaching the coverage limit. If the limits are reached, an exception is
    #  thrown.
    #
    # @exception Exception::UnableToFit is thrown if fitting cannot be performed
    libcpp_vector[libcpp_pair[double,double]] rm_outliers(
            libcpp_vector[libcpp_pair[double,double]] & pairs,
            double rsq_limit,
            double coverage_limit
            ) nogil except + # wrap-attach:MRMRTNormalizer

    # @brief This function computes Chauvenet's criterion probability for a vector
    #  and a value whose position is submitted.
    #
    # @return Chauvenet's criterion probability
    double chauvenet_probability(libcpp_vector[double] residuals, int pos
            ) nogil except + # wrap-attach:MRMRTNormalizer


    # @brief This function computes Chauvenet's criterion for a vector and a value
    #  whose position is submitted.
    #
    # @return TRUE, if Chauvenet's criterion is fullfilled and the outlier can be removed.
    bool chauvenet(libcpp_vector[double] residuals, int pos
                         ) nogil except + # wrap-attach:MRMRTNormalizer

