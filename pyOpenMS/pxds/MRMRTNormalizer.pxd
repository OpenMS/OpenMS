from libcpp.vector cimport vector as libcpp_vector
from libcpp.pair cimport pair as libcpp_pair
from libcpp cimport bool

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MRMRTNormalizer.h>" namespace "OpenMS":

    cdef cppclass MRMRTNormalizer:

        libcpp_vector[libcpp_pair[double,double]] rm_outliers(libcpp_vector[libcpp_pair[double,double]] & pairs, double rsq_limit, double coverage_limit)

