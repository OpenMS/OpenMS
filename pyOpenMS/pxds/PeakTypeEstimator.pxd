from Peak1D cimport *
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/FORMAT/PeakTypeEstimator.h>" namespace "OpenMS":

    cdef cppclass PeakTypeEstimator:
        PeakTypeEstimator()
        int estimateType(libcpp_vector[Peak1D].iterator, libcpp_vector[Peak1D].iterator) # wrap-ignore
