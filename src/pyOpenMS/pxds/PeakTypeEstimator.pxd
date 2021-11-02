from Peak1D cimport *
from libcpp.vector cimport vector as libcpp_vector

# this class has addons, see the ./addons folder

cdef extern from "<OpenMS/FORMAT/PeakTypeEstimator.h>" namespace "OpenMS":

    cdef cppclass PeakTypeEstimator:
    # wrap-doc:
    #   Estimates if the data of a spectrum is raw data or peak data
    
        PeakTypeEstimator() nogil except + 
        PeakTypeEstimator(PeakTypeEstimator &) nogil except +

        # wrpped in ../addons/PeakTypeEstimator.pyx:
        int estimateType(libcpp_vector[Peak1D].iterator,
                         libcpp_vector[Peak1D].iterator) nogil except + # wrap-ignore
