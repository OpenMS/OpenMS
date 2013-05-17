from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from AASequence cimport *
from LocalLinearMap cimport *

cdef extern from "<OpenMS/ANALYSIS/PIP/PeakIntensityPredictor.h>" namespace "OpenMS":
    
    cdef cppclass PeakIntensityPredictor "OpenMS::PeakIntensityPredictor":
        PeakIntensityPredictor() nogil except +
        PeakIntensityPredictor(PeakIntensityPredictor) nogil except + #wrap-ignore
        DoubleReal predict(AASequence & sequence) nogil except +
        DoubleReal predict(AASequence & sequence, libcpp_vector[ double ] & add_info) nogil except +
        libcpp_vector[ double ] predict(libcpp_vector[ AASequence ] & sequences) nogil except +
        libcpp_vector[ double ] predict(libcpp_vector[ AASequence ] & sequences, libcpp_vector[ libcpp_vector[ double ] ] & add_info) nogil except +

