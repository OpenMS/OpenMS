from Types cimport *
from MSSpectrum cimport *
from Peak1D cimport *

cdef extern from "<OpenMS/FILTERING/DATAREDUCTION/SplineSpectrum.h>" namespace "OpenMS":
    
    cdef cppclass SplineSpectrum "OpenMS::SplineSpectrum":

        SplineSpectrum(libcpp_vector[double] mz, libcpp_vector[double] intensity) nogil except +
        SplineSpectrum(libcpp_vector[double] mz, libcpp_vector[double] intensity, double scaling) nogil except +

        SplineSpectrum(MSSpectrum[Peak1D] raw_spectrum) nogil except +
        SplineSpectrum(MSSpectrum[Peak1D] raw_spectrum, double scaling) nogil except +

        double getMzMin() nogil except +

        double getMzMax() nogil except +

        int getSplineCount() nogil except +


