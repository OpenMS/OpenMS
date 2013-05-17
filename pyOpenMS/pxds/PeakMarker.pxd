from Types cimport *
from libcpp cimport bool
from libcpp.map cimport map as libcpp_map
from DefaultParamHandler cimport *
from MSSpectrum cimport *
from Peak1D cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/PeakMarker.h>" namespace "OpenMS":
    
    cdef cppclass PeakMarker(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        PeakMarker() nogil except +
        PeakMarker(PeakMarker) nogil except +
        # TEMPLATE # void apply(libcpp_map[ double, bool ] & , MSSpectrum[Peak1D] & ) nogil except +
        String getProductName() nogil except +

