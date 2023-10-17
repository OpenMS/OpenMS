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
        # wrap-doc:
        #  PeakMarker marks peaks that seem to fulfill some criterion

        PeakMarker() except + nogil 
        PeakMarker(PeakMarker &) except + nogil 
        # see child classes
        # void apply(libcpp_map[ double, bool ] & , MSSpectrum & ) except + nogil 
        String getProductName() except + nogil  # wrap-doc:Returns the product name

