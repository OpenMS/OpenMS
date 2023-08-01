from Types cimport *
from libcpp cimport bool
from libcpp.map cimport map as libcpp_map
from PeakMarker cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/ComplementMarker.h>" namespace "OpenMS":
    
    cdef cppclass ComplementMarker(PeakMarker) :
        # wrap-inherits:
        #  PeakMarker

        ComplementMarker() except + nogil  # wrap-doc:ComplementMarker marks peak pairs which could represent y - b ion pairs
        ComplementMarker(ComplementMarker &) except + nogil 

        void apply(libcpp_map[ double, bool ] & , MSSpectrum & ) except + nogil 
        PeakMarker * create() except + nogil  # wrap-ignore
        # String getProductName() except + nogil 
