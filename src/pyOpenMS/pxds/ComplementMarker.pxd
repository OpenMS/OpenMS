from Types cimport *
from libcpp cimport bool
from libcpp.map cimport map as libcpp_map
from PeakMarker cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/ComplementMarker.h>" namespace "OpenMS":
    
    cdef cppclass ComplementMarker(PeakMarker) :
        # wrap-inherits:
        #  PeakMarker

        ComplementMarker() nogil except + # wrap-doc:ComplementMarker marks peak pairs which could represent y - b ion pairs
        ComplementMarker(ComplementMarker &) nogil except +

        void apply(libcpp_map[ double, bool ] & , MSSpectrum & ) nogil except +
        PeakMarker * create() nogil except + # wrap-ignore
        # String getProductName() nogil except +
