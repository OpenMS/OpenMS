from Types cimport *
from libcpp cimport bool
from libcpp.map cimport map as libcpp_map
from PeakMarker cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/ComplementMarker.h>" namespace "OpenMS":
    
    cdef cppclass ComplementMarker(PeakMarker) :
        # wrap-inherits:
        #  PeakMarker
        ComplementMarker() nogil except +
        ComplementMarker(ComplementMarker) nogil except +
        void apply(libcpp_map[ double, bool ] & , MSSpectrum[Peak1D] & ) nogil except +
        PeakMarker * create() nogil except + # wrap-ignore
        # String getProductName() nogil except +

