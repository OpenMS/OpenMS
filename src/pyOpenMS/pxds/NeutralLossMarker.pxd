from Types cimport *
from libcpp cimport bool
from libcpp.map cimport map as libcpp_map
from PeakMarker cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/NeutralLossMarker.h>" namespace "OpenMS":
    
    cdef cppclass NeutralLossMarker(PeakMarker) :
        # wrap-inherits:
        #  PeakMarker
        NeutralLossMarker() nogil except +
        NeutralLossMarker(NeutralLossMarker) nogil except +
        void apply(libcpp_map[ double, bool ] & , MSSpectrum[Peak1D] & ) nogil except +
        PeakMarker * create() nogil except + # wrap-ignore
        # String getProductName() nogil except +

