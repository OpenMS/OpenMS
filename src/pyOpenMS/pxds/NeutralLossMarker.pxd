from Types cimport *
from libcpp cimport bool
from libcpp.map cimport map as libcpp_map
from PeakMarker cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/NeutralLossMarker.h>" namespace "OpenMS":
    
    cdef cppclass NeutralLossMarker(PeakMarker) :
        # wrap-inherits:
        #  PeakMarker
        # wrap-doc:
        #  NeutralLossMarker marks peak pairs which could represent an ion an its neutral loss (water, ammonia)
        
        NeutralLossMarker() nogil except +
        NeutralLossMarker(NeutralLossMarker &) nogil except +
        void apply(libcpp_map[ double, bool ] & , MSSpectrum & ) nogil except +
        PeakMarker * create() nogil except + # wrap-ignore
        # String getProductName() nogil except +

