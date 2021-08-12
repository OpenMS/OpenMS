from Types cimport *
from libcpp cimport bool
from libcpp.map cimport map as libcpp_map
from PeakMarker cimport *
from IsotopeDistribution cimport *
from MSSpectrum cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/IsotopeMarker.h>" namespace "OpenMS":
    
    cdef cppclass IsotopeMarker(PeakMarker) :
        # wrap-inherits:
        #  PeakMarker
        IsotopeMarker() nogil except + # wrap-doc:IsotopeMarker marks peak pairs which could represent an ion and its isotope
        IsotopeMarker(IsotopeMarker &) nogil except +

        void apply(libcpp_map[ double, bool ] & , MSSpectrum & ) nogil except +
        PeakMarker * create() nogil except + # wrap-ignore
        # TODO
        #String getProductName() nogil except +

