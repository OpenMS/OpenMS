from Types cimport *
from libcpp cimport bool
from PeakMarker cimport *
from IsotopeDistribution cimport *
from MSSpectrum cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/IsotopeMarker.h>" namespace "OpenMS":
    
    cdef cppclass IsotopeMarker(PeakMarker) :
        # wrap-inherits:
        #  PeakMarker
        IsotopeMarker() except + nogil  # wrap-doc:IsotopeMarker marks peak pairs which could represent an ion and its isotope
        IsotopeMarker(IsotopeMarker &) except + nogil 

        void apply(libcpp_map[ double, bool ] & , MSSpectrum & ) except + nogil 
        PeakMarker * create() except + nogil  # wrap-ignore
        # TODO
        #String getProductName() except + nogil 

