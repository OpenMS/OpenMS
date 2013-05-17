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
        IsotopeMarker() nogil except +
        IsotopeMarker(IsotopeMarker) nogil except +
        void apply(libcpp_map[ double, bool ] & , MSSpectrum[Peak1D] & ) nogil except +
        # POINTER # PeakMarker * create() nogil except +
        # String getProductName() nogil except +

