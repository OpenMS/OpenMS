from Types cimport *
from FeatureMap cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/KERNEL/PeakIndex.h>" namespace "OpenMS":
    
    cdef cppclass PeakIndex "OpenMS::PeakIndex":
        PeakIndex() nogil except +
        PeakIndex(PeakIndex) nogil except + #wrap-ignore
        Size peak
        Size spectrum
        PeakIndex(Size peak) nogil except +
        PeakIndex(Size spectrum, Size peak) nogil except +
        bool isValid() nogil except +
        void clear() nogil except +
        Feature getFeature(FeatureMap & map_) nogil except +
        Peak1D getPeak(MSExperiment & map_) nogil except +
        MSSpectrum getSpectrum(MSExperiment & map_) nogil except +
        bool operator==(PeakIndex & rhs) nogil except +
        bool operator!=(PeakIndex & rhs) nogil except +

