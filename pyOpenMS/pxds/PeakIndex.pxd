from Types cimport *
from libcpp cimport bool
from Types cimport *

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
        # TEMPLATE # # NAMESPACE # FeatureMapType::value_type  getFeature(FeatureMapType & map_) nogil except +
        # TEMPLATE # # NAMESPACE # MSExperiment[Peak1D, ChromatogramPeak]Type::PeakType  getPeak(MSExperiment[Peak1D, ChromatogramPeak]Type & map_) nogil except +
        # TEMPLATE # # NAMESPACE # MSExperiment[Peak1D, ChromatogramPeak]Type::SpectrumType  getSpectrum(MSExperiment[Peak1D, ChromatogramPeak]Type & map_) nogil except +
        bool operator==(PeakIndex & rhs) nogil except +
        bool operator!=(PeakIndex & rhs) nogil except +

