from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from MapAlignmentAlgorithm cimport *
from PeakSpectrumCompareFunctor cimport *
from String cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmSpectrumAlignment.h>" namespace "OpenMS":
    
    cdef cppclass MapAlignmentAlgorithmSpectrumAlignment(MapAlignmentAlgorithm) :
        # wrap-inherits:
        #  MapAlignmentAlgorithm
        MapAlignmentAlgorithmSpectrumAlignment() nogil except +
        MapAlignmentAlgorithmSpectrumAlignment(MapAlignmentAlgorithmSpectrumAlignment) nogil except + #wrap-ignore
        void alignMSExperiment[Peak1D, ChromatogramPeak]s(libcpp_vector[ MSExperiment[Peak1D, ChromatogramPeak] ] & , libcpp_vector[ TransformationDescription ] & ) nogil except +
        # POINTER # MapAlignmentAlgorithm * create() nogil except +
        String getProductName() nogil except +

