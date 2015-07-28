from libcpp.vector cimport vector as libcpp_vector

from ChromatogramPeak cimport *
from DefaultParamHandler cimport *
from MSExperiment cimport *
from Peak1D cimport *
from ProgressLogger cimport *
from TransformationDescription cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmSpectrumAlignment.h>" namespace "OpenMS":
    
    cdef cppclass MapAlignmentAlgorithmSpectrumAlignment(DefaultParamHandler, ProgressLogger):
        # wrap-inherits:
        #    DefaultParamHandler
        #    ProgressLogger

        MapAlignmentAlgorithmSpectrumAlignment() nogil except +
 
        void align(libcpp_vector[MSExperiment[Peak1D, ChromatogramPeak]]&, libcpp_vector[TransformationDescription]&) nogil except +

