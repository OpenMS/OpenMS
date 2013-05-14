from ProgressLogger cimport *
from libcpp.map cimport map as libcpp_map
from libcpp cimport bool
from Types cimport *

from TargetedExperiment cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/TransitionTSVReader.h>" namespace "OpenMS":

    cdef cppclass TransitionTSVReader(ProgressLogger):
        # wrap-inherits:
        #    ProgressLogger

        TransitionTSVReader()                       nogil except +
        TransitionTSVReader(TransitionTSVReader)    nogil except + # wrap-ignore

        void convertTargetedExperimentToTSV(char * filename, TargetedExperiment& targeted_exp)
    
        void convertTSVToTargetedExperiment(char * filename, TargetedExperiment& targeted_exp)
    
        void validateTargetedExperiment(TargetedExperiment targeted_exp)

