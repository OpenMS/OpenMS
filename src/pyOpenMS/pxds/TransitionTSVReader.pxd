from ProgressLogger cimport *
from libcpp.map cimport map as libcpp_map
from libcpp cimport bool
from Types cimport *
from FileTypes cimport *

from TargetedExperiment cimport *
from LightTargetedExperiment cimport LightTargetedExperiment

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/TransitionTSVReader.h>" namespace "OpenMS":

    cdef cppclass TransitionTSVReader(ProgressLogger):
        # wrap-inherits:
        #    ProgressLogger

        TransitionTSVReader()                       nogil except +
        TransitionTSVReader(TransitionTSVReader)    nogil except + # wrap-ignore

        void convertTargetedExperimentToTSV(char * filename, TargetedExperiment& targeted_exp) nogil except +
    
        void convertTSVToTargetedExperiment(char * filename, FileType filetype, TargetedExperiment& targeted_exp) nogil except +

        void convertTSVToTargetedExperiment(char * filename, FileType filetype, LightTargetedExperiment& targeted_exp) nogil except +
    
        void validateTargetedExperiment(TargetedExperiment targeted_exp) nogil except +

