from ProgressLogger cimport *
from libcpp.map cimport map as libcpp_map
from libcpp cimport bool
from Types cimport *
from FileTypes cimport *

from TargetedExperiment cimport *
from LightTargetedExperiment cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>" namespace "OpenMS":

    cdef cppclass TransitionTSVFile(ProgressLogger):
        # wrap-inherits:
        #    ProgressLogger

        TransitionTSVFile()                       nogil except +
        TransitionTSVFile(TransitionTSVFile)    nogil except + # wrap-ignore

        void convertTargetedExperimentToTSV(char * filename, TargetedExperiment& targeted_exp) nogil except +
    
        void convertTSVToTargetedExperiment(char * filename, FileType filetype, TargetedExperiment& targeted_exp) nogil except +
        void convertTSVToTargetedExperiment(char * filename, FileType filetype, LightTargetedExperiment& targeted_exp) nogil except +
    
        void validateTargetedExperiment(TargetedExperiment targeted_exp) nogil except +

