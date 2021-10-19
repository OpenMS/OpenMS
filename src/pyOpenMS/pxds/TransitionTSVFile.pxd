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

        TransitionTSVFile() nogil except +
        TransitionTSVFile(TransitionTSVFile &) nogil except + # compiler

        void convertTargetedExperimentToTSV(char * filename, TargetedExperiment& targeted_exp) nogil except + # wrap-doc:Write out a targeted experiment (TraML structure) into a tsv file
    
        void convertTSVToTargetedExperiment(char * filename, FileType filetype, TargetedExperiment& targeted_exp) nogil except + # wrap-doc:Read in a tsv/mrm file and construct a targeted experiment (TraML structure)
        void convertTSVToTargetedExperiment(char * filename, FileType filetype, LightTargetedExperiment& targeted_exp) nogil except + # wrap-doc:Read in a tsv file and construct a targeted experiment (Light transition structure)
    
        void validateTargetedExperiment(TargetedExperiment targeted_exp) nogil except + # wrap-doc:Validate a TargetedExperiment (check that all ids are unique)
