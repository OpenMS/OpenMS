from Types cimport *
from libcpp cimport bool
from TransitionTSVFile cimport *
from TargetedExperiment cimport *
from LightTargetedExperiment cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>" namespace "OpenMS":

    cdef cppclass TransitionPQPFile:

        TransitionPQPFile() nogil except +
        TransitionPQPFile(TransitionPQPFile) nogil except + #wrap-ignore

        void convertTargetedExperimentToPQP(char * filename, TargetedExperiment & targeted_exp) nogil except +
        void convertPQPToTargetedExperiment(char * filename, TargetedExperiment & targeted_exp, bool legacy_traml_id) nogil except +
        void convertPQPToTargetedExperiment(char * filename, LightTargetedExperiment & targeted_exp, bool legacy_traml_id) nogil except +

        # inherited from TransitionTSVFile
        # due to issues with Cython and overloaded inheritance
        void convertTargetedExperimentToTSV(char * filename, TargetedExperiment& targeted_exp) nogil except +

        void convertTSVToTargetedExperiment(char * filename, FileType filetype, TargetedExperiment& targeted_exp) nogil except +
        void convertTSVToTargetedExperiment(char * filename, FileType filetype, LightTargetedExperiment& targeted_exp) nogil except +

        void validateTargetedExperiment(TargetedExperiment targeted_exp) nogil except +
