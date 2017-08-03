from Types cimport *
from libcpp cimport bool
from TransitionTSVReader cimport *
from TargetedExperiment cimport *
from LightTargetedExperiment cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/TransitionPQPReader.h>" namespace "OpenMS":

    cdef cppclass TransitionPQPReader:

        TransitionPQPReader() nogil except +
        TransitionPQPReader(TransitionPQPReader) nogil except + #wrap-ignore

        void convertTargetedExperimentToPQP(char * filename, TargetedExperiment & targeted_exp) nogil except +
        void convertPQPToTargetedExperiment(char * filename, TargetedExperiment & targeted_exp, bool legacy_traml_id) nogil except +
        void convertPQPToTargetedExperiment(char * filename, LightTargetedExperiment & targeted_exp, bool legacy_traml_id) nogil except +

        # inherited from TransitionTSVReader
        # due to issues with Cython and overloaded inheritance
        void convertTargetedExperimentToTSV(char * filename, TargetedExperiment& targeted_exp) nogil except +

        void convertTSVToTargetedExperiment(char * filename, FileType filetype, TargetedExperiment& targeted_exp) nogil except +
        void convertTSVToTargetedExperiment(char * filename, FileType filetype, LightTargetedExperiment& targeted_exp) nogil except +

        void validateTargetedExperiment(TargetedExperiment targeted_exp) nogil except +
