from Types cimport *
from libcpp cimport bool
from TransitionTSVReader cimport *
from TargetedExperiment cimport *
from LightTargetedExperiment cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/TransitionPQPReader.h>" namespace "OpenMS":
    
    cdef cppclass TransitionPQPReader(TransitionTSVReader) :
        # wrap-inherits:
        #  TransitionTSVReader
        TransitionPQPReader() nogil except +
        TransitionPQPReader(TransitionPQPReader) nogil except + #wrap-ignore
        void convertTargetedExperimentToPQP(char * filename, TargetedExperiment & targeted_exp) nogil except +
        void convertPQPToTargetedExperiment(char * filename, TargetedExperiment & targeted_exp, bool legacy_traml_id) nogil except +
        void convertPQPToTargetedExperiment(char * filename, LightTargetedExperiment & targeted_exp, bool legacy_traml_id) nogil except +

