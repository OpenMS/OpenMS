from Types cimport *
from libcpp cimport bool
from TransitionTSVFile cimport *
from TargetedExperiment cimport *
from LightTargetedExperiment cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>" namespace "OpenMS":

    cdef cppclass TransitionPQPFile:

        TransitionPQPFile() nogil except +
        TransitionPQPFile(TransitionPQPFile &) nogil except + # compiler

        void convertTargetedExperimentToPQP(char * filename, TargetedExperiment & targeted_exp) nogil except +
        # wrap-doc:
                #   Write out a targeted experiment (TraML structure) into a PQP file
                #   -----
                #   :param filename: The output file
                #   :param targeted_exp: The targeted experiment

        void convertPQPToTargetedExperiment(char * filename, TargetedExperiment & targeted_exp, bool legacy_traml_id) nogil except +
        # wrap-doc:
                #   Read in a PQP file and construct a targeted experiment (TraML structure)
                #   -----
                #   :param filename: The input file
                #   :param targeted_exp: The output targeted experiment
                #   :param legacy_traml_id: Should legacy TraML IDs be used (boolean)?

        void convertPQPToTargetedExperiment(char * filename, LightTargetedExperiment & targeted_exp, bool legacy_traml_id) nogil except +
        # wrap-doc:
                #   Read in a PQP file and construct a targeted experiment (Light transition structure)
                #   -----
                #   :param filename: The input file
                #   :param targeted_exp: The output targeted experiment
                #   :param legacy_traml_id: Should legacy TraML IDs be used (boolean)?

        # inherited from TransitionTSVFile
        # due to issues with Cython and overloaded inheritance
        void convertTargetedExperimentToTSV(char * filename, TargetedExperiment& targeted_exp) nogil except +

        void convertTSVToTargetedExperiment(char * filename, FileType filetype, TargetedExperiment& targeted_exp) nogil except +
        void convertTSVToTargetedExperiment(char * filename, FileType filetype, LightTargetedExperiment& targeted_exp) nogil except +

        void validateTargetedExperiment(TargetedExperiment targeted_exp) nogil except +
