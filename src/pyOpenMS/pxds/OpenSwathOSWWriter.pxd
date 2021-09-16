from Types cimport *
from String cimport *
from LightTargetedExperiment cimport *
from FeatureMap cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWWriter.h>" namespace "OpenMS":
    
    cdef cppclass OpenSwathOSWWriter "OpenMS::OpenSwathOSWWriter":

        OpenSwathOSWWriter(String output_filename, UInt64 run_id, String input_filename, bool ms1_scores, bool sonar, bool uis_scores) nogil except +
        OpenSwathOSWWriter(OpenSwathOSWWriter &) nogil except + # compiler

        bool isActive() nogil except +
        void writeHeader() nogil except +
        String prepareLine(LightCompound & compound, LightTransition * tr, FeatureMap & output, String id_) nogil except +
        void writeLines(libcpp_vector[ String ] to_osw_output) nogil except +
