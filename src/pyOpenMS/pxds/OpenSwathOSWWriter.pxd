from Types cimport *
from String cimport *
from LightTargetedExperiment cimport *
from FeatureMap cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWWriter.h>" namespace "OpenMS":
    
    cdef cppclass OpenSwathOSWWriter "OpenMS::OpenSwathOSWWriter":

        OpenSwathOSWWriter(OpenSwathOSWWriter) nogil except + #wrap-ignore
        OpenSwathOSWWriter(String output_filename, String input_filename, bool ms1_scores, bool sonar, bool uis_scores) nogil except +

        bool isActive() nogil except +
        void writeHeader() nogil except +
        String prepareLine(LightCompound & compound, LightTransition * tr, FeatureMap & output, String id_) nogil except +
        void writeLines(libcpp_vector[ String ] to_osw_output) nogil except +

