from MSExperiment cimport *
from String cimport *
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MSPMetaboFile.h>" namespace "OpenMS":

    cdef cppclass MSPMetaboFile "OpenMS::MSPMetaboFile":

        MSPMetaboFile() nogil except +
        MSPMetaboFile(MSPMetaboFile) nogil except +
        MSPMetaboFile(const String& filename, MSExperiment& library) nogil except +

        void load(const String& filename, MSExperiment& library) nogil except +
