from MSExperiment cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/MSPMetaboFile.h>" namespace "OpenMS":

    cdef cppclass MSPMetaboFile "OpenMS::MSPMetaboFile":

        MSPMetaboFile() nogil except +
        MSPMetaboFile(MSPMetaboFile) nogil except +
        MSPMetaboFile(const String& filename, MSExperiment& library) nogil except +

        void load(const String& filename, MSExperiment& library) nogil except +
