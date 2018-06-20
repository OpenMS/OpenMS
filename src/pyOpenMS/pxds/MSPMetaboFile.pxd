from DefaultParamHandler cimport *
from MSExperiment cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/MSPMetaboFile.h>" namespace "OpenMS":
    cdef cppclass MSPMetaboFile(DefaultParamHandler):
        # wrap-inherits:
        #  DefaultParamHandler

        MSPMetaboFile() nogil except +
        MSPMetaboFile(MSPMetaboFile) nogil except +
        MSPMetaboFile(const String& filename, MSExperiment& library) nogil except +

        void load(const String& filename, MSExperiment& library) nogil except +
