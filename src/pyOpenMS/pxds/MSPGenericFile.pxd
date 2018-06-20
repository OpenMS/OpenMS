from DefaultParamHandler cimport *
from MSExperiment cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/MSPGenericFile.h>" namespace "OpenMS":
    cdef cppclass MSPGenericFile(DefaultParamHandler):
        # wrap-inherits:
        #  DefaultParamHandler

        MSPGenericFile() nogil except +
        MSPGenericFile(MSPGenericFile) nogil except +
        MSPGenericFile(const String& filename, MSExperiment& library) nogil except +

        void load(const String& filename, MSExperiment& library) nogil except +
