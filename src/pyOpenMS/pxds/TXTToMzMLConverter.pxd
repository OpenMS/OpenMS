from MSExperiment cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/TXTToMzMLConverter.h>" namespace "OpenMS":

    cdef cppclass TXTToMzMLConverter:

        TXTToMzMLConverter() nogil except +
        TXTToMzMLConverter(TXTToMzMLConverter) nogil except +
        
        MSExperiment loadInputFile(const String& filename) nogil except +
        void storeMzMLFile(const String& filename, const MSExperiment& experiment) nogil except +
