from Types cimport *
from String cimport *
from ProgressLogger cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/FORMAT/MS2File.h>" namespace "OpenMS":
    
    cdef cppclass MS2File(ProgressLogger) :
        # wrap-inherits:
        #  ProgressLogger
        MS2File() nogil except +
        MS2File(MS2File) nogil except + #wrap-ignore
        void load(const String & filename, MSExperiment & exp) nogil except +

