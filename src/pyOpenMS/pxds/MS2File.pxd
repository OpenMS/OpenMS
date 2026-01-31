from Types cimport *
from String cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/FORMAT/MS2File.h>" namespace "OpenMS":

    cdef cppclass MS2File:
        MS2File() except + nogil 
        MS2File(MS2File &) except + nogil  # compiler
        void load(const String & filename, MSExperiment & exp) except + nogil 

