from MSExperiment cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/ChromeleonFile.h>" namespace "OpenMS":

    cdef cppclass ChromeleonFile:

        ChromeleonFile() except + nogil  # wrap-doc:Load Chromeleon HPLC text file and save it into a `MSExperiment`.
        ChromeleonFile(ChromeleonFile &) except + nogil  # compiler

        void load(const String& filename, MSExperiment& experiment) except + nogil  # wrap-doc:Load the file's data and metadata, and save it into a `MSExperiment`
