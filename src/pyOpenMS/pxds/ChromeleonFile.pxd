from MSRun cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/ChromeleonFile.h>" namespace "OpenMS":

    cdef cppclass ChromeleonFile:

        ChromeleonFile() except + nogil  # wrap-doc:Load Chromeleon HPLC text file and save it into a `MSRun`.
        ChromeleonFile(ChromeleonFile &) except + nogil  # compiler

        void load(const String& filename, MSRun& experiment) except + nogil  # wrap-doc:Load the file's data and metadata, and save it into a `MSRun`
