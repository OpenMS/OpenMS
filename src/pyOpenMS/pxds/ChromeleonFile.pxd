from MSExperiment cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/ChromeleonFile.h>" namespace "OpenMS":

    cdef cppclass ChromeleonFile:

        ChromeleonFile() nogil except + # wrap-doc:Load Chromeleon HPLC text file and save it into a `MSExperiment`.
        ChromeleonFile(ChromeleonFile &) nogil except + # compiler

        void load(const String& filename, MSExperiment& experiment) nogil except + # wrap-doc:Load the file's data and metadata, and save it into a `MSExperiment`
