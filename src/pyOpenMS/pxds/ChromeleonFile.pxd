from MSExperiment cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/ChromeleonFile.h>" namespace "OpenMS":

    cdef cppclass ChromeleonFile:

        ChromeleonFile() nogil except +
        ChromeleonFile(ChromeleonFile) nogil except +

        void load(const String& filename, MSExperiment& experiment) nogil except +
