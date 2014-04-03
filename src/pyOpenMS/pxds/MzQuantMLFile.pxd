from MSQuantifications cimport *
from String cimport *
from StringList cimport *

cdef extern from "<OpenMS/FORMAT/MzQuantMLFile.h>" namespace "OpenMS":

    cdef cppclass MzQuantMLFile:
        MzQuantMLFile() nogil except +

        void load(String, MSQuantifications &) nogil except+
        void store(String, MSQuantifications &) nogil except+

        bool isSemanticallyValid(String filename,
                                 StringList & errors,
                                 StringList & warnings) nogil except +


