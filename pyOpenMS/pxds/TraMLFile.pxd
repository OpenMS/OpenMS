from String cimport *
from StringList cimport *
from TargetedExperiment cimport *

cdef extern from "<OpenMS/FORMAT/TraMLFile.h>" namespace "OpenMS":

    cdef cppclass TraMLFile:

        TraMLFile() nogil except +

        void load(String filename,
                  TargetedExperiment & id) nogil except +

        void store(String filename,
                  TargetedExperiment & id) nogil except +

        bool isSemanticallyValid(String filename, StringList & errors,
                                 StringList & warnings) nogil except +
