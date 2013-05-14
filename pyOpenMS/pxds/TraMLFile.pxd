from String cimport *
from TargetedExperiment cimport *

cdef extern from "<OpenMS/FORMAT/TraMLFile.h>" namespace "OpenMS":

    cdef cppclass TraMLFile:

        TraMLFile() nogil except +

        void load(String filename,
                  TargetedExperiment & id)

        void store(String filename,
                  TargetedExperiment & id)
