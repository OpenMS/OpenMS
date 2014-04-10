from MzTab cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/MzTabFile.h>" namespace "OpenMS":

    # TODO missing functions
    cdef cppclass MzTabFile:

        MzTabFile() nogil except +

        void store(String filename, MzTab & mz_tab) nogil except +
        void load(String filename, MzTab & mz_tab) nogil except +

