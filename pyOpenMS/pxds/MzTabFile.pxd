from MzTab cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/MzTabFile.h>" namespace "OpenMS":

    cdef cppclass MzTabFile:

        MzTabFile() nogil except +

        # void store(String filename, MzTab & mz_tab) nogil except +

