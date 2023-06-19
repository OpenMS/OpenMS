from Types cimport *
from MzTabM cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/MzTabMFile.h>" namespace "OpenMS":

    cdef cppclass MzTabMFile:

        MzTabMFile() nogil except +
        MzTabMFile(MzTabMFile &) nogil except + # compiler

        void store(String filename, MzTabM & mztab_m) nogil except + # wrap-doc:Store MzTabM file
