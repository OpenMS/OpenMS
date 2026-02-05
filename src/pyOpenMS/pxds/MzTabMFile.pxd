from Types cimport *
from MzTabM cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/MzTabMFile.h>" namespace "OpenMS":

    cdef cppclass MzTabMFile:
        # wrap-doc:
        #  File adapter for MzTab-M files

        MzTabMFile() except + nogil 
        MzTabMFile(MzTabMFile &) except + nogil  # compiler

        void store(String filename, MzTabM & mztab_m) except + nogil  # wrap-doc:Store MzTabM file
