from Param cimport *

cdef extern from "<OpenMS/FORMAT/ParamCTDFile.h>" namespace "OpenMS":

    cdef cppclass ParamCTDFile:

        ParamCTDFile() nogil except +
        void store(String, Param &) nogil except+
