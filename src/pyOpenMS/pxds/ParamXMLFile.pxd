from Param cimport *

cdef extern from "<OpenMS/FORMAT/ParamXMLFile.h>" namespace "OpenMS":

    cdef cppclass ParamXMLFile:

        ParamXMLFile() nogil except +
        ParamXMLFile(ParamXMLFile &) nogil except + # compiler

        void load(String, Param &) nogil except+
        void store(String, Param &) nogil except+
