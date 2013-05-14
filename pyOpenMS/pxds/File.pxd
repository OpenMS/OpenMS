from String cimport *
from Param cimport *
from StringList cimport *
from libcpp cimport bool

cdef extern from "<OpenMS/SYSTEM/File.h>" namespace "OpenMS":

    cdef cppclass File:
        pass


cdef extern from "<OpenMS/SYSTEM/File.h>" namespace "OpenMS::File":

    String getExecutablePath()  # wrap-attach:File
    String getOpenMSDataPath()  # wrap-attach:File
    Param  getSystemParameters()  # wrap-attach:File

