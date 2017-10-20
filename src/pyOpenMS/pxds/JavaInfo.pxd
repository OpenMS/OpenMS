from libcpp cimport bool
from Types cimport *
from String cimport *

cdef extern from "<OpenMS/SYSTEM/JavaInfo.h>" namespace "OpenMS":

    cdef cppclass JavaInfo:

        JavaInfo() nogil except +
        JavaInfo(JavaInfo) nogil except + # wrap-ignore

        bool canRun(String java_executable) nogil except +

