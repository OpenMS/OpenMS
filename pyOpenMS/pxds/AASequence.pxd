from libcpp cimport bool
from Types cimport *
from String cimport *

cdef extern from "<OpenMS/CHEMISTRY/AASequence.h>" namespace "OpenMS":

    cdef cppclass AASequence:

        AASequence() nogil except +
        AASequence(AASequence) nogil except + # wrap-ignore

        AASequence(char *) nogil except +

        String toString()  nogil except +
        String toUnmodifiedString()  nogil except +

        void setModification(Size index, String modification) nogil except +
        void setNTerminalModification(String modification) nogil except +
        void setCTerminalModification(String modification) nogil except +
        void setStringSequence(String modification) nogil except +

        String getNTerminalModification() nogil except +
        String getCTerminalModification() nogil except +

        AASequence operator+(AASequence)    nogil except +
        AASequence iadd(AASequence)   nogil except + # wrap-as:operator+=





