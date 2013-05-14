from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from String cimport *
from StringList cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/IntList.h>" namespace "OpenMS":

    cdef cppclass IntList:
        # wrap-ignore
        IntList()          nogil except +
        IntList(IntList)          nogil except +
        IntList(libcpp_vector[int])          nogil except +
        Size size()          nogil except +
        Int at(int) nogil except +

        bool contains(Int s) nogil except +
        IntList create(String & list_) nogil except +
        IntList create(StringList & list_) nogil except +
