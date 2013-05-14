from libcpp.vector cimport vector as libcpp_vector
from Types cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/IntList.h>" namespace "OpenMS":

    cdef cppclass IntList:
        # wrap-ignore
        IntList()          nogil except +
        IntList(IntList)          nogil except +
        IntList(libcpp_vector[int])          nogil except +
        size_t size()          nogil except +
        Int at(int) nogil except +
