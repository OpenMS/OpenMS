from libcpp.vector cimport vector as libcpp_vector
from Types cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/DoubleList.h>" namespace "OpenMS":

    cdef cppclass DoubleList:
        # wrap-ignore
        DoubleList() nogil except +
        DoubleList(DoubleList) nogil except +
        DoubleList(libcpp_vector[double]) nogil except +
        Size size() nogil except +
        DoubleReal at(int) nogil except +
