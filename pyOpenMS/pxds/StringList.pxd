from libcpp.string cimport string as libcpp_string
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/DATASTRUCTURES/StringList.h>" namespace "OpenMS":

    cdef cppclass StringList:
        # wrap-ignore
        StringList() nogil except +
        StringList(StringList) nogil except +
        StringList(libcpp_vector[libcpp_string]) nogil except +
        size_t size() nogil except +
        libcpp_string at(int) nogil except +
