from Types cimport *
from libcpp.string cimport string as libcpp_string
from libcpp.vector cimport vector as libcpp_vector
from libcpp.set cimport set as libcpp_set
from ParamValue cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/Param.h>" namespace "OpenMS::Param":

    cdef cppclass ParamEntry:

        libcpp_string name
        libcpp_string description
        ParamValue value
        libcpp_set[libcpp_string] tags
        libcpp_vector[libcpp_string] valid_strings
        double  max_float
        double  min_float
        Int  max_int
        Int  min_int

        ParamEntry() nogil except + # TODO
        ParamEntry(ParamEntry &) nogil except +
        ParamEntry(libcpp_string n, ParamValue v, libcpp_string d, libcpp_vector[libcpp_string] t) nogil except +
        ParamEntry(libcpp_string n, ParamValue v, libcpp_string d) nogil except +

        # TODO: wrap bool isValid(libcpp_string & message) maybe as libcpp_pair[bool, libcpp_str] isValid() if possible
        bool operator==(ParamEntry) nogil except +

