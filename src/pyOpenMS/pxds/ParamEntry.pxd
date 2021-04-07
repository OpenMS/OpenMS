from libcpp.string cimport string as libcpp_string
from libcpp.vector cimport vector as libcpp_vector
from libcpp.set cimport set as libcpp_set
from libcpp cimport bool
from ParamValue cimport *
from String cimport *
from StringList cimport *
from ParamNode cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/Param.h>" namespace "OpenMS::Param":

    cdef cppclass ParamEntry:

        String name
        String description
        ParamValue value
        libcpp_set[String] tags
        libcpp_vector[String] valid_strings
        double  max_float
        double  min_float
        Int  max_int
        Int  min_int

        ParamEntry() nogil except +
        ParamEntry(ParamEntry) nogil except +
        ParamEntry(String n, ParamValue v, String d, StringList t) nogil except +
        ParamEntry(String n, ParamValue v, String d) nogil except +

        bool isValid(String &message) nogil except +
        bool operator==(ParamEntry) nogil except +

