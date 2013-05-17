from libcpp.string cimport string as libcpp_string
from libcpp.vector cimport vector as libcpp_vector
from libcpp.set cimport set as libcpp_set
from libcpp cimport bool
from DataValue cimport *
from String cimport *
from StringList cimport *
from ParamNode cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/Param.h>" namespace "OpenMS::Param":

    cdef cppclass ParamEntry:

        String name
        String description
        DataValue value
        libcpp_set[String] tags
        libcpp_vector[String] valid_strings
        DoubleReal  max_float
        DoubleReal  min_float
        Int  max_int
        Int  min_int

        ParamEntry() nogil except +
        ParamEntry(ParamEntry) nogil except +
        ParamEntry(String n, DataValue v, String d, StringList t) nogil except +
        ParamEntry(String n, DataValue v, String d) nogil except +

        bool isValid(String &message) nogil except +
        bool operator==(ParamEntry) nogil except +

