from Param cimport *
from libcpp.string cimport string as libcpp_utf8_string
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/FORMAT/ParamCTDFile.h>" namespace "OpenMS":

    cdef cppclass ToolInfo:

        ToolInfo(ToolInfo) nogil except +

    cdef cppclass ParamCTDFile:

        ParamCTDFile() nogil except +
        void store(libcpp_utf8_string filename, Param param, ToolInfo tool_info) nogil except +
