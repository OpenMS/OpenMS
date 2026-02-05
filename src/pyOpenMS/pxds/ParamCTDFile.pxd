from Param cimport *
from libcpp.string cimport string as libcpp_utf8_string
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/FORMAT/ParamCTDFile.h>" namespace "OpenMS":

    cdef cppclass ToolInfo:

        ToolInfo(ToolInfo) except + nogil 

    cdef cppclass ParamCTDFile:
        # wrap-doc:
        #  A struct to pass information about the tool as one parameter */
        #  struct ToolInfo { std::string version_; std::string name_; std::string
        #  docurl_; std::string category_; std::string description_;
        #  std::vector<std::string> citations_; };

        ParamCTDFile() except + nogil 
        void store(libcpp_utf8_string filename, Param param, ToolInfo tool_info) except + nogil 
