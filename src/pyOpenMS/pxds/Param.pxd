from libcpp.string cimport string as libcpp_utf8_string
from libcpp.string cimport string as libcpp_string
from Types cimport *
from ParamValue cimport *
from ParamNode cimport *
from ParamEntry cimport *


# this class has addons, see the ./addons folder

cdef extern from "<OpenMS/DATASTRUCTURES/Param.h>" namespace "OpenMS":

    # pythonic helper functions in ../addons/Param.pyx !!!!

    cdef cppclass Param:

         # COMMENT: Helper functions for Python
         initPluginParam(name, version) # wrap-ignore
         asDict() # wrap-ignore
         keys() # wrap-ignore
         items() # wrap-ignore
         values() # wrap-ignore
         update(dict) # wrap-ignore
         get(bytes key, default=None) # wrap-ignore
         __getitem__(bytes key) # wrap-ignore
         __setitem__(bytes key, value) # wrap-ignore

         Param() except + nogil 
         Param(Param &) except + nogil 
         bool operator==(Param) except + nogil 

         void setValue(libcpp_utf8_string key, ParamValue val, libcpp_utf8_string desc, libcpp_vector[libcpp_utf8_string] tags) except + nogil 
         void setValue(libcpp_utf8_string key, ParamValue val, libcpp_utf8_string desc) except + nogil 
         void setValue(libcpp_utf8_string key, ParamValue val) except + nogil 
         ParamValue getValue(libcpp_utf8_string key) except + nogil 
         ValueType getValueType(libcpp_utf8_string key) except + nogil 

         ParamEntry getEntry(libcpp_utf8_string) except + nogil 
         bool exists(libcpp_utf8_string key) except + nogil 

         void addTag(libcpp_utf8_string key, libcpp_utf8_string tag) except + nogil 
         void addTags(libcpp_utf8_string key, libcpp_vector[libcpp_utf8_string] tags) except + nogil 
         int hasTag(libcpp_utf8_string key, libcpp_utf8_string tag) except + nogil 
         libcpp_vector[libcpp_string] getTags(libcpp_utf8_string key) except + nogil 
         void clearTags(libcpp_utf8_string key) except + nogil 

         libcpp_utf8_output_string getDescription(libcpp_utf8_string key) except + nogil 
         void setSectionDescription(libcpp_utf8_string key, libcpp_utf8_string desc) except + nogil 
         libcpp_utf8_output_string getSectionDescription(libcpp_utf8_string key) except + nogil 

         void addSection(libcpp_utf8_string key, libcpp_utf8_string desc) except + nogil 

         Size size() except + nogil 
         bool empty() except + nogil 

         void clear() except + nogil 
         void insert(libcpp_utf8_string prefix, Param param) except + nogil 

         void remove(libcpp_utf8_string key) except + nogil 
         void removeAll(libcpp_utf8_string prefix) except + nogil 

         Param copy(libcpp_utf8_string prefix, bool) except + nogil 
         Param copy(libcpp_utf8_string prefix) except + nogil 

         # wrapped manually for overloading with dict parameter:
         bool update(Param p_old, bool add_unknow) except + nogil  # wrap-ignore
         bool update(Param p_old) except + nogil  # wrap-ignore

         void merge(Param toMerge) except + nogil 

         void setDefaults(Param defaults, libcpp_utf8_string prefix, bool showMessage) except + nogil 
         void setDefaults(Param defaults, libcpp_utf8_string prefix) except + nogil 
         void setDefaults(Param defaults) except + nogil 

         void checkDefaults(libcpp_utf8_string name, Param defaults, libcpp_utf8_string prefix) except + nogil 
         void checkDefaults(libcpp_utf8_string name, Param defaults) except + nogil 

         libcpp_vector[libcpp_utf8_string] getValidStrings(libcpp_utf8_string key) except + nogil

         void setValidStrings(libcpp_utf8_string key, libcpp_vector[libcpp_utf8_string] strings) except + nogil 
         void setMinInt(libcpp_utf8_string key, int min) except + nogil 
         void setMaxInt(libcpp_utf8_string key, int max) except + nogil 
         void setMinFloat(libcpp_utf8_string key, double min) except + nogil 
         void setMaxFloat(libcpp_utf8_string key, double max) except + nogil 

         #void parseCommandLine(int argc, char ** argv, String prefix) # wrap-ignore
         #void parseCommandLine(int argc, char ** argv) # wrap-ignore

         ParamIterator begin() except + nogil  # wrap-ignore
         ParamIterator end()   except + nogil  # wrap-ignore

cdef extern from "<OpenMS/DATASTRUCTURES/Param.h>" namespace "OpenMS::Param":

    cdef cppclass ParamIterator:
        # wrap-ignore
        # no-pxd-import
        ParamIterator operator++() except + nogil 
        ParamIterator operator--() except + nogil 
        libcpp_utf8_string getName() except + nogil 
        int operator==(ParamIterator) except + nogil 
        int operator!=(ParamIterator) except + nogil 
        int operator<(ParamIterator) except + nogil 
        int operator>(ParamIterator) except + nogil 
        int operator<=(ParamIterator) except + nogil 
        int operator>=(ParamIterator) except + nogil 

        # Returns the traceback of the opened and closed sections
        libcpp_vector[TraceInfo] getTrace() except + nogil 

cdef extern from "<OpenMS/DATASTRUCTURES/Param.h>" namespace "OpenMS::Param::ParamIterator":

    cdef cppclass TraceInfo:

        TraceInfo(libcpp_utf8_string n, libcpp_utf8_string d, bool o) except + nogil 
        TraceInfo(TraceInfo) except + nogil 

        # name of the node
        libcpp_string name
        # description of the node
        libcpp_string description
        # If it was opened (true) or closed (false)
        bool opened
