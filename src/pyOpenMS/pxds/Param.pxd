from libcpp.string cimport string as libcpp_utf8_string
from Types cimport *
from ParamValue cimport *
from ParamNode cimport *
from ParamEntry cimport *


# this class has addons, see the ./addons folder

cdef extern from "<OpenMS/DATASTRUCTURES/Param.h>" namespace "OpenMS":

    # pythonic helper functions in ../addons/Param.pyx !!!!

    cdef cppclass Param:

         # COMMENT: Helper functions for Python
         asDict() # wrap-ignore
         keys() # wrap-ignore
         items() # wrap-ignore
         values() # wrap-ignore
         update(dict) # wrap-ignore
         get(bytes key, default=None) # wrap-ignore
         __getitem__(bytes key) # wrap-ignore
         __setitem__(bytes key, value) # wrap-ignore

         Param() nogil except +
         Param(Param &) nogil except +
         bool operator==(Param) nogil except +

         void setValue(libcpp_utf8_string key, ParamValue val, libcpp_utf8_string desc, libcpp_vector[libcpp_utf8_string] tags) nogil except +
         void setValue(libcpp_utf8_string key, ParamValue val, libcpp_utf8_string desc) nogil except +
         void setValue(libcpp_utf8_string key, ParamValue val) nogil except +
         ParamValue getValue(libcpp_utf8_string key) nogil except +
         ValueType getValueType(libcpp_utf8_string key) nogil except +

         ParamEntry getEntry(libcpp_utf8_string) nogil except +
         bool exists(libcpp_utf8_string key) nogil except +

         void addTag(libcpp_utf8_string key, libcpp_utf8_string tag) nogil except +
         void addTags(libcpp_utf8_string key, libcpp_vector[libcpp_utf8_string] tags) nogil except +
         int hasTag(libcpp_utf8_string key, libcpp_utf8_string tag) nogil except +
         libcpp_vector[libcpp_utf8_string] getTags(libcpp_utf8_string key) nogil except +
         void clearTags(libcpp_utf8_string key) nogil except +

         libcpp_utf8_output_string getDescription(libcpp_utf8_string key) nogil except +
         void setSectionDescription(libcpp_utf8_string key, libcpp_utf8_string desc) nogil except +
         libcpp_utf8_output_string getSectionDescription(libcpp_utf8_string key) nogil except +

         Size size() nogil except +
         bool empty() nogil except +

         void clear() nogil except +
         void insert(libcpp_utf8_string prefix, Param param) nogil except +

         void remove(libcpp_utf8_string key) nogil except +
         void removeAll(libcpp_utf8_string prefix) nogil except +

         Param copy(libcpp_utf8_string prefix, bool) nogil except +
         Param copy(libcpp_utf8_string prefix) nogil except +

         # wrapped manually for overloading with dict parameter:
         bool update(Param p_old, bool add_unknow) nogil except + # wrap-ignore
         bool update(Param p_old) nogil except + # wrap-ignore

         void merge(Param toMerge) nogil except +

         void setDefaults(Param defaults, libcpp_utf8_string prefix, bool showMessage) nogil except +
         void setDefaults(Param defaults, libcpp_utf8_string prefix) nogil except +
         void setDefaults(Param defaults) nogil except +

         void checkDefaults(libcpp_utf8_string name, Param defaults, libcpp_utf8_string prefix) nogil except +
         void checkDefaults(libcpp_utf8_string name, Param defaults) nogil except +

         void setValidStrings(libcpp_utf8_string key, libcpp_vector[libcpp_utf8_string] strings) nogil except +
         void setMinInt(libcpp_utf8_string key, int min) nogil except +
         void setMaxInt(libcpp_utf8_string key, int max) nogil except +
         void setMinFloat(libcpp_utf8_string key, double min) nogil except +
         void setMaxFloat(libcpp_utf8_string key, double max) nogil except +

         #void parseCommandLine(int argc, char ** argv, String prefix) # wrap-ignore
         #void parseCommandLine(int argc, char ** argv) # wrap-ignore

         ParamIterator begin() nogil except + # wrap-ignore
         ParamIterator end()   nogil except + # wrap-ignore

cdef extern from "<OpenMS/DATASTRUCTURES/Param.h>" namespace "OpenMS::Param":

    cdef cppclass ParamIterator:
        # wrap-ignore
        # no-pxd-import
        ParamIterator operator++() nogil except +
        ParamIterator operator--() nogil except +
        libcpp_utf8_string getName() nogil except +
        int operator==(ParamIterator) nogil except +
        int operator!=(ParamIterator) nogil except +
        int operator<(ParamIterator) nogil except +
        int operator>(ParamIterator) nogil except +
        int operator<=(ParamIterator) nogil except +
        int operator>=(ParamIterator) nogil except +

        # Returns the traceback of the opened and closed sections
        libcpp_vector[TraceInfo] getTrace() nogil except +

cdef extern from "<OpenMS/DATASTRUCTURES/Param.h>" namespace "OpenMS::Param::ParamIterator":

    cdef cppclass TraceInfo:

        TraceInfo(libcpp_utf8_string n, libcpp_utf8_string d, bool o) nogil except +
        TraceInfo(TraceInfo) nogil except +

        # name of the node
        libcpp_string name
        # description of the node
        libcpp_string description
        # If it was opened (true) or closed (false)
        bool opened
