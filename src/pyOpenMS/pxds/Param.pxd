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
         asDict() # wrap-ignore
         keys() # wrap-ignore
         items() # wrap-ignore
         values() # wrap-ignore
         update(dict) # wrap-ignore
         get(bytes key, default=None) # wrap-ignore
         __getitem__(bytes key) # wrap-ignore
         __setitem__(bytes key, value) # wrap-ignore

         Param() nogil except +
         Param(Param) nogil except +
         bool operator==(Param) nogil except +

         void setValue(libcpp_string key, ParamValue val, libcpp_string desc, libcpp_vector[libcpp_string] tags) nogil except +
         void setValue(libcpp_string key, ParamValue val, libcpp_string desc) nogil except +
         void setValue(libcpp_string key, ParamValue val) nogil except +
         ParamValue getValue(libcpp_string key) nogil except +
         ParamEntry getEntry(libcpp_string) nogil except +
         bool exists(libcpp_string key) nogil except +

         void addTag(libcpp_string key, libcpp_string tag) nogil except +
         void addTags(libcpp_string key, libcpp_vector[libcpp_string] tags) nogil except +
         int hasTag(libcpp_string key, libcpp_string tag) nogil except +
         libcpp_vector[libcpp_string] getTags(libcpp_string key) nogil except +
         void clearTags(libcpp_string key) nogil except +

         libcpp_string getDescription(libcpp_string key) nogil except +
         void setSectionDescription(libcpp_string key, libcpp_string desc) nogil except +
         libcpp_string getSectionDescription(libcpp_string key) nogil except +

         Size size() nogil except +
         bool empty() nogil except +

         void clear() nogil except +
         void insert(libcpp_string prefix, Param param) nogil except +

         void remove(libcpp_string key) nogil except +
         void removeAll(libcpp_string prefix) nogil except +

         Param copy(libcpp_string prefix, bool) nogil except +
         Param copy(libcpp_string prefix) nogil except +

         # wrapped manually for overloading with dict parameter:
         bool update(Param p_old, bool add_unknow) nogil except + # wrap-ignore
         bool update(Param p_old) nogil except + # wrap-ignore

         void merge(Param toMerge) nogil except +

         void setDefaults(Param defaults, libcpp_string prefix, bool showMessage) nogil except +
         void setDefaults(Param defaults, libcpp_string prefix) nogil except +
         void setDefaults(Param defaults) nogil except +

         void checkDefaults(libcpp_string name, Param defaults, libcpp_string prefix) nogil except +
         void checkDefaults(libcpp_string name, Param defaults) nogil except +

         void setValidStrings(libcpp_string key, libcpp_vector[libcpp_string] strings) nogil except +
         void setMinInt(libcpp_string key, int min) nogil except +
         void setMaxInt(libcpp_string key, int max) nogil except +
         void setMinFloat(libcpp_string key, double min) nogil except +
         void setMaxFloat(libcpp_string key, double max) nogil except +

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
        libcpp_string getName() nogil except +
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

        TraceInfo(libcpp_string n, libcpp_string d, bool o) nogil except +
        TraceInfo(TraceInfo) nogil except +

        # name of the node
        libcpp_string name
        # description of the node
        libcpp_string description
        # If it was opened (true) or closed (false)
        bool opened
