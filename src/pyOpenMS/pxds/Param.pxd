from libcpp.string cimport string as libcpp_string
from libcpp.vector cimport vector as libcpp_vector
from libcpp.set cimport set as libcpp_set
from libcpp cimport bool
from DataValue cimport *
from String cimport *
from StringList cimport *
from ParamNode cimport *
from ParamEntry cimport *

# this class has addons, see the ./addons folder

cdef extern from "<OpenMS/DATASTRUCTURES/Param.h>" namespace "OpenMS":

    # pythonic helper functions in ../addons/Param.pyx !!!!

    cdef cppclass Param:
         Param() nogil except +
         Param(Param) nogil except +
         bool operator==(Param) nogil except +

         void setValue(String key, DataValue val, String desc, StringList tags) nogil except +
         void setValue(String key, DataValue val, String desc) nogil except +
         void setValue(String key, DataValue val) nogil except +
         DataValue getValue(String key) nogil except +
         ParamEntry getEntry(String) nogil except +
         int exists(String key) nogil except +

         void addTag(String key, String tag) nogil except +
         void addTags(String key, StringList tags) nogil except +
         int hasTag(String key, String tag) nogil except +
         StringList getTags(String key) nogil except +
         void clearTags(String key) nogil except +

         libcpp_string getDescription(String key) nogil except +
         void setSectionDescription(String key, String desc) nogil except +
         libcpp_string getSectionDescription(String key) nogil except +

         Size size() nogil except +
         bool empty() nogil except +

         void clear() nogil except +
         void insert(String prefix, Param param) nogil except +

         void remove(String key) nogil except +
         void removeAll(String prefix) nogil except +

         Param copy(String prefix, bool) nogil except +
         Param copy(String prefix) nogil except +

         # wrapped manually for overloading with dict parameter:
         bool update(Param p_old, bool add_unknow) nogil except + # wrap-ignore
         bool update(Param p_old) nogil except + # wrap-ignore

         void merge(Param toMerge) nogil except +

         void setDefaults(Param defaults, String previx, bool showMessage) nogil except +
         void setDefaults(Param defaults, String previx) nogil except +
         void setDefaults(Param defaults) nogil except +

         void checkDefaults(String name, Param defaults, String prefix) nogil except +
         void checkDefaults(String name, Param defaults) nogil except +

         void setValidStrings(String key, libcpp_vector[String] strings) nogil except +
         void setMinInt(String key, int min) nogil except +
         void setMaxInt(String key, int max) nogil except +
         void setMinFloat(String key, double min) nogil except +
         void setMaxFloat(String key, double max) nogil except +

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
        String getName() nogil except +
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

        TraceInfo(String n, String d, bool o) nogil except +
        TraceInfo(TraceInfo) nogil except +

        # name of the node
        String name
        # description of the node
        String description
        # If it was opened (true) or closed (false)
        bool opened
