from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from DataValue cimport *
from String cimport *
from Types cimport *
from MetaInfoDescription cimport *

cdef extern from "<OpenMS/METADATA/DataArrays.h>" namespace "OpenMS::DataArrays":

    cdef cppclass FloatDataArray(MetaInfoDescription):
        # wrap-inherits:
        #  MetaInfoDescription

        FloatDataArray() nogil except +
        FloatDataArray(FloatDataArray) nogil except + #wrap-ignore

        Size size() nogil except +
        float operator[](int) nogil except +
        void clear() nogil except +
        void push_back(float) nogil except +

        void getKeys(libcpp_vector[String] & keys) nogil except +
        void getKeys(libcpp_vector[unsigned int] & keys) nogil except + # wrap-as:getKeysAsIntegers
        DataValue getMetaValue(unsigned int) nogil except +
        DataValue getMetaValue(String) nogil except +
        void setMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        void removeMetaValue(unsigned int) nogil except +

    cdef cppclass StringDataArray(MetaInfoDescription):
        # wrap-inherits:
        #  MetaInfoDescription

        StringDataArray() nogil except +
        StringDataArray(StringDataArray) nogil except + #wrap-ignore

        Size size() nogil except +
        String operator[](int) nogil except +
        void clear() nogil except +
        void push_back(String) nogil except +

        void getKeys(libcpp_vector[String] & keys) nogil except +
        void getKeys(libcpp_vector[unsigned int] & keys) nogil except + # wrap-as:getKeysAsIntegers
        DataValue getMetaValue(unsigned int) nogil except +
        DataValue getMetaValue(String) nogil except +
        void setMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        void removeMetaValue(unsigned int) nogil except +

    cdef cppclass IntegerDataArray(MetaInfoDescription):
        # wrap-inherits:
        #  MetaInfoDescription

        IntegerDataArray() nogil except +
        IntegerDataArray(IntegerDataArray) nogil except + #wrap-ignore

        Size size() nogil except +
        Int operator[](int) nogil except +
        void clear() nogil except +
        void push_back(Int) nogil except +

        void getKeys(libcpp_vector[String] & keys) nogil except +
        void getKeys(libcpp_vector[unsigned int] & keys) nogil except + # wrap-as:getKeysAsIntegers
        DataValue getMetaValue(unsigned int) nogil except +
        DataValue getMetaValue(String) nogil except +
        void setMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        void removeMetaValue(unsigned int) nogil except +

