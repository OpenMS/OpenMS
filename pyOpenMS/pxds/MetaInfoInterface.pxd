from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from DataValue cimport *
from String cimport *
from Types cimport *

cdef extern from "<OpenMS/METADATA/MetaInfoInterface.h>" namespace "OpenMS":

    cdef cppclass MetaInfoInterface:
        # wrap-ignore

        MetaInfoInterface() nogil except +
        MetaInfoInterface(MetaInfoInterface) nogil except +

        bool operator==(MetaInfoInterface) nogil except +
        bool operator!=(MetaInfoInterface) nogil except +
        bool isMetaEmpty() nogil except +
        void clearMetaInfo() nogil except +

        # cython has a problem with inheritance of overloaded methods,
        # so we do not declare them here, but separately in each derived
        # class which we want to be wrapped:
        #
        #void getKeys(libcpp_vector[String] & keys)
        #void getKeys(libcpp_vector[unsigned int] & keys)
        #DataValue getMetaValue(unsigned int) nogil except +
        #DataValue getMetaValue(String) nogil except +
        #void setMetaValue(unsigned int, DataValue) nogil except +
        #void setMetaValue(String, DataValue) nogil except +
        #bool metaValueExists(String) nogil except +
        #bool metaValueExists(unsigned int) nogil except +
        #void removeMetaValue(String) nogil except +
        #void removeMetaValue(unsigned int) nogil except +
