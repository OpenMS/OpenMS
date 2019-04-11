from Types cimport *
from DataValue cimport *
from String cimport *
from MetaInfoRegistry cimport *

cdef extern from "<OpenMS/METADATA/MetaInfoInterface.h>" namespace "OpenMS":

    cdef cppclass MetaInfoInterface:

        MetaInfoInterface() nogil except +
        MetaInfoInterface(MetaInfoInterface) nogil except +

        bool operator==(MetaInfoInterface) nogil except +
        bool operator!=(MetaInfoInterface) nogil except +
        bool isMetaEmpty() nogil except +
        void clearMetaInfo() nogil except +

        MetaInfoRegistry metaRegistry() nogil except +

        # Cython has a problem with inheritance of overloaded methods, so we
        # can only declare one of the two methods here (most people will want
        # to use the String-based methods).
        #
        void getKeys(libcpp_vector[String] & keys)
        #void getKeys(libcpp_vector[unsigned int] & keys)
        DataValue getMetaValue(String) nogil except +
        #DataValue getMetaValue(unsigned int) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        #void setMetaValue(unsigned int, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        #bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        #void removeMetaValue(unsigned int) nogil except +
