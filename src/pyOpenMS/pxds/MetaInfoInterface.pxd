from Types cimport *
from DataValue cimport *
from String cimport *
from MetaInfoRegistry cimport *

cdef extern from "<OpenMS/METADATA/MetaInfoInterface.h>" namespace "OpenMS":

    cdef cppclass MetaInfoInterface:
        # wrap-doc:
        #   Interface for classes that can store arbitrary meta information
        #   (Type-Name-Value tuples).
        #   -----
        #   MetaInfoInterface is a base class for all classes that use one MetaInfo
        #   object as member.  If you want to add meta information to a class, let it
        #   publicly inherit the MetaInfoInterface.  Meta information is an array of
        #   Type-Name-Value tuples.
        #   -----
        #   Usage:
        #     k = []
        #     exp.getKeys(k) # explore available key-value pairs
        #     exp.getMetaValue("someMetaName")
        #   -----

        MetaInfoInterface() nogil except +
        MetaInfoInterface(MetaInfoInterface) nogil except +

        bool operator==(MetaInfoInterface) nogil except +
        bool operator!=(MetaInfoInterface) nogil except +
        bool isMetaEmpty() nogil except +
        void clearMetaInfo() nogil except +

        MetaInfoRegistry metaRegistry() nogil except +

        # Cython has a problem with inheritance of overloaded methods, so we
        # can only declare one of the methods here (most people will want
        # to use the String-based methods).
        #
        void getKeys(libcpp_vector[String] & keys)
        #void getKeys(libcpp_vector[unsigned int] & keys)
        DataValue getMetaValue(String) nogil except +
        #DataValue getMetaValue(String, DataValue) nogil except +
        #DataValue getMetaValue(unsigned int) nogil except +
        #DataValue getMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        #void setMetaValue(unsigned int, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        #bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        #void removeMetaValue(unsigned int) nogil except +

