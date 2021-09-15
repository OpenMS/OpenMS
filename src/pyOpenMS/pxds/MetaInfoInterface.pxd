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

        MetaInfoInterface() nogil except +
        MetaInfoInterface(MetaInfoInterface &) nogil except +

        bool operator==(MetaInfoInterface) nogil except +
        bool operator!=(MetaInfoInterface) nogil except +
        bool isMetaEmpty() nogil except + # wrap-doc:Returns if the MetaInfo is empty
        void clearMetaInfo() nogil except + # wrap-doc:Removes all meta values

        MetaInfoRegistry metaRegistry() nogil except + # wrap-doc:Returns a reference to the MetaInfoRegistry

        # Cython has a problem with inheritance of overloaded methods, so we
        # can only declare one of the methods here (most people will want
        # to use the String-based methods).
        #
        void getKeys(libcpp_vector[String] & keys) nogil except + # wrap-doc:Fills the given vector with a list of all keys for which a value is set
        #void getKeys(libcpp_vector[unsigned int] & keys)
        DataValue getMetaValue(String) nogil except + # wrap-doc:Returns the value corresponding to a string, or DataValue::EMPTY if not found
        #DataValue getMetaValue(String, DataValue) nogil except +
        #DataValue getMetaValue(unsigned int) nogil except +
        #DataValue getMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except + # wrap-doc:Sets the DataValue corresponding to a name
        #void setMetaValue(unsigned int, DataValue) nogil except +
        bool metaValueExists(String) nogil except + # wrap-doc:Returns whether an entry with the given name exists
        #bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except + # wrap-doc:Removes the DataValue corresponding to `name` if it exists
        #void removeMetaValue(unsigned int) nogil except +

