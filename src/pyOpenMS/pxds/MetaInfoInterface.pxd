from Types cimport *
from DataValue cimport *
from String cimport *
from MetaInfoRegistry cimport *

cdef extern from "<OpenMS/METADATA/MetaInfoInterface.h>" namespace "OpenMS":

    cdef cppclass MetaInfoInterface:
        # wrap-doc:
        #  Interface for classes that can store arbitrary meta information
        #  (Type-Name-Value tuples).
        #  
        #  MetaInfoInterface is a base class for all classes that use one MetaInfo
        #  object as member.  If you want to add meta information to a class, let it
        #  publicly inherit the MetaInfoInterface.  Meta information is an array of
        #  Type-Name-Value tuples.
        #  
        #  Usage:
        #
        #  .. code-block:: python
        #  
        #    k = []
        #    exp.getKeys(k) # explore available key-value pairs
        #    exp.getMetaValue("someMetaName")

        MetaInfoInterface() except + nogil 
        MetaInfoInterface(MetaInfoInterface &) except + nogil 

        bool operator==(MetaInfoInterface) except + nogil 
        bool operator!=(MetaInfoInterface) except + nogil 
        bool isMetaEmpty() except + nogil  # wrap-doc:Returns if the MetaInfo is empty
        void clearMetaInfo() except + nogil  # wrap-doc:Removes all meta values

        MetaInfoRegistry metaRegistry() except + nogil  # wrap-doc:Returns a reference to the MetaInfoRegistry

        # Cython has a problem with inheritance of overloaded methods, so we
        # can only declare one of the methods here (most people will want
        # to use the String-based methods).
        #
        void getKeys(libcpp_vector[String] & keys) except + nogil  # wrap-doc:Fills the given vector with a list of all keys for which a value is set
        #void getKeys(libcpp_vector[unsigned int] & keys)
        DataValue getMetaValue(String) except + nogil  # wrap-doc:Returns the value corresponding to a string, or DataValue::EMPTY if not found
        #DataValue getMetaValue(String, DataValue) except + nogil 
        #DataValue getMetaValue(unsigned int) except + nogil 
        #DataValue getMetaValue(unsigned int, DataValue) except + nogil 
        void setMetaValue(String, DataValue) except + nogil  # wrap-doc:Sets the DataValue corresponding to a name
        #void setMetaValue(unsigned int, DataValue) except + nogil 
        bool metaValueExists(String) except + nogil  # wrap-doc:Returns whether an entry with the given name exists
        #bool metaValueExists(unsigned int) except + nogil 
        void removeMetaValue(String) except + nogil  # wrap-doc:Removes the DataValue corresponding to `name` if it exists
        #void removeMetaValue(unsigned int) except + nogil 

