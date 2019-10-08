from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from DataValue cimport *
from String cimport *
from Types cimport *
from MetaInfoRegistry cimport *

cdef extern from "<OpenMS/METADATA/MetaInfo.h>" namespace "OpenMS":

    cdef cppclass MetaInfo:

        MetaInfo() nogil except +
        MetaInfo(MetaInfo) nogil except +

        # returns the value corresponding to a string
        DataValue getValue(String name) nogil except +
        # returns the value corresponding to an index
        DataValue getValue(UInt index) nogil except +

        # returns the value corresponding to a string, with default
        DataValue getValue(String name, DataValue default_value) nogil except +
        # returns the value corresponding to an index, with default
        DataValue getValue(UInt index, DataValue default_value) nogil except +

        # returns if this MetaInfo is set
        bool exists(String name) nogil except +
        # returns if this MetaInfo is set
        bool exists(UInt index) nogil except +

        # sets the DataValue corresponding to a name
        void setValue(String name, DataValue value) nogil except +
        #  sets the DataValue corresponding to an index
        void setValue(UInt index, DataValue value) nogil except +

        # Removes the DataValue corresponding to @p name if it exists
        void removeValue(String name) nogil except +
        # Removes the DataValue corresponding to @p index if it exists
        void removeValue(UInt index) nogil except +

        # returns a reference to the MetaInfoRegistry
        # static MetaInfoRegistry registry() nogil except +

        # fills the given vector with a list of all keys for which a value is set
        void getKeys(libcpp_vector[String] & keys) nogil except +

        # fills the given vector with a list of all keys for which a value is set
        void getKeys(libcpp_vector[unsigned int] & keys) nogil except +  # wrap-as:getKeysAsIntegers

        # returns if the MetaInfo is empty
        bool empty() nogil except +

        # removes all meta values
        void clear() nogil except +

        MetaInfoRegistry  registry() nogil except +
