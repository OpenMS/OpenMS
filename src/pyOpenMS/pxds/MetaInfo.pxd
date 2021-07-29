from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from DataValue cimport *
from String cimport *
from Types cimport *
from MetaInfoRegistry cimport *

cdef extern from "<OpenMS/METADATA/MetaInfo.h>" namespace "OpenMS":

    cdef cppclass MetaInfo:
        # wrap-doc:
                #   A Type-Name-Value tuple class
                #   -----
                #   MetaInfo maps an index (an integer corresponding to a string) to
                #   DataValue objects.  The mapping of strings to the index is performed by
                #   the MetaInfoRegistry, which can be accessed by the method registry().
                #   -----
                #   There are two versions of nearly all members. One which operates with a
                #   string name and another one which operates on an index. The index version
                #   is always faster, as it does not need to look up the index corresponding
                #   to the string in the MetaInfoRegistry.
                #   -----
                #   If you wish to add a MetaInfo member to a class, consider deriving that
                #   class from MetaInfoInterface, instead of simply adding MetaInfo as
                #   member. MetaInfoInterface implements a full interface to a MetaInfo
                #   member and is more memory efficient if no meta info gets added.

        MetaInfo() nogil except +
        MetaInfo(MetaInfo) nogil except +

        # returns the value corresponding to a string
        DataValue getValue(String name) nogil except + # wrap-doc:Returns the value corresponding to a string
        # returns the value corresponding to an index
        DataValue getValue(UInt index) nogil except + # wrap-doc:Returns the value corresponding to an index

        # returns the value corresponding to a string, with default
        DataValue getValue(String name, DataValue default_value) nogil except + # wrap-doc:Returns the value corresponding to a string, with default
        # returns the value corresponding to an index, with default
        DataValue getValue(UInt index, DataValue default_value) nogil except + # wrap-doc:Returns the value corresponding to an index, with default

        # returns if this MetaInfo is set
        bool exists(String name) nogil except + # wrap-doc:Returns if this MetaInfo is set
        # returns if this MetaInfo is set
        bool exists(UInt index) nogil except + # wrap-doc:Returns if this MetaInfo is set

        # sets the DataValue corresponding to a name
        void setValue(String name, DataValue value) nogil except + # wrap-doc:Sets the DataValue corresponding to a name
        #  sets the DataValue corresponding to an index
        void setValue(UInt index, DataValue value) nogil except + # wrap-doc:Sets the DataValue corresponding to an index

        # Removes the DataValue corresponding to @p name if it exists
        void removeValue(String name) nogil except + # wrap-doc:Removes the DataValue corresponding to `name` if it exists
        # Removes the DataValue corresponding to @p index if it exists
        void removeValue(UInt index) nogil except + # wrap-doc:Removes the DataValue corresponding to `index` if it exists

        # returns a reference to the MetaInfoRegistry
        # static MetaInfoRegistry registry() nogil except +

        # fills the given vector with a list of all keys for which a value is set
        void getKeys(libcpp_vector[String] & keys) nogil except + # wrap-doc:Fills the given vector with a list of all keys for which a value is set

        # fills the given vector with a list of all keys for which a value is set
        void getKeys(libcpp_vector[unsigned int] & keys) nogil except +  # wrap-as:getKeysAsIntegers

        # returns if the MetaInfo is empty
        bool empty() nogil except + # wrap-doc:Returns if the MetaInfo is empty

        # removes all meta values
        void clear() nogil except + # wrap-doc:Removes all meta values

        MetaInfoRegistry  registry() nogil except +
