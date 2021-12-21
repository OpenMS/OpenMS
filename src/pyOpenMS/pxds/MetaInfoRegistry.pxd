from Types cimport *
from String cimport *

cdef extern from "<OpenMS/METADATA/MetaInfoRegistry.h>" namespace "OpenMS":
    
    cdef cppclass MetaInfoRegistry "OpenMS::MetaInfoRegistry":
        # wrap-doc:
                #   Registry which assigns unique integer indices to strings
                #   -----
                #   When registering a new name an index >= 1024 is assigned.
                #   Indices from 1 to 1023 are reserved for fast access and will never change:
                #   1 - isotopic_range
                #   2 - cluster_id
                #   3 - label
                #   4 - icon
                #   5 - color
                #   6 - RT
                #   7 - MZ
                #   8 - predicted_RT
                #   9 - predicted_RT_p_value
                #   10 - spectrum_reference
                #   11 - ID
                #   12 - low_quality
                #   13 - charge

        MetaInfoRegistry() nogil except +
        MetaInfoRegistry(MetaInfoRegistry &) nogil except +
        UInt registerName(const String & name, const String & description, const String & unit) nogil except + # wrap-doc:Registers a string, stores its description and unit, and returns the corresponding index. If the string is already registered, it returns the index of the string
        void setDescription(UInt index, const String & description) nogil except + # wrap-doc:Sets the description (String), corresponding to an index
        void setDescription(const String & name, const String & description) nogil except + # wrap-doc:Sets the description (String), corresponding to a name
        void setUnit(UInt index, const String & unit) nogil except + # wrap-doc:Sets the unit (String), corresponding to an index
        void setUnit(const String & name, const String & unit) nogil except + # wrap-doc:Sets the unit (String), corresponding to a name
        UInt getIndex(const String & name) nogil except + # wrap-doc:Returns the integer index corresponding to a string. If the string is not registered, returns UInt(-1) (= UINT_MAX)
        String getName(UInt index) nogil except + # wrap-doc:Returns the corresponding name to an index
        String getDescription(UInt index) nogil except + # wrap-doc:Returns the description of an index
        String getDescription(const String & name) nogil except + # wrap-doc:Returns the description of a name
        String getUnit(UInt index) nogil except + # wrap-doc:Returns the unit of an index
        String getUnit(const String & name) nogil except + # wrap-doc:Returns the unit of a name


