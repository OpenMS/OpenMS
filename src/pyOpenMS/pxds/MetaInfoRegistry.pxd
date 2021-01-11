from Types cimport *
from String cimport *

cdef extern from "<OpenMS/METADATA/MetaInfoRegistry.h>" namespace "OpenMS":
    
    cdef cppclass MetaInfoRegistry "OpenMS::MetaInfoRegistry":
        MetaInfoRegistry() nogil except +
        MetaInfoRegistry(MetaInfoRegistry) nogil except +
        UInt registerName(const String & name, const String & description, const String & unit) nogil except +
        void setDescription(UInt index, const String & description) nogil except +
        void setDescription(const String & name, const String & description) nogil except +
        void setUnit(UInt index, const String & unit) nogil except +
        void setUnit(const String & name, const String & unit) nogil except +
        UInt getIndex(const String & name) nogil except +
        String getName(UInt index) nogil except +
        String getDescription(UInt index) nogil except +
        String getDescription(const String & name) nogil except +
        String getUnit(UInt index) nogil except +
        String getUnit(const String & name) nogil except +

