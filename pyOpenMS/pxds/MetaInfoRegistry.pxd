from Types cimport *
from Types cimport *
from String cimport *

cdef extern from "<OpenMS/METADATA/MetaInfoRegistry.h>" namespace "OpenMS":
    
    cdef cppclass MetaInfoRegistry "OpenMS::MetaInfoRegistry":
        MetaInfoRegistry() nogil except +
        MetaInfoRegistry(MetaInfoRegistry) nogil except +
        UInt registerName(String & name, String & description, String & unit) nogil except +
        void setDescription(UInt index, String & description) nogil except +
        void setDescription(String & name, String & description) nogil except +
        void setUnit(UInt index, String & unit) nogil except +
        void setUnit(String & name, String & unit) nogil except +
        UInt getIndex(String & name) nogil except +
        String getName(UInt index) nogil except +
        String getDescription(UInt index) nogil except +
        String getDescription(String & name) nogil except +
        String getUnit(UInt index) nogil except +
        String getUnit(String & name) nogil except +

