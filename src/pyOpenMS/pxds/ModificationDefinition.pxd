from libcpp cimport bool
from Types cimport *
from String cimport *
from ResidueModification cimport *

cdef extern from "<OpenMS/CHEMISTRY/ModificationDefinition.h>" namespace "OpenMS":
    
    cdef cppclass ModificationDefinition "OpenMS::ModificationDefinition":
        # wrap-hash:
        #   getModificationName().c_str()

        ModificationDefinition() nogil except +
        ModificationDefinition(ModificationDefinition) nogil except +
        ModificationDefinition(const String &mod) nogil except +
        ModificationDefinition(const String &mod, bool fixed) nogil except +
        ModificationDefinition(const String &mod, bool fixed, UInt max_occur) nogil except +
        ModificationDefinition(ResidueModification &mod) nogil except +
        ModificationDefinition(ResidueModification &mod, bool fixed) nogil except +
        ModificationDefinition(ResidueModification &mod, bool fixed, UInt max_occur) nogil except +

        bool operator==(ModificationDefinition &rhs) nogil except +
        bool operator!=(ModificationDefinition &rhs) nogil except +
        bool operator<(ModificationDefinition &) nogil except +

        void setFixedModification(bool fixed) nogil except +
        bool isFixedModification() nogil except +
        void setMaxOccurrences(UInt num) nogil except +
        UInt getMaxOccurrences() nogil except +
        String getModificationName() nogil except +
        void setModification(const String &modification) nogil except +

        ResidueModification getModification() nogil except +

