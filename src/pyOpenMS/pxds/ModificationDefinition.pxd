from libcpp cimport bool
from Types cimport *
from String cimport *
from ResidueModification cimport *

cdef extern from "<OpenMS/CHEMISTRY/ModificationDefinition.h>" namespace "OpenMS":
    
    cdef cppclass ModificationDefinition "OpenMS::ModificationDefinition":
        # wrap-hash:
        #   getModificationName().c_str()

        ModificationDefinition() nogil except +
        ModificationDefinition(ModificationDefinition &) nogil except +
        ModificationDefinition(const String &mod) nogil except +
        ModificationDefinition(const String &mod, bool fixed) nogil except +
        ModificationDefinition(const String &mod, bool fixed, UInt max_occur) nogil except +
        ModificationDefinition(ResidueModification &mod) nogil except +
        ModificationDefinition(ResidueModification &mod, bool fixed) nogil except +
        ModificationDefinition(ResidueModification &mod, bool fixed, UInt max_occur) nogil except +

        bool operator==(ModificationDefinition &rhs) nogil except +
        bool operator!=(ModificationDefinition &rhs) nogil except +
        bool operator<(ModificationDefinition &) nogil except +

        void setFixedModification(bool fixed) nogil except + # wrap-doc:Sets whether this modification definition is fixed or variable (modification must occur vs. can occur)
        bool isFixedModification() nogil except + # wrap-doc:Returns if the modification if fixed true, else false
        void setMaxOccurrences(UInt num) nogil except + # wrap-doc:Sets the maximal number of occurrences per peptide (unbounded if 0)
        UInt getMaxOccurrences() nogil except + # wrap-doc:Returns the maximal number of occurrences per peptide
        String getModificationName() nogil except + # wrap-doc:Returns the name of the modification
        void setModification(const String &modification) nogil except + # wrap-doc:Sets the modification, allowed are unique names provided by ModificationsDB

        ResidueModification getModification() nogil except +

