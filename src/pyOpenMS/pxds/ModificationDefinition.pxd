from libcpp cimport bool
from Types cimport *
from String cimport *
from ResidueModification cimport *

cdef extern from "<OpenMS/CHEMISTRY/ModificationDefinition.h>" namespace "OpenMS":
    
    cdef cppclass ModificationDefinition "OpenMS::ModificationDefinition":
        # wrap-hash:
        #  getModificationName().c_str()

        ModificationDefinition() except + nogil 
        ModificationDefinition(ModificationDefinition &) except + nogil 
        ModificationDefinition(const String &mod) except + nogil 
        ModificationDefinition(const String &mod, bool fixed) except + nogil 
        ModificationDefinition(const String &mod, bool fixed, UInt max_occur) except + nogil 
        ModificationDefinition(ResidueModification &mod) except + nogil 
        ModificationDefinition(ResidueModification &mod, bool fixed) except + nogil 
        ModificationDefinition(ResidueModification &mod, bool fixed, UInt max_occur) except + nogil 

        bool operator==(ModificationDefinition &rhs) except + nogil 
        bool operator!=(ModificationDefinition &rhs) except + nogil 
        bool operator<(ModificationDefinition &) except + nogil 

        void setFixedModification(bool fixed) except + nogil  # wrap-doc:Sets whether this modification definition is fixed or variable (modification must occur vs. can occur)
        bool isFixedModification() except + nogil  # wrap-doc:Returns if the modification if fixed true, else false
        void setMaxOccurrences(UInt num) except + nogil  # wrap-doc:Sets the maximal number of occurrences per peptide (unbounded if 0)
        UInt getMaxOccurrences() except + nogil  # wrap-doc:Returns the maximal number of occurrences per peptide
        String getModificationName() except + nogil  # wrap-doc:Returns the name of the modification
        void setModification(const String &modification) except + nogil  # wrap-doc:Sets the modification, allowed are unique names provided by ModificationsDB

        ResidueModification getModification() except + nogil 

