from Types cimport *
from String cimport *
from Residue cimport *

cdef extern from "<OpenMS/CHEMISTRY/ResidueDB.h>" namespace "OpenMS":
    
    cdef cppclass ResidueDB "OpenMS::ResidueDB":
        # wrap-manual-memory:
        #   cdef AutowrapPtrHolder[_ResidueDB] inst

        ResidueDB(ResidueDB) nogil except + #wrap-ignore

        Size getNumberOfResidues() nogil except + # wrap-doc:Returns the number of residues stored
        Size getNumberOfModifiedResidues() nogil except + # wrap-doc:Returns the number of modified residues stored
        const Residue * getResidue(const String & name) nogil except + # wrap-doc:Returns a pointer to the residue with name, 3 letter code or 1 letter code name
        const Residue * getModifiedResidue(const String & name) nogil except + # wrap-doc:Returns a pointer to a modified residue given a modification name
        const Residue * getModifiedResidue(Residue * residue, const String & name) nogil except + # wrap-doc:Returns a pointer to a modified residue given a residue and a modification name
        libcpp_set[ const Residue * ] getResidues(const String & residue_set) nogil except + # wrap-doc:Returns a set of all residues stored in this residue db
        libcpp_set[ String ] getResidueSets() nogil except + # wrap-doc:Returns all residue sets that are registered which this instance
        bool hasResidue(const String & name) nogil except + # wrap-doc:Returns true if the db contains a residue with the given name
        # bool hasResidue(Residue * residue) nogil except + # does not really work as the ptr is different

# COMMENT: wrap static methods
cdef extern from "<OpenMS/CHEMISTRY/ResidueDB.h>" namespace "OpenMS::ResidueDB":
    
    ResidueDB* getInstance() nogil except + # wrap-ignore

