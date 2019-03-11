from Types cimport *
from String cimport *
from Residue cimport *

cdef extern from "<OpenMS/CHEMISTRY/ResidueDB.h>" namespace "OpenMS":
    
    cdef cppclass ResidueDB "OpenMS::ResidueDB":
        # wrap-manual-memory:
        #   cdef AutowrapPtrHolder[_ResidueDB] inst

        ResidueDB(ResidueDB) nogil except + #wrap-ignore

        Size getNumberOfResidues() nogil except +
        Size getNumberOfModifiedResidues() nogil except +
        const Residue * getResidue(const String & name) nogil except +
        const Residue * getModifiedResidue(const String & name) nogil except +
        const Residue * getModifiedResidue(Residue * residue, const String & name) nogil except +
        libcpp_set[ const Residue * ] getResidues(const String & residue_set) nogil except +
        libcpp_set[ String ] getResidueSets() nogil except +
        bool hasResidue(const String & name) nogil except +
        # bool hasResidue(Residue * residue) nogil except + # does not really work as the ptr is different

# COMMENT: wrap static methods
cdef extern from "<OpenMS/CHEMISTRY/ResidueDB.h>" namespace "OpenMS::ResidueDB":
    
    ResidueDB* getInstance() nogil except + # wrap-ignore

