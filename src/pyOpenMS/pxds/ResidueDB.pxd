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
        const Residue * getResidue(String & name) nogil except +
        const Residue * getModifiedResidue(String & name) nogil except +
        const Residue * getModifiedResidue(Residue * residue, String & name) nogil except +
        libcpp_set[ const Residue * ] getResidues(String & residue_set) nogil except +
        libcpp_set[ libcpp_string ] getResidueSets() nogil except +
        void setResidues(String & filename) nogil except +
        void addResidue(Residue & residue) nogil except +
        bool hasResidue(String & name) nogil except +
        # bool hasResidue(Residue * residue) nogil except + # does not really work as the ptr is different

# COMMENT: wrap static methods
cdef extern from "<OpenMS/CHEMISTRY/ResidueDB.h>" namespace "OpenMS::ResidueDB":
    
    ResidueDB* getInstance() nogil except + # wrap-ignore

