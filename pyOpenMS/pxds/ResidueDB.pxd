from Types cimport *
from libcpp cimport bool
from libcpp.set cimport set as libcpp_set
from Map cimport *
from String cimport *
from Residue cimport *

cdef extern from "<OpenMS/CHEMISTRY/ResidueDB.h>" namespace "OpenMS":
    
    cdef cppclass ResidueDB "OpenMS::ResidueDB":
        ResidueDB(ResidueDB) nogil except + #wrap-ignore
        Size getNumberOfResidues() nogil except +
        Size getNumberOfModifiedResidues() nogil except +
        # POINTER # Residue * getResidue(String & name) nogil except +
        # POINTER # Residue * getModifiedResidue(String & name) nogil except +
        # POINTER # Residue * getModifiedResidue(Residue * residue, String & name) nogil except +
        # POINTER # libcpp_set[ Residue * ] getResidues(String & residue_set) nogil except +
        libcpp_set[ String ]  getResidueSets() nogil except +
        void setResidues(String & filename) nogil except +
        void addResidue(Residue & residue) nogil except +
        bool hasResidue(String & name) nogil except +
        # POINTER # bool hasResidue(Residue * residue) nogil except +
        # ResidueIterator beginResidue() nogil except +
        # ResidueIterator endResidue() nogil except +
        # ResidueConstIterator beginResidue() nogil except +
        # ResidueConstIterator endResidue() nogil except +
        # POINTER # ResidueDB * getInstance() nogil except +

