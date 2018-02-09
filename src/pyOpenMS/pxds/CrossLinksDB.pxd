from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from Map cimport *
from String cimport *
from ResidueModification cimport *

# see ../addons/CrossLinksDB.pyx
cdef extern from "<OpenMS/CHEMISTRY/CrossLinksDB.h>" namespace "OpenMS":
    
    cdef cppclass CrossLinksDB:
        # wrap-manual-memory:
        #     cdef AutowrapPtrHolder[_CrossLinksDB] inst

        CrossLinksDB(CrossLinksDB) nogil except + #wrap-ignore

        Size getNumberOfModifications() nogil except +

        ResidueModification getModification(Size index) nogil except +

        void searchModifications(libcpp_set[ const ResidueModification * ] & mods,
                                 String mod_name, String residue,
                                 TermSpecificity term_spec) nogil except +

        ResidueModification getModification(String & mod_name, String & residue, TermSpecificity term_spec) nogil except +

        bool has(String modification) nogil except +

        void addModification(ResidueModification * new_mod) nogil except +

        Size findModificationIndex(String & mod_name) nogil except +

        void searchModificationsByDiffMonoMass(libcpp_vector[ String ] & mods, double mass, double max_error,
                                               String & residue, TermSpecificity term_spec) nogil except +

        const ResidueModification* getBestModificationByMonoMass(double mass, double max_error,
                                                           String residue,
                                                           TermSpecificity term_spec) nogil except +

        const ResidueModification* getBestModificationByDiffMonoMass(double mass, double max_error,
                                                               String residue, TermSpecificity term_spec) nogil except +

        void getAllSearchModifications(libcpp_vector[ String ] & modifications) nogil except +

## wrap static methods
cdef extern from "<OpenMS/CHEMISTRY/CrossLinksDB.h>" namespace "OpenMS::CrossLinksDB":
    
    CrossLinksDB* getInstance() nogil except + # wrap-ignore

