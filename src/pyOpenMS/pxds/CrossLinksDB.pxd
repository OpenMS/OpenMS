from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
from String cimport *
from ResidueModification cimport *
#from libcpp.memory cimport unique_ptr as libcpp_unique_ptr

# see ../addons/CrossLinksDB.pyx
cdef extern from "<OpenMS/CHEMISTRY/CrossLinksDB.h>" namespace "OpenMS":
    
    cdef cppclass CrossLinksDB:
        # wrap-manual-memory:
        #    cdef AutowrapPtrHolder[_CrossLinksDB] inst

        # private
        CrossLinksDB() except + nogil  #wrap-ignore
        # private
        CrossLinksDB(CrossLinksDB) except + nogil  #wrap-ignore

        Size getNumberOfModifications() except + nogil 

        void searchModifications(libcpp_set[ const ResidueModification * ] & mods,
                                 const String& mod_name,
                                 const String& residue,
                                 TermSpecificity term_spec) except + nogil 

        const ResidueModification * getModification(Size index) except + nogil 

        const ResidueModification * getModification(const String & mod_name) except + nogil 

        const ResidueModification * getModification(const String & mod_name,
                                            const String & residue,
                                            TermSpecificity term_spec) except + nogil 

        bool has(String modification) except + nogil 

        #void addModification(libcpp_unique_ptr[ResidueModification] new_mod) except + nogil 

        Size findModificationIndex(const String & mod_name) except + nogil 

        void searchModificationsByDiffMonoMass(libcpp_vector[ String ] & mods, double mass, double max_error,
                                               const String & residue, TermSpecificity term_spec) except + nogil 

        const ResidueModification* getBestModificationByDiffMonoMass(double mass, double max_error,
                                                              const String residue, TermSpecificity term_spec) except + nogil 
        void getAllSearchModifications(libcpp_vector[ String ] & modifications) except + nogil 

        void readFromOBOFile(const String & filename) except + nogil 

        bool isInstantiated() except + nogil 

## wrap static methods
cdef extern from "<OpenMS/CHEMISTRY/CrossLinksDB.h>" namespace "OpenMS::CrossLinksDB":
    
    CrossLinksDB* getInstance() except + nogil  # wrap-ignore

