from Types cimport *
from Map cimport *
from String cimport *
from ResidueModification cimport *
#from libcpp.memory cimport unique_ptr as libcpp_unique_ptr

# see ../addons/ModificationsDB.pyx
cdef extern from "<OpenMS/CHEMISTRY/ModificationsDB.h>" namespace "OpenMS":
    
    cdef cppclass ModificationsDB "OpenMS::ModificationsDB":
        # wrap-manual-memory:
        #   cdef AutowrapPtrHolder[_ModificationsDB] inst

        ModificationsDB(ModificationsDB) nogil except + #wrap-ignore

        Size getNumberOfModifications() nogil except +

        void searchModifications(libcpp_set[ const ResidueModification * ] & mods,
                                 const String& mod_name,
                                 const String& residue,
                                 TermSpecificity term_spec) nogil except +

        const ResidueModification * getModification(Size index) nogil except +

        const ResidueModification * getModification(const String & mod_name) nogil except +

        const ResidueModification * getModification(const String & mod_name,
                                            const String & residue,
                                            TermSpecificity term_spec) nogil except +

        bool has(String modification) nogil except +

        #void addModification(libcpp_unique_ptr[ResidueModification] new_mod) nogil except +

        Size findModificationIndex(const String & mod_name) nogil except +

        void searchModificationsByDiffMonoMass(libcpp_vector[ String ] & mods, double mass, double max_error,
                                               const String & residue, TermSpecificity term_spec) nogil except +


        const ResidueModification* getBestModificationByDiffMonoMass(double mass, double max_error,
                                                                     const String& residue, TermSpecificity term_spec) nogil except +

        void getAllSearchModifications(libcpp_vector[ String ] & modifications) nogil except +

        bool isInstantiated() nogil except +

## wrap static methods
cdef extern from "<OpenMS/CHEMISTRY/ModificationsDB.h>" namespace "OpenMS::ModificationsDB":
    
    ModificationsDB* getInstance() nogil except + # wrap-ignore
    
    ModificationsDB* getInstance(String unimod_file, 
                                 String psimod_file,
                                 String xlmod_file) nogil except + # wrap-ignore

