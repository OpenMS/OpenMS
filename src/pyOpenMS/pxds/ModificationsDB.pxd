from Types cimport *
from Map cimport *
from String cimport *
from ResidueModification cimport *

cdef extern from "<OpenMS/CHEMISTRY/ModificationsDB.h>" namespace "OpenMS":
    
    cdef cppclass ModificationsDB "OpenMS::ModificationsDB":
        # wrap-manual-memory

        ModificationsDB(ModificationsDB) nogil except + #wrap-ignore

        Size getNumberOfModifications() nogil except +
        ResidueModification getModification(Size index) nogil except +

        void searchTerminalModifications(libcpp_set[ const ResidueModification * ] & mods, 
                                         String & name, Term_Specificity term_spec) nogil except +

        void searchModifications(libcpp_set[ const ResidueModification * ] & mods,
                                 String & orgin, String & mod_name,
                                 Term_Specificity term_spec) nogil except +

        void searchModifications(libcpp_set[ const ResidueModification * ] & mods,
                                 String & mod_name, Term_Specificity term_spec) nogil except +

        ResidueModification getTerminalModification(String & name, Term_Specificity term_spec) nogil except +
        ResidueModification getModification(String & residue_name, String & mod_name, Term_Specificity term_spec) nogil except +
        ResidueModification getModification(String & modification) nogil except +

        Size findModificationIndex(String & mod_name) nogil except +
        void getTerminalModificationsByDiffMonoMass(libcpp_vector[ String ] & mods, double mass, double error, Term_Specificity term_spec) nogil except +
        void getModificationsByDiffMonoMass(libcpp_vector[ String ] & mods, double mass, double error) nogil except +
        void getModificationsByDiffMonoMass(libcpp_vector[ String ] & mods, String & residue, double mass, double error) nogil except +

        void readFromOBOFile(String & filename) nogil except +
        void readFromUnimodXMLFile(String & filename) nogil except +
        void getAllSearchModifications(libcpp_vector[ String ] & modifications) nogil except +
        bool has(String modification) nogil except +

        const ResidueModification * getBestModificationsByMonoMass(String residue, double mass, double max_error) nogil except +
        const ResidueModification * getBestModificationsByDiffMonoMass(String residue, double mass, double max_error) nogil except +

## wrap static methods
cdef extern from "<OpenMS/CHEMISTRY/ModificationsDB.h>" namespace "OpenMS::ModificationsDB":
    
    ModificationsDB* getInstance() nogil except + # wrap-ignore

