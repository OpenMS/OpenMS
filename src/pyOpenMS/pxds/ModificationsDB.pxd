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

        void searchModifications(libcpp_set[ const ResidueModification * ] & mods,
                                 String mod_name, String residue,
                                 Term_Specificity term_spec) nogil except +

        ResidueModification getModification(String & mod_name, String & residue, TermSpecificity term_spec) nogil except +

        bool has(String modification) nogil except +

        Size findModificationIndex(String & mod_name) nogil except +

        void searchModificationsByDiffMonoMass(libcpp_vector[ String ] & mods, double mass, double max_error,
                                               String & residue, TermSpecificity term_spec) nogil except +

        ResidueModification* getBestModificationByMonoMass(double mass, double max_error,
                                                           String residue,
                                                           TermSpecificity term_spec) nogil except +

        ResidueModification* getBestModificationByDiffMonoMass(double mass, double max_error,
                                                               String residue, TermSpecificity term_spec) nogil except +

        void readFromOBOFile(String & filename) nogil except +
        void readFromUnimodXMLFile(String & filename) nogil except +
        void getAllSearchModifications(libcpp_vector[ String ] & modifications) nogil except +

## wrap static methods
cdef extern from "<OpenMS/CHEMISTRY/ModificationsDB.h>" namespace "OpenMS::ModificationsDB":
    
    ModificationsDB* getInstance() nogil except + # wrap-ignore

