from Types cimport *
from libcpp.set cimport set as libcpp_set
from libcpp.vector cimport vector as libcpp_vector
from Map cimport *
from ResidueModification cimport *

cdef extern from "<OpenMS/CHEMISTRY/ModificationsDB.h>" namespace "OpenMS":
    
    cdef cppclass ModificationsDB "OpenMS::ModificationsDB":
        ModificationsDB(ModificationsDB) nogil except + #wrap-ignore
        # Cannot get ModificationsDB since its constructor and destructor is private ... 
        # POINTER # ModificationsDB * getInstance() nogil except +
        Size getNumberOfModifications() nogil except +
        ResidueModification getModification(Size index) nogil except +
        # NAMESPACE # # POINTER # void searchModifications(libcpp_set[ ResidueModification * ] & mods, String & mod_name, String & residue, ResidueModification::TermSpecificity term_spec) nogil except +
        ResidueModification getModification(String & mod_name, String & residue, TermSpecificity term_spec) nogil except +
        Size findModificationIndex(String & mod_name) nogil except +
        void searchModificationsByDiffMonoMass(libcpp_vector[ String ] & mods, double mass, double max_error, String & residue, TermSpecificity term_spec) nogil except +
        void readFromOBOFile(String & filename) nogil except +
        void readFromUnimodXMLFile(String & filename) nogil except +
        void getAllSearchModifications(libcpp_vector[ String ] & modifications) nogil except +
        bool has(String modification) nogil except +

        # POINTER # ResidueModification * getBestModificationByMonoMass(double mass, double max_error, String & residue, TermSpecificity term_spec) nogil except +
        # POINTER # ResidueModification * getBestModificationByDiffMonoMass(double mass, double max_error, String & residue, TermSpecificity term_spec) nogil except +

