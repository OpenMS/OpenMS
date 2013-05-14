from Types cimport *
from libcpp.set cimport set as libcpp_set
from libcpp.vector cimport vector as libcpp_vector
from Map cimport *
from ResidueModification cimport *

cdef extern from "<OpenMS/CHEMISTRY/ModificationsDB.h>" namespace "OpenMS":
    
    cdef cppclass ModificationsDB "OpenMS::ModificationsDB":
        ModificationsDB(ModificationsDB) nogil except + #wrap-ignore
        # POINTER # ModificationsDB * getInstance() nogil except +
        Size getNumberOfModifications() nogil except +
        ResidueModification  getModification(Size index) nogil except +
        # NAMESPACE # # POINTER # void searchTerminalModifications(libcpp_set[ ResidueModification * ] & mods, String & name, ResidueModification::Term_Specificity term_spec) nogil except +
        # NAMESPACE # # POINTER # void searchModifications(libcpp_set[ ResidueModification * ] & mods, String & orgin, String & mod_name, ResidueModification::Term_Specificity term_spec) nogil except +
        # NAMESPACE # # POINTER # void searchModifications(libcpp_set[ ResidueModification * ] & mods, String & mod_name, ResidueModification::Term_Specificity term_spec) nogil except +
        ResidueModification  getTerminalModification(String & name, Term_Specificity term_spec) nogil except +
        ResidueModification  getModification(String & residue_name, String & mod_name, Term_Specificity term_spec) nogil except +
        ResidueModification  getModification(String & modification) nogil except +
        Size findModificationIndex(String & mod_name) nogil except +
        void getTerminalModificationsByDiffMonoMass(libcpp_vector[ String ] & mods, DoubleReal mass, DoubleReal error, Term_Specificity term_spec) nogil except +
        void getModificationsByDiffMonoMass(libcpp_vector[ String ] & mods, DoubleReal mass, DoubleReal error) nogil except +
        void getModificationsByDiffMonoMass(libcpp_vector[ String ] & mods, String & residue, DoubleReal mass, DoubleReal error) nogil except +
        void readFromOBOFile(String & filename) nogil except +
        void readFromUnimodXMLFile(String & filename) nogil except +
        void getAllSearchModifications(libcpp_vector[ String ] & modifications) nogil except +

