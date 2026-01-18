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
        # wrap-doc:
        #  Database of cross-linker modifications. This is a singleton class.
        #  Cross-linker modifications can be loaded from OBO files.

        # private
        CrossLinksDB() except + nogil  #wrap-ignore
        # private
        CrossLinksDB(CrossLinksDB) except + nogil  #wrap-ignore

        Size getNumberOfModifications() except + nogil  # wrap-doc:Returns the number of modifications stored in the database

        void searchModifications(libcpp_set[ const ResidueModification * ] & mods,
                                 const String& mod_name,
                                 const String& residue,
                                 TermSpecificity term_spec) except + nogil  # wrap-doc:Collects all modifications which have the given name as synonym

        const ResidueModification * getModification(Size index) except + nogil  # wrap-doc:Returns the modification with the given index

        const ResidueModification * getModification(const String & mod_name) except + nogil  # wrap-doc:Returns the modification with the given name

        const ResidueModification * getModification(const String & mod_name,
                                            const String & residue,
                                            TermSpecificity term_spec) except + nogil  # wrap-doc:Returns the modification with the given name, residue, and term specificity

        bool has(String modification) except + nogil  # wrap-doc:Returns True if the modification exists in the database

        #void addModification(libcpp_unique_ptr[ResidueModification] new_mod) except + nogil

        Size findModificationIndex(const String & mod_name) except + nogil  # wrap-doc:Returns the index of the modification with the given name

        void searchModificationsByDiffMonoMass(libcpp_vector[ String ] & mods, double mass, double max_error,
                                               const String & residue, TermSpecificity term_spec) except + nogil  # wrap-doc:Collects all modifications with delta mass inside a tolerance window

        const ResidueModification* getBestModificationByDiffMonoMass(double mass, double max_error,
                                                              const String residue, TermSpecificity term_spec) except + nogil  # wrap-doc:Returns the best matching modification for the given delta mass and residue
        void getAllSearchModifications(libcpp_vector[ String ] & modifications) except + nogil  # wrap-doc:Collects all modifications that can be used for identification searches

        void readFromOBOFile(const String & filename) except + nogil  # wrap-doc:Adds modifications from a given file in OBO format

        bool isInstantiated() except + nogil  # wrap-doc:Returns True if the database has been instantiated

## wrap static methods
cdef extern from "<OpenMS/CHEMISTRY/CrossLinksDB.h>" namespace "OpenMS::CrossLinksDB":

    CrossLinksDB* getInstance() except + nogil  # wrap-ignore

