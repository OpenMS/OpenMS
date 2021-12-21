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

        # private
        ModificationsDB() nogil except + # wrap-ignore
        # private
        ModificationsDB(ModificationsDB) nogil except + # wrap-ignore

        Size getNumberOfModifications() nogil except + # wrap-doc:Returns the number of modifications read from the unimod.xml file

        void searchModifications(libcpp_set[ const ResidueModification * ] & mods,
                                 const String& mod_name,
                                 const String& residue,
                                 TermSpecificity term_spec) nogil except +
            # wrap-doc:
                #   Collects all modifications which have the given name as synonym
                #   -----
                #   If `residue` is set, only modifications with matching residue of origin are considered
                #   If `term_spec` is set, only modifications with matching term specificity are considered
                #   The resulting set of modifications will be empty if no modification exists that fulfills the criteria

        const ResidueModification * getModification(Size index) nogil except + # wrap-doc:Returns the modification with the given index

        const ResidueModification * getModification(const String & mod_name) nogil except + # wrap-doc:Returns the modification with the given name

        const ResidueModification * getModification(const String & mod_name,
                                            const String & residue,
                                            TermSpecificity term_spec) nogil except + # wrap-doc:Returns the modification with the given arguments

        bool has(String modification) nogil except + # wrap-doc:Returns true if the modification exists

        # unique_ptrs do not have a conversion provider in autowrap yet. Also, we would probably make a copy to not steal memory from the python object, defeating the purpose.
        #void addModification(libcpp_unique_ptr[ResidueModification] new_mod) nogil except +

        # TODO also do a function with bool return type to save a copy?
        const ResidueModification * addModification(const ResidueModification & new_mod) nogil except + # wrap-doc:Add a new modification to ModificationsDB. If the modification already exists (based on its fullID) it is not added. Returns the modification in the ModificationDB (which can differ from input if mod was already present).

        Size findModificationIndex(const String & mod_name) nogil except + # wrap-doc:Returns the index of the modification in the mods_ vector; a unique name must be given

        void searchModificationsByDiffMonoMass(libcpp_vector[ String ] & mods, double mass, double max_error,
                                               const String & residue, TermSpecificity term_spec) nogil except + # wrap-doc:Collects all modifications with delta mass inside a tolerance window


        const ResidueModification* getBestModificationByDiffMonoMass(double mass, double max_error,
                                                                     const String& residue, TermSpecificity term_spec) nogil except +
            # wrap-doc:
                #   Returns the best matching modification for the given delta mass and residue
                #   -----
                #   Query the modifications DB to get the best matching modification with
                #   the given delta mass at the given residue (NULL pointer means no result,
                #   maybe the maximal error tolerance needs to be increased). Possible
                #   input for CAM modification would be a delta mass of 57 and a residue
                #   of "C".
                #   -----
                #   Note: If there are multiple possible matches with equal masses, it
                #   will choose the _first_ match which defaults to the first matching
                #   UniMod entry.
                #   -----
                #   :param residue: The residue at which the modifications occurs
                #   :param mass: The monoisotopic mass of the residue including the mass of the modification
                #   :param max_error: The maximal mass error in the modification search
                #   :returns: A pointer to the best matching modification (or NULL if none was found)

        void getAllSearchModifications(libcpp_vector[ String ] & modifications) nogil except + # wrap-doc:Collects all modifications that can be used for identification searches

        bool isInstantiated() nogil except + # wrap-doc:Check whether ModificationsDB was instantiated before

## wrap static methods
cdef extern from "<OpenMS/CHEMISTRY/ModificationsDB.h>" namespace "OpenMS::ModificationsDB":
    
    ModificationsDB* getInstance() nogil except + # wrap-ignore
    
    ModificationsDB* getInstance(String unimod_file, 
                                 String psimod_file,
                                 String xlmod_file) nogil except + # wrap-ignore

