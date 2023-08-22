from libcpp.vector cimport vector as libcpp_vector
from libcpp.set cimport set as libcpp_set
from String cimport *
from StringList cimport *
from ModificationDefinition cimport *
from PeptideIdentification cimport *
from AASequence cimport *

cdef extern from "<OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>" namespace "OpenMS":

  cdef cppclass ModificationDefinitionsSet:
      # wrap-doc:
        #  Representation of a set of modification definitions
        #  
        #  This class enhances the modification definitions as defined in the
        #  class ModificationDefinition into a set of definitions. This is also
        #  e.g. used as input parameters in search engines.

    ModificationDefinitionsSet() except + nogil 
    ModificationDefinitionsSet(ModificationDefinitionsSet &) except + nogil 

    # detailed constructor with StringLists
    # The StringLists should contain UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'
    ModificationDefinitionsSet(StringList fixed_modifications, StringList variable_modifications) except + nogil 

    void setMaxModifications(Size max_mod) except + nogil  # wrap-doc:Sets the maximal number of modifications allowed per peptide
    Size getMaxModifications() except + nogil  # wrap-doc:Return the maximal number of modifications allowed per peptide
    Size getNumberOfModifications() except + nogil  # wrap-doc:Returns the number of modifications stored in this set
    Size getNumberOfFixedModifications() except + nogil  # wrap-doc:Returns the number of fixed modifications stored in this set
    Size getNumberOfVariableModifications() except + nogil  # wrap-doc:Returns the number of variable modifications stored in this set
    void addModification(ModificationDefinition &mod_def) except + nogil  # wrap-doc:Adds a modification definition to the set
    void setModifications(libcpp_set[ ModificationDefinition ] &mod_defs) except + nogil  # wrap-doc:Sets the modification definitions
    void setModifications(const String &fixed_modifications, String &variable_modifications) except + nogil 
      # wrap-doc:
          #  Set the modification definitions from a string
          #  
          #  The strings should contain a comma separated list of modifications. The names
          #  can be PSI-MOD identifier or any other unique name supported by PSI-MOD. TermSpec
          #  definitions and other specific definitions are given by the modifications themselves.

    void setModifications(StringList &fixed_modifications, StringList &variable_modifications) except + nogil  # wrap-doc:Same as above, but using StringList instead of comma separated strings
    libcpp_set[ModificationDefinition] getModifications() except + nogil  # wrap-doc:Returns the stored modification definitions
    libcpp_set[ModificationDefinition] getFixedModifications() except + nogil  # wrap-doc:Returns the stored fixed modification definitions
    libcpp_set[ModificationDefinition] getVariableModifications() except + nogil  # wrap-doc:Returns the stored variable modification definitions
    void getModificationNames(StringList &fixed_modifications, StringList &variable_modifications) except + nogil  # wrap-doc:Populates the output lists with the modification names (use e.g. for ProteinIdentification::SearchParameters)
    libcpp_set[String] getFixedModificationNames() except + nogil  # wrap-doc:Returns only the names of the fixed modifications
    libcpp_set[String] getVariableModificationNames() except + nogil  # wrap-doc:Returns only the names of the variable modifications
    libcpp_set[String] getModificationNames() except + nogil  # wrap-doc:Returns only the names of the modifications stored in the set
    bool isCompatible(AASequence &peptide) except + nogil  # wrap-doc:Returns true if the peptide is compatible with the definitions, e.g. does not contain other modifications
    void inferFromPeptides(libcpp_vector[ PeptideIdentification ] &peptides) except + nogil  # wrap-doc:Infers the sets of defined modifications from the modifications present on peptide identifications
    
