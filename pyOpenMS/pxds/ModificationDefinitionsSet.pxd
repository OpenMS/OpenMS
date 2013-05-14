from libcpp.vector cimport vector as libcpp_vector
from String cimport *
from StringList cimport *
from ModificationDefinition cimport *
from AASequence cimport *

cdef extern from "<OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>" namespace "OpenMS":

    cdef cppclass ModificationDefinitionsSet:

        ModificationDefinitionsSet() nogil except +
        ModificationDefinitionsSet(ModificationDefinitionsSet rhs) nogil except +

        # detailed constructor with comma separated list of modifications
        ModificationDefinitionsSet(String fixed_modifications, String variable_modifications) nogil except +

        # detailed constructor with StringLists
        # The StringLists should contain UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'
        ModificationDefinitionsSet(StringList fixed_modifications, StringList variable_modifications) nogil except +

        void setMaxModifications(Size max_mod)
        Size getMaxModifications()
        Size getNumberOfModifications()
        Size getNumberOfFixedModifications()
        Size getNumberOfVariableModifications()
        void addModification(ModificationDefinition &mod_def)
        void setModifications(libcpp_set[ ModificationDefinition ] &mod_defs)
        void setModifications(String &fixed_modifications, String &variable_modifications)
        void setModifications(StringList &fixed_modifications, StringList &variable_modifications)
        libcpp_set[ ModificationDefinition ] getModifications()
        libcpp_set[ ModificationDefinition ]  getFixedModifications()
        libcpp_set[ ModificationDefinition ]  getVariableModifications()
        libcpp_set[ String ] getModificationNames()
        libcpp_set[ String ] getFixedModificationNames()
        libcpp_set[ String ] getVariableModificationNames()
        bool isCompatible(AASequence &peptide)
