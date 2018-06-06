from libcpp.vector cimport vector as libcpp_vector
from String cimport *
from StringList cimport *
from ModificationDefinition cimport *
from PeptideIdentification cimport *
from AASequence cimport *

cdef extern from "<OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>" namespace "OpenMS":

    cdef cppclass ModificationDefinitionsSet:

        ModificationDefinitionsSet() nogil except +
        ModificationDefinitionsSet(ModificationDefinitionsSet rhs) nogil except +

        # detailed constructor with StringLists
        # The StringLists should contain UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'
        ModificationDefinitionsSet(StringList fixed_modifications, StringList variable_modifications) nogil except +

        void setMaxModifications(Size max_mod) nogil except +
        Size getMaxModifications() nogil except +
        Size getNumberOfModifications() nogil except +
        Size getNumberOfFixedModifications() nogil except +
        Size getNumberOfVariableModifications() nogil except +
        void addModification(ModificationDefinition &mod_def) nogil except +
        void setModifications(libcpp_set[ ModificationDefinition ] &mod_defs) nogil except +
        void setModifications(const String &fixed_modifications, String &variable_modifications) nogil except +
        void setModifications(StringList &fixed_modifications, StringList &variable_modifications) nogil except +
        libcpp_set[ ModificationDefinition ] getModifications() nogil except +
        libcpp_set[ ModificationDefinition ]  getFixedModifications() nogil except +
        libcpp_set[ ModificationDefinition ]  getVariableModifications() nogil except +
        libcpp_set[ String ] getModificationNames() nogil except +
        void getModificationNames(StringList &fixed_modifications, StringList &variable_modifications) nogil except +
        libcpp_set[ String ] getFixedModificationNames() nogil except +
        libcpp_set[ String ] getVariableModificationNames() nogil except +
        bool isCompatible(AASequence &peptide) nogil except +
        void inferFromPeptides(libcpp_vector[ PeptideIdentification ] &peptides) nogil except +
