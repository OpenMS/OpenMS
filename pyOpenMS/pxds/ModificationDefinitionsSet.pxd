from libcpp.vector cimport vector as libcpp_vector
from String cimport *
from StringList cimport *

cdef extern from "<OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>" namespace "OpenMS":
    # While this class has many functions, only the constructors are currently used
    cdef cppclass ModificationDefinitionsSet:

        ModificationDefinitionsSet() nogil except +
        ModificationDefinitionsSet(ModificationDefinitionsSet rhs) nogil except +

        # /// detailed constructor with comma separated list of modifications
        ModificationDefinitionsSet(String fixed_modifications, String variable_modifications) nogil except +

        # /// detailed constructor with StringLists
        # The StringLists should contain UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'
        ModificationDefinitionsSet(StringList fixed_modifications, StringList variable_modifications) nogil except +

