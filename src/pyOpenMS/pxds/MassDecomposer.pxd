from Types cimport *
from libcpp cimport bool

cdef extern from "<OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/MassDecomposer.h>" namespace "OpenMS::ims":
    
    cdef cppclass MassDecomposer[ValueType,DecompositionValueType]:
        # wrap-ignore
        # ABSTRACT class
        # no-pxd-import
        MassDecomposer(MassDecomposer &) except + nogil  # compiler

        bool exist(ValueType mass) except + nogil 
            # wrap-doc:
                #  Returns true if the decomposition for the given `mass` exists, otherwise - false
                #  
                #  
                #  :param mass: Mass to be checked on decomposing
                #  :return: True, if the decomposition for `mass` exist, otherwise - false

        # decomposition_type getDecomposition(ValueType mass) except + nogil 
        # decompositions_type getAllDecompositions(ValueType mass) except + nogil 
        DecompositionValueType getNumberOfDecompositions(ValueType mass) except + nogil 

