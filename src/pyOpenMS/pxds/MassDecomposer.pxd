from Types cimport *
from libcpp cimport bool

cdef extern from "<OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/MassDecomposer.h>" namespace "OpenMS::ims":
    
    cdef cppclass MassDecomposer[ValueType,DecompositionValueType]:
        # wrap-ignore
        # ABSTRACT class
        # no-pxd-import
        MassDecomposer(MassDecomposer &) nogil except + # compiler

        bool exist(ValueType mass) nogil except +
            # wrap-doc:
                #   Returns true if the decomposition for the given `mass` exists, otherwise - false
                #   -----
                #   :param mass: Mass to be checked on decomposing
                #   :returns: True, if the decomposition for `mass` exist, otherwise - false

        # decomposition_type getDecomposition(ValueType mass) nogil except +
        # decompositions_type getAllDecompositions(ValueType mass) nogil except +
        DecompositionValueType getNumberOfDecompositions(ValueType mass) nogil except +

