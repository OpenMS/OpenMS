from Types cimport *
from libcpp cimport bool

cdef extern from "<OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/MassDecomposer.h>" namespace "OpenMS::ims":
    
    cdef cppclass MassDecomposer[ValueType,DecompositionValueType]:
        # wrap-ignore
        # ABSTRACT class
        MassDecomposer(MassDecomposer) nogil except + #wrap-ignore
        bool exist(ValueType mass) nogil except +
        # decomposition_type getDecomposition(ValueType mass) nogil except +
        # decompositions_type getAllDecompositions(ValueType mass) nogil except +
        DecompositionValueType getNumberOfDecompositions(ValueType mass) nogil except +

