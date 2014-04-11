from Types cimport *
from libcpp cimport bool
from Weights cimport *
from MassDecomposer cimport *

cdef extern from "<OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IntegerMassDecomposer.h>" namespace "OpenMS::ims":
    
    cdef cppclass IntegerMassDecomposer[ValueType,DecompositionValueType]:
        # wrap-ignore

        # We could get this to work but only with int,int and not with long
        # unsigned int,unsigned int and then we may run into all kinds of
        # problems.

        # wrap-instances:
        #   IntegerMassDecomposer := IntegerMassDecomposer[int, int]
        IntegerMassDecomposer() nogil except + # wrap-ignore
        IntegerMassDecomposer(IntegerMassDecomposer) nogil except + #wrap-ignore
        IntegerMassDecomposer(IMSWeights & alphabet) nogil except +
        bool exist(ValueType mass) nogil except +

        # Works in theory but may produce duplicate 
        # libcpp_vector[libcpp_vector[int]] getDecomposition(ValueType mass) nogil except +
        # DecompositionValueType getNumberOfDecompositions(ValueType mass) nogil except +

        # Does not work
        # libcpp_vector[libcpp_vector[int]] getAllDecompositions(ValueType mass) nogil except +


