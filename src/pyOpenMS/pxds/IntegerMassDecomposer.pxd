from Types cimport *
from libcpp cimport bool
from Weights cimport *
from MassDecomposer cimport *

cdef extern from "<OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IntegerMassDecomposer.h>" namespace "OpenMS::ims":
    
    cdef cppclass IntegerMassDecomposer[ValueType,DecompositionValueType]:
        # wrap-ignore
        # no-pxd-import

        # We could get this to work but only with int,int and not with long
        # unsigned int,unsigned int and then we may run into all kinds of
        # problems.

        # wrap-instances:
        #  IntegerMassDecomposer := IntegerMassDecomposer[int, int]
        IntegerMassDecomposer(IMSWeights & alphabet) except + nogil 
            # wrap-doc:
            #  Constructor with weights
            #  
            #  
            #  :param alphabet: Weights over which masses to be decomposed

        IntegerMassDecomposer(IntegerMassDecomposer &) except + nogil 

        bool exist(ValueType mass) except + nogil 
            # wrap-doc:
            #  Returns true if decomposition over the 'mass' exists, otherwise - false
            #  
            #  
            #  :param mass: Mass to be decomposed
            #  :return: true if decomposition over a given mass exists, otherwise - false

        # Works in theory but may produce duplicate 
        # libcpp_vector[libcpp_vector[int]] getDecomposition(ValueType mass) except + nogil 
        # DecompositionValueType getNumberOfDecompositions(ValueType mass) except + nogil 

        # Does not work
        # libcpp_vector[libcpp_vector[int]] getAllDecompositions(ValueType mass) except + nogil 


