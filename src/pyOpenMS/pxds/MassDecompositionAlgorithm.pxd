from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from DefaultParamHandler cimport *
from MassDecomposition cimport *
# from RealMassDecomposer cimport *
# from IMSAlphabet cimport *
# from Weights cimport *

cdef extern from "<OpenMS/CHEMISTRY/MASSDECOMPOSITION/MassDecompositionAlgorithm.h>" namespace "OpenMS":
    
    cdef cppclass MassDecompositionAlgorithm(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        MassDecompositionAlgorithm() except + nogil 
        # private
        MassDecompositionAlgorithm(MassDecompositionAlgorithm) except + nogil  #wrap-ignore
        void getDecompositions(libcpp_vector[ MassDecomposition ] & decomps, double weight) except + nogil  # wrap-doc:Returns the possible decompositions given the weight

