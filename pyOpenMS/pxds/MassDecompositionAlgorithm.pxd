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
        MassDecompositionAlgorithm() nogil except +
        MassDecompositionAlgorithm(MassDecompositionAlgorithm) nogil except + #wrap-ignore
        void getDecompositions(libcpp_vector[ MassDecomposition ] & decomps, DoubleReal weight) nogil except +

