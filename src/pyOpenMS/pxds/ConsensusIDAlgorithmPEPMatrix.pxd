from Types cimport *
from ConsensusIDAlgorithmSimilarity cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmPEPMatrix.h>" namespace "OpenMS":
    
    cdef cppclass ConsensusIDAlgorithmPEPMatrix(ConsensusIDAlgorithmSimilarity) :
        # wrap-inherits:
        #  ConsensusIDAlgorithmSimilarity
        ConsensusIDAlgorithmPEPMatrix() except + nogil 
        # private
        ConsensusIDAlgorithmPEPMatrix(ConsensusIDAlgorithmPEPMatrix) except + nogil  #wrap-ignore

