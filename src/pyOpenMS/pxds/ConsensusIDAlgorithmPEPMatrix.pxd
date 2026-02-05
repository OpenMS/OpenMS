from Types cimport *
from ConsensusIDAlgorithmSimilarity cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmPEPMatrix.h>" namespace "OpenMS":
    
    cdef cppclass ConsensusIDAlgorithmPEPMatrix(ConsensusIDAlgorithmSimilarity) :
        # wrap-doc:
        #  Calculates a consensus from multiple ID runs based on PEPs and
        #  sequence similarities
        # wrap-inherits:
        #  ConsensusIDAlgorithmSimilarity
        ConsensusIDAlgorithmPEPMatrix() except + nogil 
        # private
        ConsensusIDAlgorithmPEPMatrix(ConsensusIDAlgorithmPEPMatrix) except + nogil  #wrap-ignore

