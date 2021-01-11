from Types cimport *
from ConsensusIDAlgorithmSimilarity cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmPEPMatrix.h>" namespace "OpenMS":
    
    cdef cppclass ConsensusIDAlgorithmPEPMatrix(ConsensusIDAlgorithmSimilarity) :
        # wrap-inherits:
        #  ConsensusIDAlgorithmSimilarity
        ConsensusIDAlgorithmPEPMatrix() nogil except +
        ConsensusIDAlgorithmPEPMatrix(ConsensusIDAlgorithmPEPMatrix) nogil except + #wrap-ignore

