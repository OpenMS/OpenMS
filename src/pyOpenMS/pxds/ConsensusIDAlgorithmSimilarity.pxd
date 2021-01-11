from Types cimport *
from ConsensusIDAlgorithm cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmSimilarity.h>" namespace "OpenMS":
    
    cdef cppclass ConsensusIDAlgorithmSimilarity(ConsensusIDAlgorithm) :
        # wrap-inherits:
        #  ConsensusIDAlgorithm
        # wrap-ignore
        # ABSTRACT class
        # no-pxd-import
        ConsensusIDAlgorithmSimilarity() nogil except + #wrap-ignore
        ConsensusIDAlgorithmSimilarity(ConsensusIDAlgorithmSimilarity) nogil except + #wrap-ignore

