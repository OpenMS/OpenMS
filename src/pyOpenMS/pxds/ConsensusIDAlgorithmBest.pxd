from Types cimport *
from ConsensusIDAlgorithmIdentity cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmBest.h>" namespace "OpenMS":
    
    cdef cppclass ConsensusIDAlgorithmBest(ConsensusIDAlgorithmIdentity) :
        # wrap-doc:
        #  Calculates a consensus from multiple ID runs by taking the best
        #  search score
        # wrap-inherits:
        #  ConsensusIDAlgorithmIdentity

        ConsensusIDAlgorithmBest() except + nogil 
        # private
        ConsensusIDAlgorithmBest(ConsensusIDAlgorithmBest) except + nogil  # wrap-ignore

