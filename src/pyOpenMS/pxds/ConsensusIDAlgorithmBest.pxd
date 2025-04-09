from Types cimport *
from ConsensusIDAlgorithmIdentity cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmBest.h>" namespace "OpenMS":
    
    cdef cppclass ConsensusIDAlgorithmBest(ConsensusIDAlgorithmIdentity) :
        # wrap-inherits:
        #  ConsensusIDAlgorithmIdentity

        ConsensusIDAlgorithmBest() except + nogil 
        # private
        ConsensusIDAlgorithmBest(ConsensusIDAlgorithmBest) except + nogil  # wrap-ignore

