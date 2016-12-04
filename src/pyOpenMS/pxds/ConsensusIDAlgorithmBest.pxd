from Types cimport *
from ConsensusIDAlgorithmIdentity cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmBest.h>" namespace "OpenMS":
    
    cdef cppclass ConsensusIDAlgorithmBest(ConsensusIDAlgorithmIdentity) :
        # wrap-inherits:
        #  ConsensusIDAlgorithmIdentity
        ConsensusIDAlgorithmBest() nogil except +
        ConsensusIDAlgorithmBest(ConsensusIDAlgorithmBest) nogil except + #wrap-ignore

