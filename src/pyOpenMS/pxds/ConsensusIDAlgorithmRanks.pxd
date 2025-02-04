from Types cimport *
from ConsensusIDAlgorithmIdentity cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmRanks.h>" namespace "OpenMS":
    
    cdef cppclass ConsensusIDAlgorithmRanks(ConsensusIDAlgorithmIdentity) :
        # wrap-inherits:
        #  ConsensusIDAlgorithmIdentity
        ConsensusIDAlgorithmRanks() except + nogil 
        # private
        ConsensusIDAlgorithmRanks(ConsensusIDAlgorithmRanks) except + nogil  # wrap-ignore

