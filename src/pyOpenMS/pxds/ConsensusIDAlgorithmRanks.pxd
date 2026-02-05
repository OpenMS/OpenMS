from Types cimport *
from ConsensusIDAlgorithmIdentity cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmRanks.h>" namespace "OpenMS":
    
    cdef cppclass ConsensusIDAlgorithmRanks(ConsensusIDAlgorithmIdentity) :
        # wrap-doc:
        #  Calculates a consensus from multiple ID runs based on the ranks of
        #  the search hits
        # wrap-inherits:
        #  ConsensusIDAlgorithmIdentity
        ConsensusIDAlgorithmRanks() except + nogil 
        # private
        ConsensusIDAlgorithmRanks(ConsensusIDAlgorithmRanks) except + nogil  # wrap-ignore

