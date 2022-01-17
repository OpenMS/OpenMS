from Types cimport *
from ConsensusIDAlgorithmIdentity cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmRanks.h>" namespace "OpenMS":
    
    cdef cppclass ConsensusIDAlgorithmRanks(ConsensusIDAlgorithmIdentity) :
        # wrap-inherits:
        #  ConsensusIDAlgorithmIdentity
        ConsensusIDAlgorithmRanks() nogil except +
        # private
        ConsensusIDAlgorithmRanks(ConsensusIDAlgorithmRanks) nogil except + # wrap-ignore

