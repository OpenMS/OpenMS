from Types cimport *
from ConsensusIDAlgorithmIdentity cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmWorst.h>" namespace "OpenMS":
    
    cdef cppclass ConsensusIDAlgorithmWorst(ConsensusIDAlgorithmIdentity) :
        # wrap-doc:
        #  Calculates a consensus from multiple ID runs by taking the worst
        #  search score (conservative approach)
        # wrap-inherits:
        #  ConsensusIDAlgorithmIdentity
        ConsensusIDAlgorithmWorst() except + nogil 

        # private
        ConsensusIDAlgorithmWorst(ConsensusIDAlgorithmWorst &) except + nogil  # wrap-ignore

