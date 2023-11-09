from Types cimport *
from ConsensusIDAlgorithmIdentity cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmWorst.h>" namespace "OpenMS":
    
    cdef cppclass ConsensusIDAlgorithmWorst(ConsensusIDAlgorithmIdentity) :
        # wrap-inherits:
        #  ConsensusIDAlgorithmIdentity
        ConsensusIDAlgorithmWorst() except + nogil 

        # private
        ConsensusIDAlgorithmWorst(ConsensusIDAlgorithmWorst &) except + nogil  # wrap-ignore

