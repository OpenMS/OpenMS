from Types cimport *
from ConsensusIDAlgorithmIdentity cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmWorst.h>" namespace "OpenMS":
    
    cdef cppclass ConsensusIDAlgorithmWorst(ConsensusIDAlgorithmIdentity) :
        # wrap-inherits:
        #  ConsensusIDAlgorithmIdentity
        ConsensusIDAlgorithmWorst() nogil except +

        # private
        ConsensusIDAlgorithmWorst(ConsensusIDAlgorithmWorst &) nogil except + # wrap-ignore

