from Types cimport *
from ConsensusIDAlgorithmIdentity cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmAverage.h>" namespace "OpenMS":
    
    cdef cppclass ConsensusIDAlgorithmAverage(ConsensusIDAlgorithmIdentity) :
        # wrap-inherits:
        #  ConsensusIDAlgorithmIdentity
        ConsensusIDAlgorithmAverage() nogil except +
        # private
        ConsensusIDAlgorithmAverage(ConsensusIDAlgorithmAverage &) nogil except + # wrap-ignore

