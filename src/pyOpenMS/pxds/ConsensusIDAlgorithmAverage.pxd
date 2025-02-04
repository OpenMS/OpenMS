from Types cimport *
from ConsensusIDAlgorithmIdentity cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmAverage.h>" namespace "OpenMS":
    
    cdef cppclass ConsensusIDAlgorithmAverage(ConsensusIDAlgorithmIdentity) :
        # wrap-inherits:
        #  ConsensusIDAlgorithmIdentity
        ConsensusIDAlgorithmAverage() except + nogil 
        # private
        ConsensusIDAlgorithmAverage(ConsensusIDAlgorithmAverage &) except + nogil  # wrap-ignore

