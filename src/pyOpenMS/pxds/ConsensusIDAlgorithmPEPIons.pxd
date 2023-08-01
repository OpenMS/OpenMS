from Types cimport *
from ConsensusIDAlgorithmSimilarity cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmPEPIons.h>" namespace "OpenMS":
    
    cdef cppclass ConsensusIDAlgorithmPEPIons(ConsensusIDAlgorithmSimilarity) :
        # wrap-inherits:
        #  ConsensusIDAlgorithmSimilarity
        ConsensusIDAlgorithmPEPIons() except + nogil 
        # private
        ConsensusIDAlgorithmPEPIons(ConsensusIDAlgorithmPEPIons &) except + nogil  #wrap-ignore

