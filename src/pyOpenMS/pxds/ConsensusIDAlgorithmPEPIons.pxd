from Types cimport *
from ConsensusIDAlgorithmSimilarity cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmPEPIons.h>" namespace "OpenMS":
    
    cdef cppclass ConsensusIDAlgorithmPEPIons(ConsensusIDAlgorithmSimilarity) :
        # wrap-inherits:
        #  ConsensusIDAlgorithmSimilarity
        ConsensusIDAlgorithmPEPIons() nogil except +
        # private
        ConsensusIDAlgorithmPEPIons(ConsensusIDAlgorithmPEPIons &) nogil except + #wrap-ignore

