from Types cimport *
from ConsensusIDAlgorithm cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmIdentity.h>" namespace "OpenMS":
    
    cdef cppclass ConsensusIDAlgorithmIdentity(ConsensusIDAlgorithm) :
        # wrap-inherits:
        #  ConsensusIDAlgorithm
        # wrap-ignore
        # ABSTRACT class
        # no-pxd-import
        ConsensusIDAlgorithmIdentity() nogil except + #wrap-ignore
        ConsensusIDAlgorithmIdentity(ConsensusIDAlgorithmIdentity) nogil except + #wrap-ignore

