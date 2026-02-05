from Types cimport *
from ConsensusIDAlgorithm cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmIdentity.h>" namespace "OpenMS":
    
    cdef cppclass ConsensusIDAlgorithmIdentity(ConsensusIDAlgorithm) :
        # wrap-doc:
        #  Abstract base class for ConsensusID algorithms that compare only
        #  identical sequences
        # wrap-inherits:
        #  ConsensusIDAlgorithm
        # wrap-ignore
        # ABSTRACT class
        # no-pxd-import

        # protected
        ConsensusIDAlgorithmIdentity() except + nogil  #wrap-ignore
        # private
        ConsensusIDAlgorithmIdentity(ConsensusIDAlgorithmIdentity) except + nogil  #wrap-ignore

