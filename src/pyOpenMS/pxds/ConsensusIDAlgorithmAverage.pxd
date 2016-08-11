from Types cimport *
from ConsensusIDAlgorithmIdentity cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmAverage.h>" namespace "OpenMS":
    
    cdef cppclass ConsensusIDAlgorithmAverage(ConsensusIDAlgorithmIdentity) :
        # wrap-inherits:
        #  ConsensusIDAlgorithmIdentity
        ConsensusIDAlgorithmAverage() nogil except +
        ConsensusIDAlgorithmAverage(ConsensusIDAlgorithmAverage) nogil except + #wrap-ignore

        ## void apply(libcpp_vector[ PeptideIdentification ] & ids, Size number_of_runs) nogil except +

