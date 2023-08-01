from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from DefaultParamHandler cimport *
from PeptideIdentification cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/ConsensusIDAlgorithm.h>" namespace "OpenMS":
    
    cdef cppclass ConsensusIDAlgorithm(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        # wrap-ignore
        # ABSTRACT class
        # no-pxd-import

        # protected
        ConsensusIDAlgorithm() except + nogil  # wrap-ignore
        # private
        ConsensusIDAlgorithm(ConsensusIDAlgorithm) except + nogil  # wrap-ignore

        void apply(libcpp_vector[ PeptideIdentification ] & ids, Size number_of_runs) except + nogil  # wrap-doc:Calculates the consensus ID for a set of peptide identifications of one spectrum or (consensus) feature
