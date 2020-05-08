from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
from DefaultParamHandler cimport *
from PeptideIdentification cimport *
from String cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/ConsensusIDAlgorithm.h>" namespace "OpenMS":
    
    cdef cppclass ConsensusIDAlgorithm(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        # wrap-ignore
        # ABSTRACT class
        # no-pxd-import
        ConsensusIDAlgorithm() nogil except + #wrap-ignore
        ConsensusIDAlgorithm(ConsensusIDAlgorithm) nogil except + #wrap-ignore
        void apply(libcpp_vector[ PeptideIdentification ] & ids, libcpp_map[String, String] & se_info, Size number_of_runs) nogil except +

