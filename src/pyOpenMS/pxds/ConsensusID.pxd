from libcpp.vector cimport vector as libcpp_vector
from PeptideIdentification cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/ConsensusID.h>" namespace "OpenMS":

    cdef cppclass ConsensusID(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler

        ConsensusID()                  nogil except +
        ConsensusID(ConsensusID)   nogil except + #wrap-ignore

        void apply(libcpp_vector[PeptideIdentification] & ids) nogil except + 

