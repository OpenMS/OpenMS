from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool

from DefaultParamHandler cimport *
from PeptideIdentification cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/IDDecoyProbability.h>" namespace "OpenMS":

    cdef cppclass IDDecoyProbability(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler

        IDDecoyProbability()                    nogil except +
        IDDecoyProbability(IDDecoyProbability)  nogil except +

        void apply(libcpp_vector[PeptideIdentification] & prob_ids, libcpp_vector[PeptideIdentification] & fwd_ids, libcpp_vector[PeptideIdentification] & rev_ids)  nogil except +
        void apply(libcpp_vector[PeptideIdentification] & ids)  nogil except +

