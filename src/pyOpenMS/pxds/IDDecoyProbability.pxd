from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool

from DefaultParamHandler cimport *
from PeptideIdentification cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/IDDecoyProbability.h>" namespace "OpenMS":

    cdef cppclass IDDecoyProbability(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler

        IDDecoyProbability() nogil except + # wrap-doc:IDDecoyProbability calculates probabilities using decoy approach
        IDDecoyProbability(IDDecoyProbability) nogil except +

        void apply(libcpp_vector[PeptideIdentification] & prob_ids, libcpp_vector[PeptideIdentification] & fwd_ids, libcpp_vector[PeptideIdentification] & rev_ids)  nogil except +
            # wrap-doc:
            #   Converts the forward and reverse identification into probabilities
            #   -----
            #   :param prob_ids: Output of the algorithm which includes identifications with probability based scores
            #   :param fwd_ids: Input parameter which represents the identifications of the forward search
            #   :param rev_ids: Input parameter which represents the identifications of the reversed search

        void apply(libcpp_vector[PeptideIdentification] & ids)  nogil except +
