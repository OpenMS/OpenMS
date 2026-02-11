from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool

from DefaultParamHandler cimport *
from PeptideIdentification cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/IDDecoyProbability.h>" namespace "OpenMS":

    cdef cppclass IDDecoyProbability(DefaultParamHandler):
        # wrap-inherits:
        #   DefaultParamHandler

        IDDecoyProbability() except + nogil  # wrap-doc:IDDecoyProbability calculates probabilities using decoy approach
        IDDecoyProbability(IDDecoyProbability) except + nogil 

        void apply(PeptideIdentificationList & prob_ids, PeptideIdentificationList & fwd_ids, PeptideIdentificationList & rev_ids)  except + nogil 
            # wrap-doc:
            #  Converts the forward and reverse identification into probabilities
            #  
            #  
            #  :param prob_ids: Output of the algorithm which includes identifications with probability based scores
            #  :param fwd_ids: Input parameter which represents the identifications of the forward search
            #  :param rev_ids: Input parameter which represents the identifications of the reversed search

        void apply(PeptideIdentificationList & ids)  except + nogil 
