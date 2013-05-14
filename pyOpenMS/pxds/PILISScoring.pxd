from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from PeptideIdentification cimport *
from DefaultParamHandler cimport *
from Map cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/PILISScoring.h>" namespace "OpenMS":
    
    cdef cppclass PILISScoring(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        PILISScoring() nogil except +
        PILISScoring(PILISScoring) nogil except +
        void getScores(libcpp_vector[ PeptideIdentification ] & ids) nogil except +
        void getScore(PeptideIdentification & id_) nogil except +

