from libcpp cimport vector as libcpp_vector
from MSExperiment cimport *

cdef extern from "<OpenMS/IONMOBILITY/MSRunIMSplitter.h>" namespace "OpenMS":

    cdef cppclass MSRunIMSplitter:

        MSRunIMSplitter() nogil except +
        MSRunIMSplitter(MSRunIMSplitter) nogil except + # wrap-ignore


# COMMENT: wrap static methods
cdef extern from "<OpenMS/IONMOBILITY/MSRunIMSplitter.h>" namespace "OpenMS::MSRunIMSplitter":
        
        # static members
        libcpp_vector[ MSExperiment ] splitByFAIMSCV(MSExperiment exp) nogil except +  # wrap-attach:MSRunIMSplitter



