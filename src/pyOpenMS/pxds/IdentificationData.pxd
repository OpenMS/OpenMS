from Types cimport *
from libcpp cimport bool
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector
from MetaInfoInterface cimport *

cdef extern from "<OpenMS/METADATA/ID/IdentificationData.h>" namespace "OpenMS":
    
    cdef cppclass IdentificationData(MetaInfoInterface) :
        # wrap-inherits:
        #  MetaInfoInterface

        IdentificationData() nogil except +
