from libcpp cimport vector as libcpp_vector
from ConsensusMap cimport *
from String cimport *
from StringList cimport *

cdef extern from "<OpenMS/METADATA/MSQuantifications.h>" namespace "OpenMS":

    cdef cppclass MSQuantifications:
        MSQuantifications()  nogil
        MSQuantifications(MSQuantifications &) nogil
        bool operator==(MSQuantifications &) nogil
        bool operator!=(MSQuantifications &) nogil

        libcpp_vector[ConsensusMap] getConsensusMaps() nogil
        void setConsensusMaps(libcpp_vector[ConsensusMap]) nogil
