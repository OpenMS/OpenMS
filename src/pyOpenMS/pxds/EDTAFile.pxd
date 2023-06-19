from libcpp.vector cimport vector as libcpp_vector
from String cimport *

from FeatureMap cimport *
from ConsensusMap cimport *
from Feature cimport *

from libcpp cimport bool

cdef extern from "<OpenMS/FORMAT/EDTAFile.h>" namespace "OpenMS":

    cdef cppclass EDTAFile:

        EDTAFile() nogil except +
        EDTAFile(EDTAFile &) nogil except + # compiler

        void store(String filename, FeatureMap & map) nogil except +
        void store(String filename, ConsensusMap & map)  nogil except +
        void load(String filename, ConsensusMap & consensus_map) nogil except +

