from libcpp.vector cimport vector as libcpp_vector
from String cimport *

from FeatureMap cimport *
from ConsensusMap cimport *
from Feature cimport *

from libcpp cimport bool

cdef extern from "<OpenMS/FORMAT/EDTAFile.h>" namespace "OpenMS":

    cdef cppclass EDTAFile:

        EDTAFile() except + nogil 
        EDTAFile(EDTAFile &) except + nogil  # compiler

        void store(String filename, FeatureMap & map) except + nogil 
        void store(String filename, ConsensusMap & map)  except + nogil 
        void load(String filename, ConsensusMap & consensus_map) except + nogil 

