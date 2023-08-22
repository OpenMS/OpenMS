from ConsensusMap cimport *
from String cimport *
from PeakFileOptions cimport *

cdef extern from "<OpenMS/FORMAT/ConsensusXMLFile.h>" namespace "OpenMS":

    cdef cppclass ConsensusXMLFile:
        ConsensusXMLFile() except + nogil 

        void load(String, ConsensusMap &) except + nogil  # wrap-doc:Loads a consensus map from file and calls updateRanges
        void store(String, ConsensusMap &) except + nogil  # wrap-doc:Stores a consensus map to file

        PeakFileOptions getOptions() except + nogil  # wrap-doc:Mutable access to the options for loading/storing
