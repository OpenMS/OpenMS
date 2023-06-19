from ConsensusMap cimport *
from String cimport *
from PeakFileOptions cimport *

cdef extern from "<OpenMS/FORMAT/ConsensusXMLFile.h>" namespace "OpenMS":

    cdef cppclass ConsensusXMLFile:
        ConsensusXMLFile() nogil except +

        void load(String, ConsensusMap &) nogil except + # wrap-doc:Loads a consensus map from file and calls updateRanges
        void store(String, ConsensusMap &) nogil except + # wrap-doc:Stores a consensus map to file

        PeakFileOptions getOptions() nogil except + # wrap-doc:Mutable access to the options for loading/storing
