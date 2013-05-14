from ConsensusMap cimport *
from String cimport *
from PeakFileOptions cimport *

cdef extern from "<OpenMS/FORMAT/ConsensusXMLFile.h>" namespace "OpenMS":

    cdef cppclass ConsensusXMLFile:
        ConsensusXMLFile() nogil except +

        void load(String, ConsensusMap &) nogil except+
        void store(String, ConsensusMap &) nogil except+

        PeakFileOptions getOptions() nogil except +
