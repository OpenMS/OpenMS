from Types cimport *
from ConsensusMap cimport *
from Peak2D cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ProteinInference.h>" namespace "OpenMS":
    
    cdef cppclass ProteinInference "OpenMS::ProteinInference":
        ProteinInference() nogil except +
        ProteinInference(ProteinInference) nogil except +
        void infer(ConsensusMap & consensus_map, UInt reference_map) nogil except +

