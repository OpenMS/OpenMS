from Types cimport *
from String cimport *
from ConsensusMap cimport *

cdef extern from "<OpenMS/FORMAT/IBSpectraFile.h>" namespace "OpenMS":
    
    cdef cppclass IBSpectraFile "OpenMS::IBSpectraFile":
        IBSpectraFile() nogil except +
        IBSpectraFile(IBSpectraFile) nogil except +
        void store(const String & filename, ConsensusMap & cm) nogil except +
