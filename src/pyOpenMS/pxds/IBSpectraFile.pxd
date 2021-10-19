from Types cimport *
from String cimport *
from ConsensusMap cimport *

cdef extern from "<OpenMS/FORMAT/IBSpectraFile.h>" namespace "OpenMS":
    
    cdef cppclass IBSpectraFile "OpenMS::IBSpectraFile":

        IBSpectraFile() nogil except + # wrap-doc:Implements the export of consensusmaps into the IBSpectra format used by isobar to load quantification results
        IBSpectraFile(IBSpectraFile &) nogil except +

        void store(const String & filename, ConsensusMap & cm) nogil except +
            # wrap-doc:
            #   Writes the contents of the ConsensusMap cm into the file named by filename
            #   -----
            #   :param filename: The name of the file where the contents of cm should be stored
            #   :param cm: The ConsensusMap that should be exported to filename
            #   :raises:
            #     Exception: InvalidParameter if the ConsensusMap does not hold the result of an isobaric quantification experiment (e.g., itraq)
