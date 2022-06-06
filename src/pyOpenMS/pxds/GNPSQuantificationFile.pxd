from Types cimport *
from String cimport *
from ConsensusMap cimport *

cdef extern from "<OpenMS/FORMAT/GNPSQuantificationFile.h>" namespace "OpenMS":

    cdef cppclass GNPSQuantificationFile:    

        GNPSQuantificationFile() nogil except +
        GNPSQuantificationFile(GNPSQuantificationFile &) nogil except +

        void store(const ConsensusMap& consensus_map, const String& output_file) nogil except +
        # wrap-doc:
        #  Write feature quantification table (txt file) from a ConsensusMap. Required for GNPS FBMN.
        #  
        #  The table contains map information on the featureXML files from which the ConsensusMap was generated as well as
        #  a row for every consensus feature with information on rt, mz, intensity, width and quality. The same information is
        #  added for each original feature in the consensus feature.
        #  
        #  Parameters
        #  ----------
        #  consensus_map : ConsensusMap
        #    Input ConsensusMap annotated with IonIdentityMolecularNetworking.annotateConsensusMap.
        #  
        #  output_file : str
        #    Output file path for the feature quantification table.  
        
