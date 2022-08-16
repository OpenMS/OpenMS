from Types cimport *
from String cimport *
from ConsensusMap cimport *

cdef extern from "<OpenMS/FORMAT/GNPSMetaValueFile.h>" namespace "OpenMS":

    cdef cppclass GNPSMetaValueFile:    

        GNPSMetaValueFile() nogil except +
        GNPSMetaValueFile(GNPSMetaValueFile &) nogil except +

        void store(const ConsensusMap& consensus_map, const String& output_file) nogil except +
        # wrap-doc:
        #  Write meta value table (tsv file) from a list of mzML files. Required for GNPS FBMN.
        #  
        #  This will produce the minimal required meta values and can be extended manually.
        #  
        #  Parameters
        #  ----------
        #  consensus_map: ConsensusMap
        #    Input ConsensusMap from which the input mzML files will be determined.
        #  
        #  output_file : str
        #    Output file path for the meta value table.  
        
