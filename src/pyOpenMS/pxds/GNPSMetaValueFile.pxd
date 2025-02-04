from Types cimport *
from String cimport *
from ConsensusMap cimport *

cdef extern from "<OpenMS/FORMAT/GNPSMetaValueFile.h>" namespace "OpenMS":

    cdef cppclass GNPSMetaValueFile:    

        GNPSMetaValueFile() except + nogil 
        GNPSMetaValueFile(GNPSMetaValueFile &) except + nogil 

        void store(const ConsensusMap& consensus_map, const String& output_file) except + nogil 
        # wrap-doc:
        #  Write meta value table (tsv file) from a list of mzML files. Required for GNPS FBMN.
        #  
        #  This will produce the minimal required meta values and can be extended manually.
        #  
        #  :param consensus_map: Input ConsensusMap from which the input mzML files will be determined.
        #  :param output_file: Output file path for the meta value table.  
        
