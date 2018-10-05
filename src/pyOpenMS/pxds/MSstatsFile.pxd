from libcpp cimport bool
from Types cimport *
from String cimport *
from ConsensusMap cimport *
from ExperimentalDesign cimport *
from StringList cimport *

cdef extern from "<OpenMS/FORMAT/MSstatsFile.h>" namespace "OpenMS":

    cdef cppclass MSstatsFile:

        MSstatsFile() nogil except +
        
        void store(String & filename, ConsensusMap & consensus_map, ExperimentalDesign & design, 
          StringList & reannotate_filenames, bool is_isotope_label_type, 
          String & bioreplicate, String & condition, String & retention_time_summarization_method) nogil except +
