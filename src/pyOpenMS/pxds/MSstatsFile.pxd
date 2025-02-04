from libcpp cimport bool
from Types cimport *
from String cimport *
from ConsensusMap cimport *
from ExperimentalDesign cimport *
from StringList cimport *


cdef extern from "<OpenMS/FORMAT/MSstatsFile.h>" namespace "OpenMS":

    cdef cppclass MSstatsFile:

        MSstatsFile() except + nogil 
        MSstatsFile(MSstatsFile &) except + nogil  # compiler

        void storeLFQ(String & filename,
                      ConsensusMap & consensus_map,
                      ExperimentalDesign & design,
                      StringList & reannotate_filenames,
                      bool is_isotope_label_type,
                      String & bioreplicate,
                      String & condition,
                      String & retention_time_summarization_method) except + nogil  # wrap-doc:Store label free experiment (MSstats)

        void storeISO(String & filename,
                      ConsensusMap & consensus_map,
                      ExperimentalDesign & design,
                      StringList & reannotate_filenames,
                      String & bioreplicate,
                      String & condition,
                      String & mixture,
                      String & retention_time_summarization_method) except + nogil  # wrap-doc:Store isobaric experiment (MSstatsTMT)
