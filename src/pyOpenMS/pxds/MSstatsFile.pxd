
cdef extern from "<OpenMS/FORMAT/MSstatsFile.h>" namespace "OpenMS":

    cdef cppclass MSstatsFile:

        MSstatsFile() nogil except +
        
        store(const String& filename, ConsensusMap& consensus_map, const ExperimentalDesign& design, 
          const StringList& reannotate_filenames, const bool is_isotope_label_type, 
          const String& bioreplicate, const String& condition, const String& retention_time_summarization_method)) nogil except+
