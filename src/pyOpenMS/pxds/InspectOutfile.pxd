from Types cimport *
from libcpp cimport bool
from libcpp.map cimport map as libcpp_map
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector
from String cimport *
from MSExperiment cimport *
from PeptideIdentification cimport *
from FileHandler cimport *

cdef extern from "<OpenMS/FORMAT/InspectOutfile.h>" namespace "OpenMS":
    
    cdef cppclass InspectOutfile "OpenMS::InspectOutfile":
        InspectOutfile() nogil except +
        InspectOutfile(InspectOutfile) nogil except +

        bool operator==(InspectOutfile & inspect_outfile) nogil except +

        libcpp_vector[ size_t ] load(const String & result_filename, 
                                     libcpp_vector[ PeptideIdentification ] & peptide_identifications,
                                     ProteinIdentification & protein_identification, double
                                     p_value_threshold, const String & database_filename) nogil except +

        libcpp_vector[ size_t ] getWantedRecords(const String & result_filename,
                                                 double p_value_threshold) nogil except +

        void compressTrieDB(const String & database_filename, const String &
                            index_filename, libcpp_vector[ size_t ] &
                            wanted_records, const String & snd_database_filename,
                            const String & snd_index_filename, bool append) nogil except +

        void generateTrieDB(const String & source_database_filename, const String &
                            database_filename, const String & index_filename, bool
                            append, const String species) nogil except +

        void getACAndACType(String line,
                            String & accession,
                            String & accession_type) nogil except +

        # TODO mixed, nested STL
        void getPrecursorRTandMZ(
            libcpp_vector[ libcpp_pair[ String, libcpp_vector[ libcpp_pair[ size_t, size_t ] ] ] ] & files_and_peptide_identification_with_scan_number,
            libcpp_vector[ PeptideIdentification ] & ids) nogil except + # wrap-ignore

        void getLabels(const String & source_database_filename,
                       String & ac_label,
                       String & sequence_start_label,
                       String & sequence_end_label,
                       String & comment_label,
                       String & species_label) nogil except +

        libcpp_vector[ size_t ] getSequences(const String & database_filename,
                                             libcpp_map[ size_t, size_t ] & wanted_records,
                                             libcpp_vector[ String ] & sequences) nogil except +

        void getExperiment(MSExperiment & exp, String & type_, const String & in_filename) nogil except +

        bool getSearchEngineAndVersion(const String & cmd_output, ProteinIdentification & protein_identification) nogil except +

        void readOutHeader(const String & filename,
                           const String & header_line,
                           Int & spectrum_file_column,
                           Int & scan_column,
                           Int & peptide_column,
                           Int & protein_column,
                           Int & charge_column,
                           Int & MQ_score_column,
                           Int & p_value_column,
                           Int & record_number_column,
                           Int & DB_file_pos_column,
                           Int & spec_file_pos_column,
                           Size & number_of_columns) nogil except +

