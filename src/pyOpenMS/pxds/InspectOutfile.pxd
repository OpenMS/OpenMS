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
      
        InspectOutfile() nogil except + # wrap-doc:This class serves to read in an Inspect outfile and write an idXML file
        InspectOutfile(InspectOutfile &) nogil except +


        bool operator==(InspectOutfile & inspect_outfile) nogil except +

        libcpp_vector[ size_t ] load(const String & result_filename, 
                                     libcpp_vector[ PeptideIdentification ] & peptide_identifications,
                                     ProteinIdentification & protein_identification, double
                                     p_value_threshold, const String & database_filename) nogil except +
            # wrap-doc:
                #   Load the results of an Inspect search
                #   -----
                #   :param result_filename: Input parameter which is the file name of the input file
                #   :param peptide_identifications: Output parameter which holds the peptide identifications from the given file
                #   :param protein_identification: Output parameter which holds the protein identifications from the given file
                #   :param p_value_threshold
                #   :param database_filename
                #   :raises:
                #     Exception: FileNotFound is thrown if the given file could not be found
                #   :raises:
                #     Exception: ParseError is thrown if the given file could not be parsed
                #   :raises:
                #     Exception: FileEmpty is thrown if the given file is empty

        libcpp_vector[ size_t ] getWantedRecords(const String & result_filename,
                                                 double p_value_threshold) nogil except +
            # wrap-doc:
                #   Loads only results which exceeds a given p-value threshold
                #   -----
                #   :param result_filename: The filename of the results file
                #   :param p_value_threshold: Only identifications exceeding this threshold are read
                #   :raises:
                #     Exception: FileNotFound is thrown if the given file could not be found
                #   :raises:
                #     Exception: FileEmpty is thrown if the given file is empty

        void compressTrieDB(const String & database_filename, const String &
                            index_filename, libcpp_vector[ size_t ] &
                            wanted_records, const String & snd_database_filename,
                            const String & snd_index_filename, bool append) nogil except + # wrap-doc:Generates a trie database from another one, using the wanted records only

        void generateTrieDB(const String & source_database_filename, const String &
                            database_filename, const String & index_filename, bool
                            append, const String species) nogil except + # wrap-doc:Generates a trie database from a given one (the type of database is determined by getLabels)

        void getACAndACType(String line,
                            String & accession,
                            String & accession_type) nogil except + # wrap-doc:Retrieve the accession type and accession number from a protein description line

        # TODO mixed, nested STL
        void getPrecursorRTandMZ(
            libcpp_vector[ libcpp_pair[ String, libcpp_vector[ libcpp_pair[ size_t, size_t ] ] ] ] & files_and_peptide_identification_with_scan_number,
            libcpp_vector[ PeptideIdentification ] & ids) nogil except + # wrap-ignore

        void getLabels(const String & source_database_filename,
                       String & ac_label,
                       String & sequence_start_label,
                       String & sequence_end_label,
                       String & comment_label,
                       String & species_label) nogil except + # wrap-doc:Retrieve the labels of a given database (at the moment FASTA and Swissprot)

        libcpp_vector[ size_t ] getSequences(const String & database_filename,
                                             libcpp_map[ size_t, size_t ] & wanted_records,
                                             libcpp_vector[ String ] & sequences) nogil except + # wrap-doc:Retrieve sequences from a trie database

        void getExperiment(MSExperiment & exp, String & type_, const String & in_filename) nogil except + # wrap-doc:Get the experiment from a file

        bool getSearchEngineAndVersion(const String & cmd_output, ProteinIdentification & protein_identification) nogil except + # wrap-doc:Get the search engine and its version from the output of the InsPecT executable without parameters. Returns true on success, false otherwise

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
                           Size & number_of_columns) nogil except + # wrap-doc:Read the header of an inspect output file and retrieve various information

