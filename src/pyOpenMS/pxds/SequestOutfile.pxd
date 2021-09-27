from libcpp cimport bool
from libcpp.map cimport map as libcpp_map
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from String cimport *
from PeptideIdentification cimport *
from ProteinIdentification cimport *
from DateTime cimport *

cdef extern from "<OpenMS/FORMAT/SequestOutfile.h>" namespace "OpenMS":
    
    cdef cppclass SequestOutfile "OpenMS::SequestOutfile":
        # wrap-doc:
        #   Representation of a Sequest output file
        
        SequestOutfile() nogil except + # wrap-doc:Representation of a Sequest output file
        SequestOutfile(SequestOutfile &) nogil except +

        bool operator==(SequestOutfile &sequest_outfile) nogil except +
        void load(const String &result_filename,
                  libcpp_vector[ PeptideIdentification ] &peptide_identifications,
                  ProteinIdentification &protein_identification,
                  double p_value_threshold,
                  libcpp_vector[ double ] &pvalues,
                  const String &database,
                  bool ignore_proteins_per_peptide) nogil except +
            # wrap-doc:
                #   Loads data from a Sequest outfile
                #   -----
                #   :param result_filename: The file to be loaded
                #   :param peptide_identifications: The identifications
                #   :param protein_identification: The protein identifications
                #   :param p_value_threshold: The significance level (for the peptide hit scores)
                #   :param pvalues: A list with the pvalues of the peptides (pvalues computed with peptide prophet)
                #   :param database: The database used for the search
                #   :param ignore_proteins_per_peptide: This is a hack to deal with files that use a suffix like "+1" in column "Reference", but do not actually list extra protein references in subsequent lines

        bool getColumns(const String &line, libcpp_vector[ String ] &substrings, Size number_of_columns, Size reference_column) nogil except + # wrap-doc:Retrieves columns from a Sequest outfile line
        void getSequences(const String &database_filename,
                          libcpp_map[ String, size_t ] &ac_position_map,
                          libcpp_vector[ String ] &sequences,
                          libcpp_vector[ libcpp_pair[ String, size_t ] ] &found,
                          libcpp_map[ String, size_t ] &not_found) nogil except + # wrap-ignore
        void getACAndACType(String line, String &accession, String &accession_type) nogil except + # wrap-doc:Retrieves the accession type and accession number from a protein description line

        # TODO immutable types by reference
        # 
        # void readOutHeader(String &result_filename, DateTime &datetime,
        # double &precursor_mz_value, Int &charge, Size &precursor_mass_type,
        # Size &ion_mass_type, Size &displayed_peptides, String &sequest,
        # String &sequest_version, String &database_type, Int &number_column,
        # Int &rank_sp_column, Int &id_column, Int &mh_column, Int
        # &delta_cn_column, Int &xcorr_column, Int &sp_column, Int &sf_column,
        # Int &ions_column, Int &reference_column, Int &peptide_column, Int
        # &score_column, Size &number_of_columns) nogil except +

