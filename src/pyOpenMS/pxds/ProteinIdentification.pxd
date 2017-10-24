from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from DataValue cimport *
from String cimport *
from Types cimport *
from MetaInfoInterface cimport *
from ProteinHit cimport *
from DigestionEnzymeProtein cimport *
from PeptideIdentification cimport *
from DateTime cimport *

cdef extern from "<OpenMS/METADATA/ProteinIdentification.h>" namespace "OpenMS":

    cdef cppclass ProteinIdentification(MetaInfoInterface):
        # wrap-inherits:
        #    MetaInfoInterface

        ProteinIdentification() nogil except +
        ProteinIdentification(ProteinIdentification) nogil except +

        bool operator==(ProteinIdentification) nogil except +
        bool operator!=(ProteinIdentification) nogil except +

        # cython has a problem with inheritance of overloaded methods,
        # so we do not declare them here, but separately in each derived
        # class which we want to be wrapped:
        void getKeys(libcpp_vector[String] & keys) nogil except +
        void getKeys(libcpp_vector[unsigned int] & keys) nogil except + # wrap-as:getKeysAsIntegers
        DataValue getMetaValue(unsigned int) nogil except +
        DataValue getMetaValue(String) nogil except +
        void setMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        void removeMetaValue(unsigned int) nogil except +

        # Returns the protein hits (mutable)
        libcpp_vector[ProteinHit] getHits() nogil except +
        # Appends a protein hit
        void insertHit(ProteinHit input) nogil except +
        # Sets the protein hits
        void setHits(libcpp_vector[ProteinHit] hits) nogil except +
        # Finds a protein hit by accession (returns past-the-end iterator if not found)
        # libcpp_vector[ProteinHit].iterator findHit(String accession) nogil except +

        # Returns the protein groups (mutable)
        libcpp_vector[ProteinGroup] getProteinGroups() nogil except +
        # Appends a new protein group
        void insertProteinGroup(ProteinGroup group) nogil except +

        # Returns the indistinguishable proteins (mutable)
        libcpp_vector[ProteinGroup] getIndistinguishableProteins() nogil except +
        # Appends new indistinguishable proteins
        void insertIndistinguishableProteins(ProteinGroup group) nogil except +

        # Returns the protein significance threshold value
        double getSignificanceThreshold() nogil except +
        # Sets the protein significance threshold value
        void setSignificanceThreshold(double value) nogil except +
        # Returns the protein score type
        String getScoreType() nogil except +
        # Sets the protein score type
        void setScoreType(String type) nogil except +
        # Returns true if a higher score represents a better score
        bool isHigherScoreBetter() nogil except +
        # Sets the orientation of the score (is higher better?)
        void setHigherScoreBetter(bool higher_is_better) nogil except +
        # Sorts the protein hits according to their score
        void sort() nogil except +
        # Sorts the protein hits by score and assigns ranks (best score has rank 1)
        void assignRanks() nogil except +

        # Compute the coverage (in percent) of all ProteinHits given PeptideHits
        void computeCoverage(libcpp_vector[PeptideIdentification] pep_ids) nogil except +

        # Returns the date of the protein identification run
        DateTime getDateTime() nogil except +
        # Sets the date of the protein identification run
        void setDateTime(DateTime date) nogil except +
        # Sets the search engine type
        void setSearchEngine(String search_engine) nogil except +
        # Returns the type of search engine used
        String getSearchEngine() nogil except +
        # Sets the search engine version
        void setSearchEngineVersion(String search_engine_version) nogil except +
        # Returns the search engine version
        String getSearchEngineVersion() nogil except +
        # Sets the search parameters
        void setSearchParameters(SearchParameters search_parameters) nogil except +
        # Returns the search parameters
        SearchParameters getSearchParameters() nogil except +
        # Returns the identifier
        String getIdentifier() nogil except +
        # Sets the identifier
        void setIdentifier(String id_) nogil except +

        void setPrimaryMSRunPath(StringList& s) nogil except +
        void getPrimaryMSRunPath(StringList& toFill) nogil except +

cdef extern from "<OpenMS/METADATA/ProteinIdentification.h>" namespace "OpenMS::ProteinIdentification":

    cdef enum PeakMassType:
        # wrap-attach:
        #    ProteinIdentification
        MONOISOTOPIC, AVERAGE, SIZE_OF_PEAKMASSTYPE

    cdef cppclass ProteinGroup:

      ProteinGroup()  nogil except +
      ProteinGroup(ProteinGroup)  nogil except +

      # Probability of this group
      double probability
      # Accessions of (indistinguishable) proteins that belong to the same group
      StringList accessions


    # Search parameters of the DB search
    cdef cppclass SearchParameters(MetaInfoInterface):
        # wrap-inherits:
        #    MetaInfoInterface

      SearchParameters()  nogil except +
      SearchParameters(SearchParameters) nogil except +

      String db            #< The used database
      String db_version            #< The database version
      String taxonomy            #< The taxonomy restriction
      String charges            #< The allowed charges for the search
      PeakMassType mass_type            #< Mass type of the peaks
      libcpp_vector[String] fixed_modifications            #< Used fixed modifications
      libcpp_vector[String] variable_modifications            #< Allowed variable modifications
      UInt missed_cleavages            #< The number of allowed missed cleavages
      double fragment_mass_tolerance            #< Mass tolerance of fragment ions (Dalton)
      bool fragment_mass_tolerance_ppm
      double precursor_mass_tolerance            #< Mass tolerance of precursor ions (Dalton)
      bool precursor_mass_tolerance_ppm
      DigestionEnzymeProtein digestion_enzyme            #< The enzyme for cleavage
