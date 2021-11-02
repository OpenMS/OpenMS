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
# from MSExperiment cimport *

cdef extern from "<OpenMS/METADATA/ProteinIdentification.h>" namespace "OpenMS":

    cdef cppclass ProteinIdentification(MetaInfoInterface):
        # wrap-inherits:
        #    MetaInfoInterface

        ProteinIdentification() nogil except +
        ProteinIdentification(ProteinIdentification &) nogil except +

        bool operator==(ProteinIdentification) nogil except +
        bool operator!=(ProteinIdentification) nogil except +

        
        libcpp_vector[ProteinHit] getHits() nogil except + # wrap-doc:Returns the protein hits
        
        void insertHit(ProteinHit input) nogil except + # wrap-doc:Appends a protein hit
        
        void setHits(libcpp_vector[ProteinHit] hits) nogil except + # wrap-doc:Sets the protein hits
        # Finds a protein hit by accession (returns past-the-end iterator if not found)
        # libcpp_vector[ProteinHit].iterator findHit(String accession) nogil except +

        
        libcpp_vector[ProteinGroup] getProteinGroups() nogil except + # wrap-doc:Returns the protein groups
        
        void insertProteinGroup(ProteinGroup group) nogil except + # wrap-doc:Appends a new protein group

        
        libcpp_vector[ProteinGroup] getIndistinguishableProteins() nogil except + # wrap-doc:Returns the indistinguishable proteins
        
        void insertIndistinguishableProteins(ProteinGroup group) nogil except + # wrap-doc:Appends new indistinguishable proteins

        
        double getSignificanceThreshold() nogil except + # wrap-doc:Returns the protein significance threshold value
        
        void setSignificanceThreshold(double value) nogil except + # wrap-doc:Sets the protein significance threshold value
        
        String getScoreType() nogil except + # wrap-doc:Returns the protein score type
        
        void setScoreType(String type) nogil except + # wrap-doc:Sets the protein score type
        
        bool isHigherScoreBetter() nogil except + # wrap-doc:Returns true if a higher score represents a better score
        
        void setHigherScoreBetter(bool higher_is_better) nogil except + # wrap-doc:Sets the orientation of the score (is higher better?)
        
        void sort() nogil except + # wrap-doc:Sorts the protein hits according to their score
        
        void assignRanks() nogil except + # wrap-doc:Sorts the protein hits by score and assigns ranks (best score has rank 1)

        
        void computeCoverage(libcpp_vector[PeptideIdentification] pep_ids) nogil except + # wrap-doc:Compute the coverage (in percent) of all ProteinHits given PeptideHits

        
        DateTime getDateTime() nogil except + # wrap-doc:Returns the date of the protein identification run
        
        void setDateTime(DateTime date) nogil except + # wrap-doc:Sets the date of the protein identification run
        
        void setSearchEngine(String search_engine) nogil except + # wrap-doc:Sets the search engine type
        
        String getSearchEngine() nogil except + # wrap-doc:Returns the type of search engine used
        
        void setSearchEngineVersion(String search_engine_version) nogil except + # wrap-doc:Sets the search engine version
        
        String getSearchEngineVersion() nogil except + # wrap-doc:Returns the search engine version
        
        void setSearchParameters(SearchParameters search_parameters) nogil except + # wrap-doc:Sets the search parameters
        
        SearchParameters getSearchParameters() nogil except + # wrap-doc:Returns the search parameters
        
        String getIdentifier() nogil except + # wrap-doc:Returns the identifier
        
        void setIdentifier(String id_) nogil except + # wrap-doc:Sets the identifier

        void setPrimaryMSRunPath(StringList& s) nogil except +
          # wrap-doc:
            #   Set the file paths to the primary MS runs (usually the mzML files obtained after data conversion from raw files)
            #   -----
            #   :param raw: Store paths to the raw files (or equivalent) rather than mzMLs

        void addPrimaryMSRunPath(StringList& s) nogil except +
        void getPrimaryMSRunPath(StringList& output) nogil except +

        void setPrimaryMSRunPath(StringList& s, bool raw) nogil except +
        void addPrimaryMSRunPath(StringList& s, bool raw) nogil except +
        void getPrimaryMSRunPath(StringList& output, bool raw) nogil except +

        # This causes as problem with circular dependencies when trying to use
        # ExperimentalSettings in MSExperiment
        # TODO: use addons if we really need this
        # void setPrimaryMSRunPath(StringList& s, MSExperiment& e) nogil except +

cdef extern from "<OpenMS/METADATA/ProteinIdentification.h>" namespace "OpenMS::ProteinIdentification":

    cdef enum PeakMassType:
        # wrap-attach:
        #    ProteinIdentification
        MONOISOTOPIC, AVERAGE, SIZE_OF_PEAKMASSTYPE

    cdef cppclass ProteinGroup:

      ProteinGroup()  nogil except +
      ProteinGroup(ProteinGroup &)  nogil except +

      # Probability of this group
      double probability
      # Accessions of (indistinguishable) proteins that belong to the same group
      StringList accessions


    # Search parameters of the DB search
    cdef cppclass SearchParameters(MetaInfoInterface):
        # wrap-inherits:
        #    MetaInfoInterface

      SearchParameters()  nogil except +
      SearchParameters(SearchParameters &) nogil except +

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

