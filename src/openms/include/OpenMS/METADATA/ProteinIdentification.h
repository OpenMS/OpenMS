// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Nico Pfeifer, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/CHEMISTRY/DigestionEnzymeProtein.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/METADATA/DataArrays.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <set>

namespace OpenMS
{
  class MSExperiment;
  class PeptideIdentification;
  class PeptideEvidence;
  class ConsensusMap;

  /**
    @brief Representation of a protein identification run

    This class stores the general information and the protein hits of a protein identification run.

    The actual peptide hits are stored in PeptideIdentification instances that are part of spectra or features.

    In order to be able to connect the ProteinIdentification and the
    corresponding peptide identifications, both classes have a string
    identifier. We recommend using the search engine name and the date as
    identifier.
    Setting this identifier is especially important when there are several
    protein identification runs for a map, i.e. several ProteinIdentification
    instances.

    @todo Add MetaInfoInterface to modifications => update IdXMLFile and ProteinIdentificationVisualizer (Andreas)

    @ingroup Metadata
  */
  class OPENMS_DLLAPI ProteinIdentification :
    public MetaInfoInterface
  {
public:
    /// Hit type definition
    typedef ProteinHit HitType;

    /// two way mapping from ms-run-path to protID|pepID-identifier
    struct Mapping
    {
      std::map<String, StringList> identifier_to_msrunpath;
      std::map<StringList, String> runpath_to_identifier;

      Mapping() = default;

      explicit Mapping(const std::vector<ProteinIdentification>& prot_ids)
      {
        create(prot_ids);
      }

      void create(const std::vector<ProteinIdentification>& prot_ids)
      {
        identifier_to_msrunpath.clear();
        runpath_to_identifier.clear();
        StringList filenames;
        for (const ProteinIdentification& prot_id : prot_ids)
        {
          prot_id.getPrimaryMSRunPath(filenames);
          if (filenames.empty())
          {
            throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No MS run path annotated in ProteinIdentification.");
          }
          identifier_to_msrunpath[prot_id.getIdentifier()] = filenames;
          const auto& it = runpath_to_identifier.find(filenames);
          if (it != runpath_to_identifier.end())
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Multiple protein identifications with the same ms-run-path in Consensus/FeatureXML. Check input!\n",
                                          ListUtils::concatenate(filenames, ","));
          }
          runpath_to_identifier[filenames] = prot_id.getIdentifier();
        }
      }

      String getPrimaryMSRunPath(const PeptideIdentification& pepid) const
      { 
        // if a merge index n is annotated, we use the filename annotated at index n in the protein identification, otherwise the one at index 0        
        size_t merge_index = pepid.getMetaValue(Constants::UserParam::ID_MERGE_INDEX, 0);
        const auto& filenames = identifier_to_msrunpath.at(pepid.getIdentifier());        
        return (merge_index < filenames.size()) ?  filenames[merge_index] : ""; // return filename or empty string if missing
      }
    };

    /**
        @brief Bundles multiple (e.g. indistinguishable) proteins in a group
    */
    class OPENMS_DLLAPI ProteinGroup
    {
    public:
      /// Float data array vector type
      typedef OpenMS::DataArrays::FloatDataArray FloatDataArray ;
      typedef std::vector<FloatDataArray> FloatDataArrays;
      /// String data array vector type
      typedef OpenMS::DataArrays::StringDataArray StringDataArray ;
      typedef std::vector<StringDataArray> StringDataArrays;
      /// Integer data array vector type
      typedef OpenMS::DataArrays::IntegerDataArray IntegerDataArray ;
      typedef std::vector<IntegerDataArray> IntegerDataArrays;

      /// Probability of this group
      double probability;

      /// Accessions of (indistinguishable) proteins that belong to the same group
      std::vector<String> accessions;

      ProteinGroup();

      /// Equality operator
      bool operator==(const ProteinGroup& rhs) const;

      /*
        @brief Comparison operator (for sorting)

        This operator is intended for sorting protein groups in a "best first"
        manner. That means higher probabilities are "less" than lower
        probabilities (!); smaller groups are "less" than larger groups;
        everything else being equal, accessions are compared lexicographically.
      */
      bool operator<(const ProteinGroup& rhs) const;

      /// Float data arrays
      /**
      @name data array methods

      These methods are used to annotate protein group meta information.

      These statements should help you chose which approach to use
        - Access to meta info arrays is slower than to a member variable
        - Access to meta info arrays is faster than to a %MetaInfoInterface
        - Meta info arrays are stored when using mzML format for storing
    */
    //@{
      /// Returns a const reference to the float meta data arrays
      const FloatDataArrays& getFloatDataArrays() const;

      /// Returns a mutable reference to the float meta data arrays
      FloatDataArrays& getFloatDataArrays()
      {
        return float_data_arrays_;
      }

      /// Sets the float meta data arrays
      void setFloatDataArrays(const FloatDataArrays& fda);

      /// Returns a const reference to the string meta data arrays
      const StringDataArrays& getStringDataArrays() const;

      /// Returns a mutable reference to the string meta data arrays
      StringDataArrays& getStringDataArrays();

      /// Sets the string meta data arrays
      void setStringDataArrays(const StringDataArrays& sda);

      /// Returns a const reference to the integer meta data arrays
      const IntegerDataArrays& getIntegerDataArrays() const;

      /// Returns a mutable reference to the integer meta data arrays
      IntegerDataArrays& getIntegerDataArrays();

      /// Sets the integer meta data arrays
      void setIntegerDataArrays(const IntegerDataArrays& ida);

      /// Returns a mutable reference to the first integer meta data array with the given name
      inline IntegerDataArray& getIntegerDataArrayByName(String name)
      {
        return *std::find_if(integer_data_arrays_.begin(), integer_data_arrays_.end(),
          [&name](const IntegerDataArray& da) { return da.getName() == name; } );
      }

      /// Returns a mutable reference to the first string meta data array with the given name
      inline StringDataArray& getStringDataArrayByName(String name)
      {
        return *std::find_if(string_data_arrays_.begin(), string_data_arrays_.end(),
          [&name](const StringDataArray& da) { return da.getName() == name; } );
      }

      /// Returns a mutable reference to the first float meta data array with the given name
      inline FloatDataArray& getFloatDataArrayByName(String name)
      {
        return *std::find_if(float_data_arrays_.begin(), float_data_arrays_.end(),
          [&name](const FloatDataArray& da) { return da.getName() == name; } );
      }

      /// Returns a const reference to the first integer meta data array with the given name
      inline const IntegerDataArray& getIntegerDataArrayByName(String name) const
      {
        return *std::find_if(integer_data_arrays_.begin(), integer_data_arrays_.end(),
          [&name](const IntegerDataArray& da) { return da.getName() == name; } );
      }

      /// Returns a const reference to the first string meta data array with the given name
      inline const StringDataArray& getStringDataArrayByName(String name) const
      {
        return *std::find_if(string_data_arrays_.begin(), string_data_arrays_.end(),
          [&name](const StringDataArray& da) { return da.getName() == name; } );
      }

      /// Returns a const reference to the first float meta data array with the given name
      inline const FloatDataArray& getFloatDataArrayByName(String name) const
      {
        return *std::find_if(float_data_arrays_.begin(), float_data_arrays_.end(),
          [&name](const FloatDataArray& da) { return da.getName() == name; } );
      }

    private:
      /// Float data arrays
      FloatDataArrays float_data_arrays_;

      /// String data arrays
      StringDataArrays string_data_arrays_;

      /// Integer data arrays
      IntegerDataArrays integer_data_arrays_;
    };

    /// Peak mass type
    enum PeakMassType
    {
      MONOISOTOPIC,
      AVERAGE,
      SIZE_OF_PEAKMASSTYPE
    };

    /// Names corresponding to peak mass types
    static const std::string NamesOfPeakMassType[SIZE_OF_PEAKMASSTYPE];

    /// Search parameters of the DB search
    struct OPENMS_DLLAPI SearchParameters :
      public MetaInfoInterface
    {
      String db; ///< The used database
      String db_version; ///< The database version
      String taxonomy; ///< The taxonomy restriction
      String charges; ///< The allowed charges for the search
      PeakMassType mass_type; ///< Mass type of the peaks
      std::vector<String> fixed_modifications; ///< Used fixed modifications
      std::vector<String> variable_modifications; ///< Allowed variable modifications
      UInt missed_cleavages; ///< The number of allowed missed cleavages
      double fragment_mass_tolerance; ///< Mass tolerance of fragment ions (Dalton or ppm)
      bool fragment_mass_tolerance_ppm; ///< Mass tolerance unit of fragment ions (true: ppm, false: Dalton)
      double precursor_mass_tolerance; ///< Mass tolerance of precursor ions (Dalton or ppm)
      bool precursor_mass_tolerance_ppm; ///< Mass tolerance unit of precursor ions (true: ppm, false: Dalton)
      Protease digestion_enzyme; ///< The cleavage site information in details (from ProteaseDB)
      EnzymaticDigestion::Specificity enzyme_term_specificity; ///< The number of required cutting-rule matching termini during search (none=0, semi=1, or full=2)

      SearchParameters();
      /// Copy constructor
      SearchParameters(const SearchParameters&) = default;
      /// Move constructor
      SearchParameters(SearchParameters&&) = default;
      /// Destructor
      ~SearchParameters() = default;

      /// Assignment operator
      SearchParameters& operator=(const SearchParameters&) = default;
      /// Move assignment operator
      SearchParameters& operator=(SearchParameters&&)& = default;

      bool operator==(const SearchParameters& rhs) const;

      bool operator!=(const SearchParameters& rhs) const;

      /// returns the charge range from the search engine settings as a pair of ints
      std::pair<int,int> getChargeRange() const;

      /// Tests if these search engine settings are mergeable with @p sp
      /// depending on the given @p experiment_type.
      /// Modifications are compared as sets. Databases based on filename.
      /// "labeled_MS1" experiments additionally allow different modifications.
      bool mergeable(const ProteinIdentification::SearchParameters& sp, const String& experiment_type) const;

      private:
      int getChargeValue_(String& charge_str) const;
    };

    /** @name Constructors, destructors, assignment operator <br> */
    //@{
    /// Default constructor
    ProteinIdentification();
    /// Copy constructor
    ProteinIdentification(const ProteinIdentification&) = default;
    /// Move constructor
    ProteinIdentification(ProteinIdentification&&) = default;
    /// Destructor
    virtual ~ProteinIdentification();

    /// Assignment operator
    ProteinIdentification& operator=(const ProteinIdentification&) = default;
    /// Move assignment operator
    ProteinIdentification& operator=(ProteinIdentification&&) = default;

    /// Equality operator
    bool operator==(const ProteinIdentification& rhs) const;
    /// Inequality operator
    bool operator!=(const ProteinIdentification& rhs) const;
    //@}

    ///@name Protein hit information (public members)
    //@{
    /// Returns the protein hits
    const std::vector<ProteinHit>& getHits() const;
    /// Returns the protein hits (mutable)
    std::vector<ProteinHit>& getHits();
    /// Appends a protein hit
    void insertHit(const ProteinHit& input);
    /// Appends a protein hit
    void insertHit(ProteinHit&& input);

    /**
        @brief Sets the protein hits

        @note This may invalidate (indistinguishable) protein groups! If necessary, use e.g. @p IDFilter::updateProteinGroups to update the groupings.
     */
    void setHits(const std::vector<ProteinHit>& hits);

    /// Finds a protein hit by accession (returns past-the-end iterator if not found)
    std::vector<ProteinHit>::iterator findHit(const String& accession);

    /// Returns the protein groups
    const std::vector<ProteinGroup>& getProteinGroups() const;
    /// Returns the protein groups (mutable)
    std::vector<ProteinGroup>& getProteinGroups();
    /// Appends a new protein group
    void insertProteinGroup(const ProteinGroup& group);

    /// Returns the indistinguishable proteins
    const std::vector<ProteinGroup>& getIndistinguishableProteins() const;
    /// Returns the indistinguishable proteins (mutable)
    std::vector<ProteinGroup>& getIndistinguishableProteins();
    /// Appends new indistinguishable proteins
    void insertIndistinguishableProteins(const ProteinGroup& group);
    /// Appends singleton groups (with the current score) for every yet ungrouped protein hit
    void fillIndistinguishableGroupsWithSingletons();

    /// Returns the protein significance threshold value
    double getSignificanceThreshold() const;
    /// Sets the protein significance threshold value
    void setSignificanceThreshold(double value);
    /// Returns the protein score type
    const String& getScoreType() const;
    /// Sets the protein score type
    void setScoreType(const String& type);
    /// Returns true if a higher score represents a better score
    bool isHigherScoreBetter() const;
    /// Sets the orientation of the score (is higher better?)
    void setHigherScoreBetter(bool higher_is_better);
    /// Sorts the protein hits according to their score
    void sort();
    /// Sorts the protein hits by score and assigns ranks (best score has rank 1)
    void assignRanks();
    /**
       @brief Compute the coverage (in percent) of all ProteinHits given PeptideHits

       @throws Exception::MissingInformation if ProteinsHits do not have sequence information

       Does not return anything but stores the coverage inside the ProteinHit objects
    */
    void computeCoverage(const std::vector<PeptideIdentification>& pep_ids);
    void computeCoverage(const ConsensusMap& cmap, bool use_unassigned_ids);
    //@}

    /**
       @brief Compute the modifications of all ProteinHits given PeptideHits

      For every protein accession, the pair of position and modification is returned.
      Because fixed modifications might not be of interest, a list can be provided to skip those.
    */
    void computeModifications(
      const std::vector<PeptideIdentification>& pep_ids,
      const StringList& skip_modifications);
    void computeModifications(
      const ConsensusMap& cmap,
      const StringList& skip_modifications,
      bool use_unassigned_ids);


    ///@name General information
    //@{
    /// Returns the date of the protein identification run
    const DateTime& getDateTime() const;
    /// Sets the date of the protein identification run
    void setDateTime(const DateTime& date);
    /// Sets the search engine type
    void setSearchEngine(const String& search_engine);
    /// Returns the type of search engine used
    const String& getSearchEngine() const;
    /// Return the type of search engine that was first applied (e.g., before percolator or consensusID) or "Unknown"
    const String getOriginalSearchEngineName() const;
    /// Sets the search engine version
    void setSearchEngineVersion(const String& search_engine_version);
    /// Returns the search engine version
    const String& getSearchEngineVersion() const;
    /// Sets the inference engine type
    void setInferenceEngine(const String& search_engine);
    /// Returns the type of search engine used
    const String getInferenceEngine() const;
    /// Sets the search engine version
    void setInferenceEngineVersion(const String& inference_engine_version);
    /// Returns the search engine version
    const String getInferenceEngineVersion() const;
    /// Sets the search parameters
    void setSearchParameters(const SearchParameters& search_parameters);
    /// Sets the search parameters (move)
    void setSearchParameters(SearchParameters&& search_parameters);
    /// Returns the search parameters
    const SearchParameters& getSearchParameters() const;
    /// Returns the search parameters (mutable)
    SearchParameters& getSearchParameters();
    /// Returns the identifier
    const String& getIdentifier() const;
    /// Sets the identifier
    void setIdentifier(const String& id);
    /**
       Set the file paths to the primary MS runs (usually the mzML files obtained after data conversion from raw files)

       @param s The file paths
       @param raw Store paths to the raw files (or equivalent) rather than mzMLs
    */
    void setPrimaryMSRunPath(const StringList& s, bool raw = false);

    /// set the file path to the primary MS run but try to use the mzML annotated in the MSExperiment.
    void setPrimaryMSRunPath(const StringList& s, MSExperiment& e);
    void addPrimaryMSRunPath(const String& s, bool raw = false);
    void addPrimaryMSRunPath(const StringList& s, bool raw = false);

    /**
       Get the file paths to the primary MS runs

       @param[out] output The file paths
       @param raw Get raw files (or equivalent) instead of mzMLs
    */
    void getPrimaryMSRunPath(StringList& output, bool raw = false) const;

    /// get the number of primary MS runs involve in this ID run
    Size nrPrimaryMSRunPaths(bool raw = false) const;

    /// Checks if this object has inference data. Looks for "InferenceEngine" metavalue.
    /// If not, falls back to old behaviour of reading the search engine name.
    bool hasInferenceData() const;

    /// Checks if the search engine name matches an inference engine known to OpenMS.
    bool hasInferenceEngineAsSearchEngine() const;

    /// Checks if the peptide IDs of this IDRun are mergeable with another @p id_run
    /// given an @p experiment_type .
    /// Checks search engine and search engine settings.
    bool peptideIDsMergeable(const ProteinIdentification& id_run, const String& experiment_type) const;

    /// Collects all search engine settings registered for the given search engine @p se.
    /// If @p se is empty, the main search engine is used, otherwise it will also search the metavalues.
    std::vector<std::pair<String,String>> getSearchEngineSettingsAsPairs(const String& se = "") const;

    //@}

    /// Copies only metadata (no protein hits or protein groups)
    void copyMetaDataOnly(const ProteinIdentification&);
protected:
    ///@name General information (search engine, parameters and database)
    //@{
    String id_;
    String search_engine_;
    String search_engine_version_;
    SearchParameters search_parameters_;
    DateTime date_;
    //@}

    ///@name Protein hit information (protected members)
    //@{
    String protein_score_type_;
    bool higher_score_better_;
    std::vector<ProteinHit> protein_hits_;
    std::vector<ProteinGroup> protein_groups_;
    /// Indistinguishable proteins: @p accessions[0] is "group leader", @p probability is meaningless
    std::vector<ProteinGroup> indistinguishable_proteins_;
    double protein_significance_threshold_;
    //@}

  private:
    void computeCoverageFromEvidenceMapping_(const std::unordered_map<String, std::set<PeptideEvidence>>& map);
    void fillEvidenceMapping_(std::unordered_map<String, std::set<PeptideEvidence> >& map_acc_2_evidence,
                              const std::vector<PeptideIdentification>& pep_ids) const;

    void fillModMapping_(const std::vector<PeptideIdentification>& pep_ids, const StringList& skip_modifications,
                         std::unordered_map<String, std::set<std::pair<Size, ResidueModification>>>& prot2mod) const;
  };


} //namespace OpenMS
