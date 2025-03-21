// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Erhan Kenar, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/FORMAT/MzTabM.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CHEMISTRY/AdductInfo.h>


#include <iosfwd>
#include <vector>

namespace OpenMS
{
  class OPENMS_DLLAPI AccurateMassSearchResult
  {
  public:
    /// Default constructor
    AccurateMassSearchResult();

    /// Default destructor
    ~AccurateMassSearchResult();

    /// copy constructor
    AccurateMassSearchResult(const AccurateMassSearchResult&);

    /// assignment operator
    AccurateMassSearchResult& operator=(const AccurateMassSearchResult&);

    /// get the m/z of the small molecule + adduct
    double getObservedMZ() const;

    /// set the m/z of the small molecule + adduct
    void setObservedMZ(const double&);

    /// get the theoretical m/z of the small molecule + adduct
    double getCalculatedMZ() const;

    /// set the theoretical m/z of the small molecule + adduct
    void setCalculatedMZ(const double&);

    /// get the mass used to query the database (uncharged small molecule)
    double getQueryMass() const;

    /// set the mass used to query the database (uncharged small molecule)
    void setQueryMass(const double&);

    /// get the mass returned by the query (uncharged small molecule)
    double getFoundMass() const;

    /// set the mass returned by the query (uncharged small molecule)
    void setFoundMass(const double&);

    /// get the charge
    Int getCharge() const;

    /// set the charge
    void setCharge(const Int&);

    /// get the error between observed and theoretical m/z in ppm
    double getMZErrorPPM() const;

    /// set the error between observed and theoretical m/z in ppm
    void setMZErrorPPM(const double);

    /// get the observed rt
    double getObservedRT() const;

    /// set the observed rt
    void setObservedRT(const double& rt);

    /// get the observed intensity
    double getObservedIntensity() const;

    /// set the observed intensity
    void setObservedIntensity(const double&);

    /// get the observed intensities
    std::vector<double> getIndividualIntensities() const;

    /// set the observed intensities
    void setIndividualIntensities(const std::vector<double>&);

    Size getMatchingIndex() const;
    void setMatchingIndex(const Size&);

    Size getSourceFeatureIndex() const;
    void setSourceFeatureIndex(const Size&);

    const String& getFoundAdduct() const;
    void setFoundAdduct(const String&);

    const String& getFormulaString() const;
    void setEmpiricalFormula(const String&);

    const std::vector<String>& getMatchingHMDBids() const;
    void setMatchingHMDBids(const std::vector<String>&);

    /// return trace intensities of the underlying feature;
    const std::vector<double>& getMasstraceIntensities() const;
    void setMasstraceIntensities(const std::vector<double>&);

    double getIsotopesSimScore() const;
    void setIsotopesSimScore(const double&);

    // debug/output functions
    friend OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const AccurateMassSearchResult& amsr);

private:
    /// Stored information/results of DB query
    double observed_mz_;
    double theoretical_mz_;
    double searched_mass_;
    double db_mass_;
    Int charge_;
    double mz_error_ppm_;
    double observed_rt_;
    double observed_intensity_;
    std::vector<double> individual_intensities_;
    Size matching_index_;
    Size source_feature_index_;

    String found_adduct_;
    String empirical_formula_;
    std::vector<String> matching_hmdb_ids_;

    std::vector<double> mass_trace_intensities_;
    double isotopes_sim_score_;
  };

  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const AccurateMassSearchResult& amsr);

  /**
    @brief An algorithm to search for exact mass matches from a spectrum against a database (e.g. HMDB).

    For each peak, neutral masses are reconstructed from observed (spectrum) m/z values by enumerating all possible adducts with matching charge.
    The resulting neutral masses (can be more than one, depending on list of possible adducts) are matched against masses from
    a database within a certain mass error (Da or ppm).

    Supports any database which contains an identifier, chemical sum formula and (optional) mass.
    If masses in the database are not given (= set to 0), they are computed from sum formulas.

    Both positive and negative ion mode is supported. Charge for (Consensus-)Features can be either positive or negative, but
    only the absolute value is used since many FeatureFinders will only report positive charges even in negative ion mode.
    Entities with charge=0 are treated as "unknown charge" and are tested with all potential adducts and subsequently matched against the database.

    A file with a list of potential adducts can be given for each mode separately.
    Each line contains a chemical formula (plus quantor) and a charge (separated by semicolon), e.g.
    M+H;1+
    The M can be preceded by a quantor (e.g.2M, 3M), implicitly assumed as 1.
    The chemical formula can contain multiple segments, separated by + or - operators, e.g. M+H-H2O;+1 (water loss in positive mode).
    Brackets are implicit per segment, i.e. M+H-H2O is parsed as M + (H) - (H2O).
    Each segment can also be preceded by a quantor, e.g. M+H-H2O would parse as
    M + (H) - 2x(H2O).
    If debug mode is enabled, the masses of each segment are printed for verification.
    In particular, typing H20 (twenty H) is different from H2O (water).

    Ionization mode of the observed m/z values can be determined automatically if the input map (either FeatureMap or ConsensusMap) is annotated
    with a meta value, as done by @ref TOPP_FeatureFinderMetabo.


    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI AccurateMassSearchEngine :
    public DefaultParamHandler,
    public ProgressLogger
  {
public:

    /// uses 'AccurateMassSearchEngine' as search engine id for protein and peptide ids which are generated by AMS
    static constexpr char search_engine_identifier[] = "AccurateMassSearchEngine";

    /// Default constructor
    AccurateMassSearchEngine();

    /// Default destructor
    ~AccurateMassSearchEngine() override;

    /**
      @brief search for a specific observed mass by enumerating all possible adducts and search M+X against database.
      If use_feature_adducts is activated, queryByMZ uses annotated, observed adducts as EmpiricalFormulas, restricting M+X candidates.

       */
    void queryByMZ(const double& observed_mz, const Int& observed_charge, const String& ion_mode, std::vector<AccurateMassSearchResult>& results, const EmpiricalFormula& observed_adduct = EmpiricalFormula()) const;
    void queryByFeature(const Feature& feature, const Size& feature_index, const String& ion_mode, std::vector<AccurateMassSearchResult>& results) const;
    void queryByConsensusFeature(const ConsensusFeature& cfeat, const Size& cf_index, const Size& number_of_maps, const String& ion_mode, std::vector<AccurateMassSearchResult>& results) const;

    /// main method of AccurateMassSearchEngine
    /// input map is not const, since it will get annotated with results
    void run(FeatureMap&, MzTab&) const;

    void run(FeatureMap&, MzTabM&) const;

    /// main method of AccurateMassSearchEngine
    /// input map is not const, since it will get annotated with results
    /// @note Call init() before calling run!
    void run(ConsensusMap&, MzTab&) const;

    /// parse database and adduct files
    void init();

protected:
    void updateMembers_() override;

private:
    /// private member functions

    /// if ion-mode is auto, this will set the internal mode according to input data
    /// @throw InvalidParameter if ion mode cannot be resolved
    template <typename MAPTYPE> String resolveAutoMode_(const MAPTYPE& map) const
    {
      String ion_mode_internal;
      String ion_mode_detect_msg = "";
      if (map.size() > 0)
      {
        if (map[0].metaValueExists("scan_polarity"))
        {
          StringList pols = ListUtils::create<String>(String(map[0].getMetaValue("scan_polarity")), ';');
          if (pols.size() == 1 && !pols[0].empty())
          {
            pols[0].toLower();
            if (pols[0] == "positive" || pols[0] == "negative")
            {
              ion_mode_internal = pols[0];
              OPENMS_LOG_INFO << "Setting auto ion-mode to '" << ion_mode_internal << "' for file " << File::basename(map.getLoadedFilePath()) << std::endl;
            }
            else ion_mode_detect_msg = String("Meta value 'scan_polarity' does not contain unknown ion mode") + String(map[0].getMetaValue("scan_polarity"));
          }
          else
          {
            ion_mode_detect_msg = String("ambiguous ion mode: ") + String(map[0].getMetaValue("scan_polarity"));
          }
        }
        else
        {
          ion_mode_detect_msg = String("Meta value 'scan_polarity' not found in (Consensus-)Feature map");
        }
      }
      else
      { // do nothing, since map is
        OPENMS_LOG_INFO << "Meta value 'scan_polarity' cannot be determined since (Consensus-)Feature map is empty!" << std::endl;
      }

      if (!ion_mode_detect_msg.empty())
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Auto ionization mode could not resolve ion mode of data (") + ion_mode_detect_msg + "!");
      }

      return ion_mode_internal;
    }

    void parseMappingFile_(const StringList&);
    void parseStructMappingFile_(const StringList&);
    void parseAdductsFile_(const String& filename, std::vector<AdductInfo>& result);
    void searchMass_(double neutral_query_mass, double diff_mass, std::pair<Size, Size>& hit_indices) const;

    /// Add search results to a Consensus/Feature
    void annotate_(const std::vector<AccurateMassSearchResult>&, BaseFeature&) const;

    /// Extract query results from feature
    std::vector<AccurateMassSearchResult> extractQueryResults_(const Feature& feature, const Size& feature_index, const String& ion_mode_internal, Size& dummy_count) const;

    /// Add resulting matches to IdentificationData
    void addMatchesToID_(
      IdentificationData& id,
      const std::vector<AccurateMassSearchResult>& amr, 
      const IdentificationData::InputFileRef& file_ref,
      const IdentificationData::ScoreTypeRef& mass_error_ppm_score_ref,
      const IdentificationData::ScoreTypeRef& mass_error_Da_score_ref,
      const IdentificationData::ProcessingStepRef& step_ref,
      BaseFeature& f) const;

    /// For two vectors of identical length, compute the cosine of the angle between them.
    /// Since we look at the angle, scaling of the vectors does not change the result (when ignoring numerical instability).
    double computeCosineSim_(const std::vector<double>& x, const std::vector<double>& y) const;

    double computeIsotopePatternSimilarity_(const Feature& feat, const EmpiricalFormula& form) const;

    typedef std::vector<std::vector<AccurateMassSearchResult> > QueryResultsTable;

    void exportMzTab_(const QueryResultsTable& overall_results, const Size number_of_maps, MzTab& mztab_out, const std::vector<String>& file_locations) const;

    void exportMzTabM_(const FeatureMap& fmap, MzTabM& mztabm_out) const;

    /// private member variables
    typedef std::vector<std::vector<String> > MassIDMapping;
    typedef std::map<String, std::vector<String> > HMDBPropsMapping;

    struct MappingEntry_
    {
      double mass;
      std::vector<String> massIDs;
      String formula;
    };
    std::vector<MappingEntry_> mass_mappings_;

    struct CompareEntryAndMass_ // defined here to allow for inlining by compiler
    {
      double asMass(const MappingEntry_& v) const
      {
        return v.mass;
      }

      double asMass(double t) const
      {
        return t;
      }

      template <typename T1, typename T2>
      bool operator()(T1 const& t1, T2 const& t2) const
      {
        return asMass(t1) < asMass(t2);
      }

    };

    HMDBPropsMapping hmdb_properties_mapping_;

    bool is_initialized_; ///< true if init_() was called without any subsequent param changes

    bool legacyID_ = true;

    /// parameter stuff
    double mass_error_value_;
    String mass_error_unit_;
    String ion_mode_;
    bool iso_similarity_;

    String pos_adducts_fname_;
    String neg_adducts_fname_;

    StringList db_mapping_file_;
    StringList db_struct_file_;

    std::vector<AdductInfo> pos_adducts_;
    std::vector<AdductInfo> neg_adducts_;

    String database_name_;
    String database_version_;
    String database_location_;

    bool keep_unidentified_masses_;
  };

}
