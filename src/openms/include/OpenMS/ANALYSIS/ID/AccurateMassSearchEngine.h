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
  /**
    @brief One small-molecule hit produced by AccurateMassSearchEngine.

    Represents a single (observed feature/peak) — (database compound +
    adduct) match. AccurateMassSearchEngine emits one of these per accepted
    (feature, adduct, candidate-compound) combination, so the same input
    feature can appear in several results if it is consistent with more
    than one adduct or matches multiple isobaric compounds in the database.

    @section AccurateMassSearchResult_lifecycle Lifecycle

    Each result is built up incrementally by the search engine via the
    setters and exposed read-only to consumers via the getters. The class
    holds no logic of its own — it is a plain value object grouping the
    observable, the database hit, and the diagnostics that explain why the
    two were matched.

    @section AccurateMassSearchResult_fields Field semantics

    Observed quantities (from the spectrum / feature):
      - @c observed_mz       : observed @em m/z of the feature.
      - @c observed_rt       : retention time of the feature.
      - @c observed_intensity: intensity of the feature.
      - @c individual_intensities : per-FeatureHandle intensities, populated
                                    only on the @c ConsensusMap overload of
                                    @ref AccurateMassSearchEngine::run; stays
                                    empty when the engine is invoked on a plain
                                    @ref FeatureMap.
      - @c mass_trace_intensities : per-isotopologue intensities of the
                                    underlying mass traces (used for the
                                    isotope-pattern similarity score).

    Match-derived quantities:
      - @c searched_mass : neutral mass back-calculated from
                           @c observed_mz under the assumption of
                           @c found_adduct and @c charge. Used as the
                           database query.
      - @c db_mass       : neutral mass of the matched database compound
                           (the "ground truth" for this hit).
      - @c theoretical_mz: @em m/z that the matched compound + adduct
                           would produce; differs from @c observed_mz by
                           @c mz_error_ppm.
      - @c mz_error_ppm  : signed @em m/z error in ppm, observed vs.
                           theoretical.
      - @c charge        : ion charge under which the match was attempted
                           (positive for cations, negative for anions).

    Compound annotation:
      - @c found_adduct      : adduct that produced the match, in the
                               notation used by AdductInfo
                               (e.g. @c "M+H;1+", @c "M-H2O+H;1+").
      - @c empirical_formula : empirical formula of the matched compound.
      - @c matching_hmdb_ids : zero or more database IDs that share the
                               formula and mass (multiple isobars are
                               common in HMDB-style databases).

    Back-references:
      - @c matching_index      : index of the hit within the engine's
                                 internal candidate list (useful when
                                 reconstructing scoring details after
                                 export).
      - @c source_feature_index: index of the source feature/peak in the
                                 input map (lets callers traceback to the
                                 raw observation).

    Quality score:
      - @c isotopes_sim_score : isotope-pattern similarity between the
                                feature's observed mass-trace intensities
                                and the theoretical pattern of the
                                candidate compound. Higher = better.

    @ingroup Analysis_ID
    @see AccurateMassSearchEngine
    @see AdductInfo
  */
  class OPENMS_DLLAPI AccurateMassSearchResult
  {
  public:
    /// @brief Default constructor; all numeric fields zero, strings and
    /// vectors empty.
    AccurateMassSearchResult();

    /// @brief Default destructor.
    ~AccurateMassSearchResult();

    /// @brief Copy constructor.
    /// @param[in] other Source result to copy from.
    AccurateMassSearchResult(const AccurateMassSearchResult& other);

    /// @brief Copy-assignment operator.
    /// @param[in] other Source result to copy from.
    AccurateMassSearchResult& operator=(const AccurateMassSearchResult& other);

    /// @brief Observed @em m/z of the source feature (Th).
    double getObservedMZ() const;

    /**
      @brief Set the observed @em m/z of the source feature.
      @param[in] mz Observed @em m/z (Th).
    */
    void setObservedMZ(const double& mz);

    /// @brief Theoretical @em m/z of the matched compound + adduct (Th).
    double getCalculatedMZ() const;

    /**
      @brief Set the theoretical @em m/z of the matched compound + adduct.
      @param[in] mz Theoretical @em m/z of compound + adduct (Th).
    */
    void setCalculatedMZ(const double& mz);

    /// @brief Neutral mass used to query the database (Dalton), back-
    /// calculated from #getObservedMZ assuming #getFoundAdduct and
    /// #getCharge.
    double getQueryMass() const;

    /**
      @brief Set the neutral mass used to query the database.
      @param[in] mass Neutral mass used to query the database (Da).
    */
    void setQueryMass(const double& mass);

    /// @brief Neutral mass of the matched database compound (Dalton).
    double getFoundMass() const;

    /**
      @brief Set the neutral mass of the matched database compound.
      @param[in] mass Neutral mass of the matched compound (Da).
    */
    void setFoundMass(const double& mass);

    /// @brief Ion charge under which this match was attempted (positive
    /// for cations, negative for anions).
    Int getCharge() const;

    /**
      @brief Set the ion charge under which this match was attempted.
      @param[in] ch Signed ion charge (positive for cations).
    */
    void setCharge(const Int& ch);

    /// @brief Signed @em m/z error in ppm: (observed − theoretical) /
    /// theoretical · 1e6.
    double getMZErrorPPM() const;

    /**
      @brief Set the signed @em m/z error in ppm.
      @param[in] ppm Signed @em m/z error (observed − theoretical) in ppm.
    */
    void setMZErrorPPM(const double ppm);

    /// @brief Retention time of the source feature (seconds).
    double getObservedRT() const;

    /**
      @brief Set the retention time of the source feature (seconds).
      @param[in] rt Retention time of the source feature (s).
    */
    void setObservedRT(const double& rt);

    /// @brief Intensity of the source feature (consensus-aggregated when
    /// the input was a ConsensusMap).
    double getObservedIntensity() const;

    /**
      @brief Set the intensity of the source feature.
      @param[in] intensity Intensity of the source feature.
    */
    void setObservedIntensity(const double& intensity);

    /// @brief Per-sample intensities of the source feature when the input
    /// is a ConsensusMap; empty otherwise.
    std::vector<double> getIndividualIntensities() const;

    /**
      @brief Set the per-sample intensities of the source feature.
      @param[in] indiv_ints Per-sample intensities (ConsensusMap input).
    */
    void setIndividualIntensities(const std::vector<double>& indiv_ints);

    /// @brief Index of this hit in the engine's internal candidate list
    /// for the source feature. Useful when reconstructing scoring details
    /// after export.
    Size getMatchingIndex() const;

    /**
      @brief Set the candidate-list index.
      @param[in] idx Index in the engine's internal candidate list.
    */
    void setMatchingIndex(const Size& idx);

    /// @brief Index of the source feature in the input map (back-
    /// reference to the raw observation).
    Size getSourceFeatureIndex() const;

    /**
      @brief Set the source feature index.
      @param[in] idx Index in the input feature/peak map.
    */
    void setSourceFeatureIndex(const Size& idx);

    /// @brief Adduct that produced the match, in AdductInfo notation
    /// (e.g. @c "M+H;1+", @c "M-H2O+H;1+").
    const String& getFoundAdduct() const;

    /**
      @brief Set the adduct that produced the match.
      @param[in] add Adduct in AdductInfo notation (e.g. @c "M+H;1+").
    */
    void setFoundAdduct(const String& add);

    /// @brief Empirical formula of the matched database compound.
    const String& getFormulaString() const;

    /**
      @brief Set the empirical formula of the matched database compound.
      @param[in] ep Empirical formula of the matched compound.
    */
    void setEmpiricalFormula(const String& ep);

    /// @brief All database identifiers (HMDB-style) that share the
    /// matched compound's formula and mass; can be empty when no
    /// identifier mapping is available.
    const std::vector<String>& getMatchingHMDBids() const;

    /**
      @brief Set the matched database identifiers.
      @param[in] ids Database identifiers (HMDB-style).
    */
    void setMatchingHMDBids(const std::vector<String>& ids);

    /// @brief Per-isotopologue intensities of the underlying feature's
    /// mass traces; used as the empirical isotope pattern when computing
    /// #getIsotopesSimScore.
    const std::vector<double>& getMasstraceIntensities() const;

    /**
      @brief Set the underlying feature's per-isotopologue intensities.
      @param[in] intensities Per-isotopologue intensities from the feature.
    */
    void setMasstraceIntensities(const std::vector<double>& intensities);

    /// @brief Similarity between the observed isotope pattern (from
    /// #getMasstraceIntensities) and the theoretical pattern of the
    /// matched compound. Higher values indicate a better match.
    double getIsotopesSimScore() const;

    /**
      @brief Set the isotope-pattern similarity score.
      @param[in] score Isotope-pattern similarity score (higher = better.
    */
    void setIsotopesSimScore(const double& score);

    /// @brief Diagnostic stream output of all fields (not a parseable
    /// serialization).
    friend OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const AccurateMassSearchResult& amsr);

private:
    /// Stored information/results of DB query
    double observed_mz_;                   ///< observed m/z of source feature (Th)
    double theoretical_mz_;                ///< theoretical m/z of matched compound + adduct (Th)
    double searched_mass_;                 ///< neutral mass used to query the database (Da)
    double db_mass_;                       ///< neutral mass of the matched database compound (Da)
    Int charge_;                           ///< charge assumed for this match (positive for cations)
    double mz_error_ppm_;                  ///< signed m/z error (observed − theoretical), ppm
    double observed_rt_;                   ///< RT of source feature (s)
    double observed_intensity_;            ///< intensity of source feature
    std::vector<double> individual_intensities_; ///< per-sample intensities (ConsensusMap input)
    Size matching_index_;                  ///< index in the engine's internal candidate list
    Size source_feature_index_;            ///< index in the input feature/peak map

    String found_adduct_;                  ///< adduct in AdductInfo notation (e.g. "M+H;1+")
    String empirical_formula_;             ///< empirical formula of matched compound
    std::vector<String> matching_hmdb_ids_;///< zero or more DB IDs (multiple isobars share a formula)

    std::vector<double> mass_trace_intensities_; ///< per-isotopologue intensities from the feature
    double isotopes_sim_score_;            ///< observed vs. theoretical isotope-pattern similarity
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
