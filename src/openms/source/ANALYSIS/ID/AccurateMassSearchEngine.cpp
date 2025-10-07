// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Erhan Kenar, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/AccurateMassSearchEngine.h>
#include <OpenMS/CHEMISTRY/AdductInfo.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
#include <OpenMS/SYSTEM/File.h>

#include <numeric>

namespace OpenMS
{
  /// default constructor
  AccurateMassSearchResult::AccurateMassSearchResult() :
  observed_mz_(),
  theoretical_mz_(),
  searched_mass_(),
  db_mass_(),
  charge_(),
  mz_error_ppm_(),
  observed_rt_(),
  observed_intensity_(),
  individual_intensities_(),
  matching_index_(),
  source_feature_index_(),
  found_adduct_(),
  empirical_formula_(),
  matching_hmdb_ids_(),
  mass_trace_intensities_(),
  isotopes_sim_score_(-1.0)
  {
  }

  /// default destructor
  AccurateMassSearchResult::~AccurateMassSearchResult() = default;

  /// copy constructor
  AccurateMassSearchResult::AccurateMassSearchResult(const AccurateMassSearchResult& source) = default;

  /// assignment operator
  AccurateMassSearchResult& AccurateMassSearchResult::operator=(const AccurateMassSearchResult& rhs)
  {
    if (this == &rhs) return *this;

    observed_mz_ = rhs.observed_mz_;
    theoretical_mz_ = rhs.theoretical_mz_;
    searched_mass_ = rhs.searched_mass_;
    db_mass_ = rhs.db_mass_;
    charge_ = rhs.charge_;
    mz_error_ppm_ = rhs.mz_error_ppm_;
    observed_rt_ = rhs.observed_rt_;
    observed_intensity_ = rhs.observed_intensity_;
    individual_intensities_ = rhs.individual_intensities_;
    matching_index_ = rhs.matching_index_;
    source_feature_index_ = rhs.source_feature_index_;
    found_adduct_ = rhs.found_adduct_;
    empirical_formula_ = rhs.empirical_formula_;
    matching_hmdb_ids_ = rhs.matching_hmdb_ids_;
    mass_trace_intensities_ = rhs.mass_trace_intensities_;
    isotopes_sim_score_ = rhs.isotopes_sim_score_;

    return *this;
  }

  double AccurateMassSearchResult::getObservedMZ() const
  {
    return observed_mz_;
  }

  void AccurateMassSearchResult::setObservedMZ(const double& m)
  {
    observed_mz_ = m;
  }

  double AccurateMassSearchResult::getCalculatedMZ() const
  {
    return theoretical_mz_;
  }

  void AccurateMassSearchResult::setCalculatedMZ(const double& m)
  {
    theoretical_mz_ = m;
  }

  double AccurateMassSearchResult::getQueryMass() const
  {
    return searched_mass_;
  }

  void AccurateMassSearchResult::setQueryMass(const double& m)
  {
    searched_mass_ = m;
  }

  double AccurateMassSearchResult::getFoundMass() const
  {
    return db_mass_;
  }

  void AccurateMassSearchResult::setFoundMass(const double& m)
  {
    db_mass_ = m;
  }

  Int AccurateMassSearchResult::getCharge() const
  {
    return charge_;
  }

  void AccurateMassSearchResult::setCharge(const Int& ch)
  {
    charge_ = ch;
  }

  double AccurateMassSearchResult::getMZErrorPPM() const
  {
    return mz_error_ppm_;
  }

  void AccurateMassSearchResult::setMZErrorPPM(const double ppm)
  {
    mz_error_ppm_ = ppm;
  }

  double AccurateMassSearchResult::getObservedRT() const
  {
    return observed_rt_;
  }

  void AccurateMassSearchResult::setObservedRT(const double& rt)
  {
    observed_rt_ = rt;
  }

  double AccurateMassSearchResult::getObservedIntensity() const
  {
    return observed_intensity_;
  }

  void AccurateMassSearchResult::setObservedIntensity(const double& intensity)
  {
    observed_intensity_ = intensity;
  }

  std::vector<double> AccurateMassSearchResult::getIndividualIntensities() const
  {
    return individual_intensities_;
  }

  void AccurateMassSearchResult::setIndividualIntensities(const std::vector<double>& indiv_ints)
  {
    individual_intensities_ = indiv_ints;
  }

  Size AccurateMassSearchResult::getMatchingIndex() const
  {
    return matching_index_;
  }

  void AccurateMassSearchResult::setMatchingIndex(const Size& idx)
  {
    matching_index_ = idx;
  }

  Size AccurateMassSearchResult::getSourceFeatureIndex() const
  {
    return source_feature_index_;
  }

  void AccurateMassSearchResult::setSourceFeatureIndex(const Size& idx)
  {
    source_feature_index_ = idx;
  }

  const String& AccurateMassSearchResult::getFoundAdduct() const
  {
    return found_adduct_;
  }

  void AccurateMassSearchResult::setFoundAdduct(const String& add)
  {
    found_adduct_ = add;
  }

  const String& AccurateMassSearchResult::getFormulaString() const
  {
    return empirical_formula_;
  }

  void AccurateMassSearchResult::setEmpiricalFormula(const String& ep)
  {
    empirical_formula_ = ep;
  }

  const std::vector<String>& AccurateMassSearchResult::getMatchingHMDBids() const
  {
    return matching_hmdb_ids_;
  }

  void AccurateMassSearchResult::setMatchingHMDBids(const std::vector<String>& match_ids)
  {
    matching_hmdb_ids_ = match_ids;
  }

  const std::vector<double>& AccurateMassSearchResult::getMasstraceIntensities() const
  {
    return mass_trace_intensities_;
  }

  void AccurateMassSearchResult::setMasstraceIntensities(const std::vector<double>& mti)
  {
    mass_trace_intensities_ = mti;
  }

  double AccurateMassSearchResult::getIsotopesSimScore() const
  {
    return isotopes_sim_score_;
  }

  void AccurateMassSearchResult::setIsotopesSimScore(const double& sim_score)
  {
    isotopes_sim_score_ = sim_score;
  }

  std::ostream& operator<<(std::ostream& os, const AccurateMassSearchResult& amsr)
  {
    // set maximum precision
    std::streamsize old_precision = os.precision(std::numeric_limits<double>::digits10 + 2);
    os << "observed RT: " << amsr.observed_rt_ << "\n";
    os << "observed intensity: " << amsr.observed_intensity_ << "\n";
    os << "observed m/z: " <<  amsr.observed_mz_ << "\n";
    os << "m/z error ppm: " << amsr.mz_error_ppm_ << "\n";
    os << "charge: " << amsr.charge_ << "\n";
    os << "query mass (searched): " << amsr.searched_mass_ << "\n";
    os << "theoretical (neutral) mass: " << amsr.db_mass_ << "\n";
    os << "matching idx: " << amsr.matching_index_ << "\n";
    os << "emp. formula: " << amsr.empirical_formula_ << "\n";
    os << "adduct: " << amsr.found_adduct_ << "\n";
    os << "matching HMDB ids:";
    for (Size i = 0; i < amsr.matching_hmdb_ids_.size(); ++i)
    {
      os << " " << amsr.matching_hmdb_ids_[i];
    }
    os << "\n";
    os << "isotope similarity score: " << amsr.isotopes_sim_score_ << "\n";

    // restore precision
    os.precision(old_precision);
    return os;
  }

  AccurateMassSearchEngine::AccurateMassSearchEngine() :
    DefaultParamHandler("AccurateMassSearchEngine"),
    ProgressLogger(),
    is_initialized_(false)
  {
    defaults_.setValue("mass_error_value", 5.0, "Tolerance allowed for accurate mass search.");

    defaults_.setValue("mass_error_unit", "ppm", "Unit of mass error (ppm or Da)");
    defaults_.setValidStrings("mass_error_unit", {"ppm", "Da"});

    defaults_.setValue("ionization_mode", "positive", "Positive or negative ionization mode? If 'auto' is used, the first feature of the input map must contain the meta-value 'scan_polarity'. If its missing, the tool will exit with error.");
    defaults_.setValidStrings("ionization_mode", {"positive", "negative", "auto"});

    defaults_.setValue("isotopic_similarity", "false", "Computes a similarity score for each hit (only if the feature exhibits at least two isotopic mass traces).");
    defaults_.setValidStrings("isotopic_similarity", {"false", "true"});

    defaults_.setValue("db:mapping", std::vector<std::string>{"CHEMISTRY/HMDBMappingFile.tsv"}, "Database input file(s), containing three tab-separated columns of mass, formula, identifier. "
                                                                      "If 'mass' is 0, it is re-computed from the molecular sum formula. "
                                                                      "By default CHEMISTRY/HMDBMappingFile.tsv in OpenMS/share is used! If empty, the default will be used.");
    defaults_.setValue("db:struct", std::vector<std::string>{"CHEMISTRY/HMDB2StructMapping.tsv"}, "Database input file(s), containing four tab-separated columns of identifier, name, SMILES, INCHI."
                                                                        "The identifier should match with mapping file. SMILES and INCHI are reported in the output, but not used otherwise. "
                                                                        "By default CHEMISTRY/HMDB2StructMapping.tsv in OpenMS/share is used! If empty, the default will be used.");
    defaults_.setValue("positive_adducts", "CHEMISTRY/PositiveAdducts.tsv", "This file contains the list of potential positive adducts that will be looked for in the database. "
                                                                                 "Edit the list if you wish to exclude/include adducts. "
                                                                                 "By default CHEMISTRY/PositiveAdducts.tsv in OpenMS/share is used.", {"advanced"});
    defaults_.setValue("negative_adducts", "CHEMISTRY/NegativeAdducts.tsv", "This file contains the list of potential negative adducts that will be looked for in the database. "
                                                                                 "Edit the list if you wish to exclude/include adducts. "
                                                                                 "By default CHEMISTRY/NegativeAdducts.tsv in OpenMS/share is used.", {"advanced"});

    defaults_.setValue("use_feature_adducts", "false", "Whether to filter AMS candidates mismatching available feature adduct annotation.");
    defaults_.setValidStrings("use_feature_adducts", {"false", "true"});

    defaults_.setValue("keep_unidentified_masses", "true", "Keep features that did not yield any DB hit.");
    defaults_.setValidStrings("keep_unidentified_masses", {"true", "false"});

    defaults_.setValue("mzTab:exportIsotopeIntensities", "false", "[featureXML input only] Export column with available isotope trace intensities (opt_global_MTint)");
    defaults_.setValidStrings("mzTab:exportIsotopeIntensities", {"false", "true"});

    defaults_.setValue("id_format", "legacy", "Use legacy (ProteinID/PeptideID based storage of metabolomics data) with mzTab-v1.0.0 as output format or novel Identification Data (ID) with mzTab-v2.0.0-M as output format (ID and its MzTab-M output is currently only support for featureXML files).");
    defaults_.setValidStrings("id_format", {"legacy", "ID"});

    defaultsToParam_();
  }

  AccurateMassSearchEngine::~AccurateMassSearchEngine() = default;

/// public methods

  void AccurateMassSearchEngine::queryByMZ(const double& observed_mz, const Int& observed_charge, const String& ion_mode, std::vector<AccurateMassSearchResult>& results, const EmpiricalFormula& observed_adduct) const
  {
    if (!is_initialized_)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "AccurateMassSearchEngine::init() was not called!");
    }

    // Depending on ion_mode_internal_, either positive or negative adducts are used
    std::vector<AdductInfo>::const_iterator it_s, it_e;
    if (ion_mode == "positive")
    {
      it_s = pos_adducts_.begin();
      it_e = pos_adducts_.end();
    }
    else if (ion_mode == "negative")
    {
      it_s = neg_adducts_.begin();
      it_e = neg_adducts_.end();
    }
    else
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Ion mode cannot be set to '") + ion_mode + "'. Must be 'positive' or 'negative'!");
    }

    std::pair<Size, Size> hit_idx;
    for (std::vector<AdductInfo>::const_iterator it = it_s; it != it_e; ++it)
    {
      if (observed_charge != 0 && (std::abs(observed_charge) != std::abs(it->getCharge())))
      { // charge of evidence and adduct must match in absolute terms (absolute, since any FeatureFinder gives only positive charges, even for negative-mode spectra)
        // observed_charge==0 will pass, since we basically do not know its real charge (apparently, no isotopes were found)
        continue;
      }

      if ((observed_adduct != EmpiricalFormula()) && (observed_adduct != it->getEmpiricalFormula()))
      { // If feature has no adduct annotation, method call defaults to empty EF(). If feature is annotated with an adduct, it must match.
        continue;
      }

      // get potential hits as indices in masskey_table
      double neutral_mass = it->getNeutralMass(observed_mz); // calculate mass of uncharged small molecule without adduct mass

      // Our database is just a set of neutral masses (i.e., without adducts)
      // However, given is either an absolute m/z tolerance or a ppm tolerance for the observed m/z
      // We now need an upper bound on the absolute allowed mass difference, given the above tolerance in m/z.
      // The selected candidates then have an mass tolerance which corresponds to the user's m/z tolerance.
      // (the other approach is to pre-compute m/z values for all combinations of adducts, charges and DB entries -- too much)
      double diff_mz;
      // check if mass error window is given in ppm or Da
      if (mass_error_unit_ == "ppm")
      {
        // convert ppm to absolute m/z tolerance for the current candidate
        diff_mz = (observed_mz / 1e6) * mass_error_value_;
      }
      else
      {
        diff_mz = mass_error_value_;
      }
      // convert absolute m/z diff to absolute mass diff
      // What about the adduct?
      // absolute mass error: the adduct itself is irrelevant here since its a constant for both the theoretical and observed mass
      //       ppm tolerance: the diff_mz accounts for it already (heavy adducts lead to larger m/z tolerance)

      // The adduct mass multiplier has to be taken into account when calculating the diff_mass (observed = 228 Da; Multiplier = 2M; theoretical mass = 114 Da)
      // if not the allowed mass error will be the one from 228 Da instead of 114 Da (in this example twice as high).

      double diff_mass = (diff_mz * std::abs(it->getCharge())) / it->getMolMultiplier(); // do not use observed charge (could be 0=unknown)

      searchMass_(neutral_mass, diff_mass, hit_idx);

      //std::cerr << ion_mode_internal_ << " adduct: " << adduct_name << ", " << adduct_mass << " Da, " << query_mass << " qm(against DB), " << charge << " q\n";

      // store information from query hits in AccurateMassSearchResult objects
      for (Size i = hit_idx.first; i < hit_idx.second; ++i)
      {
        // check if DB entry is compatible to the adduct
        if (!it->isCompatible(EmpiricalFormula(mass_mappings_[i].formula)))
        {
          // only written if TOPP tool has --debug
          OPENMS_LOG_DEBUG << "'" << mass_mappings_[i].formula << "' cannot have adduct '" << it->getName() << "'. Omitting.\n";
          continue;
        }

        // compute ppm errors
        double db_mass = mass_mappings_[i].mass;
        double theoretical_mz = it->getMZ(db_mass);
        double error_ppm_mz = Math::getPPM(observed_mz, theoretical_mz); // negative values are allowed!

        AccurateMassSearchResult ams_result;
        ams_result.setObservedMZ(observed_mz);
        ams_result.setCalculatedMZ(theoretical_mz);
        ams_result.setQueryMass(neutral_mass);
        ams_result.setFoundMass(db_mass);
        ams_result.setCharge(std::abs(it->getCharge())); // use theoretical adducts charge (is always valid); native charge might be zero
        ams_result.setMZErrorPPM(error_ppm_mz);
        ams_result.setMatchingIndex(i);
        ams_result.setFoundAdduct(it->getName());
        ams_result.setEmpiricalFormula(mass_mappings_[i].formula);
        ams_result.setMatchingHMDBids(mass_mappings_[i].massIDs);

        results.push_back(ams_result);
      }

    }

    // if result is empty, add a 'not-found' indicator if empty hits should be stored
    if (results.empty() && keep_unidentified_masses_)
    {
      AccurateMassSearchResult ams_result;
      ams_result.setObservedMZ(observed_mz);
      ams_result.setCalculatedMZ(std::numeric_limits<double>::quiet_NaN());
      ams_result.setQueryMass(std::numeric_limits<double>::quiet_NaN());
      ams_result.setFoundMass(std::numeric_limits<double>::quiet_NaN());
      ams_result.setCharge(observed_charge);
      ams_result.setMZErrorPPM(std::numeric_limits<double>::quiet_NaN());
      ams_result.setMatchingIndex(-1); // this is checked to identify 'not-found'
      ams_result.setFoundAdduct("null");
      ams_result.setEmpiricalFormula("");
      ams_result.setMatchingHMDBids(std::vector<String>(1, "null"));
      results.push_back(ams_result);
    }

    return;
  }

  void AccurateMassSearchEngine::queryByFeature(const Feature& feature, const Size& feature_index, const String& ion_mode, std::vector<AccurateMassSearchResult>& results) const
  {
    if (!is_initialized_)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "AccurateMassSearchEngine::init() was not called!");
    }

    std::vector<AccurateMassSearchResult> results_part;

    bool use_feature_adducts = param_.getValue("use_feature_adducts").toString() == "true";
    if (use_feature_adducts && feature.metaValueExists(Constants::UserParam::DC_CHARGE_ADDUCTS))
    {
      queryByMZ(feature.getMZ(), feature.getCharge(), ion_mode, results_part, EmpiricalFormula(feature.getMetaValue(Constants::UserParam::DC_CHARGE_ADDUCTS)));
    }
    else
    {
      queryByMZ(feature.getMZ(), feature.getCharge(), ion_mode, results_part);
    }

    bool isotope_export = param_.getValue("mzTab:exportIsotopeIntensities").toString() == "true";

    for (Size hit_idx = 0; hit_idx < results_part.size(); ++hit_idx)
    {
      results_part[hit_idx].setObservedRT(feature.getRT());
      results_part[hit_idx].setSourceFeatureIndex(feature_index);
      results_part[hit_idx].setObservedIntensity(feature.getIntensity());

      std::vector<double> mti;
      if (isotope_export)
      {
          if (feature.metaValueExists("masstrace_intensity"))
          {
            mti = feature.getMetaValue("masstrace_intensity");
          }
        results_part[hit_idx].setMasstraceIntensities(mti);
      }

      // append
      results.push_back(results_part[hit_idx]);
    }
  }

  void AccurateMassSearchEngine::queryByConsensusFeature(const ConsensusFeature& cfeat, const Size& cf_index, const Size& number_of_maps, const String& ion_mode, std::vector<AccurateMassSearchResult>& results) const
  {
    if (!is_initialized_)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "AccurateMassSearchEngine::init() was not called!");
    }
    results.clear();
    // get hits
    queryByMZ(cfeat.getMZ(), cfeat.getCharge(), ion_mode, results);

    // collect meta data:
    // intensities for all maps as given in handles; 0 if no handle is present for a map
    const ConsensusFeature::HandleSetType& ind_feats(cfeat.getFeatures()); // sorted by MapIndices
    ConsensusFeature::const_iterator f_it = ind_feats.begin();
    std::vector<double> tmp_f_ints;
    for (Size map_idx = 0; map_idx < number_of_maps; ++map_idx)
    {
      if (f_it != ind_feats.end() && map_idx == f_it->getMapIndex())
      {
        tmp_f_ints.push_back(f_it->getIntensity());
        ++f_it;
      }
      else
      {
        tmp_f_ints.push_back(0.0);
      }
    }

    // augment all hits with meta data
    for (Size hit_idx = 0; hit_idx < results.size(); ++hit_idx)
    {
      results[hit_idx].setObservedRT(cfeat.getRT());
      results[hit_idx].setSourceFeatureIndex(cf_index);
      // results_part[hit_idx].setObservedIntensity(cfeat.getIntensity());
      results[hit_idx].setIndividualIntensities(tmp_f_ints);
    }
  }

  void AccurateMassSearchEngine::init()
  {
    // Loads the default mapping file (chemical formulas -> HMDB IDs)
    parseMappingFile_(db_mapping_file_);
    // This loads additional properties like common name, smiles, and inchi key for each HMDB id
    parseStructMappingFile_(db_struct_file_);

    parseAdductsFile_(pos_adducts_fname_, pos_adducts_);
    parseAdductsFile_(neg_adducts_fname_, neg_adducts_);

    is_initialized_ = true;
  }

  void AccurateMassSearchEngine::run(FeatureMap& fmap, MzTabM& mztabm_out) const
  {
    if (!is_initialized_)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "AccurateMassSearchEngine::init() was not called!");
    }

    IdentificationData& id = fmap.getIdentificationData();
    IdentificationData::InputFileRef file_ref;
    IdentificationData::ScoreTypeRef mass_error_ppm_score_ref;
    IdentificationData::ScoreTypeRef mass_error_Da_score_ref;
    IdentificationData::ProcessingStepRef step_ref;

    StringList ms_run_paths;
    fmap.getPrimaryMSRunPath(ms_run_paths);

    // set identifier for FeatureMap if missing (mandatory for OMS output)
    if (fmap.getIdentifier().empty())
    {
      fmap.setIdentifier(File::basename(ms_run_paths[0]));
    }

    // check ion_mode
    String ion_mode_internal(ion_mode_);
    if (ion_mode_ == "auto")
    {
      ion_mode_internal = resolveAutoMode_(fmap);
    }

    // register input file
    IdentificationData::InputFile file(ms_run_paths[0]);
    file_ref = id.registerInputFile(file);
    std::vector<IdentificationData::InputFileRef> file_refs;
    file_refs.emplace_back(file_ref);

    // add previous DataProcessingStep(s) from FeatureMap
    auto data_processing = fmap.getDataProcessing();
    for (const auto& it : data_processing)
    {
      // software
      IdentificationData::ProcessingSoftware sw(it.getSoftware().getName(), it.getSoftware().getVersion());
      // transfer previous metadata
      sw.addMetaValues(it);
      IdentificationDataInternal::ProcessingSoftwareRef sw_ref = id.registerProcessingSoftware(sw);
      // ProcessingStep: software, input_file_refs, data_time, actions
      IdentificationData::ProcessingStep step(sw_ref, file_refs, it.getCompletionTime(), it.getProcessingActions());
      step_ref = id.registerProcessingStep(step);
      id.setCurrentProcessingStep(step_ref);
    }

    // add information about current tool
    // register a score type
    IdentificationData::ScoreType mass_error_ppm_score("MassErrorPPMScore", false);
    mass_error_ppm_score_ref = id.registerScoreType(mass_error_ppm_score);
    IdentificationData::ScoreType mass_error_Da_score("MassErrorDaScore", false);
    mass_error_Da_score_ref = id.registerScoreType(mass_error_Da_score);

    // add the same score_refs to the ProcessingSoftware - to reference the Software with the
    // ObservationMatch - the order is important - the most important score first.
    std::vector<IdentificationDataInternal::ScoreTypeRef> assigned_scores{mass_error_ppm_score_ref, mass_error_Da_score_ref};

    // register software (connected to score)
    // CVTerm will be set in mztab-m based on the name
    // if the name is not available in PSI-OBO "analysis software" will be used.
    IdentificationData::ProcessingSoftware sw("AccurateMassSearch", VersionInfo::getVersion(), assigned_scores);
    sw.setMetaValue("reliability", "2");
    IdentificationData::ProcessingSoftwareRef sw_ref = id.registerProcessingSoftware(sw);

    // all supported search settings
    IdentificationData::DBSearchParam search_param;
    search_param.database = database_name_;
    search_param.database_version = database_version_;
    search_param.setMetaValue("database_location", database_location_);

    search_param.precursor_mass_tolerance = this->mass_error_value_;
    search_param.precursor_tolerance_ppm = this->mass_error_unit_ == "ppm" ? true : false;
    IdentificationData::SearchParamRef search_param_ref = id.registerDBSearchParam(search_param);

    // file has been processed by software performing a specific processing action.
    std::set<DataProcessing::ProcessingAction> actions;
    actions.insert(DataProcessing::IDENTIFICATION);
    IdentificationData::ProcessingStep step(sw_ref, file_refs, DateTime::now(), actions);
    step_ref = id.registerProcessingStep(step, search_param_ref);
    id.setCurrentProcessingStep(step_ref); // add the new step

    // map for storing overall results
    QueryResultsTable overall_results;
    Size dummy_count(0);
    for (Size i = 0; i < fmap.size(); ++i)
    {
      std::vector<AccurateMassSearchResult> query_results = extractQueryResults_(fmap[i], i, ion_mode_internal, dummy_count);
      if (query_results.empty())
      {
        continue;
      }
      overall_results.push_back(query_results);

      addMatchesToID_(id, query_results, file_ref, mass_error_ppm_score_ref, mass_error_Da_score_ref, step_ref, fmap[i]); // MztabM
    }

    // filter FeatureMap to only have entries with an PrimaryID attached
    if (!keep_unidentified_masses_)
    {
      fmap.erase(std::remove_if(fmap.begin(), fmap.end(), [](const Feature& f){ return !f.hasPrimaryID(); }), fmap.end());
    }

    // add the identification data to the featureXML
    // to allow featureXML export (without the use of legacy_ID)
    // been transferred from the previous data stored within
    // the feature.
    IdentificationDataConverter::exportFeatureIDs(fmap, false);

    if (fmap.empty())
    {
      OPENMS_LOG_INFO << "FeatureMap was empty! No hits found!" << std::endl;
    }
    else
    { // division by 0 if used on empty fmap
      OPENMS_LOG_INFO << "\nFound " << (overall_results.size() - dummy_count) << " matched masses (with at least one hit each)\nfrom " << fmap.size() << " features\n  --> " << (overall_results.size()-dummy_count)*100/fmap.size() << "% explained" << std::endl;
    }

    exportMzTabM_(fmap, mztabm_out);

    return;
  }

  void AccurateMassSearchEngine::run(FeatureMap& fmap, MzTab& mztab_out) const
  {
    if (!is_initialized_)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "AccurateMassSearchEngine::init() was not called!");
    }

    StringList ms_run_paths;
    fmap.getPrimaryMSRunPath(ms_run_paths);

    // check ion_mode
    String ion_mode_internal(ion_mode_);
    if (ion_mode_ == "auto")
    {
      ion_mode_internal = resolveAutoMode_(fmap);
    }

    // corresponding file locations
    std::vector<String> file_locations;
    if (!ms_run_paths.empty()) // if the file location is not available it will be set to UNKNOWN by MzTab
    {
      file_locations.emplace_back(ms_run_paths[0]);
    }

    // map for storing overall results
    QueryResultsTable overall_results;
    Size dummy_count(0);
    for (Size i = 0; i < fmap.size(); ++i)
    {
      std::vector<AccurateMassSearchResult> query_results = extractQueryResults_(fmap[i], i, ion_mode_internal, dummy_count);
      if (query_results.empty())
      {
        continue;
      }
      overall_results.push_back(query_results);

      annotate_(query_results, fmap[i]);
    }

    // filter FeatureMap to only have entries with an identification
    if (!keep_unidentified_masses_)
    {
      fmap.erase(std::remove_if(fmap.begin(), fmap.end(), [](Feature f){ return f.getPeptideIdentifications().size() == 0; }), fmap.end());
    }

    // add dummy ProteinIdentification which is required to keep PeptideHits alive during store()
    fmap.getProteinIdentifications().resize(fmap.getProteinIdentifications().size() + 1);
    fmap.getProteinIdentifications().back().setIdentifier("AccurateMassSearchEngine");
    fmap.getProteinIdentifications().back().setSearchEngine("AccurateMassSearch");
    fmap.getProteinIdentifications().back().setDateTime(DateTime().now());

    if (fmap.empty())
    {
      OPENMS_LOG_INFO << "FeatureMap was empty! No hits found!" << std::endl;
    }
    else
    { // division by 0 if used on empty fmap
      OPENMS_LOG_INFO << "\nFound " << (overall_results.size() - dummy_count) << " matched masses (with at least one hit each)\nfrom " << fmap.size() << " features\n  --> " << (overall_results.size()-dummy_count)*100/fmap.size() << "% explained" << std::endl;
    }

    exportMzTab_(overall_results, 1, mztab_out, file_locations);

    return;
  }

  void AccurateMassSearchEngine::addMatchesToID_(
    IdentificationData& id,
    const std::vector<AccurateMassSearchResult>& amr,
    const IdentificationData::InputFileRef& file_ref,
    const IdentificationData::ScoreTypeRef& mass_error_ppm_score_ref,
    const IdentificationData::ScoreTypeRef& mass_error_Da_score_ref,
    const IdentificationData::ProcessingStepRef& step_ref,
    BaseFeature& f) const
  {
    // register feature as search item associated with input file
    IdentificationData::Observation obs(String(f.getUniqueId()), file_ref, f.getRT(), f.getMZ());
    auto obs_ref = id.registerObservation(obs);

    for (const AccurateMassSearchResult& r : amr)
    {
      for (Size i = 0; i < r.getMatchingHMDBids().size(); ++i)
      {
        if (!hmdb_properties_mapping_.count(r.getMatchingHMDBids()[i]))
        {
          throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("DB entry '") + r.getMatchingHMDBids()[i] + "' not found in struct file!");
        }
        // get name from index 0 (2nd column in structMapping file)
        HMDBPropsMapping::const_iterator entry = hmdb_properties_mapping_.find(r.getMatchingHMDBids()[i]);
        if  (entry == hmdb_properties_mapping_.end())
        {
          throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("DB entry '") + r.getMatchingHMDBids()[i] + "' found in struct file but missing in mapping file!");
        }

        double mass_error_Da = r.getObservedMZ() - r.getCalculatedMZ();
        double mass_error_ppm =  r.getMZErrorPPM();

        std::map<IdentificationDataInternal::ScoreTypeRef, double> scores{{mass_error_ppm_score_ref, mass_error_ppm},
                                                                          {mass_error_Da_score_ref, mass_error_Da}
                                                                         };
        IdentificationDataInternal::AppliedProcessingStep applied_processing_step(step_ref, scores);
        IdentificationDataInternal::AppliedProcessingSteps applied_processing_steps;
        applied_processing_steps.emplace_back(applied_processing_step);

        // register compound
        const String& name = entry->second[0];
        const String& smiles = entry->second[1];
        const String& inchi_key = entry->second[2];
        std::vector<String> names = {name}; // to fit legacy format - MetaValue
        std::vector<String> identifiers = {r.getMatchingHMDBids()[i]}; // to fit legacy format - MetaValue
        IdentificationData::IdentifiedCompound compound(r.getMatchingHMDBids()[i],
                                                        EmpiricalFormula(r.getFormulaString()),
                                                        name,
                                                        smiles,
                                                        inchi_key,
                                                        applied_processing_steps);

        auto compound_ref = id.registerIdentifiedCompound(compound); // if already in DB -> NOP

        // compound-feature match
        IdentificationData::ObservationMatch match(compound_ref, obs_ref, r.getCharge());
        match.addScore(mass_error_ppm_score_ref, mass_error_ppm, step_ref);
        match.addScore(mass_error_Da_score_ref, mass_error_Da, step_ref);
        match.setMetaValue("identifier", identifiers);
        match.setMetaValue("description", names);
        match.setMetaValue("modifications", r.getFoundAdduct());
        match.setMetaValue("chemical_formula", r.getFormulaString());
        match.setMetaValue("mz_error_ppm", mass_error_ppm);
        match.setMetaValue("mz_error_Da", mass_error_Da);

        // add adduct to the ObservationMatch
        String adduct = r.getFoundAdduct(); // M+Na;1+
        if (!adduct.empty() && adduct != "null")
        {
          AdductInfo ainfo = AdductInfo::parseAdductString(adduct);
          auto adduct_ref = id.registerAdduct(ainfo);
          match.adduct_opt = adduct_ref;
        }

        // register ObservationMatch
        auto obs_match_ref = id.registerObservationMatch(match);
        IdentificationData::IdentifiedMolecule molecule(compound_ref);
        // add to Feature (set PrimaryID to add a reference to a specific molecule)
        f.setPrimaryID(molecule);
        f.addIDMatch(obs_match_ref);
      }
    }
  }

  void AccurateMassSearchEngine::annotate_(const std::vector<AccurateMassSearchResult>& amr, BaseFeature& f) const
  {
    f.getPeptideIdentifications().resize(f.getPeptideIdentifications().size() + 1);
    f.getPeptideIdentifications().back().setIdentifier(search_engine_identifier);
    for (const AccurateMassSearchResult& result : amr)
    {
      PeptideHit hit;
      hit.setMetaValue("identifier", result.getMatchingHMDBids());
      StringList names;
      for (Size i = 0; i < result.getMatchingHMDBids().size(); ++i)
      { // mapping ok?
        if (!hmdb_properties_mapping_.count(result.getMatchingHMDBids()[i]))
        {
          throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("DB entry '") + result.getMatchingHMDBids()[i] + "' not found in struct file!");
        }
        // get name from index 0 (2nd column in structMapping file)
        HMDBPropsMapping::const_iterator entry = hmdb_properties_mapping_.find(result.getMatchingHMDBids()[i]);
        if  (entry == hmdb_properties_mapping_.end())
        {
          throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("DB entry '") + result.getMatchingHMDBids()[i] + "' found in struct file but missing in mapping file!");
        }
        names.push_back(entry->second[0]);
      }
      hit.setCharge(result.getCharge());
      hit.setMetaValue("description", names);
      hit.setMetaValue("modifications", result.getFoundAdduct());
      hit.setMetaValue("chemical_formula", result.getFormulaString());
      hit.setMetaValue("mz_error_ppm", result.getMZErrorPPM());
      hit.setMetaValue("mz_error_Da", result.getObservedMZ() - result.getCalculatedMZ());
      f.getPeptideIdentifications().back().insertHit(hit);
    }
  }

  void AccurateMassSearchEngine::run(ConsensusMap& cmap, MzTab& mztab_out)  const
  {
    if (!is_initialized_)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "AccurateMassSearchEngine::init() was not called!");
    }

    String ion_mode_internal(ion_mode_);
    if (ion_mode_ == "auto")
    {
      ion_mode_internal = resolveAutoMode_(cmap);
    }

    ConsensusMap::ColumnHeaders fd_map = cmap.getColumnHeaders();
    Size num_of_maps = fd_map.size();

    // corresponding file locations
    std::vector<String> file_locations;
    for (const auto& fd : fd_map)
    {
      file_locations.emplace_back(fd.second.filename);
    }

    // map for storing overall results
    QueryResultsTable overall_results;
    for (Size i = 0; i < cmap.size(); ++i)
    {
      std::vector<AccurateMassSearchResult> query_results;
      queryByConsensusFeature(cmap[i], i, num_of_maps, ion_mode_internal, query_results);
      annotate_(query_results, cmap[i]);
      overall_results.push_back(query_results);
    }
    // add dummy protein identification which is required to keep peptidehits alive during store()
    cmap.getProteinIdentifications().resize(cmap.getProteinIdentifications().size() + 1);
    cmap.getProteinIdentifications().back().setIdentifier(search_engine_identifier);
    cmap.getProteinIdentifications().back().setSearchEngine(search_engine_identifier);
    cmap.getProteinIdentifications().back().setDateTime(DateTime().now());

    exportMzTab_(overall_results, num_of_maps, mztab_out, file_locations);
    return;
  }

  // FeatureMap with IdentificationData attached!
  void AccurateMassSearchEngine::exportMzTabM_(const FeatureMap& fmap, MzTabM& mztabm_out) const
  {
    mztabm_out = MzTabM::exportFeatureMapToMzTabM(fmap);
  }

  void AccurateMassSearchEngine::exportMzTab_(const QueryResultsTable& overall_results, const Size number_of_maps, MzTab& mztab_out, const std::vector<String>& file_locations) const
  {
    if (overall_results.empty())
    {
      return;
    }

    MzTabMetaData md = mztab_out.getMetaData();

    // may contain quantification data so we choose quantification
    md.mz_tab_type.fromCellString("Quantification");

    // we don't report assay so we mark this as a summary file
    md.mz_tab_mode.fromCellString("Summary");

    md.description.fromCellString("Result summary from accurate mass search.");

    // Set meta data
    MzTabParameter search_engine_score;
    search_engine_score.fromCellString("[,,MassErrorPPMScore,]");

    md.smallmolecule_search_engine_score[1] = search_engine_score;

    for (size_t i = 0; i != file_locations.size(); ++i)
    {
      MzTabMSRunMetaData run_md;
      run_md.location.set(file_locations[i]);
      md.ms_run[i + 1] = run_md; // +1 due to start at 1 in MzTab
    }

    // do not use overall_results.begin()->at(0).getIndividualIntensities().size(); since the first entry might be empty (no hit)
    Size n_study_variables = number_of_maps;

    for (Size i = 0; i != n_study_variables; ++i)
    {
      MzTabStudyVariableMetaData sv_md;
      sv_md.description.fromCellString("Accurate mass search result file.");
      md.study_variable[i + 1] = sv_md;
    }

    mztab_out.setMetaData(md);

    // iterate the overall results table
    MzTabSmallMoleculeSectionRows all_sm_rows;

    Size id_group(1);

    std::map<String, UInt> adduct_stats; // adduct --> # occurrences
    std::map<String, std::set<Size> > adduct_stats_unique; // adduct --> # occurrences (count each feature only once)

    bool isotope_export = param_.getValue("mzTab:exportIsotopeIntensities").toString() == "true";

    for (QueryResultsTable::const_iterator tab_it = overall_results.begin(); tab_it != overall_results.end(); ++tab_it)
    {
      for (Size hit_idx = 0; hit_idx < tab_it->size(); ++hit_idx)
      {
        std::vector<String> matching_ids = (*tab_it)[hit_idx].getMatchingHMDBids();
        // iterate over multiple IDs, generate a new row for each one
        for (Size id_idx = 0; id_idx < matching_ids.size(); ++id_idx)
        {
          MzTabSmallMoleculeSectionRow mztab_row_record;
          // set the identifier field
          String hid_temp = matching_ids[id_idx];
          bool db_hit = (hid_temp != "null");
          if (db_hit)
          {
            MzTabString hmdb_id;
            hmdb_id.set(hid_temp);
            std::vector<MzTabString> hmdb_id_dummy;
            hmdb_id_dummy.push_back(hmdb_id);
            MzTabStringList string_dummy_list;
            string_dummy_list.set(hmdb_id_dummy);
            mztab_row_record.identifier = string_dummy_list;

            // set the chemical formula field
            MzTabString chem_form;
            String form_temp = (*tab_it)[hit_idx].getFormulaString();
            chem_form.set(form_temp);

            mztab_row_record.chemical_formula = chem_form;

            HMDBPropsMapping::const_iterator entry = hmdb_properties_mapping_.find(hid_temp);

            // set the smiles field
            String smi_temp = entry->second[1]; // extract SMILES from struct mapping file
            MzTabString smi_string;
            smi_string.set(smi_temp);

            mztab_row_record.smiles = smi_string;

            // set the inchi_key field
            String inchi_temp = entry->second[2]; // extract INCHIKEY from struct mapping file
            MzTabString inchi_key;
            inchi_key.set(inchi_temp);

            mztab_row_record.inchi_key = inchi_key;

            // set description field (we use it for the common name of the compound)
            MzTabString common_name;
            common_name.set(entry->second[0]);
            mztab_row_record.description = common_name;

            // set the calc_mass_to_charge field (theoretical mass)
            MzTabDouble mass_to_charge;
            mass_to_charge.set((*tab_it)[hit_idx].getCalculatedMZ());
            mztab_row_record.calc_mass_to_charge = mass_to_charge;

            // set charge field
            MzTabInteger mcharge;
            mcharge.set((*tab_it)[hit_idx].getCharge());
            mztab_row_record.charge = mcharge;
          }

          // experimental RT, m/z, database field and version, search engine and (null) score is also set if no db entry was matched
          // set RT field
          MzTabDouble rt_temp;
          rt_temp.set((*tab_it)[hit_idx].getObservedRT());
          std::vector<MzTabDouble> rt_temp3(1, rt_temp);
          MzTabDoubleList observed_rt;
          observed_rt.set(rt_temp3);
          mztab_row_record.retention_time = observed_rt;

          MzTabDouble exp_mass_to_charge;
          exp_mass_to_charge.set((*tab_it)[hit_idx].getObservedMZ());
          mztab_row_record.exp_mass_to_charge = exp_mass_to_charge;

          // set database field
          String dbname_temp = database_name_;
          MzTabString dbname;
          dbname.set(dbname_temp);
          mztab_row_record.database = dbname;

          // set database_version field
          String dbver_temp = database_version_;
          MzTabString dbversion;
          dbversion.set(dbver_temp);
          mztab_row_record.database_version = dbversion;

          MzTabParameterList search_engines;
          search_engines.fromCellString("[,,AccurateMassSearch,]");
          mztab_row_record.search_engine = search_engines;

          // same score for all files since it used the mass-to-charge of the ConsensusFeature
          // for identification -> set as best_search_engine_score
          mztab_row_record.best_search_engine_score[1] = MzTabDouble((*tab_it)[hit_idx].getMZErrorPPM());

          // set search_engine_score per hit -> null
          MzTabDouble null_score;
          for (size_t i = 0; i != number_of_maps; ++i) // start from one since it is already filled.
          {
            mztab_row_record.search_engine_score_ms_run[1][i] = null_score;
          }

          // check if we deal with a feature or consensus feature
          std::vector<double> indiv_ints(tab_it->at(hit_idx).getIndividualIntensities());
          std::vector<MzTabDouble> int_temp3;

          bool single_intensity = (indiv_ints.empty());
          if (single_intensity)
          {
            double int_temp((*tab_it)[hit_idx].getObservedIntensity());
            MzTabDouble int_temp2;
            int_temp2.set(int_temp);
            int_temp3.push_back(int_temp2);
          }
          else
          {
            for (Size ii = 0; ii < indiv_ints.size(); ++ii)
            {
              double int_temp(indiv_ints[ii]);
              MzTabDouble int_temp2;
              int_temp2.set(int_temp);
              int_temp3.push_back(int_temp2);
            }
          }

          for (Size i = 0; i != int_temp3.size(); ++i)
          {
            mztab_row_record.smallmolecule_abundance_study_variable[i + 1] = int_temp3[i];
          }

          // set smallmolecule_abundance_stdev_sub; not applicable for a single feature intensity, however must be filled. Otherwise, the mzTab export fails.
          MzTabDouble stdev_temp;
          stdev_temp.set(0.0);
          std::vector<MzTabDouble> stdev_temp3;

          if (indiv_ints.empty())
          {
            stdev_temp3.push_back(stdev_temp);
          }
          else
          {
            for (Size ii = 0; ii < indiv_ints.size(); ++ii)
            {
              stdev_temp3.push_back(stdev_temp);
            }
          }

          for (Size i = 0; i != stdev_temp3.size(); ++i)
          {
            mztab_row_record.smallmolecule_abundance_stdev_study_variable[i + 1] = stdev_temp3[i];
          }

          // set smallmolecule_abundance_std_error_sub; not applicable for a single feature intensity, however must be filled. Otherwise, the mzTab export fails.
          MzTabDouble stderr_temp2;
          stderr_temp2.set(0.0);
          std::vector<MzTabDouble> stderr_temp3;

          if (indiv_ints.empty())
          {
            stderr_temp3.push_back(stderr_temp2);
          }
          else
          {
            for (Size ii = 0; ii < indiv_ints.size(); ++ii)
            {
              stderr_temp3.push_back(stderr_temp2);
            }
          }

          for (Size i = 0; i != stderr_temp3.size(); ++i)
          {
            mztab_row_record.smallmolecule_abundance_std_error_study_variable[i + 1] = stderr_temp3[i];
          }

          // optional columns:
          std::vector<MzTabOptionalColumnEntry> optionals;

          // ppm error
          MzTabString ppmerr;
          if (db_hit)
          {
            ppmerr.set(String((*tab_it)[hit_idx].getMZErrorPPM()));
          }
          MzTabOptionalColumnEntry col0;
          col0.first = "opt_global_mz_ppm_error";
          col0.second = ppmerr;
          optionals.push_back(col0);

          // set found adduct ion
          MzTabString addion;
          if (db_hit)
          {
            String addion_temp((*tab_it)[hit_idx].getFoundAdduct());
            addion.set(addion_temp);
            ++adduct_stats[addion_temp]; // just some stats
            adduct_stats_unique[addion_temp].insert(id_group); // stats ...
          }
          MzTabOptionalColumnEntry col1;
          col1.first = "opt_global_adduct_ion";
          col1.second = addion;
          optionals.push_back(col1);

          // set isotope similarity score
          MzTabString sim_score;
          if (db_hit)
          {
            double sim_score_temp((*tab_it)[hit_idx].getIsotopesSimScore());
            std::stringstream read_in;
            read_in << sim_score_temp;
            String sim_score_temp2(read_in.str());
            sim_score.set(sim_score_temp2);
          }

          MzTabOptionalColumnEntry col2;
          col2.first = "opt_global_isosim_score";
          col2.second = sim_score;
          optionals.push_back(col2);

          // mass trace intensities (use NULL if not present)
          if (isotope_export)
          {
              MzTabString trace_int; // implicitly NULL

              std::vector<double> mt_int = (*tab_it)[hit_idx].getMasstraceIntensities();
              std::vector<std::string> mt_int_strlist;
              std::transform(std::begin(mt_int),
                             std::end(mt_int),
                             std::back_inserter(mt_int_strlist),
                             [](double d) { return std::to_string(d); }
              );

              String mt_int_str = ListUtils::concatenate(mt_int_strlist, ",");

              MzTabOptionalColumnEntry col_mt;
              col_mt.first = String("opt_global_MTint");
              col_mt.second = MzTabString(mt_int_str);
              optionals.push_back(col_mt);
          }

          // set neutral mass
          MzTabString neutral_mass_string;
          if (db_hit)
          {
            String neutral_mass((*tab_it)[hit_idx].getQueryMass());
            neutral_mass_string.fromCellString(neutral_mass);
          }

          MzTabOptionalColumnEntry col3;
          col3.first = "opt_global_neutral_mass";
          col3.second = neutral_mass_string;
          optionals.push_back(col3);

          // set id group; rows with the same id group number originated from the same feature
          String id_group_temp(id_group);
          MzTabString id_group_str;
          id_group_str.set(id_group_temp);
          MzTabOptionalColumnEntry col4;
          col4.first = "opt_global_id_group";
          col4.second = id_group_str;
          optionals.push_back(col4);
          mztab_row_record.opt_ = optionals;
          all_sm_rows.push_back(mztab_row_record);
        }
      }
      ++id_group;
    }

    mztab_out.setSmallMoleculeSectionRows(all_sm_rows);

    // print some adduct stats:
    OPENMS_LOG_INFO << "Hits by adduct: #peaks explained (# matching db entries)'\n";
    for (std::map<String, UInt>::const_iterator it = adduct_stats.begin(); it != adduct_stats.end(); ++it)
    {
      OPENMS_LOG_INFO << "  '" << it->first << "' : " << adduct_stats_unique[it->first].size() << " (" << it->second << ")\n";
    }
    OPENMS_LOG_INFO << std::endl;

  }

  /// protected methods

  void AccurateMassSearchEngine::updateMembers_()
  {
    mass_error_value_ = (double)param_.getValue("mass_error_value");
    mass_error_unit_ = param_.getValue("mass_error_unit").toString();
    ion_mode_ = param_.getValue("ionization_mode").toString();

    iso_similarity_ = param_.getValue("isotopic_similarity").toBool();

    // use defaults if empty for all .tsv files
    db_mapping_file_ = ListUtils::toStringList<std::string>(param_.getValue("db:mapping"));
    if (db_mapping_file_.empty()) db_mapping_file_ = ListUtils::toStringList<std::string>(defaults_.getValue("db:mapping"));
    db_struct_file_ = ListUtils::toStringList<std::string>(param_.getValue("db:struct"));
    if (db_struct_file_.empty()) db_struct_file_ = ListUtils::toStringList<std::string>(defaults_.getValue("db:struct"));

    pos_adducts_fname_ = param_.getValue("positive_adducts").toString();
    neg_adducts_fname_ = param_.getValue("negative_adducts").toString();

    keep_unidentified_masses_ = param_.getValue("keep_unidentified_masses").toBool();
    // database names might have changed, so parse files again before next query
    is_initialized_ = false;

    legacyID_ = param_.getValue("id_format") == "legacy";
  }

/// private methods

  void AccurateMassSearchEngine::parseMappingFile_(const StringList& db_mapping_file)
  {
    mass_mappings_.clear();

    database_location_ = ListUtils::concatenate(db_mapping_file, '|');

    // load map_fname mapping file
    for (StringList::const_iterator it_f = db_mapping_file.begin(); it_f != db_mapping_file.end(); ++it_f)
    {
      String filename = *it_f;
      // load map_fname mapping file
      if (!File::readable(filename))
      {
        // throws Exception::FileNotFound if not found
        filename = File::find(filename);
      }

      String line;
      Size line_count(0);
      std::stringstream str_buf;
      std::istream_iterator<String> eol;

      // OPENMS_LOG_DEBUG << "parsing " << fname << " file..." << std::endl;

      std::ifstream ifs(filename.c_str());
      while (getline(ifs, line))
      {
        line.trim();
        // skip empty lines
        if (line.empty()) continue;
        ++line_count;

        if (line_count == 1)
        {
          std::vector<String> fields;
          line.trim().split('\t', fields);
          if (fields[0] == "database_name")
          {
            database_name_ = fields[1];
            continue;
          }
          else
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Mapping file (") + filename + "') must contain \"database_name\t{NAME}\" as first line.!", line);
          }
        }
        else if (line_count == 2)
        {
          std::vector<String> fields;
          line.trim().split('\t', fields);
          if (fields[0] == "database_version")
          {
            database_version_ = fields[1];
            continue;
          }
          else
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Mapping file (") + filename + "') must contain \"database_version\t{VERSION}\" as second line.!", line);
          }
        }

        str_buf.clear();
        str_buf << line;
        std::istream_iterator<String> istr_it(str_buf);

        Size word_count(0);
        MappingEntry_ entry;

        while (istr_it != eol)
        {
          // OPENMS_LOG_DEBUG << *istr_it << " ";
          if (word_count == 0)
          {
            entry.mass = istr_it->toDouble();
          }
          else if (word_count == 1)
          {
            entry.formula = *istr_it;
            if (entry.mass == 0)
            { // recompute mass from formula
              entry.mass = EmpiricalFormula(entry.formula).getMonoWeight();
              //std::cerr << "mass of " << entry.formula << " is " << entry.mass << "\n";
            }
          }
          else // one or more IDs can follow
          {
            entry.massIDs.push_back(*istr_it);
          }

          ++word_count;
          ++istr_it;
        }
        // OPENMS_LOG_DEBUG << std::endl;

        if (entry.massIDs.empty())
        {
          throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("File '") + filename + "' in line " + line_count + " as '" + line + "' cannot be parsed. Found " + word_count + " entries, expected at least three!");
        }
        mass_mappings_.push_back(entry);
      }
    }
    std::sort(mass_mappings_.begin(), mass_mappings_.end(), CompareEntryAndMass_());

    OPENMS_LOG_INFO << "Read " << mass_mappings_.size() << " entries from mapping file!" << std::endl;

    return;
  }

  void AccurateMassSearchEngine::parseStructMappingFile_(const StringList& db_struct_file)
  {
    hmdb_properties_mapping_.clear();

    for (StringList::const_iterator it_f = db_struct_file.begin(); it_f != db_struct_file.end(); ++it_f)
    {
      String filename = *it_f;

      // load map_fname mapping file
      if (!File::readable(filename))
      {
        // throws Exception::FileNotFound if not found
        filename = File::find(filename);
      }

      std::ifstream ifs(filename.c_str());
      String line;
      // OPENMS_LOG_DEBUG << "parsing " << fname << " file..." << std::endl;

      std::vector<String> parts;
      while (getline(ifs, line))
      {
        line.trim();
        line.split("\t", parts);

        if (parts.size() == 4)
        {
          String hmdb_id_key(parts[0]);

          if (hmdb_properties_mapping_.count(hmdb_id_key))
          {
            throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("File '") + filename + "' in line '" + line + "' cannot be parsed. The ID entry was already used (see above)!");
          }
          std::copy(parts.begin() + 1, parts.end(), std::back_inserter(hmdb_properties_mapping_[hmdb_id_key]));
        }
        else
        {
          throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("File '") + filename + "' in line '" + line + "' cannot be parsed. Expected four entries separated by tab. Found " + parts.size() + " entries!");
        }
      }
    }

    // add a null entry, so mzTab annotation does not discard 'not-found' features
    std::vector<String> dummy_data(3, "null");
    hmdb_properties_mapping_["null"] = dummy_data;

    return;
  }

  void AccurateMassSearchEngine::parseAdductsFile_(const String& filename, std::vector<AdductInfo>& result)
  {
    result.clear();

    String fname = filename;
    // search for mapping file
    if (!File::readable(fname))
    { // throws Exception::FileNotFound if not found
      fname = File::find(filename);
    }
    TextFile tf(fname, true, -1, true); // trim & skip_empty
    for (TextFile::ConstIterator it = tf.begin(); it != tf.end(); ++it)
    {
      result.push_back(AdductInfo::parseAdductString(*it));
    }

    OPENMS_LOG_INFO << "Read " << result.size() << " entries from adduct file '" << fname << "'." << std::endl;

    return;
  }

  void AccurateMassSearchEngine::searchMass_(double neutral_query_mass, double diff_mass, std::pair<Size, Size>& hit_indices) const
  {
    //OPENMS_LOG_INFO << "searchMass: neutral_query_mass=" << neutral_query_mass << " diff_mz=" << diff_mz << " ppm allowed:" << mass_error_value_ << std::endl;

    // binary search for formulas which are within diff_mz distance
    if (mass_mappings_.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "There are no entries found in mass-to-ids mapping file! Aborting... ", "0");
    }

    std::vector<MappingEntry_>::const_iterator lower_it = std::lower_bound(mass_mappings_.begin(), mass_mappings_.end(), neutral_query_mass - diff_mass, CompareEntryAndMass_()); // first element equal or larger
    std::vector<MappingEntry_>::const_iterator upper_it = std::upper_bound(mass_mappings_.begin(), mass_mappings_.end(), neutral_query_mass + diff_mass, CompareEntryAndMass_()); // first element greater than

    Size start_idx = std::distance(mass_mappings_.begin(), lower_it);
    Size end_idx = std::distance(mass_mappings_.begin(), upper_it);

    hit_indices.first = start_idx;
    hit_indices.second = end_idx;

    return;
  }

  double AccurateMassSearchEngine::computeCosineSim_( const std::vector<double>& x, const std::vector<double>& y ) const
  {
    if (x.size() != y.size())
    {
      return 0.0;
    }

    double mixed_sum(0.0);
    double x_squared_sum(0.0);
    double y_squared_sum(0.0);


    for (Size i = 0; i < x.size(); ++i)
    {
      mixed_sum += x[i] * y[i];
      x_squared_sum += x[i] * x[i];
      y_squared_sum += y[i] * y[i];
    }

    double denom(std::sqrt(x_squared_sum) * std::sqrt(y_squared_sum));

    return (denom > 0.0) ? mixed_sum / denom : 0.0;
  }

  double AccurateMassSearchEngine::computeIsotopePatternSimilarity_(const Feature& feat, const EmpiricalFormula& form) const
  {
    Size num_traces = (Size)feat.getMetaValue(Constants::UserParam::NUM_OF_MASSTRACES);
    const Size MAX_THEORET_ISOS(5);

    Size common_size = std::min(num_traces, MAX_THEORET_ISOS);

    // compute theoretical isotope distribution
    IsotopeDistribution iso_dist(form.getIsotopeDistribution(CoarseIsotopePatternGenerator((UInt)common_size)));
    std::vector<double> theoretical_iso_dist;
    std::transform(
      iso_dist.begin(),
      iso_dist.end(),
      back_inserter(theoretical_iso_dist),
      [](const IsotopeDistribution::MassAbundance& p)
      {
        return p.getIntensity();
      });

    // same for observed isotope distribution
    std::vector<double> observed_iso_dist;
    if (num_traces > 0)
    {
      observed_iso_dist = feat.getMetaValue("masstrace_intensity");
    }

    return computeCosineSim_(theoretical_iso_dist, observed_iso_dist);
  }

  std::vector<AccurateMassSearchResult> AccurateMassSearchEngine::extractQueryResults_(const Feature& feature, const Size& feature_index, const String& ion_mode_internal, Size& dummy_count) const
  {
    std::vector<AccurateMassSearchResult> query_results;

    queryByFeature(feature, feature_index, ion_mode_internal, query_results);

    if (query_results.empty())
    {
      return query_results;
    }

    bool is_dummy = (query_results[0].getMatchingIndex() == (Size) - 1);
    if (is_dummy)
      ++dummy_count;

    if (iso_similarity_ && !is_dummy)
    {
      if (!feature.metaValueExists(Constants::UserParam::NUM_OF_MASSTRACES))
      {
        OPENMS_LOG_WARN
        << "Feature does not contain meta value '" << Constants::UserParam::NUM_OF_MASSTRACES << "'. Cannot compute isotope similarity.";
      }
      else if ((Size) feature.getMetaValue(Constants::UserParam::NUM_OF_MASSTRACES) > 1)
      { // compute isotope pattern similarities (do not take the best-scoring one, since it might have really bad ppm or other properties --
        // it is impossible to decide here which one is best
        for (Size hit_idx = 0; hit_idx < query_results.size(); ++hit_idx)
        {
          String emp_formula(query_results[hit_idx].getFormulaString());
          double iso_sim(computeIsotopePatternSimilarity_(feature, EmpiricalFormula(emp_formula)));
          query_results[hit_idx].setIsotopesSimScore(iso_sim);
        }
      }
    }
    return query_results;
  }

} // closing namespace OpenMS
