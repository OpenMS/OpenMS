// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>
///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

typedef std::map<std::string, MRMTransitionGroup< MSChromatogram, OpenSwath::LightTransition > > TransitionGroupMapType;

namespace OpenMS
{
  // Exercise the exact INTERNAL production entry without exposing a public API.
  class MRMFeatureFinderScoringTestAccess
  {
  public:
    static std::vector<std::exception_ptr> scoreBatch(const MRMFeatureFinderScoring& scorer,
                          const std::vector<std::pair<MRMFeatureFinderScoring::MRMTransitionGroupType*, FeatureMap*>>& groups,
                          const TransformationDescription& trafo, const std::vector<OpenSwath::SwathMap>& maps,
                          bool ms1only = false, Size* peak_live_contexts = nullptr)
    {
      return scorer.scorePeakgroupsBatch_(groups, trafo, maps, ms1only, nullptr, peak_live_contexts);
    }
  };
}

namespace
{
  class RecordingSpectrumAccess : public SpectrumAccessOpenMS
  {
  public:
    explicit RecordingSpectrumAccess(const std::shared_ptr<PeakMap>& spectra) : SpectrumAccessOpenMS(spectra) {}
    mutable std::vector<double> requested_rts;

    std::vector<std::size_t> getSpectraByRT(double rt, double delta_rt) const override
    {
      requested_rts.push_back(rt);
      if (!std::isfinite(rt))
      {
        throw std::invalid_argument("Nonfinite RT in recording spectrum access");
      }
      return SpectrumAccessOpenMS::getSpectraByRT(rt, delta_rt);
    }
  };
}

START_TEST(MRMFeatureFinderScoring, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MRMFeatureFinderScoring* ptr = nullptr;
MRMFeatureFinderScoring* nullPointer = nullptr;

START_SECTION(MRMFeatureFinderScoring())
{
  ptr = new MRMFeatureFinderScoring();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~MRMFeatureFinderScoring())
{
  delete ptr;
}
END_SECTION

START_SECTION(void pickExperiment(OpenSwath::SpectrumAccessPtr input, FeatureMap< Feature > &output, OpenSwath::LightTargetedExperiment &transition_exp, TransformationDescription trafo, OpenSwath::SpectrumAccessPtr swath_map, TransitionGroupMapType &transition_group_map))
{
  MRMFeatureFinderScoring ff;
  // Param picker_param = ff.getDefaults();
  // picker_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:method", "legacy"); // old parameters
  // picker_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:peak_width", 40.0); // old parameters
  // ff.setParameters(picker_param);
      
  MRMFeature feature;
  FeatureMap featureFile;
  TransformationDescription trafo;
  std::shared_ptr<PeakMap> swath_map (new PeakMap);
  TransitionGroupMapType transition_group_map;
  MRMFeatureFinderScoring::MRMTransitionGroupType transition_group;

  // Load the chromatograms (mzML) and the meta-information (TraML)
  std::shared_ptr<PeakMap> exp (new PeakMap);
  OpenSwath::LightTargetedExperiment transitions;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("OpenSwath_generic_input.mzML"), *exp);
  {
    TargetedExperiment transition_exp_;
    TraMLFile().load(OPENMS_GET_TEST_DATA_PATH("OpenSwath_generic_input.TraML"), transition_exp_);
    OpenSwathDataAccessHelper::convertTargetedExp(transition_exp_, transitions);
  }

  // Pick features in the experiment
#ifdef USE_SP_INTERFACE
  OpenSwath::SpectrumAccessPtr swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);
  OpenSwath::SpectrumAccessPtr chromatogram_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);
  std::vector< OpenSwath::SwathMap > swath_maps(1);
  swath_maps[0].sptr = swath_ptr;
  ff.pickExperiment(chromatogram_ptr, featureFile, transitions, trafo, swath_maps, transition_group_map);
#else
  ff.pickExperiment(exp, featureFile, transitions, trafo, *swath_map, transition_group_map);
#endif

  // Test the number of features found 
  TEST_EQUAL(transition_group_map.size(), 2)

  ///////////////////////////////////////////////////////////////////////////
  //// Scores for the first group
  transition_group = transition_group_map["tr_gr1"];
  TEST_EQUAL(transition_group.size(), 2)
  TEST_EQUAL(transition_group.getFeatures().size(), 1)
  
  // Look closely at the feature we found in the first group
  feature = transition_group.getFeatures()[0];
  TOLERANCE_ABSOLUTE(0.1);
  TEST_REAL_SIMILAR (feature.getRT(), 3119.092);
  TEST_REAL_SIMILAR (feature.getIntensity(), 3599.98);

  // feature attributes
  TEST_REAL_SIMILAR(feature.getMetaValue("leftWidth" ), 3092.85009765625);
  TEST_REAL_SIMILAR(feature.getMetaValue("rightWidth"), 3151.10009765625);
  TEST_REAL_SIMILAR(feature.getMetaValue("total_xic"), 3680.16);

  // feature scores
  TEST_REAL_SIMILAR(feature.getMetaValue("var_xcorr_coelution"), 0);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_xcorr_shape"), 0.9981834605);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_library_rmsd"), 0.108663236);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_library_corr"), 1);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_elution_model_fit_score"), 0.9854);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_intensity_score"), 0.971);
  TEST_REAL_SIMILAR(feature.getMetaValue("sn_ratio"), 65.2167155356372);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_log_sn_score"), 4.17771580932057);

  TOLERANCE_RELATIVE(1.001);
  TEST_REAL_SIMILAR(feature.getMetaValue("rt_score"), 3118.651968);
  TOLERANCE_ABSOLUTE(0.1);

  ///////////////////////////////////////////////////////////////////////////
  //// Scores for the second group
  transition_group = transition_group_map["tr_gr2"];
  TEST_EQUAL(transition_group.size(), 3)
  TEST_EQUAL(transition_group.getFeatures().size(), 2)
  TEST_EQUAL(featureFile.size(), 3)

  // Look closely at the feature we found in the second group
  feature = transition_group.getFeatures()[0];
  TOLERANCE_ABSOLUTE(0.1);
  TEST_REAL_SIMILAR(feature.getRT(), 3119.092);
  TEST_REAL_SIMILAR(feature.getIntensity(), 1049.65);

  // feature attributes
  TEST_REAL_SIMILAR(feature.getMetaValue("leftWidth" ), 3092.85009765625);
  TEST_REAL_SIMILAR(feature.getMetaValue("rightWidth"), 3147.67993164062);
  TEST_REAL_SIMILAR(feature.getMetaValue("total_xic"), 1610.27);

  // feature scores
  TEST_REAL_SIMILAR(feature.getMetaValue("var_xcorr_coelution"), 7);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_xcorr_shape"), 0.7245);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_library_rmsd"), 0.43566);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_library_corr"), -0.784);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_elution_model_fit_score"), 0.902);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_intensity_score"), 0.642);
  TEST_REAL_SIMILAR(feature.getMetaValue("sn_ratio"), 22.9654496908248);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_log_sn_score"), 3.13399089915908);

  // test legacy parameters
  {
    
    Param picker_param = ff.getDefaults();
    picker_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:method", "legacy"); // old parameters
    picker_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:peak_width", 40.0); // old parameters
    ff.setParameters(picker_param);

    transition_group_map.clear();
    featureFile.clear();
        
    // Pick features in the experiment
#ifdef USE_SP_INTERFACE
    OpenSwath::SpectrumAccessPtr swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);
    OpenSwath::SpectrumAccessPtr chromatogram_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);
    std::vector< OpenSwath::SwathMap > swath_maps(1);
    swath_maps[0].sptr = swath_ptr;
    ff.pickExperiment(chromatogram_ptr, featureFile, transitions, trafo, swath_maps, transition_group_map);
#else
    ff.pickExperiment(exp, featureFile, transitions, trafo, *swath_map, transition_group_map);
#endif

    ///////////////////////////////////////////////////////////////////////////
    //// Scores for the first group
    transition_group = transition_group_map["tr_gr1"];
    TEST_EQUAL(transition_group.size(), 2)
    TEST_EQUAL(transition_group.getFeatures().size(), 1)
    // Look closely at the feature we found in the first group
    feature = transition_group.getFeatures()[0];
    TOLERANCE_ABSOLUTE(0.1);
    TEST_REAL_SIMILAR (feature.getRT(), 3119.092);
    TEST_REAL_SIMILAR (feature.getIntensity(), 3574.23);

    // feature attributes
    TEST_REAL_SIMILAR(feature.getMetaValue("leftWidth" ), 3099.69995117188);
    TEST_REAL_SIMILAR(feature.getMetaValue("rightWidth"), 3147.68);
    TEST_REAL_SIMILAR(feature.getMetaValue("total_xic"), 3680.16);

    ///////////////////////////////////////////////////////////////////////////
    //// Scores for the second group
    // Note: the more accurate (interpolated) S/N estimate now resolves a small,
    // previously-merged nearby peak as its own separate feature (RT ~3148.6,
    // intensity ~59.4), so this group now yields 3 features instead of 2.
    transition_group = transition_group_map["tr_gr2"];
    TEST_EQUAL(transition_group.size(), 3)
    TEST_EQUAL(transition_group.getFeatures().size(), 3)
    TEST_EQUAL(featureFile.size(), 4)
    // Look closely at the feature we found in the second group
    feature = transition_group.getFeatures()[0];
    TOLERANCE_ABSOLUTE(0.1);
    TEST_REAL_SIMILAR(feature.getRT(), 3119.092);
    TEST_REAL_SIMILAR(feature.getIntensity(), 1034.55);

    // feature attributes
    TEST_REAL_SIMILAR(feature.getMetaValue("leftWidth" ), 3099.7);
    TEST_REAL_SIMILAR(feature.getMetaValue("rightWidth"), 3144.25);
    TEST_REAL_SIMILAR(feature.getMetaValue("total_xic"), 1610.27);
    TEST_REAL_SIMILAR(feature.getMetaValue("var_xcorr_coelution"), 2.265);
  }

}
END_SECTION

START_SECTION(void pickExperiment(OpenSwath::SpectrumAccessPtr input, FeatureMap< Feature > &output, OpenSwath::LightTargetedExperiment &transition_exp, TransformationDescription trafo, OpenSwath::SpectrumAccessPtr swath_map, TransitionGroupMapType &transition_group_map))
{ 
  MRMFeatureFinderScoring ff;
  Param ff_param = MRMFeatureFinderScoring().getDefaults();
  Param scores_to_use;
  scores_to_use.setValue("use_uis_scores", "true", "Use UIS scores for peptidoform identification ", {"advanced"});
  scores_to_use.setValidStrings("use_uis_scores", {"true","false"});
  ff_param.insert("Scores:", scores_to_use);
  ff_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:method", "legacy"); // old parameters
  ff_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:peak_width", 40.0); // old parameters
  ff.setParameters(ff_param);

  MRMFeature feature;
  FeatureMap featureFile;
  TransformationDescription trafo;
  std::shared_ptr<PeakMap> swath_map (new PeakMap);
  TransitionGroupMapType transition_group_map;

  // Load the chromatograms (mzML) and the meta-information (TraML)
  std::shared_ptr<PeakMap> exp (new PeakMap);
  OpenSwath::LightTargetedExperiment transitions;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("OpenSwath_generic_input.mzML"), *exp);
  { 
    TargetedExperiment transition_exp_;
    TraMLFile().load(OPENMS_GET_TEST_DATA_PATH("OpenSwath_identification_input.TraML"), transition_exp_);
    OpenSwathDataAccessHelper::convertTargetedExp(transition_exp_, transitions);
  }

  // Pick features in the experiment
#ifdef USE_SP_INTERFACE
  OpenSwath::SpectrumAccessPtr swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);
  OpenSwath::SpectrumAccessPtr chromatogram_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);
  std::vector< OpenSwath::SwathMap > swath_maps(1);
  swath_maps[0].sptr = swath_ptr;
  ff.pickExperiment(chromatogram_ptr, featureFile, transitions, trafo, swath_maps, transition_group_map);
#else
  ff.pickExperiment(exp, featureFile, transitions, trafo, *swath_map, transition_group_map);
#endif

  // Test the number of features found 
  TEST_EQUAL(transition_group_map.size(), 2)

  ///////////////////////////////////////////////////////////////////////////
  //// Scores for the second group
  // Note: the more accurate (interpolated) S/N estimate now resolves a small,
  // previously-merged nearby peak as its own separate feature (RT ~3148.6,
  // intensity ~59.4), so this group now yields 3 features instead of 2.
  MRMFeatureFinderScoring::MRMTransitionGroupType transition_group;
  transition_group = transition_group_map["tr_gr2"];
  TEST_EQUAL(transition_group.size(), 3)
  TEST_EQUAL(transition_group.getFeatures().size(), 3)
  TEST_EQUAL(featureFile.size(), 4)

  // Look closely at the feature we found in the second group
  feature = transition_group.getFeatures()[0];
  TOLERANCE_ABSOLUTE(0.1);
  TEST_REAL_SIMILAR(feature.getRT(), 3119.092);
  TEST_REAL_SIMILAR(feature.getIntensity(), 1034.55);

  // feature attributes
  TEST_REAL_SIMILAR(feature.getMetaValue("leftWidth" ), 3099.7);
  TEST_REAL_SIMILAR(feature.getMetaValue("rightWidth"), 3144.25);
  TEST_REAL_SIMILAR(feature.getMetaValue("total_xic"), 1610.27);

  // feature scores
  TEST_REAL_SIMILAR(feature.getMetaValue("var_xcorr_coelution"), 2.265);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_xcorr_shape"), 0.7245);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_library_rmsd"), 0.43566);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_library_corr"), -0.784);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_elution_model_fit_score"), 0.902);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_intensity_score"), 2.36342573991536);
  TEST_REAL_SIMILAR(feature.getMetaValue("sn_ratio"), 22.9654496908248);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_log_sn_score"), 3.13399089915908);

  // feature identification scores
  TEST_EQUAL(feature.getMetaValue("id_target_transition_names").toStringList()[0], "tr5");
  TEST_EQUAL(feature.getMetaValue("id_target_transition_names").toStringList()[1], "tr2");
  TEST_REAL_SIMILAR(feature.getMetaValue("id_target_ind_log_intensity").toDoubleList()[0], 5.03352);
  TEST_REAL_SIMILAR(feature.getMetaValue("id_target_ind_log_intensity").toDoubleList()[1], 7.92704);
  TEST_REAL_SIMILAR(feature.getMetaValue("id_target_num_transitions"), 2);
  TEST_REAL_SIMILAR(feature.getMetaValue("id_target_ind_xcorr_coelution").toDoubleList()[0], 1);
  TEST_REAL_SIMILAR(feature.getMetaValue("id_target_ind_xcorr_coelution").toDoubleList()[1], 1.66667);
  TEST_REAL_SIMILAR(feature.getMetaValue("id_target_ind_xcorr_shape").toDoubleList()[0], 0.68631);
  TEST_REAL_SIMILAR(feature.getMetaValue("id_target_ind_xcorr_shape").toDoubleList()[1], 0.690494);
  TEST_REAL_SIMILAR(feature.getMetaValue("id_target_ind_log_sn_score").toDoubleList()[0], 1.16692);
  TEST_REAL_SIMILAR(feature.getMetaValue("id_target_ind_log_sn_score").toDoubleList()[1], 4.18437274576997);
  TEST_EQUAL(feature.getMetaValue("id_target_ind_isotope_correlation").toDoubleList().size(), 0);
  TEST_EQUAL(feature.getMetaValue("id_target_ind_isotope_overlap").toDoubleList().size(), 0);
  TEST_EQUAL(feature.getMetaValue("id_target_ind_massdev_score").toDoubleList().size(), 0);

}
END_SECTION

START_SECTION(void mapExperimentToTransitionList(OpenSwath::SpectrumAccessPtr input, OpenSwath::LightTargetedExperiment &transition_exp, TransitionGroupMapType &transition_group_map, TransformationDescription trafo, double rt_extraction_window))
{
  MRMFeatureFinderScoring ff;
  MRMFeature feature;
  FeatureMap featureFile;
  TransformationDescription trafo;
  PeakMap swath_map;
  TransitionGroupMapType transition_group_map;
  MRMFeatureFinderScoring::MRMTransitionGroupType transition_group;

  // Load the chromatograms (mzML) and the meta-information (TraML)
  std::shared_ptr<PeakMap> exp (new PeakMap);
  OpenSwath::LightTargetedExperiment transitions;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("OpenSwath_generic_input.mzML"), *exp);
  {
    TargetedExperiment transition_exp_;
    TraMLFile().load(OPENMS_GET_TEST_DATA_PATH("OpenSwath_generic_input.TraML"), transition_exp_);
    OpenSwathDataAccessHelper::convertTargetedExp(transition_exp_, transitions);
  }

  // Pick features in the experiment
  ff.prepareProteinPeptideMaps_(transitions);
#ifdef USE_SP_INTERFACE
  OpenSwath::SpectrumAccessPtr chromatogram_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);
  ff.mapExperimentToTransitionList(chromatogram_ptr, transitions, transition_group_map, trafo, -1);
#else
  ff.mapExperimentToTransitionList(exp, transitions, transition_group_map, trafo, -1);
#endif

  // Test the number of features found 
  TEST_EQUAL(transition_group_map.size(), 2)

  ///////////////////////////////////////////////////////////////////////////
  //// The first group
  transition_group = transition_group_map["tr_gr1"];
  TEST_EQUAL(transition_group.size(), 2)
  TEST_EQUAL(transition_group.getTransitions().size(), 2)
  TEST_EQUAL(transition_group.getChromatograms().size(), 2)

  TEST_EQUAL(transition_group.hasChromatogram("tr1"), true)
  TEST_EQUAL(transition_group.hasChromatogram("tr2"), true)

  TEST_EQUAL(transition_group.getChromatogram("tr2").getNativeID(), "tr2")
  TEST_EQUAL(transition_group.getTransition("tr2").getNativeID(), "tr2")

  ///////////////////////////////////////////////////////////////////////////
  //// The second group
  transition_group = transition_group_map["tr_gr2"];
  TEST_EQUAL(transition_group.size(), 3)
  TEST_EQUAL(transition_group.getTransitions().size(), 3)
  TEST_EQUAL(transition_group.getChromatograms().size(), 3)

  TEST_EQUAL(transition_group.hasChromatogram("tr3"), true)
  TEST_EQUAL(transition_group.hasChromatogram("tr4"), true)
  TEST_EQUAL(transition_group.hasChromatogram("tr5"), true)
 
  TEST_EQUAL(transition_group.getChromatogram("tr5").getNativeID(), "tr5")
  TEST_EQUAL(transition_group.getTransition("tr5").getNativeID(), "tr5")

}
END_SECTION

START_SECTION( void scorePeakgroups(MRMTransitionGroupType& transition_group, TransformationDescription & trafo, OpenSwath::SpectrumAccessPtr swath_map, FeatureMap& output) ) 
{
  NOT_TESTABLE // tested above
}
END_SECTION

START_SECTION(void prepareProteinPeptideMaps_(OpenSwath::LightTargetedExperiment& transition_exp) )
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION(INTERNAL batch scoring preserves candidate scores and input-group output order)
{
  #ifdef _OPENMP
  const int previous_threads = omp_get_max_threads();
  omp_set_num_threads(1); // deterministic original-order reference, including quality ties
  #endif
  auto chromatograms = std::make_shared<PeakMap>();
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("OpenSwath_generic_input.mzML"), *chromatograms);
  TransformationDescription trafo;
  TEST_TRUE(MRMFeatureFinderScoringTestAccess::scoreBatch(MRMFeatureFinderScoring(), {}, trafo, {}).empty())
  for (int scoring_mode : {0, 1, 2})
  {
    const bool uis = scoring_mode == 1;
    const bool im = scoring_mode == 2;
    TargetedExperiment library;
    TraMLFile().load(uis ? OPENMS_GET_TEST_DATA_PATH("OpenSwath_identification_input.TraML") :
                          OPENMS_GET_TEST_DATA_PATH("OpenSwath_generic_input.TraML"), library);
    OpenSwath::LightTargetedExperiment transitions;
    OpenSwathDataAccessHelper::convertTargetedExp(library, transitions);
    MRMFeatureFinderScoring scorer;
    Param parameters = scorer.getDefaults();
    parameters.setValue("Scores:use_uis_scores", uis ? "true" : "false");
    parameters.setValue("Scores:use_ion_mobility_scores", im ? "true" : "false");
    parameters.setValue("TransitionGroupPicker:PeakPickerChromatogram:method", "legacy");
    parameters.setValue("TransitionGroupPicker:PeakPickerChromatogram:peak_width", 40.0);
    scorer.setParameters(parameters);
    TransitionGroupMapType seed;
    scorer.prepareProteinPeptideMaps_(transitions);
    scorer.mapExperimentToTransitionList(SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(chromatograms),
                                        transitions, seed, trafo, -1);
    MRMTransitionGroupPicker picker;
    Param picker_parameters = parameters.copy("TransitionGroupPicker:", true);
    picker_parameters.setValue("compute_total_mi", "true");
    picker.setParameters(picker_parameters);
    for (auto& [id, group] : seed)
    {
      picker.pickTransitionGroup(group);
    }
    TEST_EQUAL(seed.size(), 2)

    // Deliberately interleave candidate apex RTs across groups and duplicate an
    // entire candidate for exact quality/RT ties. Preserve its full traces.
    auto& first = seed.at("tr_gr1").getFeaturesMuteable();
    auto& second = seed.at("tr_gr2").getFeaturesMuteable();
    MRMFeature first_candidate = first.front();
    MRMFeature second_candidate = second.front();
    first.assign(3, first_candidate);
    second.assign(3, second_candidate);
    first[0].setRT(3130.0);
    first[1].setRT(3105.0);
    first[2].setRT(3130.0);
    second[0].setRT(3120.0);
    second[1].setRT(3100.0);
    second[2].setRT(3130.0);
    for (auto& [id, group] : seed)
    {
      Size candidate_index = 0;
      for (auto& feature : group.getFeaturesMuteable())
      {
        feature.setUniqueId();
        feature.setMetaValue("batch_candidate_index", candidate_index++);
        std::vector<std::string> feature_ids;
        feature.getFeatureIDs(feature_ids);
        for (const auto& feature_id : feature_ids)
        {
          feature.getFeature(feature_id).ensureUniqueId();
        }
      }
    }

    // Small real spectrum-access input exercises DIA scoring and deferred
    // RANKING/REMAINDER scoring, without any cache/prefetch substitute.
    auto spectra = std::make_shared<PeakMap>();
    for (double rt : {3090.0, 3100.0, 3110.0, 3120.0, 3130.0, 3140.0, 3150.0})
    {
      MSSpectrum spectrum;
      spectrum.setRT(rt);
      spectrum.setMSLevel(2);
      for (const auto& transition : transitions.transitions)
      {
        Peak1D peak;
        peak.setMZ(transition.getProductMZ());
        peak.setIntensity(100.0);
        spectrum.push_back(peak);
      }
      spectrum.sortByPosition();
      spectra->addSpectrum(spectrum);
    }
    auto access = std::make_shared<RecordingSpectrumAccess>(spectra);
    std::vector<OpenSwath::SwathMap> maps(1);
    maps[0].sptr = access;
    maps[0].lower = 0;
    maps[0].upper = 2000;

    for (bool separated : {false, true})
    {
      auto scheduled_seed = seed;
      if (separated)
      {
        // Reverse group output order relative to execution, with disjoint RT
        // ranges. Both UIS subsets and ordinary scoring use the real scorer.
        auto& early = scheduled_seed.at("tr_gr1").getFeaturesMuteable();
        auto& late = scheduled_seed.at("tr_gr2").getFeaturesMuteable();
        early[0].setRT(3100.0);
        early[1].setRT(3105.0);
        early[2].setRT(3100.0);
        late[0].setRT(3130.0);
        late[1].setRT(3120.0);
        late[2].setRT(3130.0);
      }
      for (int report_limit : {-1, 0, 1, 2, 10})
      {
        parameters.setValue("stop_report_after_feature", report_limit);
        scorer.setParameters(parameters);
        auto original = scheduled_seed;
        auto batch = scheduled_seed;
        // Include an empty candidate group in the same batch.
        auto empty = seed.at("tr_gr1");
        empty.getFeaturesMuteable().clear();
        FeatureMap original_output, batch_output;
        std::vector<double> expected_remainder_rts;
        std::vector<std::pair<MRMFeatureFinderScoring::MRMTransitionGroupType*, FeatureMap*>> work;
        for (const auto& id : {"tr_gr2", "tr_gr1"})
        {
          access->requested_rts.clear();
          scorer.scorePeakgroups(original.at(id), trafo, maps, original_output);
          if (!uis && !im && report_limit >= 0 && report_limit < 3)
          {
            TEST_EQUAL(access->requested_rts.size(), 3 + report_limit)
            expected_remainder_rts.insert(expected_remainder_rts.end(),
              access->requested_rts.end() - report_limit, access->requested_rts.end());
          }
          work.emplace_back(&batch.at(id), &batch_output);
        }
        work.emplace_back(&empty, &batch_output);
        access->requested_rts.clear();
        Size peak_live_contexts = 0;
        const auto errors = MRMFeatureFinderScoringTestAccess::scoreBatch(scorer, work, trafo, maps, false, &peak_live_contexts);
        if (uis || im || report_limit == -1 || report_limit == 0 || report_limit >= 3)
        {
          // The upfront implementation retained all three contexts, including
          // the empty group. Now disjoint groups never coexist in memory.
          TEST_EQUAL(peak_live_contexts, separated ? 1 : 2)
        }
        else
        {
          // Nonzero REMAINDER work must wait for all candidate raw requests.
          TEST_EQUAL(peak_live_contexts, 2)
        }
        TEST_EQUAL(errors.size(), work.size())
        for (const auto& error : errors)
        {
          TEST_FALSE(static_cast<bool>(error))
          if (error) { std::rethrow_exception(error); }
        }
        TEST_EQUAL(batch_output.size(), original_output.size())
        for (Size i = 0; i < original_output.size(); ++i)
        {
          auto actual = batch_output[i];
          auto expected = original_output[i];
          // Mixed detecting subsets intentionally recreate MRMFeatures, including
          // their generated top-level IDs, in both scoring entrypoints.
          actual.clearUniqueId();
          expected.clearUniqueId();
          TEST_TRUE(actual == expected)
        }
        for (const auto& [id, group] : original)
        {
          const auto& result = batch.at(id);
          TEST_EQUAL(result.getTransitions().size(), group.getTransitions().size())
          TEST_EQUAL(result.getChromatograms().size(), group.getChromatograms().size())
          TEST_TRUE(result.getChromatograms() == group.getChromatograms())
          TEST_EQUAL(result.getFeatures().size(), group.getFeatures().size())
          for (Size i = 0; i < group.getFeatures().size(); ++i)
          {
            auto expected = group.getFeatures()[i];
            auto actual = result.getFeatures()[i];
            actual.clearUniqueId();
            expected.clearUniqueId();
            TEST_TRUE(actual == expected)
            TEST_TRUE(actual.getFeatures() == expected.getFeatures())
            TEST_REAL_SIMILAR(actual.getScores().raw_rt_score, expected.getScores().raw_rt_score)
            TEST_REAL_SIMILAR(actual.getScores().massdev_score, expected.getScores().massdev_score)
            TEST_REAL_SIMILAR(actual.getScores().isotope_correlation, expected.getScores().isotope_correlation)
          }
        }
        TEST_TRUE(empty.getFeatures().empty())
        TEST_FALSE(access->requested_rts.empty())
        if (!uis)
        {
          std::vector<double> expected_rts;
          for (const auto& id : {"tr_gr2", "tr_gr1"})
          {
            for (const auto& feature : scheduled_seed.at(id).getFeatures())
            {
              expected_rts.push_back(feature.getRT());
            }
          }
          std::stable_sort(expected_rts.begin(), expected_rts.end());
          expected_rts.insert(expected_rts.end(), expected_remainder_rts.begin(), expected_remainder_rts.end());
          // REMAINDER requests must follow every candidate and retain input-group
          // order, even when disjoint groups finish in the opposite order.
          TEST_TRUE(access->requested_rts == expected_rts)
        }
        // Unlimited scoring has no deferred remainder queries. Every actual DIA
        // request, including UIS requests at the same apex, must be RT ordered.
        if (report_limit == -1)
        {
          TEST_TRUE(std::is_sorted(access->requested_rts.begin(), access->requested_rts.end()))
          TEST_REAL_SIMILAR(access->requested_rts.front(), 3100.0)
          TEST_REAL_SIMILAR(access->requested_rts.back(), 3130.0)
        }
      }
    }

    if (!uis)
    {
      // A nonfinite apex disables sorting for the whole batch. Record and reject
      // that query before the spectrum reader sees it, then finish the other group.
      parameters.setValue("stop_report_after_feature", -1);
      scorer.setParameters(parameters);
      {
        auto invalid_nonfinite = seed.at("tr_gr1");
        invalid_nonfinite.setTransitionGroupID("");
        invalid_nonfinite.getFeaturesMuteable()[0].setRT(std::numeric_limits<double>::quiet_NaN());
        auto valid = seed;
        FeatureMap failed_output, valid_output;
        access->requested_rts.clear();
        const auto errors = MRMFeatureFinderScoringTestAccess::scoreBatch(scorer,
          {{&invalid_nonfinite, &failed_output}, {&valid.at("tr_gr1"), &valid_output},
           {&valid.at("tr_gr2"), &valid_output}}, trafo, maps);
        TEST_TRUE(static_cast<bool>(errors[0]))
        TEST_FALSE(static_cast<bool>(errors[1]))
        TEST_FALSE(static_cast<bool>(errors[2]))
        TEST_TRUE(failed_output.empty())
        TEST_EQUAL(valid_output.size(), 6)
        TEST_TRUE(std::is_sorted(access->requested_rts.begin(), access->requested_rts.end()))
      }
      for (double nonfinite : {std::numeric_limits<double>::quiet_NaN(),
                               std::numeric_limits<double>::infinity(),
                               -std::numeric_limits<double>::infinity()})
      {
        auto batch = seed;
        batch.at("tr_gr1").getFeaturesMuteable()[2].setRT(nonfinite);
        FeatureMap failed_output, valid_output;
        access->requested_rts.clear();
        const auto errors = MRMFeatureFinderScoringTestAccess::scoreBatch(scorer,
          {{&batch.at("tr_gr1"), &failed_output}, {&batch.at("tr_gr2"), &valid_output}}, trafo, maps);
        TEST_EQUAL(errors.size(), 2)
        TEST_TRUE(static_cast<bool>(errors[0]))
        TEST_FALSE(static_cast<bool>(errors[1]))
        TEST_TRUE(failed_output.empty())
        TEST_EQUAL(valid_output.size(), 3)
        TEST_EQUAL(access->requested_rts.size(), 6)
        TEST_REAL_SIMILAR(access->requested_rts.at(0), 3130.0)
        TEST_REAL_SIMILAR(access->requested_rts.at(1), 3105.0)
        TEST_FALSE(std::isfinite(access->requested_rts.at(2)))
        TEST_REAL_SIMILAR(access->requested_rts.at(3), 3120.0)
        TEST_REAL_SIMILAR(access->requested_rts.at(4), 3100.0)
        TEST_REAL_SIMILAR(access->requested_rts.at(5), 3130.0)
      }

      // Reuse the picked traces as precursor traces for an MS1-only assay.
      // This checks the separate precursor S/N, RT, DIA and output formulas.
      TransitionGroupMapType ms1_seed;
      for (const auto& [id, source] : seed)
      {
        auto& group = ms1_seed[id];
        group.setTransitionGroupID(id);
        auto chromatogram = source.getChromatograms().front();
        const std::string precursor_id = id + "_Precursor_i0";
        chromatogram.setNativeID(precursor_id);
        group.addPrecursorChromatogram(chromatogram, precursor_id);
        for (const auto& source_feature : source.getFeatures())
        {
          auto feature = source_feature;
          auto precursor = source_feature.getFeatures().front();
          precursor.setMetaValue("native_id", precursor_id);
          feature.addPrecursorFeature(precursor, precursor_id);
          feature.setMZ(source.getTransitions().front().getPrecursorMZ());
          group.addFeature(feature);
        }
      }
      scorer.setMS1Map(access);
      for (int report_limit : {-1, 0, 1})
      {
        parameters.setValue("stop_report_after_feature", report_limit);
        scorer.setParameters(parameters);
        auto original = ms1_seed;
        auto batch = ms1_seed;
        FeatureMap expected, actual;
        std::vector<std::pair<MRMFeatureFinderScoring::MRMTransitionGroupType*, FeatureMap*>> work;
        for (const auto& id : {"tr_gr2", "tr_gr1"})
        {
          scorer.scorePeakgroups(original.at(id), trafo, {}, expected, true);
          work.emplace_back(&batch.at(id), &actual);
        }
        const auto errors = MRMFeatureFinderScoringTestAccess::scoreBatch(scorer, work, trafo, {}, true);
        for (const auto& error : errors)
        {
          TEST_FALSE(static_cast<bool>(error))
          if (error) { std::rethrow_exception(error); }
        }
        TEST_EQUAL(actual.size(), expected.size())
        for (Size i = 0; i < expected.size(); ++i)
        {
          TEST_TRUE(actual[i] == expected[i])
        }
        for (const auto& [id, group] : original)
        {
          const auto& result = batch.at(id);
          TEST_EQUAL(result.getPrecursorChromatograms().size(), group.getPrecursorChromatograms().size())
          for (Size i = 0; i < group.getFeatures().size(); ++i)
          {
            TEST_TRUE(result.getFeatures()[i] == group.getFeatures()[i])
            TEST_REAL_SIMILAR(result.getFeatures()[i].getScores().ms1_ppm_score, group.getFeatures()[i].getScores().ms1_ppm_score)
          }
        }
      }
    }

    // A failed preparation must not prevent another group's scoring/output.
    auto invalid = seed.at("tr_gr1");
    invalid.setTransitionGroupID("");
    auto valid = seed.at("tr_gr2");
    parameters.setValue("stop_report_after_feature", 1);
    scorer.setParameters(parameters);
    FeatureMap failed_output, valid_output;
    const auto errors = MRMFeatureFinderScoringTestAccess::scoreBatch(scorer,
      {{&invalid, &failed_output}, {&valid, &valid_output}}, trafo, maps);
    TEST_EQUAL(errors.size(), 2)
    TEST_TRUE(static_cast<bool>(errors[0]))
    TEST_FALSE(static_cast<bool>(errors[1]))
    TEST_TRUE(failed_output.empty())
    TEST_EQUAL(valid_output.size(), 1)
  }
  #ifdef _OPENMP
  omp_set_num_threads(previous_threads);
  #endif
}
END_SECTION

START_SECTION(void setStrictFlag(bool f))
  NOT_TESTABLE
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
