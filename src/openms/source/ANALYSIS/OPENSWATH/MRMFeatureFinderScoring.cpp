// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>

// data access
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/MRMFeatureAccessOpenMS.h>

// peak picking & noise estimation
#include <OpenMS/ANALYSIS/OPENSWATH/MRMScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h>

// Helpers
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>

#include <boost/range/adaptor/map.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>

#define run_identifier "unique_run_identifier"

bool SortDoubleDoublePairFirst(const std::pair<double, double>& left, const std::pair<double, double>& right)
{
  return left.first < right.first;
}


void processFeatureForOutput(OpenMS::Feature& curr_feature, bool write_convex_hull_, double
                             quantification_cutoff_, double& total_intensity, double& total_peak_apices, const std::string& ms_level)
{
  // Save some space when writing out the featureXML
  if (!write_convex_hull_)
  {
    curr_feature.getConvexHulls().clear();
  }

  // Ensure a unique id is present
  curr_feature.ensureUniqueId();

  // Sum up intensities of the features
  if (curr_feature.getMZ() > quantification_cutoff_)
  {
    total_intensity += curr_feature.getIntensity();
    total_peak_apices += (double)curr_feature.getMetaValue("peak_apex_int");
  }

  curr_feature.setMetaValue("FeatureLevel", ms_level);
}

namespace OpenMS
{

  MRMFeatureFinderScoring::MRMFeatureFinderScoring() :
    DefaultParamHandler("MRMFeatureFinderScoring"),
    ProgressLogger()
  {
    defaults_.setValue("stop_report_after_feature", -1, "Stop reporting after feature (ordered by quality; -1 means do not stop).");
    defaults_.setValue("rt_extraction_window", -1.0, "Only extract RT around this value (-1 means extract over the whole range, a value of 500 means to extract around +/- 500 s of the expected elution). For this to work, the TraML input file needs to contain normalized RT values.");
    defaults_.setValue("rt_normalization_factor", 1.0, "The normalized RT is expected to be between 0 and 1. If your normalized RT has a different range, pass this here (e.g. it goes from 0 to 100, set this value to 100)");
    defaults_.setValue("quantification_cutoff", 0.0, "Cutoff in m/z below which peaks should not be used for quantification any more", {"advanced"});
    defaults_.setMinFloat("quantification_cutoff", 0.0);
    defaults_.setValue("write_convex_hull", "false", "Whether to write out all points of all features into the featureXML", {"advanced"});
    defaults_.setValidStrings("write_convex_hull", {"true","false"});
    defaults_.setValue("spectrum_addition_method", "simple", "For spectrum addition, either use simple concatenation or use peak resampling", {"advanced"});
    defaults_.setValidStrings("spectrum_addition_method", {"simple", "resample"});
    defaults_.setValue("add_up_spectra", 1, "Add up spectra around the peak apex (needs to be a non-even integer)", {"advanced"});
    defaults_.setMinInt("add_up_spectra", 1);
    defaults_.setValue("spacing_for_spectra_resampling", 0.005, "If spectra are to be added, use this spacing to add them up", {"advanced"});
    defaults_.setMinFloat("spacing_for_spectra_resampling", 0.0);
    defaults_.setValue("uis_threshold_sn", -1, "S/N threshold to consider identification transition (set to -1 to consider all)");
    defaults_.setValue("uis_threshold_peak_area", 0, "Peak area threshold to consider identification transition (set to -1 to consider all)");
    defaults_.setValue("scoring_model", "default", "Scoring model to use", {"advanced"});
    defaults_.setValidStrings("scoring_model", {"default","single_transition"});
    defaults_.setValue("im_extra_drift", 0.0, "Extra drift time to extract for IM scoring (as a fraction, e.g. 0.25 means 25% extra on each side)", {"advanced"});
    defaults_.setMinFloat("im_extra_drift", 0.0);
    defaults_.setValue("strict", "true", "Whether to error (true) or skip (false) if a transition in a transition group does not have a corresponding chromatogram.", {"advanced"});
    defaults_.setValidStrings("strict", {"true","false"});
    defaults_.setValue("use_ms1_ion_mobility", "true", "Performs ion mobility extraction in MS1. Set to false if MS1 spectra do not contain ion mobility", {"advanced"});

    defaults_.insert("TransitionGroupPicker:", MRMTransitionGroupPicker().getDefaults());

    defaults_.insert("DIAScoring:", DIAScoring().getDefaults());

    defaults_.insert("EMGScoring:", EmgScoring().getDefaults());

    // One can turn on / off each score individually
    Param scores_to_use;
    scores_to_use.setValue("use_shape_score", "true", "Use the shape score (this score measures the similarity in shape of the transitions using a cross-correlation)", {"advanced"});
    scores_to_use.setValidStrings("use_shape_score", {"true","false"});
    scores_to_use.setValue("use_coelution_score", "true", "Use the coelution score (this score measures the similarity in coelution of the transitions using a cross-correlation)", {"advanced"});
    scores_to_use.setValidStrings("use_coelution_score", {"true","false"});
    scores_to_use.setValue("use_rt_score", "true", "Use the retention time score (this score measure the difference in retention time)", {"advanced"});
    scores_to_use.setValidStrings("use_rt_score", {"true","false"});
    scores_to_use.setValue("use_library_score", "true", "Use the library score", {"advanced"});
    scores_to_use.setValidStrings("use_library_score", {"true","false"});
    scores_to_use.setValue("use_elution_model_score", "true", "Use the elution model (EMG) score (this score fits a gaussian model to the peak and checks the fit)", {"advanced"});
    scores_to_use.setValidStrings("use_elution_model_score", {"true","false"});
    scores_to_use.setValue("use_intensity_score", "true", "Use the intensity score", {"advanced"});
    scores_to_use.setValidStrings("use_intensity_score", {"true","false"});
    scores_to_use.setValue("use_nr_peaks_score", "true", "Use the number of peaks score", {"advanced"});
    scores_to_use.setValidStrings("use_nr_peaks_score", {"true","false"});
    scores_to_use.setValue("use_total_xic_score", "true", "Use the total XIC score", {"advanced"});
    scores_to_use.setValidStrings("use_total_xic_score", {"true","false"});
    scores_to_use.setValue("use_total_mi_score", "false", "Use the total MI score", {"advanced"});
    scores_to_use.setValidStrings("use_total_mi_score", {"true","false"});
    scores_to_use.setValue("use_sn_score", "true", "Use the SN (signal to noise) score", {"advanced"});
    scores_to_use.setValidStrings("use_sn_score", {"true","false"});
    scores_to_use.setValue("use_mi_score", "false", "Use the MI (mutual information) score", {"advanced"});
    scores_to_use.setValidStrings("use_mi_score", {"true","false"});
    scores_to_use.setValue("use_dia_scores", "true", "Use the DIA (SWATH) scores. If turned off, will not use fragment ion spectra for scoring.", {"advanced"});
    scores_to_use.setValidStrings("use_dia_scores", {"true","false"});
    scores_to_use.setValue("use_ms1_correlation", "false", "Use the correlation scores with the MS1 elution profiles", {"advanced"});
    scores_to_use.setValidStrings("use_ms1_correlation", {"true","false"});
    scores_to_use.setValue("use_ion_mobility_scores", "false", "Use the scores for Ion Mobility scans", {"advanced"});
    scores_to_use.setValidStrings("use_ion_mobility_scores", {"true","false"});
    scores_to_use.setValue("use_ms1_fullscan", "false", "Use the full MS1 scan at the peak apex for scoring (ppm accuracy of precursor and isotopic pattern)", {"advanced"});
    scores_to_use.setValidStrings("use_ms1_fullscan", {"true","false"});
    scores_to_use.setValue("use_ms1_mi", "false", "Use the MS1 MI score", {"advanced"});
    scores_to_use.setValidStrings("use_ms1_mi", {"true","false"});
    scores_to_use.setValue("use_uis_scores", "false", "Use UIS scores for peptidoform identification", {"advanced"});
    scores_to_use.setValidStrings("use_uis_scores", {"true","false"});
    scores_to_use.setValue("use_peak_shape_metrics", "false", "Use peak shape metrics for scoring", {"advanced"});
    scores_to_use.setValue("use_ionseries_scores", "true", "Use MS2-level b/y ion-series scores for peptidoform identification", {"advanced"});
    scores_to_use.setValidStrings("use_ionseries_scores", {"true","false"});
    scores_to_use.setValue("use_ms2_isotope_scores", "true", "Use MS2-level isotope scores (pearson & manhattan) across product transitions (based on ID if annotated or averagine)", {"advanced"});
    scores_to_use.setValidStrings("use_ms2_isotope_scores", {"true","false"});
    defaults_.insert("Scores:", scores_to_use);

    // write defaults into Param object param_
    defaultsToParam_();
  }

  MRMFeatureFinderScoring::~MRMFeatureFinderScoring() = default;

  void MRMFeatureFinderScoring::pickExperiment(const PeakMap& chromatograms,
                                               FeatureMap& output,
                                               const TargetedExperiment& transition_exp_,
                                               const TransformationDescription& trafo,
                                               const PeakMap& swath_map)
  {
    OpenSwath::LightTargetedExperiment transition_exp;
    OpenSwathDataAccessHelper::convertTargetedExp(transition_exp_, transition_exp);
    TransitionGroupMapType transition_group_map;

    boost::shared_ptr<PeakMap > sh_chromatograms = boost::make_shared<PeakMap >(chromatograms);
    boost::shared_ptr<PeakMap > sh_swath_map = boost::make_shared<PeakMap >(swath_map);

    OpenSwath::SpectrumAccessPtr chromatogram_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(sh_chromatograms);
    OpenSwath::SpectrumAccessPtr swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(sh_swath_map);

    OpenSwath::SwathMap m;
    m.sptr = swath_ptr;
    std::vector<OpenSwath::SwathMap> swath_ptrs;
    swath_ptrs.push_back(m);

    pickExperiment(chromatogram_ptr, output, transition_exp, trafo, swath_ptrs, transition_group_map);
  }

  void MRMFeatureFinderScoring::pickExperiment(const OpenSwath::SpectrumAccessPtr& input,
                                               FeatureMap& output,
                                               const OpenSwath::LightTargetedExperiment& transition_exp,
                                               const TransformationDescription& trafo,
                                               const std::vector<OpenSwath::SwathMap>& swath_maps,
                                               TransitionGroupMapType& transition_group_map)
  {
    //
    // Step 1
    //
    // Store the peptide retention times in an intermediate map
    prepareProteinPeptideMaps_(transition_exp);

    // Store the proteins from the input in the output feature map
    std::vector<ProteinHit> protein_hits;
    for (const ProteinType& prot : transition_exp.getProteins())
    {
      ProteinHit prot_hit = ProteinHit();
      prot_hit.setSequence(prot.sequence);
      prot_hit.setAccession(prot.id);
      protein_hits.push_back(prot_hit);
    }

    ProteinIdentification prot_id = ProteinIdentification();
    prot_id.setHits(protein_hits);
    prot_id.setIdentifier(run_identifier);
    output.getProteinIdentifications().push_back(prot_id);

    //
    // Step 2
    //
    // Create all MRM transition groups from the individual transitions.
    mapExperimentToTransitionList(input, transition_exp, transition_group_map, trafo, rt_extraction_window_);
    int counter = 0;
    for (const auto& trgroup : transition_group_map)
    {
      if (!trgroup.second.getChromatograms().empty()) {counter++; }
    }
    OPENMS_LOG_INFO << "Will analyse " << counter << " peptides with a total of " << transition_exp.getTransitions().size() << " transitions " << std::endl;

    //
    // Step 3
    //
    // Go through all transition groups: first create consensus features, then score them
    MRMTransitionGroupPicker trgroup_picker;
    Param trgroup_picker_param = param_.copy("TransitionGroupPicker:", true);
    // If use_total_mi_score is defined, we need to instruct MRMTransitionGroupPicker to compute the score
    if (su_.use_total_mi_score_)
    {
      trgroup_picker_param.setValue("compute_total_mi", "true");
    }
    trgroup_picker.setParameters(trgroup_picker_param);

    Size progress = 0;
    startProgress(0, transition_group_map.size(), "picking peaks");
    for (TransitionGroupMapType::iterator trgroup_it = transition_group_map.begin(); trgroup_it != transition_group_map.end(); ++trgroup_it)
    {

      setProgress(++progress);
      MRMTransitionGroupType& transition_group = trgroup_it->second;
      if (transition_group.getChromatograms().empty() || transition_group.getTransitions().empty())
      {
        continue;
      }

      trgroup_picker.pickTransitionGroup(transition_group);
      scorePeakgroups(trgroup_it->second, trafo, swath_maps, output);
    }
    endProgress();

    //output.sortByPosition(); // if the exact same order is needed
    return;
  }

  void MRMFeatureFinderScoring::prepareProteinPeptideMaps_(const OpenSwath::LightTargetedExperiment& transition_exp)
  {
    for (Size i = 0; i < transition_exp.getCompounds().size(); i++)
    {
      PeptideRefMap_[transition_exp.getCompounds()[i].id] = &transition_exp.getCompounds()[i];
    }
  }

  void MRMFeatureFinderScoring::splitTransitionGroupsDetection_(const MRMTransitionGroupType& transition_group,
                                                                MRMTransitionGroupType& transition_group_detection) const
  {
    std::vector<TransitionType> tr = transition_group.getTransitions();
    std::vector<std::string> detecting_transitions;
    for (std::vector<TransitionType>::const_iterator tr_it = tr.begin(); tr_it != tr.end(); ++tr_it)
    {
      if (tr_it->isDetectingTransition())
      {
        detecting_transitions.push_back(tr_it->getNativeID());
      }
    }

    if (detecting_transitions.size() == transition_group.getTransitions().size())
    {
      transition_group_detection = transition_group;
    }
    else
    {
      transition_group_detection = transition_group.subset(detecting_transitions);
    }
  }

  void MRMFeatureFinderScoring::splitTransitionGroupsIdentification_(const MRMTransitionGroupType& transition_group,
                                                                     MRMTransitionGroupType& transition_group_identification,
                                                                     MRMTransitionGroupType& transition_group_identification_decoy) const
  {
    std::vector<TransitionType> tr = transition_group.getTransitions();
    std::vector<std::string> identifying_transitions, identifying_transitions_decoy;
    for (std::vector<TransitionType>::iterator tr_it = tr.begin(); tr_it != tr.end(); ++tr_it)
    {
      if (tr_it->isIdentifyingTransition())
      {
        if (tr_it->decoy)
        {
          identifying_transitions_decoy.push_back(tr_it->getNativeID());
        }
        else
        {
          identifying_transitions.push_back(tr_it->getNativeID());
        }
      }
    }

    transition_group_identification = transition_group.subsetDependent(identifying_transitions);
    transition_group_identification_decoy = transition_group.subsetDependent(identifying_transitions_decoy);
  }

  OpenSwath_Ind_Scores MRMFeatureFinderScoring::scoreIdentification_(MRMTransitionGroupType& trgr_ident,
                                                                     OpenSwathScoring& scorer,
                                                                     const size_t feature_idx,
                                                                     const std::vector<std::string>& native_ids_detection,
                                                                     const double det_intensity_ratio_score,
                                                                     const double det_mi_ratio_score,
                                                                     const std::vector<OpenSwath::SwathMap>& swath_maps) const
  {
    MRMFeature idmrmfeature = trgr_ident.getFeaturesMuteable()[feature_idx];
    OpenSwath::IMRMFeature* idimrmfeature;
    idimrmfeature = new MRMFeatureOpenMS(idmrmfeature);

    // get drift time upper/lower offset (this assumes that all chromatograms
    // are derived from the same precursor with the same drift time)
    RangeMobility im_range;

    if ( (!trgr_ident.getChromatograms().empty()) || (!trgr_ident.getPrecursorChromatograms().empty()) )
    {
      auto & prec = trgr_ident.getChromatograms()[0].getPrecursor();
      im_range.setMin(prec.getDriftTime()); // sets the minimum and maximum
      im_range.minSpanIfSingular(prec.getDriftTimeWindowLowerOffset());
    }

    std::vector<std::string> native_ids_identification;
    std::vector<OpenSwath::ISignalToNoisePtr> signal_noise_estimators_identification;

    for (Size i = 0; i < trgr_ident.size(); i++)
    {
      OpenSwath::ISignalToNoisePtr snptr(new OpenMS::SignalToNoiseOpenMS< MSChromatogram >(
            trgr_ident.getChromatogram(trgr_ident.getTransitions()[i].getNativeID()),
            sn_win_len_, sn_bin_count_, write_log_messages_));
      if (  (snptr->getValueAtRT(idmrmfeature.getRT()) > uis_threshold_sn_)
            && (idmrmfeature.getFeature(trgr_ident.getTransitions()[i].getNativeID()).getIntensity() > uis_threshold_peak_area_))
      {
        signal_noise_estimators_identification.push_back(snptr);
        native_ids_identification.push_back(trgr_ident.getTransitions()[i].getNativeID());
      }
    }

    OpenSwath_Ind_Scores idscores;
    if (!native_ids_identification.empty())
    {
      scorer.calculateChromatographicIdScores(idimrmfeature,
                                              native_ids_identification,
                                              native_ids_detection,
                                              signal_noise_estimators_identification,
                                              idscores);

      std::vector<double> ind_mi_score;
      if (su_.use_mi_score_)
      {
        ind_mi_score = idscores.ind_mi_score;
      }

      for (size_t i = 0; i < native_ids_identification.size(); i++)
      {
        idscores.ind_transition_names.emplace_back(native_ids_identification[i]);
        if (idmrmfeature.getFeature(native_ids_identification[i]).getIntensity() > 0)
        {
          double intensity_score = double(idmrmfeature.getFeature(native_ids_identification[i]).getIntensity()) / double(idmrmfeature.getFeature(native_ids_identification[i]).getMetaValue("total_xic"));

          double intensity_ratio = 0;
          if (det_intensity_ratio_score > 0) { intensity_ratio = intensity_score / det_intensity_ratio_score; }
          if (intensity_ratio > 1) { intensity_ratio = 1 / intensity_ratio; }

          double total_mi = 0;
          if (su_.use_total_mi_score_)
          {
            total_mi = double(idmrmfeature.getFeature(native_ids_identification[i]).getMetaValue("total_mi"));
          }

          double mi_ratio = 0;
          if (su_.use_mi_score_ && su_.use_total_mi_score_)
          {
            if (det_mi_ratio_score > 0) { mi_ratio = (ind_mi_score[i] / total_mi) / det_mi_ratio_score; }
            if (mi_ratio > 1) { mi_ratio = 1 / mi_ratio; }
          }

          idscores.ind_area_intensity.push_back(idmrmfeature.getFeature(native_ids_identification[i]).getIntensity());
          idscores.ind_total_area_intensity.push_back(idmrmfeature.getFeature(native_ids_identification[i]).getMetaValue("total_xic"));
          idscores.ind_intensity_score.push_back(intensity_score);
          idscores.ind_apex_intensity.push_back(idmrmfeature.getFeature(native_ids_identification[i]).getMetaValue("peak_apex_int"));
          idscores.ind_apex_position.push_back(idmrmfeature.getFeature(native_ids_identification[i]).getMetaValue("peak_apex_position"));
          idscores.ind_fwhm.push_back(idmrmfeature.getFeature(native_ids_identification[i]).getMetaValue("width_at_50"));
          idscores.ind_total_mi .push_back(total_mi);
          idscores.ind_log_intensity.push_back(std::log(idmrmfeature.getFeature(native_ids_identification[i]).getIntensity()));
          idscores.ind_intensity_ratio.push_back(intensity_ratio);
          idscores.ind_mi_ratio.push_back(mi_ratio);

          if (su_.use_peak_shape_metrics)
          {
            idscores.ind_start_position_at_5.push_back(idmrmfeature.getFeature(native_ids_identification[i]).getMetaValue("start_position_at_5"));
            idscores.ind_end_position_at_5.push_back(idmrmfeature.getFeature(native_ids_identification[i]).getMetaValue("end_position_at_5"));
            idscores.ind_start_position_at_10.push_back(idmrmfeature.getFeature(native_ids_identification[i]).getMetaValue("start_position_at_10"));
            idscores.ind_end_position_at_10.push_back(idmrmfeature.getFeature(native_ids_identification[i]).getMetaValue("end_position_at_10"));
            idscores.ind_start_position_at_50.push_back(idmrmfeature.getFeature(native_ids_identification[i]).getMetaValue("start_position_at_50"));
            idscores.ind_end_position_at_50.push_back(idmrmfeature.getFeature(native_ids_identification[i]).getMetaValue("end_position_at_50"));
            idscores.ind_total_width.push_back(idmrmfeature.getFeature(native_ids_identification[i]).getMetaValue("total_width"));
            idscores.ind_tailing_factor.push_back(idmrmfeature.getFeature(native_ids_identification[i]).getMetaValue("tailing_factor"));
            idscores.ind_asymmetry_factor.push_back(idmrmfeature.getFeature(native_ids_identification[i]).getMetaValue("asymmetry_factor"));
            idscores.ind_slope_of_baseline.push_back(idmrmfeature.getFeature(native_ids_identification[i]).getMetaValue("slope_of_baseline"));
            idscores.ind_baseline_delta_2_height.push_back(idmrmfeature.getFeature(native_ids_identification[i]).getMetaValue("baseline_delta_2_height"));
            idscores.ind_points_across_baseline.push_back(idmrmfeature.getFeature(native_ids_identification[i]).getMetaValue("points_across_baseline"));
            idscores.ind_points_across_half_height.push_back(idmrmfeature.getFeature(native_ids_identification[i]).getMetaValue("points_across_half_height"));
          }
        }
        else
        {
          idscores.ind_area_intensity.push_back(0);
          idscores.ind_total_area_intensity.push_back(0);
          idscores.ind_intensity_score.push_back(0);
          idscores.ind_apex_intensity.push_back(0);
          idscores.ind_apex_position.push_back(0);
          idscores.ind_fwhm.push_back(0);
          idscores.ind_total_mi.push_back(0);
          idscores.ind_log_intensity.push_back(0);
          idscores.ind_intensity_ratio.push_back(0);
          idscores.ind_mi_ratio.push_back(0);

          if (su_.use_peak_shape_metrics)
          {
            idscores.ind_start_position_at_5.push_back(0);
            idscores.ind_end_position_at_5.push_back(0);
            idscores.ind_start_position_at_10.push_back(0);
            idscores.ind_end_position_at_10.push_back(0);
            idscores.ind_start_position_at_50.push_back(0);
            idscores.ind_end_position_at_50.push_back(0);
            idscores.ind_total_width.push_back(0);
            idscores.ind_tailing_factor.push_back(0);
            idscores.ind_asymmetry_factor.push_back(0);
            idscores.ind_slope_of_baseline.push_back(0);
            idscores.ind_baseline_delta_2_height.push_back(0);
            idscores.ind_points_across_baseline.push_back(0);
            idscores.ind_points_across_half_height.push_back(0);
          }
        }
      idscores.ind_num_transitions = native_ids_identification.size();

      }
    }

    // Compute DIA scores only on the identification transitions
    bool swath_present = (!swath_maps.empty() && swath_maps[0].sptr->getNrSpectra() > 0);
    if (swath_present && su_.use_dia_scores_ && !native_ids_identification.empty())
    {
      std::vector<double> ind_isotope_correlation, ind_isotope_overlap, ind_massdev_score;
      for (size_t i = 0; i < native_ids_identification.size(); i++)
      {
        OpenSwath_Scores tmp_scores;

        scorer.calculateDIAIdScores(idimrmfeature,
                                    trgr_ident.getTransition(native_ids_identification[i]),
                                    swath_maps, im_range, diascoring_, tmp_scores);

        ind_isotope_correlation.push_back(tmp_scores.isotope_correlation);
        ind_isotope_overlap.push_back(tmp_scores.isotope_overlap);
        ind_massdev_score.push_back(tmp_scores.massdev_score);
      }
      idscores.ind_isotope_correlation = ind_isotope_correlation;
      idscores.ind_isotope_overlap = ind_isotope_overlap;
      idscores.ind_massdev_score = ind_massdev_score;
    }

    delete idimrmfeature;
    return idscores;
  }

  void MRMFeatureFinderScoring::scorePeakgroups(MRMTransitionGroupType& transition_group,
                                                const TransformationDescription& trafo,
                                                const std::vector<OpenSwath::SwathMap>& swath_maps,
                                                FeatureMap& output,
                                                bool ms1only) const
  {
    if (PeptideRefMap_.empty())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "Error: Peptide reference map is empty, please call prepareProteinPeptideMaps_ first.");
    }
    if (transition_group.getTransitionGroupID().empty())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "Error: Transition group id is empty, please set it.");
    }

    MRMTransitionGroupType transition_group_detection, transition_group_identification, transition_group_identification_decoy;
    splitTransitionGroupsDetection_(transition_group, transition_group_detection);
    if (su_.use_uis_scores)
    {
      splitTransitionGroupsIdentification_(transition_group, transition_group_identification, transition_group_identification_decoy);
    }

    std::vector<OpenSwath::ISignalToNoisePtr> signal_noise_estimators;
    std::vector<MRMFeature> feature_list;

    // get drift time upper/lower offset (this assumes that all chromatograms
    // are derived from the same precursor with the same drift time)
    RangeMobility im_range;
    double drift_target(0);

    auto setDriftTarget = [](auto& prec){
      double lower_bound = prec.getDriftTime() - prec.getDriftTimeWindowLowerOffset();
      double upper_bound = prec.getDriftTime() + prec.getDriftTimeWindowUpperOffset();
      return RangeMobility(lower_bound, upper_bound);
    };

    if ( !transition_group_detection.getChromatograms().empty() )
    {
      auto & prec = transition_group_detection.getChromatograms()[0].getPrecursor();
      drift_target = prec.getDriftTime();

      if (drift_target > 0) im_range = setDriftTarget(prec);
    }

    else if ( !transition_group_detection.getPrecursorChromatograms().empty() )
    {
      auto & prec = transition_group_detection.getPrecursorChromatograms()[0].getPrecursor();
      drift_target = prec.getDriftTime();

      if (drift_target > 0) im_range = setDriftTarget(prec);
    }

    // currently we cannot do much about the log messages and they mostly occur in decoy transition signals
    for (Size k = 0; k < transition_group_detection.getChromatograms().size(); k++)
    {
      OpenSwath::ISignalToNoisePtr snptr(new OpenMS::SignalToNoiseOpenMS< MSChromatogram >(
            transition_group_detection.getChromatograms()[k], sn_win_len_, sn_bin_count_, write_log_messages_));
      signal_noise_estimators.push_back(snptr);
    }

    // skip MS1 noise estimator if we perform fragment ion analysis
    std::vector<OpenSwath::ISignalToNoisePtr> ms1_signal_noise_estimators;
    if (ms1only)
    {
      for (Size k = 0; k < transition_group_detection.getPrecursorChromatograms().size(); k++)
      {
        OpenSwath::ISignalToNoisePtr snptr(new OpenMS::SignalToNoiseOpenMS< MSChromatogram >(
              transition_group_detection.getPrecursorChromatograms()[k], sn_win_len_, sn_bin_count_, write_log_messages_));
        ms1_signal_noise_estimators.push_back(snptr);
      }
    }

    // get the expected rt value for this compound
    const PeptideType* pep = PeptideRefMap_.at(transition_group_detection.getTransitionGroupID());
    double expected_rt = pep->rt;
    TransformationDescription newtr = trafo;
    newtr.invert();
    expected_rt = newtr.apply(expected_rt);

    OpenSwathScoring scorer;
    scorer.initialize(rt_normalization_factor_, add_up_spectra_,
                      spacing_for_spectra_resampling_,
                      im_extra_drift_,
                      su_,
                      spectrum_addition_method_,
                      use_ms1_ion_mobility_);

    ProteaseDigestion pd;
    pd.setEnzyme("Trypsin");

    auto& mrmfeatures = transition_group_detection.getFeaturesMuteable();

    // Go through all peak groups (found MRM features) and score them
    #ifdef _OPENMP
    int in_parallel = omp_in_parallel();
    #endif
    #pragma omp parallel for if (in_parallel == 0)
    for (SignedSize feature_idx = 0; feature_idx < (SignedSize) mrmfeatures.size(); ++feature_idx)
    {
      auto& mrmfeature = mrmfeatures[feature_idx];
      OpenSwath::IMRMFeature* imrmfeature;
      imrmfeature = new MRMFeatureOpenMS(mrmfeature);

      OPENMS_LOG_DEBUG << "Scoring feature " << (mrmfeature) << " == " << mrmfeature.getMetaValue("PeptideRef") <<
        " [ expected RT " << PeptideRefMap_.at(mrmfeature.getMetaValue("PeptideRef"))->rt << " / " << expected_rt << " ]" <<
        " with " << transition_group_detection.size()  << " transitions and " <<
        transition_group_detection.getChromatograms().size() << " chromatograms" << std::endl;

      int group_size = boost::numeric_cast<int>(transition_group_detection.size());
      if (group_size == 0 && !ms1only)
      {
        delete imrmfeature; // free resources before continuing
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "Error: Transition group " + transition_group_detection.getTransitionGroupID() +
                                         " has no chromatograms.");
      }
      bool swath_present = (!swath_maps.empty() && swath_maps[0].sptr->getNrSpectra() > 0);
      double xx_lda_prescore;
      double precursor_mz(-1);

      if (ms1only)
      {
        ///////////////////////////////////
        // Call the scoring for MS1 only
        ///////////////////////////////////

        OpenSwath_Scores& scores = mrmfeature.getScores();
        precursor_mz = mrmfeature.getMZ();

        // S/N scores
        OpenSwath::MRMScoring mrmscore_;
        scores.sn_ratio = mrmscore_.calcSNScore(imrmfeature, ms1_signal_noise_estimators);
        // everything below S/N 1 can be set to zero (and the log safely applied)
        if (scores.sn_ratio < 1)
        {
          scores.log_sn_score = 0;
        }
        else
        {
          scores.log_sn_score = std::log(scores.sn_ratio);
        }
        if (su_.use_sn_score_)
        {
          mrmfeature.addScore("sn_ratio", scores.sn_ratio);
          mrmfeature.addScore("var_log_sn_score", scores.log_sn_score);
          // compute subfeature log-SN values
          for (Size k = 0; k < transition_group_detection.getPrecursorChromatograms().size(); k++)
          {
            Feature & f = mrmfeature.getPrecursorFeature(transition_group_detection.getPrecursorChromatograms()[k].getNativeID());
            double sn_value = ms1_signal_noise_estimators[k]->getValueAtRT(mrmfeature.getRT());
            if (sn_value < 1) {sn_value = 1.0;}
            f.setMetaValue("logSN", std::log(sn_value));
          }
        }

        // RT scores
        double normalized_experimental_rt = trafo.apply(imrmfeature->getRT());
        {
          // rt score is delta iRT
          double rt_score = mrmscore_.calcRTScore(*pep, normalized_experimental_rt);

          scores.normalized_experimental_rt = normalized_experimental_rt;
          scores.raw_rt_score = rt_score;
          scores.norm_rt_score = rt_score / rt_normalization_factor_;
        }
        if (su_.use_rt_score_)
        {
          mrmfeature.addScore("delta_rt", mrmfeature.getRT() - expected_rt);
          mrmfeature.addScore("assay_rt", expected_rt);
          mrmfeature.addScore("norm_RT", scores.normalized_experimental_rt);
          mrmfeature.addScore("rt_score", scores.raw_rt_score);
          mrmfeature.addScore("var_norm_rt_score", scores.norm_rt_score);
        }

        // full spectra scores
        if (ms1_map_ && ms1_map_->getNrSpectra() > 0 && mrmfeature.getMZ() > 0)
        {
          scorer.calculatePrecursorDIAScores(ms1_map_, diascoring_, precursor_mz, imrmfeature->getRT(), *pep, im_range, scores);
        }
        if (su_.use_ms1_fullscan)
        {
          mrmfeature.addScore("var_ms1_ppm_diff", scores.ms1_ppm_score);
          mrmfeature.addScore("var_ms1_isotope_correlation", scores.ms1_isotope_correlation);
          mrmfeature.addScore("var_ms1_isotope_overlap", scores.ms1_isotope_overlap);
        }
        xx_lda_prescore = -scores.calculate_lda_prescore(scores);
        if (scoring_model_ == "single_transition")
        {
          xx_lda_prescore = -scores.calculate_lda_single_transition(scores);
        }
        mrmfeature.addScore("main_var_xx_lda_prelim_score", xx_lda_prescore);
        mrmfeature.addScore("xx_lda_prelim_score", xx_lda_prescore);
        mrmfeature.setOverallQuality(xx_lda_prescore);
      }
      else //!ms1only
      {

        ///////////////////////////////////
        // Call the scoring for fragment ions
        ///////////////////////////////////

        std::vector<double> normalized_library_intensity;
        transition_group_detection.getLibraryIntensity(normalized_library_intensity);
        OpenSwath::Scoring::normalize_sum(&normalized_library_intensity[0], boost::numeric_cast<int>(normalized_library_intensity.size()));

        std::vector<std::string> native_ids_detection;
        for (Size i = 0; i < transition_group_detection.size(); i++)
        {
          std::string native_id = transition_group_detection.getTransitions()[i].getNativeID();
          native_ids_detection.push_back(native_id);
        }

        std::vector<std::string> precursor_ids;
        for (Size i = 0; i < transition_group_detection.getPrecursorChromatograms().size(); i++)
        {
          std::string precursor_id = transition_group_detection.getPrecursorChromatograms()[i].getNativeID();
          precursor_ids.push_back(precursor_id);
        }

        ///////////////////////////////////
        // Library and chromatographic scores
        OpenSwath_Scores& scores = mrmfeature.getScores();
        scorer.calculateChromatographicScores(imrmfeature, native_ids_detection, precursor_ids, normalized_library_intensity,
                                              signal_noise_estimators, scores);

        double normalized_experimental_rt = trafo.apply(imrmfeature->getRT());
        scorer.calculateLibraryScores(imrmfeature, transition_group_detection.getTransitions(), *pep, normalized_experimental_rt, scores);

        ///////////////////////////////////
        // DIA scores
        if (swath_present && su_.use_dia_scores_)
        {
          std::vector<double> masserror_ppm;
          scorer.calculateDIAScores(imrmfeature,
                                    transition_group_detection.getTransitions(),
                                    swath_maps, ms1_map_, diascoring_, *pep, scores, masserror_ppm,
                                    drift_target, im_range);
          mrmfeature.setMetaValue("masserror_ppm", masserror_ppm);
        }
        
        double det_intensity_ratio_score = 0;
        if ((double)mrmfeature.getMetaValue("total_xic") > 0)
        {
          det_intensity_ratio_score = mrmfeature.getIntensity() / (double)mrmfeature.getMetaValue("total_xic");
        }

        ///////////////////////////////////
        // Mutual Information scores
        double det_mi_ratio_score = 0;
        if (su_.use_mi_score_ && su_.use_total_mi_score_)
        {
          if ((double)mrmfeature.getMetaValue("total_mi") > 0)
          {
            det_mi_ratio_score = scores.mi_score / (double)mrmfeature.getMetaValue("total_mi");
          }
        }

        ///////////////////////////////////
        // Unique Ion Signature (UIS) scores
        if (su_.use_uis_scores && !transition_group_identification.getTransitions().empty())
        {
          OpenSwath_Ind_Scores idscores = scoreIdentification_(transition_group_identification, scorer, feature_idx,
                                                               native_ids_detection, det_intensity_ratio_score,
                                                               det_mi_ratio_score, swath_maps);
          mrmfeature.IDScoresAsMetaValue(false, idscores);
        }
        if (su_.use_uis_scores && !transition_group_identification_decoy.getTransitions().empty())
        {
          OpenSwath_Ind_Scores idscores = scoreIdentification_(transition_group_identification_decoy, scorer, feature_idx,
                                                               native_ids_detection, det_intensity_ratio_score,
                                                               det_mi_ratio_score, swath_maps);
          mrmfeature.IDScoresAsMetaValue(true, idscores);
        }

        if (su_.use_coelution_score_)
        {
          mrmfeature.addScore("var_xcorr_coelution", scores.xcorr_coelution_score);
          mrmfeature.addScore("var_xcorr_coelution_weighted", scores.weighted_coelution_score);
        }
        if (su_.use_shape_score_)
        {
          mrmfeature.addScore("var_xcorr_shape", scores.xcorr_shape_score);
          mrmfeature.addScore("var_xcorr_shape_weighted", scores.weighted_xcorr_shape);
        }
        if (su_.use_library_score_)
        {
          mrmfeature.addScore("var_library_corr", scores.library_corr);
          mrmfeature.addScore("var_library_rmsd", scores.library_norm_manhattan);
          mrmfeature.addScore("var_library_sangle", scores.library_sangle);
          mrmfeature.addScore("var_library_rootmeansquare", scores.library_rootmeansquare);
          mrmfeature.addScore("var_library_manhattan", scores.library_manhattan);
          mrmfeature.addScore("var_library_dotprod", scores.library_dotprod);
        }
        if (su_.use_rt_score_)
        {
          mrmfeature.addScore("delta_rt", mrmfeature.getRT() - expected_rt);
          mrmfeature.addScore("assay_rt", expected_rt);
          mrmfeature.addScore("norm_RT", scores.normalized_experimental_rt);
          mrmfeature.addScore("rt_score", scores.raw_rt_score);
          mrmfeature.addScore("var_norm_rt_score", scores.norm_rt_score);
        }
        // TODO do we really want these intensity scores ?
        if (su_.use_intensity_score_)
        {
          if ((double)mrmfeature.getMetaValue("total_xic") > 0)
          {
            mrmfeature.addScore("var_intensity_score", mrmfeature.getIntensity() / (double)mrmfeature.getMetaValue("total_xic"));
          }
          else
          {
            mrmfeature.addScore("var_intensity_score", 0);
          }
        }
        if (su_.use_total_xic_score_) { mrmfeature.addScore("total_xic", (double)mrmfeature.getMetaValue("total_xic")); }
        if (su_.use_total_mi_score_) { mrmfeature.addScore("total_mi", (double)mrmfeature.getMetaValue("total_mi")); }

        if (su_.use_nr_peaks_score_) { mrmfeature.addScore("nr_peaks", scores.nr_peaks); }
        if (su_.use_sn_score_)
        {
          mrmfeature.addScore("sn_ratio", scores.sn_ratio);
          mrmfeature.addScore("var_log_sn_score", scores.log_sn_score);
          // compute subfeature log-SN values
          for (Size k = 0; k < transition_group_detection.getChromatograms().size(); k++)
          {
            Feature & f = mrmfeature.getFeature(transition_group_detection.getChromatograms()[k].getNativeID());
            double sn_value = signal_noise_estimators[k]->getValueAtRT(mrmfeature.getRT());
            if (sn_value < 1) {sn_value = 1.0;}
            f.setMetaValue("logSN", std::log(sn_value));
          }
        }

        if (su_.use_mi_score_)
        {
          mrmfeature.addScore("var_mi_score", scores.mi_score);
          mrmfeature.addScore("var_mi_weighted_score", scores.weighted_mi_score);
          if (su_.use_total_mi_score_)
          {
            if (((double)mrmfeature.getMetaValue("total_mi")) > 0)
            {
              mrmfeature.addScore("var_mi_ratio_score", scores.mi_score  / (double)mrmfeature.getMetaValue("total_mi"));
            }
            else
            {
              mrmfeature.addScore("var_mi_ratio_score", 0);
            }
          }
        }

        // TODO get it working with imrmfeature
        if (su_.use_elution_model_score_)
        {
          //TODO wouldn't a weighted elution model score be much better? lower intensity traces usually will not have
          // a nice profile
          scores.elution_model_fit_score = emgscoring_.calcElutionFitScore(mrmfeature, transition_group_detection);
          mrmfeature.addScore("var_elution_model_fit_score", scores.elution_model_fit_score);
        }

        xx_lda_prescore = -scores.calculate_lda_prescore(scores);
        if (scoring_model_ == "single_transition")
        {
          xx_lda_prescore = -scores.calculate_lda_single_transition(scores);
        }
        if (!swath_present)
        {
          mrmfeature.addScore("main_var_xx_lda_prelim_score", xx_lda_prescore);
        }
        mrmfeature.setOverallQuality(xx_lda_prescore);
        mrmfeature.addScore("xx_lda_prelim_score", xx_lda_prescore);

        // Add the DIA scores and ion mobility scores 
        if (swath_present && su_.use_dia_scores_)
        {
          if (su_.use_ms2_isotope_scores)
          {
            mrmfeature.addScore("var_isotope_correlation_score", scores.isotope_correlation);
            mrmfeature.addScore("var_isotope_overlap_score", scores.isotope_overlap);
          }

          mrmfeature.addScore("var_massdev_score", scores.massdev_score);
          mrmfeature.addScore("var_massdev_score_weighted", scores.weighted_massdev_score);

          if (su_.use_ionseries_scores)
          {
            mrmfeature.addScore("var_bseries_score", scores.bseries_score);
            mrmfeature.addScore("var_yseries_score", scores.yseries_score);
          }

          if (su_.use_ms2_isotope_scores)
          {
            mrmfeature.addScore("var_dotprod_score", scores.dotprod_score_dia);
            mrmfeature.addScore("var_manhatt_score", scores.manhatt_score_dia);
          }

          if (su_.use_ms1_correlation)
          {
            if (scores.ms1_xcorr_shape_score > -1)
            {
              mrmfeature.addScore("var_ms1_xcorr_shape", scores.ms1_xcorr_shape_score);
            }
            if (scores.ms1_xcorr_coelution_score > -1)
            {
              mrmfeature.addScore("var_ms1_xcorr_coelution", scores.ms1_xcorr_coelution_score);
            }
            mrmfeature.addScore("var_ms1_xcorr_shape_contrast", scores.ms1_xcorr_shape_contrast_score);
            mrmfeature.addScore("var_ms1_xcorr_shape_combined", scores.ms1_xcorr_shape_combined_score);
            mrmfeature.addScore("var_ms1_xcorr_coelution_contrast", scores.ms1_xcorr_coelution_contrast_score);
            mrmfeature.addScore("var_ms1_xcorr_coelution_combined", scores.ms1_xcorr_coelution_combined_score);
          }
          if (su_.use_ms1_mi)
          {
            if (scores.ms1_mi_score > -1)
            {
              mrmfeature.addScore("var_ms1_mi_score", scores.ms1_mi_score);
            }
            mrmfeature.addScore("var_ms1_mi_contrast_score", scores.ms1_mi_contrast_score);
            mrmfeature.addScore("var_ms1_mi_combined_score", scores.ms1_mi_combined_score);
          }
          if (su_.use_ms1_fullscan)
          {
            mrmfeature.addScore("var_ms1_ppm_diff", scores.ms1_ppm_score);
            mrmfeature.addScore("var_ms1_isotope_correlation", scores.ms1_isotope_correlation);
            mrmfeature.addScore("var_ms1_isotope_overlap", scores.ms1_isotope_overlap);
          }

          double xx_swath_prescore = -scores.calculate_swath_lda_prescore(scores);
          mrmfeature.addScore("main_var_xx_swath_prelim_score", xx_swath_prescore);
          mrmfeature.setOverallQuality(xx_swath_prescore);
        }

        if (swath_present && su_.use_im_scores)
        {
          mrmfeature.addScore("var_im_xcorr_shape", scores.im_xcorr_shape_score);
          mrmfeature.addScore("var_im_xcorr_coelution", scores.im_xcorr_coelution_score);
          mrmfeature.addScore("var_im_delta_score", scores.im_delta_score);
          mrmfeature.addScore("var_im_ms1_delta_score", scores.im_ms1_delta_score);
          mrmfeature.addScore("im_drift", scores.im_drift); // MS2 level
          mrmfeature.addScore("im_drift_weighted", scores.im_drift_weighted); // MS2 level
          mrmfeature.addScore("im_ms1_drift", scores.im_ms1_drift); // MS1 level
          mrmfeature.addScore("im_ms1_delta", scores.im_ms1_delta); // MS1 level
          mrmfeature.addScore("im_delta", scores.im_delta); // MS2 level
        }

        precursor_mz = transition_group_detection.getTransitions()[0].getPrecursorMZ();

      }

      ///////////////////////////////////////////////////////////////////////////
      // add the peptide hit information to the feature
      ///////////////////////////////////////////////////////////////////////////
      PeptideIdentification pep_id_;
      PeptideHit pep_hit_;

      if (pep->getChargeState() != 0)
      {
        pep_hit_.setCharge(pep->getChargeState());
      }
      pep_hit_.setScore(xx_lda_prescore);
      if (swath_present && mrmfeature.metaValueExists("xx_swath_prelim_score"))
      {
        pep_hit_.setScore(mrmfeature.getMetaValue("xx_swath_prelim_score"));
      }

      if (pep->isPeptide() && !pep->sequence.empty())
      {
        pep_hit_.setSequence(AASequence::fromString(pep->sequence));
        mrmfeature.setMetaValue("missedCleavages", pd.peptideCount(pep_hit_.getSequence()) - 1);
      }

      // set protein accession numbers
      for (Size k = 0; k < pep->protein_refs.size(); k++)
      {
        PeptideEvidence pe;
        pe.setProteinAccession(pep->protein_refs[k]);
        pep_hit_.addPeptideEvidence(pe);
      }
      pep_id_.insertHit(pep_hit_);
      pep_id_.setIdentifier(run_identifier);

      mrmfeature.getPeptideIdentifications().push_back(pep_id_);
      mrmfeature.ensureUniqueId();

      mrmfeature.setMetaValue("PrecursorMZ", precursor_mz);
      prepareFeatureOutput_(mrmfeature, ms1only, pep->getChargeState());
      mrmfeature.setMetaValue("xx_swath_prelim_score", 0.0);
      #pragma omp critical
      feature_list.push_back(mrmfeature);

      delete imrmfeature;
    }

    // Order by quality (high to low, via reverse iterator)
    std::sort(feature_list.rbegin(), feature_list.rend(), OpenMS::Feature::OverallQualityLess());

    for (Size i = 0; i < feature_list.size(); i++)
    {
      if (stop_report_after_feature_ >= 0 && i >= (Size)stop_report_after_feature_) {break;}
      output.push_back(feature_list[i]);
    }

    // store all data manipulation performed on the features of the transition group
    transition_group = transition_group_detection;
  }

  void MRMFeatureFinderScoring::prepareFeatureOutput_(OpenMS::MRMFeature& mrmfeature, bool ms1only, int charge) const
  {
    // Prepare the subordinates for the mrmfeature (process all current
    // features and then append all precursor subordinate features)
    std::vector<Feature> allFeatures = mrmfeature.getFeatures();
    double total_intensity = 0, total_peak_apices = 0;
    double ms1_total_intensity = 0, ms1_total_peak_apices = 0;

    for (std::vector<Feature>::iterator f_it = allFeatures.begin(); f_it != allFeatures.end(); ++f_it)
    {
      processFeatureForOutput(*f_it, write_convex_hull_, quantification_cutoff_, total_intensity, total_peak_apices, "MS2");
    }
    // Also append data for MS1 precursors
    std::vector<String> precursors_ids;
    mrmfeature.getPrecursorFeatureIDs(precursors_ids);
    for (std::vector<String>::iterator id_it = precursors_ids.begin(); id_it != precursors_ids.end(); ++id_it)
    {
      Feature curr_feature = mrmfeature.getPrecursorFeature(*id_it);
      if (charge != 0)
      {
        curr_feature.setCharge(charge);
      }
      processFeatureForOutput(curr_feature, write_convex_hull_, quantification_cutoff_, ms1_total_intensity, ms1_total_peak_apices, "MS1");
      if (ms1only)
      {
        total_intensity += curr_feature.getIntensity();
        total_peak_apices += (double)curr_feature.getMetaValue("peak_apex_int");
      }
      allFeatures.push_back(curr_feature);
    }
    mrmfeature.setSubordinates(allFeatures); // add all the subfeatures as subordinates

    // overwrite the reported intensities with those above the m/z cutoff
    mrmfeature.setIntensity(total_intensity);
    mrmfeature.setMetaValue("peak_apices_sum", total_peak_apices);
    mrmfeature.setMetaValue("ms1_area_intensity", ms1_total_intensity);
    mrmfeature.setMetaValue("ms1_apex_intensity", ms1_total_peak_apices);
  }

  void MRMFeatureFinderScoring::updateMembers_()
  {
    stop_report_after_feature_ = (int)param_.getValue("stop_report_after_feature");
    rt_extraction_window_ = (double)param_.getValue("rt_extraction_window");
    rt_normalization_factor_ = (double)param_.getValue("rt_normalization_factor");
    quantification_cutoff_ = (double)param_.getValue("quantification_cutoff");
    write_convex_hull_ = param_.getValue("write_convex_hull").toBool();
    add_up_spectra_ = param_.getValue("add_up_spectra");
    spectrum_addition_method_ = param_.getValue("spectrum_addition_method").toString();
    spacing_for_spectra_resampling_ = param_.getValue("spacing_for_spectra_resampling");
    im_extra_drift_ = (double)param_.getValue("im_extra_drift");
    uis_threshold_sn_ = param_.getValue("uis_threshold_sn");
    uis_threshold_peak_area_ = param_.getValue("uis_threshold_peak_area");
    scoring_model_ = param_.getValue("scoring_model").toString();

    sn_win_len_ = (double)param_.getValue("TransitionGroupPicker:PeakPickerChromatogram:sn_win_len");
    sn_bin_count_ = (unsigned int)param_.getValue("TransitionGroupPicker:PeakPickerChromatogram:sn_bin_count");
    write_log_messages_ = (bool)param_.getValue("TransitionGroupPicker:PeakPickerChromatogram:write_sn_log_messages").toBool();

    diascoring_.setParameters(param_.copy("DIAScoring:", true));

    emgscoring_.setFitterParam(param_.copy("EMGScoring:", true));
    strict_ = (bool)param_.getValue("strict").toBool();
    use_ms1_ion_mobility_ = (bool)param_.getValue("use_ms1_ion_mobility").toBool();

    su_.use_coelution_score_     = param_.getValue("Scores:use_coelution_score").toBool();
    su_.use_shape_score_         = param_.getValue("Scores:use_shape_score").toBool();
    su_.use_rt_score_            = param_.getValue("Scores:use_rt_score").toBool();
    su_.use_library_score_       = param_.getValue("Scores:use_library_score").toBool();
    su_.use_elution_model_score_ = param_.getValue("Scores:use_elution_model_score").toBool();
    su_.use_intensity_score_     = param_.getValue("Scores:use_intensity_score").toBool();
    su_.use_total_xic_score_     = param_.getValue("Scores:use_total_xic_score").toBool();
    su_.use_total_mi_score_      = param_.getValue("Scores:use_total_mi_score").toBool();
    su_.use_nr_peaks_score_      = param_.getValue("Scores:use_nr_peaks_score").toBool();
    su_.use_sn_score_            = param_.getValue("Scores:use_sn_score").toBool();
    su_.use_mi_score_            = param_.getValue("Scores:use_mi_score").toBool();

    su_.use_dia_scores_          = param_.getValue("Scores:use_dia_scores").toBool();
    su_.use_im_scores            = param_.getValue("Scores:use_ion_mobility_scores").toBool();
    su_.use_ms1_correlation      = param_.getValue("Scores:use_ms1_correlation").toBool();
    su_.use_ms1_fullscan         = param_.getValue("Scores:use_ms1_fullscan").toBool();
    su_.use_ms1_mi               = param_.getValue("Scores:use_ms1_mi").toBool();
    su_.use_uis_scores           = param_.getValue("Scores:use_uis_scores").toBool();
    su_.use_ionseries_scores     = param_.getValue("Scores:use_ionseries_scores").toBool();
    su_.use_ms2_isotope_scores   = param_.getValue("Scores:use_ms2_isotope_scores").toBool();
    su_.use_peak_shape_metrics   = param_.getValue("Scores:use_peak_shape_metrics").toBool();
  }

  void MRMFeatureFinderScoring::mapExperimentToTransitionList(const OpenSwath::SpectrumAccessPtr& input,
                                                              const TargetedExpType& transition_exp,
                                                              TransitionGroupMapType& transition_group_map,
                                                              TransformationDescription trafo,
                                                              double rt_extraction_window)
  {
    double rt_min, rt_max, expected_rt;
    trafo.invert();

    std::map<String, int> chromatogram_map;
    Size nr_chromatograms = input->getNrChromatograms();
    for (Size i = 0; i < input->getNrChromatograms(); i++)
    {
      chromatogram_map[input->getChromatogramNativeID(i)] = boost::numeric_cast<int>(i);
    }

    // Iterate through all transitions and store the transition with the
    // corresponding chromatogram in the corresponding transition group
    Size progress = 0;
    startProgress(0, nr_chromatograms, "Mapping transitions to chromatograms ");
    for (Size i = 0; i < transition_exp.getTransitions().size(); i++)
    {
      // get the current transition and try to find the corresponding chromatogram
      const TransitionType* transition = &transition_exp.getTransitions()[i];
      if (chromatogram_map.find(transition->getNativeID()) == chromatogram_map.end())
      {
        OPENMS_LOG_DEBUG << "Error: Transition " + transition->getNativeID() + " from group " +
          transition->getPeptideRef() + " does not have a corresponding chromatogram" << std::endl;
        if (strict_)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                           "Error: Transition " + transition->getNativeID() + " from group " +
                                           transition->getPeptideRef() + " does not have a corresponding chromatogram");
        }
        continue;
      }

      //-----------------------------------
      // Retrieve chromatogram and filter it by the desired RT
      //-----------------------------------
      OpenSwath::ChromatogramPtr cptr = input->getChromatogramById(chromatogram_map[transition->getNativeID()]);
      MSChromatogram chromatogram;

      // Get the expected retention time, apply the RT-transformation
      // (which describes the normalization) and then take the difference.
      // Note that we inverted the transformation in the beginning because
      // we want to transform from normalized to real RTs here and not the
      // other way round.
      if (rt_extraction_window > 0)
      {
        expected_rt = PeptideRefMap_[transition->getPeptideRef()]->rt;
        double de_normalized_experimental_rt = trafo.apply(expected_rt);
        rt_max = de_normalized_experimental_rt + rt_extraction_window;
        rt_min = de_normalized_experimental_rt - rt_extraction_window;
        OpenSwathDataAccessHelper::convertToOpenMSChromatogramFilter(chromatogram, cptr, rt_min, rt_max);
      }
      else
      {
        OpenSwathDataAccessHelper::convertToOpenMSChromatogram(cptr, chromatogram);
      }

      // Check for empty chromatograms (e.g. RT transformation is off)
      if (chromatogram.empty())
      {
        std::cerr << "Error: Could not find any points for chromatogram " + transition->getNativeID() + \
          ". Maybe your retention time transformation is off?" << std::endl;
        if (strict_)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                           "Error: Could not find any points for chromatogram " + transition->getNativeID() + \
                                           ". Maybe your retention time transformation is off?");
        }
      }
      chromatogram.setMetaValue("product_mz", transition->getProductMZ());
      chromatogram.setMetaValue("precursor_mz", transition->getPrecursorMZ());
      Precursor prec; prec.setPosition(transition->getPrecursorMZ());
      Product prod; prod.setMZ(transition->getProductMZ());
      chromatogram.setPrecursor(prec);
      chromatogram.setProduct(prod);
      chromatogram.setNativeID(transition->getNativeID());

      // Create new transition group if there is none for this peptide
      if (transition_group_map.find(transition->getPeptideRef()) == transition_group_map.end())
      {
        MRMTransitionGroupType transition_group;
        transition_group.setTransitionGroupID(transition->getPeptideRef());
        transition_group_map[transition->getPeptideRef()] = transition_group;
      }

      // Now add the transition and the chromatogram to the group
      MRMTransitionGroupType& transition_group = transition_group_map[transition->getPeptideRef()];
      transition_group.addTransition(*transition, transition->getNativeID());
      transition_group.addChromatogram(chromatogram, chromatogram.getNativeID());

      setProgress(++progress);
    }
    endProgress();

    // The assumption is that for each transition that is in the TargetedExperiment we have exactly one chromatogram
    for (TransitionGroupMapType::iterator trgroup_it = transition_group_map.begin(); trgroup_it != transition_group_map.end(); ++trgroup_it)
    {
      if (!trgroup_it->second.isInternallyConsistent())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Error: Could not match all transition to all chromatograms:\nFor chromatogram " + \
                                         trgroup_it->second.getTransitionGroupID() + " I found " + String(trgroup_it->second.getChromatograms().size()) + \
                                         " chromatograms but " + String(trgroup_it->second.getTransitions().size()) + " transitions.");
      }
    }
  }

}

