// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SONARScoring.h>

// peak picking & noise estimation
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/MRMScoring.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h>

#include <boost/range/adaptor/map.hpp>
#include <boost/foreach.hpp>

#define run_identifier "unique_run_identifier"

bool SortDoubleDoublePairFirst(const std::pair<double, double>& left, const std::pair<double, double>& right)
{
  return left.first < right.first;
}

void processFeatureForOutput(OpenMS::Feature& curr_feature, bool write_convex_hull_, double
                             quantification_cutoff_, double& total_intensity, double& total_peak_apices, std::string ms_level)
{
  // Save some space when writing out the featureXML
  if (!write_convex_hull_)
  {
    curr_feature.getConvexHulls().clear();
  }

  // Ensure a unique id is present
  curr_feature.ensureUniqueId();

  // Sum up intensities of the MS2 features
  if (curr_feature.getMZ() > quantification_cutoff_ && ms_level == "MS2")
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
    defaults_.setValue("quantification_cutoff", 0.0, "Cutoff in m/z below which peaks should not be used for quantification any more", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("quantification_cutoff", 0.0);
    defaults_.setValue("write_convex_hull", "false", "Whether to write out all points of all features into the featureXML", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("write_convex_hull", ListUtils::create<String>("true,false"));
    defaults_.setValue("add_up_spectra", 1, "Add up spectra around the peak apex (needs to be a non-even integer)", ListUtils::create<String>("advanced"));
    defaults_.setMinInt("add_up_spectra", 1);
    defaults_.setValue("spacing_for_spectra_resampling", 0.005, "If spectra are to be added, use this spacing to add them up", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("spacing_for_spectra_resampling", 0.0);
    defaults_.setValue("uis_threshold_sn", -1, "S/N threshold to consider identification transition (set to -1 to consider all)");
    defaults_.setValue("uis_threshold_peak_area", 0, "Peak area threshold to consider identification transition (set to -1 to consider all)");

    defaults_.insert("TransitionGroupPicker:", MRMTransitionGroupPicker().getDefaults());

    defaults_.insert("DIAScoring:", DIAScoring().getDefaults());

    defaults_.insert("EMGScoring:", EmgScoring().getDefaults());

    // One can turn on / off each score individually
    Param scores_to_use;
    scores_to_use.setValue("use_shape_score", "true", "Use the shape score (this score measures the similarity in shape of the transitions using a cross-correlation)", ListUtils::create<String>("advanced"));
    scores_to_use.setValidStrings("use_shape_score", ListUtils::create<String>("true,false"));
    scores_to_use.setValue("use_coelution_score", "true", "Use the coelution score (this score measures the similarity in coelution of the transitions using a cross-correlation)", ListUtils::create<String>("advanced"));
    scores_to_use.setValidStrings("use_coelution_score", ListUtils::create<String>("true,false"));
    scores_to_use.setValue("use_rt_score", "true", "Use the retention time score (this score measure the difference in retention time)", ListUtils::create<String>("advanced"));
    scores_to_use.setValidStrings("use_rt_score", ListUtils::create<String>("true,false"));
    scores_to_use.setValue("use_library_score", "true", "Use the library score", ListUtils::create<String>("advanced"));
    scores_to_use.setValidStrings("use_library_score", ListUtils::create<String>("true,false"));
    scores_to_use.setValue("use_elution_model_score", "true", "Use the elution model (EMG) score (this score fits a gaussian model to the peak and checks the fit)", ListUtils::create<String>("advanced"));
    scores_to_use.setValidStrings("use_elution_model_score", ListUtils::create<String>("true,false"));
    scores_to_use.setValue("use_intensity_score", "true", "Use the intensity score", ListUtils::create<String>("advanced"));
    scores_to_use.setValidStrings("use_intensity_score", ListUtils::create<String>("true,false"));
    scores_to_use.setValue("use_nr_peaks_score", "true", "Use the number of peaks score", ListUtils::create<String>("advanced"));
    scores_to_use.setValidStrings("use_nr_peaks_score", ListUtils::create<String>("true,false"));
    scores_to_use.setValue("use_total_xic_score", "true", "Use the total XIC score", ListUtils::create<String>("advanced"));
    scores_to_use.setValidStrings("use_total_xic_score", ListUtils::create<String>("true,false"));
    scores_to_use.setValue("use_sn_score", "true", "Use the SN (signal to noise) score", ListUtils::create<String>("advanced"));
    scores_to_use.setValidStrings("use_sn_score", ListUtils::create<String>("true,false"));
    scores_to_use.setValue("use_dia_scores", "true", "Use the DIA (SWATH) scores. If turned off, will not use fragment ion spectra for scoring.", ListUtils::create<String>("advanced"));
    scores_to_use.setValidStrings("use_dia_scores", ListUtils::create<String>("true,false"));
    scores_to_use.setValue("use_ms1_correlation", "false", "Use the correlation scores with the MS1 elution profiles", ListUtils::create<String>("advanced"));
    scores_to_use.setValidStrings("use_ms1_correlation", ListUtils::create<String>("true,false"));
    scores_to_use.setValue("use_sonar_scores", "false", "Use the scores for SONAR scans (scanning swath)", ListUtils::create<String>("advanced"));
    scores_to_use.setValidStrings("use_sonar_scores", ListUtils::create<String>("true,false"));
    scores_to_use.setValue("use_ms1_fullscan", "false", "Use the full MS1 scan at the peak apex for scoring (ppm accuracy of precursor and isotopic pattern)", ListUtils::create<String>("advanced"));
    scores_to_use.setValidStrings("use_ms1_fullscan", ListUtils::create<String>("true,false"));
    scores_to_use.setValue("use_uis_scores", "false", "Use UIS scores for peptidoform identification ", ListUtils::create<String>("advanced"));
    scores_to_use.setValidStrings("use_uis_scores", ListUtils::create<String>("true,false"));
    defaults_.insert("Scores:", scores_to_use);

    // write defaults into Param object param_
    defaultsToParam_();

    strict_ = true;
  }

  MRMFeatureFinderScoring::~MRMFeatureFinderScoring()
  {
  }

  void MRMFeatureFinderScoring::pickExperiment(PeakMap& chromatograms,
                                               FeatureMap& output, TargetedExperiment& transition_exp_,
                                               TransformationDescription trafo, PeakMap& swath_map)
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

  void MRMFeatureFinderScoring::pickExperiment(OpenSwath::SpectrumAccessPtr input,
                                               FeatureMap& output, OpenSwath::LightTargetedExperiment& transition_exp,
                                               TransformationDescription trafo, 
                                               std::vector<OpenSwath::SwathMap> swath_maps,
                                               TransitionGroupMapType& transition_group_map)
  {
    updateMembers_();

    //
    // Step 1
    //
    // Store the peptide retention times in an intermediate map
    prepareProteinPeptideMaps_(transition_exp);

    // Store the proteins from the input in the output feature map
    std::vector<ProteinHit> protein_hits;
    for (Size i = 0; i < transition_exp.getProteins().size(); i++)
    {
      const ProteinType& prot = transition_exp.getProteins()[i];
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
    for (TransitionGroupMapType::iterator trgroup_it = transition_group_map.begin(); trgroup_it != transition_group_map.end(); ++trgroup_it)
    {
      if (trgroup_it->second.getChromatograms().size() > 0) {counter++; }
    }
    LOG_INFO << "Will analyse " << counter << " peptides with a total of " << transition_exp.getTransitions().size() << " transitions " << std::endl;

    //
    // Step 3
    //
    // Go through all transition groups: first create consensus features, then score them
    Size progress = 0;
    startProgress(0, transition_group_map.size(), "picking peaks");
    for (TransitionGroupMapType::iterator trgroup_it = transition_group_map.begin(); trgroup_it != transition_group_map.end(); ++trgroup_it)
    {

      setProgress(++progress);
      MRMTransitionGroupType& transition_group = trgroup_it->second;
      if (transition_group.getChromatograms().size() == 0 || transition_group.getTransitions().size() == 0)
      {
        continue;
      }

      MRMTransitionGroupPicker trgroup_picker;
      trgroup_picker.setParameters(param_.copy("TransitionGroupPicker:", true));
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

  void MRMFeatureFinderScoring::splitTransitionGroupsDetection_(MRMTransitionGroupType& transition_group, MRMTransitionGroupType& transition_group_detection)
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

  void MRMFeatureFinderScoring::splitTransitionGroupsIdentification_(MRMTransitionGroupType& transition_group, MRMTransitionGroupType& transition_group_identification, MRMTransitionGroupType& transition_group_identification_decoy)
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

  OpenSwath_Scores MRMFeatureFinderScoring::scoreIdentification_(MRMTransitionGroupType& trgr_ident,
                                                                 OpenSwathScoring& scorer,
                                                                 const size_t feature_idx,
                                                                 const std::vector<std::string>& native_ids_detection,
                                                                 const double sn_win_len_,
                                                                 const unsigned int sn_bin_count_,
                                                                 bool write_log_messages,
                                                                 std::vector<OpenSwath::SwathMap> swath_maps)
  {
    MRMFeature idmrmfeature = trgr_ident.getFeaturesMuteable()[feature_idx];
    OpenSwath::IMRMFeature* idimrmfeature;
    idimrmfeature = new MRMFeatureOpenMS(idmrmfeature);  

    std::vector<std::string> native_ids_identification;
    std::vector<OpenSwath::ISignalToNoisePtr> signal_noise_estimators_identification;

    for (Size i = 0; i < trgr_ident.size(); i++)
    {
      OpenSwath::ISignalToNoisePtr snptr(new OpenMS::SignalToNoiseOpenMS< MSChromatogram >(
            trgr_ident.getChromatogram(trgr_ident.getTransitions()[i].getNativeID()),
            sn_win_len_, sn_bin_count_, write_log_messages));
      if (  (snptr->getValueAtRT(idmrmfeature.getRT()) > uis_threshold_sn_) 
            && (idmrmfeature.getFeature(trgr_ident.getTransitions()[i].getNativeID()).getIntensity() > uis_threshold_peak_area_))
      {
        signal_noise_estimators_identification.push_back(snptr);
        native_ids_identification.push_back(trgr_ident.getTransitions()[i].getNativeID());
      }
    }

    OpenSwath_Scores idscores;
    if (native_ids_identification.size() > 0)
    {
      scorer.calculateChromatographicIdScores(idimrmfeature,
                                              native_ids_identification, 
                                              native_ids_detection,
                                              signal_noise_estimators_identification,
                                              idscores);

      std::stringstream ind_transition_names;
      std::stringstream ind_area_intensity;
      std::stringstream ind_apex_intensity;
      std::stringstream ind_log_intensity;
      for (size_t i = 0; i < native_ids_identification.size(); i++)
      {
        if (i != 0)
        {
          ind_transition_names << ";";
          ind_area_intensity << ";";
          ind_apex_intensity << ";";
          ind_log_intensity << ";";
        }
        ind_transition_names << native_ids_identification[i];
        if (idmrmfeature.getFeature(native_ids_identification[i]).getIntensity() > 0)
        {
          ind_area_intensity << idmrmfeature.getFeature(native_ids_identification[i]).getIntensity();
          ind_apex_intensity << idmrmfeature.getFeature(native_ids_identification[i]).getMetaValue("peak_apex_int");
          ind_log_intensity << std::log(idmrmfeature.getFeature(native_ids_identification[i]).getIntensity());
        }
        else
        {
          ind_area_intensity << 0;
          ind_apex_intensity << 0;
          ind_log_intensity << 0;
        }
      }
      idscores.ind_transition_names = ind_transition_names.str();
      idscores.ind_area_intensity = ind_area_intensity.str();
      idscores.ind_apex_intensity = ind_apex_intensity.str();
      idscores.ind_log_intensity = ind_log_intensity.str();
      idscores.ind_num_transitions = native_ids_identification.size();
    }

    // Compute DIA scores only on the identification transitions
    bool swath_present = (!swath_maps.empty() && swath_maps[0].sptr->getNrSpectra() > 0);
    if (swath_present && su_.use_dia_scores_ && native_ids_identification.size() > 0)
    {
      std::stringstream ind_isotope_correlation, ind_isotope_overlap, ind_massdev_score;
      for (size_t i = 0; i < native_ids_identification.size(); i++)
      {
        OpenSwath_Scores tmp_scores;

        scorer.calculateDIAIdScores(idimrmfeature, 
                                    trgr_ident.getTransition(native_ids_identification[i]),
                                    swath_maps, diascoring_, tmp_scores);

        if (i != 0)
        {
          ind_isotope_correlation << ";";
          ind_isotope_overlap << ";";
          ind_massdev_score << ";";
        }
        ind_isotope_correlation << tmp_scores.isotope_correlation;
        ind_isotope_overlap << tmp_scores.isotope_overlap;
        ind_massdev_score << tmp_scores.massdev_score;
      }
      idscores.ind_isotope_correlation = ind_isotope_correlation.str();
      idscores.ind_isotope_overlap = ind_isotope_overlap.str();
      idscores.ind_massdev_score = ind_massdev_score.str();
    }

    delete idimrmfeature;
    return idscores;
  }

  void MRMFeatureFinderScoring::scorePeakgroups(MRMTransitionGroupType& transition_group,
                                                TransformationDescription& trafo, 
                                                std::vector<OpenSwath::SwathMap> swath_maps,
                                                FeatureMap& output, 
                                                bool ms1only)
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

    double sn_win_len_ = (double)param_.getValue("TransitionGroupPicker:PeakPickerMRM:sn_win_len");
    unsigned int sn_bin_count_ = (unsigned int)param_.getValue("TransitionGroupPicker:PeakPickerMRM:sn_bin_count");
    bool write_log_messages = (bool)param_.getValue("TransitionGroupPicker:PeakPickerMRM:write_sn_log_messages").toBool();
    // currently we cannot do much about the log messages and they mostly occur in decoy transition signals
    for (Size k = 0; k < transition_group_detection.getChromatograms().size(); k++)
    {
      OpenSwath::ISignalToNoisePtr snptr(new OpenMS::SignalToNoiseOpenMS< MSChromatogram >(
            transition_group_detection.getChromatograms()[k], sn_win_len_, sn_bin_count_, write_log_messages));
      signal_noise_estimators.push_back(snptr);
    }

    // skip MS1 noise estimator if we perform fragment ion analysis
    std::vector<OpenSwath::ISignalToNoisePtr> ms1_signal_noise_estimators;
    if (ms1only)
    {
      for (Size k = 0; k < transition_group_detection.getPrecursorChromatograms().size(); k++)
      {
        OpenSwath::ISignalToNoisePtr snptr(new OpenMS::SignalToNoiseOpenMS< MSChromatogram >(
              transition_group_detection.getPrecursorChromatograms()[k], sn_win_len_, sn_bin_count_, write_log_messages));
        ms1_signal_noise_estimators.push_back(snptr);
      }
    }

    // get the expected rt value for this compound
    const PeptideType* pep = PeptideRefMap_[transition_group_detection.getTransitionGroupID()];
    double expected_rt = pep->rt;
    TransformationDescription newtr = trafo;
    newtr.invert();
    expected_rt = newtr.apply(expected_rt);

    OpenSwathScoring scorer;
    scorer.initialize(rt_normalization_factor_, add_up_spectra_, spacing_for_spectra_resampling_, su_);

    size_t feature_idx = 0;
    // Go through all peak groups (found MRM features) and score them
    for (std::vector<MRMFeature>::iterator mrmfeature = transition_group_detection.getFeaturesMuteable().begin();
         mrmfeature != transition_group_detection.getFeaturesMuteable().end(); ++mrmfeature)
    {
      OpenSwath::IMRMFeature* imrmfeature;
      imrmfeature = new MRMFeatureOpenMS(*mrmfeature);

      LOG_DEBUG << "scoring feature " << (*mrmfeature) << " == " << mrmfeature->getMetaValue("PeptideRef") <<
        " [ expected RT " << PeptideRefMap_[mrmfeature->getMetaValue("PeptideRef")]->rt << " / " << expected_rt << " ]" <<
        " with " << transition_group_detection.size()  << " nr transitions and nr chromats " << 
        transition_group_detection.getChromatograms().size() << std::endl;

      int group_size = boost::numeric_cast<int>(transition_group_detection.size());
      if (group_size == 0 && !ms1only)
      {
        delete imrmfeature; // free resources before continuing
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "Error: Transition group " + transition_group_detection.getTransitionGroupID() + 
                                         " has no chromatograms.");
      }
      if (group_size < 2 && !ms1only)
      {
        LOG_ERROR << "Error: Transition group " << transition_group_detection.getTransitionGroupID()
                  << " has only one chromatogram." << std::endl;
        delete imrmfeature; // free resources before continuing
        continue;
      }

      bool swath_present = (!swath_maps.empty() && swath_maps[0].sptr->getNrSpectra() > 0);
      bool sonar_present = (swath_maps.size() > 1);
      double xx_lda_prescore;
      double precursor_mz(-1);

      if (ms1only)
      {
        ///////////////////////////////////
        // Call the scoring for MS1 only
        ///////////////////////////////////

        OpenSwath_Scores scores;
        precursor_mz = mrmfeature->getMZ();

        // S/N scores
        OpenSwath::MRMScoring mrmscore_;
        scores.sn_ratio = mrmscore_.calcSNScore(imrmfeature, ms1_signal_noise_estimators);
        // everything below S/N 1 can be set to zero (and the log safely applied)
        if (scores.sn_ratio < 1) { scores.log_sn_score = 0; }
        else { scores.log_sn_score = std::log(scores.sn_ratio); }
        if (su_.use_sn_score_) 
        { 
          mrmfeature->addScore("sn_ratio", scores.sn_ratio);
          mrmfeature->addScore("var_log_sn_score", scores.log_sn_score); 
          // compute subfeature log-SN values
          for (Size k = 0; k < transition_group_detection.getPrecursorChromatograms().size(); k++)
          {
            Feature & f = mrmfeature->getPrecursorFeature(transition_group_detection.getPrecursorChromatograms()[k].getNativeID());
            double sn_value = ms1_signal_noise_estimators[k]->getValueAtRT(mrmfeature->getRT());
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
          mrmfeature->addScore("delta_rt", mrmfeature->getRT() - expected_rt);
          mrmfeature->addScore("assay_rt", expected_rt);
          mrmfeature->addScore("norm_RT", scores.normalized_experimental_rt);
          mrmfeature->addScore("rt_score", scores.raw_rt_score);
          mrmfeature->addScore("var_norm_rt_score", scores.norm_rt_score);
        }

        // full spectra scores 
        if (ms1_map_ && ms1_map_->getNrSpectra() > 0 && mrmfeature->getMZ() > 0) 
        {
          scorer.calculatePrecursorDIAScores(ms1_map_, diascoring_, precursor_mz, imrmfeature->getRT(), *pep, scores);
        }
        if (su_.use_ms1_fullscan)
        {
          mrmfeature->addScore("var_ms1_ppm_diff", scores.ms1_ppm_score);
          mrmfeature->addScore("var_ms1_isotope_correlation", scores.ms1_isotope_correlation);
          mrmfeature->addScore("var_ms1_isotope_overlap", scores.ms1_isotope_overlap);
        }
        xx_lda_prescore = -scores.calculate_lda_prescore(scores);
        mrmfeature->addScore("main_var_xx_lda_prelim_score", xx_lda_prescore);
        mrmfeature->setOverallQuality(xx_lda_prescore);
      }
      else
      {

        ///////////////////////////////////
        // Call the scoring for fragment ions
        ///////////////////////////////////

        std::vector<double> normalized_library_intensity;
        transition_group_detection.getLibraryIntensity(normalized_library_intensity);
        OpenSwath::Scoring::normalize_sum(&normalized_library_intensity[0], boost::numeric_cast<int>(normalized_library_intensity.size()));
        std::vector<std::string> native_ids_detection;
        std::string precursor_id;
        for (Size i = 0; i < transition_group_detection.size(); i++)
        {
          native_ids_detection.push_back(transition_group_detection.getTransitions()[i].getNativeID());
        }
        for (Size i = 0; i < transition_group_detection.getPrecursorChromatograms().size(); i++)
        {
          // try to identify the correct precursor native id
          String precursor_chrom_id = transition_group_detection.getPrecursorChromatograms()[i].getNativeID();
          if (precursor_chrom_id.hasSuffix("_Precursor_i0")) {precursor_id = precursor_chrom_id;}
        }

        OpenSwath_Scores scores;
        scorer.calculateChromatographicScores(imrmfeature, native_ids_detection, precursor_id, normalized_library_intensity,
                                              signal_noise_estimators, scores);

        double normalized_experimental_rt = trafo.apply(imrmfeature->getRT());
        scorer.calculateLibraryScores(imrmfeature, transition_group_detection.getTransitions(), *pep, normalized_experimental_rt, scores);
        if (swath_present && su_.use_dia_scores_)
        {
          scorer.calculateDIAScores(imrmfeature, transition_group_detection.getTransitions(),
                                    swath_maps, ms1_map_, diascoring_, *pep, scores);
        }
        if (sonar_present && su_.use_sonar_scores)
        {
          sonarscoring_.computeSonarScores(imrmfeature, transition_group_detection.getTransitions(), swath_maps, scores);
        }

        if (su_.use_uis_scores && transition_group_identification.getTransitions().size() > 0)
        {
          OpenSwath_Scores idscores = scoreIdentification_(transition_group_identification, 
                                                           scorer, feature_idx,
                                                           native_ids_detection,
                                                           sn_win_len_,
                                                           sn_bin_count_,
                                                           write_log_messages,
                                                           swath_maps);

          mrmfeature->setMetaValue("id_target_transition_names", idscores.ind_transition_names);
          mrmfeature->addScore("id_target_num_transitions", idscores.ind_num_transitions);
          mrmfeature->setMetaValue("id_target_area_intensity", idscores.ind_area_intensity);
          mrmfeature->setMetaValue("id_target_apex_intensity", idscores.ind_apex_intensity);
          mrmfeature->setMetaValue("id_target_transition_names", idscores.ind_transition_names);
          mrmfeature->setMetaValue("id_target_ind_log_intensity", idscores.ind_log_intensity);
          mrmfeature->setMetaValue("id_target_ind_xcorr_coelution", idscores.ind_xcorr_coelution_score);
          mrmfeature->setMetaValue("id_target_ind_xcorr_shape", idscores.ind_xcorr_shape_score);
          mrmfeature->setMetaValue("id_target_ind_log_sn_score", idscores.ind_log_sn_score);
          mrmfeature->setMetaValue("id_target_ind_isotope_correlation", idscores.ind_isotope_correlation);
          mrmfeature->setMetaValue("id_target_ind_isotope_overlap", idscores.ind_isotope_overlap);
          mrmfeature->setMetaValue("id_target_ind_massdev_score", idscores.ind_massdev_score);
        }

        if (su_.use_uis_scores && transition_group_identification_decoy.getTransitions().size() > 0)
        {
          OpenSwath_Scores idscores = scoreIdentification_(transition_group_identification_decoy, 
                                                           scorer, feature_idx,
                                                           native_ids_detection,
                                                           sn_win_len_,
                                                           sn_bin_count_,
                                                           write_log_messages,
                                                           swath_maps);

          mrmfeature->setMetaValue("id_decoy_transition_names", idscores.ind_transition_names);
          mrmfeature->addScore("id_decoy_num_transitions", idscores.ind_num_transitions);
          mrmfeature->setMetaValue("id_decoy_area_intensity", idscores.ind_area_intensity);
          mrmfeature->setMetaValue("id_decoy_apex_intensity", idscores.ind_apex_intensity);
          mrmfeature->setMetaValue("id_decoy_ind_log_intensity", idscores.ind_log_intensity);
          mrmfeature->setMetaValue("id_decoy_ind_xcorr_coelution", idscores.ind_xcorr_coelution_score);
          mrmfeature->setMetaValue("id_decoy_ind_xcorr_shape", idscores.ind_xcorr_shape_score);
          mrmfeature->setMetaValue("id_decoy_ind_log_sn_score", idscores.ind_log_sn_score);
          mrmfeature->setMetaValue("id_decoy_ind_isotope_correlation", idscores.ind_isotope_correlation);
          mrmfeature->setMetaValue("id_decoy_ind_isotope_overlap", idscores.ind_isotope_overlap);
          mrmfeature->setMetaValue("id_decoy_ind_massdev_score", idscores.ind_massdev_score);
        }

        if (su_.use_coelution_score_)
        {
          mrmfeature->addScore("var_xcorr_coelution", scores.xcorr_coelution_score);
          mrmfeature->addScore("var_xcorr_coelution_weighted", scores.weighted_coelution_score);
        }
        if (su_.use_shape_score_)
        {
          mrmfeature->addScore("var_xcorr_shape", scores.xcorr_shape_score);
          mrmfeature->addScore("var_xcorr_shape_weighted", scores.weighted_xcorr_shape);
        }
        if (su_.use_library_score_)
        {
          mrmfeature->addScore("var_library_corr", scores.library_corr);
          mrmfeature->addScore("var_library_rmsd", scores.library_norm_manhattan);
          mrmfeature->addScore("var_library_sangle", scores.library_sangle);
          mrmfeature->addScore("var_library_rootmeansquare", scores.library_rootmeansquare);
          mrmfeature->addScore("var_library_manhattan", scores.library_manhattan);
          mrmfeature->addScore("var_library_dotprod", scores.library_dotprod);
        }
        if (su_.use_rt_score_)
        {
          mrmfeature->addScore("delta_rt", mrmfeature->getRT() - expected_rt);
          mrmfeature->addScore("assay_rt", expected_rt);
          mrmfeature->addScore("norm_RT", scores.normalized_experimental_rt);
          mrmfeature->addScore("rt_score", scores.raw_rt_score);
          mrmfeature->addScore("var_norm_rt_score", scores.norm_rt_score);
        }
        // TODO do we really want these intensity scores ?
        if (su_.use_intensity_score_) { mrmfeature->addScore("var_intensity_score", mrmfeature->getIntensity() / (double)mrmfeature->getMetaValue("total_xic")); }
        if (su_.use_total_xic_score_) { mrmfeature->addScore("total_xic", (double)mrmfeature->getMetaValue("total_xic")); }
        if (su_.use_nr_peaks_score_) { mrmfeature->addScore("nr_peaks", scores.nr_peaks); }
        if (su_.use_sn_score_)
        {
          mrmfeature->addScore("sn_ratio", scores.sn_ratio);
          mrmfeature->addScore("var_log_sn_score", scores.log_sn_score);
          // compute subfeature log-SN values
          for (Size k = 0; k < transition_group_detection.getChromatograms().size(); k++)
          {
            Feature & f = mrmfeature->getFeature(transition_group_detection.getChromatograms()[k].getNativeID());
            double sn_value = signal_noise_estimators[k]->getValueAtRT(mrmfeature->getRT());
            if (sn_value < 1) {sn_value = 1.0;}
            f.setMetaValue("logSN", std::log(sn_value));
          }
        }

        // TODO get it working with imrmfeature
        if (su_.use_elution_model_score_)
        {
          scores.elution_model_fit_score = emgscoring_.calcElutionFitScore((*mrmfeature), transition_group_detection);
          mrmfeature->addScore("var_elution_model_fit_score", scores.elution_model_fit_score);
        }

        xx_lda_prescore = -scores.calculate_lda_prescore(scores);
        if (!swath_present)
        {
          mrmfeature->addScore("main_var_xx_lda_prelim_score", xx_lda_prescore);
          mrmfeature->setOverallQuality(xx_lda_prescore);
        }
        else
        {
          mrmfeature->addScore("xx_lda_prelim_score", xx_lda_prescore);
        }

        // Add the DIA / SWATH scores
        if (swath_present && su_.use_dia_scores_)
        {
          mrmfeature->addScore("var_isotope_correlation_score", scores.isotope_correlation);
          mrmfeature->addScore("var_isotope_overlap_score", scores.isotope_overlap);
          mrmfeature->addScore("var_massdev_score", scores.massdev_score);
          mrmfeature->addScore("var_massdev_score_weighted", scores.weighted_massdev_score);
          mrmfeature->addScore("var_bseries_score", scores.bseries_score);
          mrmfeature->addScore("var_yseries_score", scores.yseries_score);
          mrmfeature->addScore("var_dotprod_score", scores.dotprod_score_dia);
          mrmfeature->addScore("var_manhatt_score", scores.manhatt_score_dia);
          if (su_.use_ms1_correlation)
          {
            mrmfeature->addScore("var_ms1_xcorr_shape", scores.xcorr_ms1_shape_score);
            mrmfeature->addScore("var_ms1_xcorr_coelution", scores.xcorr_ms1_coelution_score);
          }
          if (su_.use_ms1_fullscan)
          {
            mrmfeature->addScore("var_ms1_ppm_diff", scores.ms1_ppm_score);
            mrmfeature->addScore("var_ms1_isotope_correlation", scores.ms1_isotope_correlation);
            mrmfeature->addScore("var_ms1_isotope_overlap", scores.ms1_isotope_overlap);
          }

          double xx_swath_prescore = -scores.calculate_swath_lda_prescore(scores);
          mrmfeature->addScore("main_var_xx_swath_prelim_score", xx_swath_prescore);
          mrmfeature->setOverallQuality(xx_swath_prescore);
        }

        precursor_mz = transition_group_detection.getTransitions()[0].getPrecursorMZ();

        if (sonar_present && su_.use_sonar_scores)
        {

          // set all scores less than 1 to zero (do not over-punish large negative scores)
          double log_sn = 0;
          if (scores.sonar_sn > 1) log_sn = std::log(scores.sonar_sn);
          double log_trend = 0;
          if (scores.sonar_trend > 1) log_trend = std::log(scores.sonar_trend);
          double log_diff = 0;
          if (scores.sonar_diff > 1) log_diff = std::log(scores.sonar_diff);

          mrmfeature->addScore("var_sonar_lag", scores.sonar_lag);
          mrmfeature->addScore("var_sonar_shape", scores.sonar_shape);
          mrmfeature->addScore("var_sonar_log_sn", log_sn);
          mrmfeature->addScore("var_sonar_log_diff", log_diff);
          mrmfeature->addScore("var_sonar_log_trend", log_trend);
          mrmfeature->addScore("var_sonar_rsq", scores.sonar_rsq);
        }
      }

      ///////////////////////////////////////////////////////////////////////////
      // add the peptide hit information to the feature
      ///////////////////////////////////////////////////////////////////////////
      PeptideIdentification pep_id_ = PeptideIdentification();
      PeptideHit pep_hit_ = PeptideHit();

      if (pep->getChargeState() != 0)
      {
        pep_hit_.setCharge(pep->getChargeState());
      }
      pep_hit_.setScore(xx_lda_prescore);
      if (swath_present)
      {
        pep_hit_.setScore(mrmfeature->getScore("xx_swath_prelim_score"));
      }

      if (pep->isPeptide())
      {
        pep_hit_.setSequence(AASequence::fromString(pep->sequence));
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

      mrmfeature->getPeptideIdentifications().push_back(pep_id_);
      mrmfeature->ensureUniqueId();

      mrmfeature->setMetaValue("PrecursorMZ", precursor_mz);

      // Prepare the subordinates for the mrmfeature (process all current
      // features and then append all precursor subordinate features)
      std::vector<Feature> allFeatures = mrmfeature->getFeatures();
      double total_intensity = 0, total_peak_apices = 0;
      for (std::vector<Feature>::iterator f_it = allFeatures.begin(); f_it != allFeatures.end(); ++f_it)
      {
        processFeatureForOutput(*f_it, write_convex_hull_, quantification_cutoff_, total_intensity, total_peak_apices, "MS2");
      }
      // Also append data for MS1 precursors
      std::vector<String> precursors_ids;
      mrmfeature->getPrecursorFeatureIDs(precursors_ids);
      for (std::vector<String>::iterator id_it = precursors_ids.begin(); id_it != precursors_ids.end(); ++id_it)
      {
        Feature curr_feature = mrmfeature->getPrecursorFeature(*id_it);
        if (pep->getChargeState() != 0)
        {
          curr_feature.setCharge(pep->getChargeState());
        }
        processFeatureForOutput(curr_feature, write_convex_hull_, quantification_cutoff_, total_intensity, total_peak_apices, "MS1");
        allFeatures.push_back(curr_feature);
      }
      mrmfeature->setSubordinates(allFeatures); // add all the subfeatures as subordinates

      // overwrite the reported intensities with those above the m/z cutoff
      mrmfeature->setIntensity(total_intensity);
      mrmfeature->setMetaValue("peak_apices_sum", total_peak_apices);
      feature_list.push_back((*mrmfeature));

      delete imrmfeature;
      feature_idx++;
    }

    // Order by quality
    std::sort(feature_list.begin(), feature_list.end(), OpenMS::Feature::OverallQualityLess());
    std::reverse(feature_list.begin(), feature_list.end());

    for (Size i = 0; i < feature_list.size(); i++)
    {
      if (stop_report_after_feature_ >= 0 && i >= (Size)stop_report_after_feature_) {break;}
      output.push_back(feature_list[i]);
    }
    transition_group = transition_group_detection;
  }

  void MRMFeatureFinderScoring::updateMembers_()
  {
    stop_report_after_feature_ = (int)param_.getValue("stop_report_after_feature");
    rt_extraction_window_ = (double)param_.getValue("rt_extraction_window");
    rt_normalization_factor_ = (double)param_.getValue("rt_normalization_factor");
    quantification_cutoff_ = (double)param_.getValue("quantification_cutoff");
    write_convex_hull_ = param_.getValue("write_convex_hull").toBool();
    add_up_spectra_ = param_.getValue("add_up_spectra");
    spacing_for_spectra_resampling_ = param_.getValue("spacing_for_spectra_resampling");
    uis_threshold_sn_ = param_.getValue("uis_threshold_sn");
    uis_threshold_peak_area_ = param_.getValue("uis_threshold_peak_area");

    // set SONAR values
    Param p = sonarscoring_.getDefaults();
    p.setValue("dia_extraction_window", param_.getValue("DIAScoring:dia_extraction_window"));
    p.setValue("dia_centroided", param_.getValue("DIAScoring:dia_centroided"));
    sonarscoring_.setParameters(p);

    diascoring_.setParameters(param_.copy("DIAScoring:", true));
    emgscoring_.setFitterParam(param_.copy("EmgScoring:", true));

    su_.use_coelution_score_     = param_.getValue("Scores:use_coelution_score").toBool();
    su_.use_shape_score_         = param_.getValue("Scores:use_shape_score").toBool();
    su_.use_rt_score_            = param_.getValue("Scores:use_rt_score").toBool();
    su_.use_library_score_       = param_.getValue("Scores:use_library_score").toBool();
    su_.use_elution_model_score_ = param_.getValue("Scores:use_elution_model_score").toBool();
    su_.use_intensity_score_     = param_.getValue("Scores:use_intensity_score").toBool();
    su_.use_total_xic_score_     = param_.getValue("Scores:use_total_xic_score").toBool();
    su_.use_nr_peaks_score_      = param_.getValue("Scores:use_nr_peaks_score").toBool();
    su_.use_sn_score_            = param_.getValue("Scores:use_sn_score").toBool();
    su_.use_dia_scores_          = param_.getValue("Scores:use_dia_scores").toBool();
    su_.use_sonar_scores         = param_.getValue("Scores:use_sonar_scores").toBool();
    su_.use_ms1_correlation      = param_.getValue("Scores:use_ms1_correlation").toBool();
    su_.use_ms1_fullscan         = param_.getValue("Scores:use_ms1_fullscan").toBool();
    su_.use_uis_scores           = param_.getValue("Scores:use_uis_scores").toBool();
  }

  void MRMFeatureFinderScoring::mapExperimentToTransitionList(OpenSwath::SpectrumAccessPtr input,
                                                              TargetedExpType& transition_exp, TransitionGroupMapType& transition_group_map,
                                                              TransformationDescription trafo, double rt_extraction_window)
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
        std::cerr << "Error: Transition " + transition->getNativeID() + " from group " +
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
      if (!trgroup_it->second.isInternallyConsistent() || trgroup_it->second.getChromatograms().size() != trgroup_it->second.getTransitions().size())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Error: Could not match all transition to all chromatograms:\nFor chromatogram " + \
                                         trgroup_it->second.getTransitionGroupID() + " I found " + String(trgroup_it->second.getChromatograms().size()) + \
                                         " chromatograms but " + String(trgroup_it->second.getTransitions().size()) + " transitions.");
      }
    }
  }

}

