// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/TargetedSpectraExtractor.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>
#include <OpenMS/PROCESSING/SMOOTHING/GaussFilter.h>
#include <OpenMS/PROCESSING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/PROCESSING/DEISOTOPING/Deisotoper.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/ANALYSIS/ID/AccurateMassSearchEngine.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>

namespace OpenMS
{
  TargetedSpectraExtractor::TargetedSpectraExtractor() :
    DefaultParamHandler("TargetedSpectraExtractor")
  {
    getDefaultParameters(defaults_);

    subsections_.emplace_back("SavitzkyGolayFilter");
    defaults_.setValue("SavitzkyGolayFilter:frame_length", 15);
    defaults_.setValue("SavitzkyGolayFilter:polynomial_order", 3);

    subsections_.emplace_back("GaussFilter");
    defaults_.setValue("GaussFilter:gaussian_width", 0.2);

    subsections_.emplace_back("PeakPickerHiRes");
    defaults_.setValue("PeakPickerHiRes:signal_to_noise", 1.0);

    defaults_.insert("AccurateMassSearchEngine:", AccurateMassSearchEngine().getDefaults());
    defaults_.setValue("AccurateMassSearchEngine:keep_unidentified_masses", "false");
    defaults_.setValidStrings("AccurateMassSearchEngine:keep_unidentified_masses", {"true","false"});

    // write defaults into Param object param_
    defaultsToParam_();
    updateMembers_();
  }

  void TargetedSpectraExtractor::updateMembers_()
  {
    rt_window_ = (double)param_.getValue("rt_window");
    min_select_score_ = (double)param_.getValue("min_select_score");
    mz_tolerance_ = (double)param_.getValue("mz_tolerance");
    mz_unit_is_Da_ = param_.getValue("mz_unit_is_Da").toBool();
    use_gauss_ = param_.getValue("use_gauss").toBool();
    peak_height_min_ = (double)param_.getValue("peak_height_min");
    peak_height_max_ = (double)param_.getValue("peak_height_max");
    fwhm_threshold_ = (double)param_.getValue("fwhm_threshold");
    tic_weight_ = (double)param_.getValue("tic_weight");
    fwhm_weight_ = (double)param_.getValue("fwhm_weight");
    snr_weight_ = (double)param_.getValue("snr_weight");
    top_matches_to_report_ = (Size)param_.getValue("top_matches_to_report");
    min_match_score_ = (double)param_.getValue("min_match_score");
    
    min_fragment_mz_ = (double) param_.getValue("min_fragment_mz");
    max_fragment_mz_ = (double) param_.getValue("max_fragment_mz");
    relative_allowable_product_mass_ = (double) param_.getValue("relative_allowable_product_mass");

    deisotoping_use_deisotoper_ = param_.getValue("deisotoping:use_deisotoper").toBool();
    deisotoping_fragment_tolerance_ = (double) param_.getValue("deisotoping:fragment_tolerance");
    deisotoping_fragment_unit_ = param_.getValue("deisotoping:fragment_unit").toString();
    deisotoping_min_charge_ = param_.getValue("deisotoping:min_charge");
    deisotoping_max_charge_ = param_.getValue("deisotoping:max_charge");
    deisotoping_min_isopeaks_ = param_.getValue("deisotoping:min_isopeaks");
    deisotoping_max_isopeaks_ = param_.getValue("deisotoping:max_isopeaks");
    deisotoping_keep_only_deisotoped_ = param_.getValue("deisotoping:keep_only_deisotoped").toBool();
    deisotoping_annotate_charge_ = param_.getValue("deisotoping:annotate_charge").toBool();

    max_precursor_mass_threashold_ = (double) param_.getValue("max_precursor_mass_threashold");
  }

  void TargetedSpectraExtractor::getDefaultParameters(Param& params) const
  {
    params.clear();

    params.setValue(
      "rt_window",
      30.0,
      "Precursor Retention Time window used during the annotation phase.\n"
      "For each transition in the target list, annotateSpectra() looks for "
      "the first spectrum whose RT time falls within the RT Window, whose "
      "left and right limits are computed at each analyzed spectrum.\n"
      "Also the spectrum's percursor MZ is checked against the transition MZ."
    );

    params.setValue(
      "min_select_score",
      0.7,
      "Used in selectSpectra(), after the spectra have been assigned a score.\n"
      "Remained transitions will have at least one spectrum assigned.\n"
      "Each spectrum needs to have a score >= min_select_score_ to be valid, "
      "otherwise it gets filtered out."
    );
    params.setMinFloat("min_select_score", 0.0);

    params.setValue(
      "mz_tolerance",
      0.1,
      "Precursor MZ tolerance used during the annotation phase.\n"
      "For each transition in the target list, annotateSpectra() looks for "
      "the first spectrum whose precursor MZ is close enough (+-mz_tolerance_) "
      "to the transition's MZ.\n"
      "Also the spectrum's precursor RT is checked against the transition RT."
    );

    params.setValue("mz_unit_is_Da", "true", "Unit to use for mz_tolerance_ and fwhm_threshold_: true for Da, false for ppm.");
    params.setValidStrings("mz_unit_is_Da", {"false","true"});

    params.setValue("use_gauss", "true", "Use Gaussian filter for smoothing (alternative is Savitzky-Golay filter)");
    params.setValidStrings("use_gauss", {"false","true"});

    params.setValue("peak_height_min", 0.0, "Used in pickSpectrum(), a peak's intensity needs to be >= peak_height_min_ for it to be picked.");
    params.setMinFloat("peak_height_min", 0.0);
    params.setValue("peak_height_max", std::numeric_limits<double>::max(), "Used in pickSpectrum(), a peak's intensity needs to be <= peak_height_max_ for it to be picked.");
    params.setMinFloat("peak_height_max", 0.0);
    params.setValue("fwhm_threshold", 0.0, "Used in pickSpectrum(), a peak's FWHM needs to be >= fwhm_threshold_ for it to be picked.");
    params.setMinFloat("fwhm_threshold", 0.0);

    params.setValue("tic_weight", 1.0, "TIC weight when scoring spectra.");
    params.setMinFloat("tic_weight", 0.0);
    params.setValue("fwhm_weight", 1.0, "FWHM weight when scoring spectra.");
    params.setMinFloat("fwhm_weight", 0.0);
    params.setValue("snr_weight", 1.0, "SNR weight when scoring spectra.");
    params.setMinFloat("snr_weight", 0.0);

    params.setValue(
      "top_matches_to_report",
      5,
      "The number of matches to output from `matchSpectrum()`. "
      "These will be the matches of highest scores, sorted in descending order."
    );
    params.setMinInt("top_matches_to_report", 1);

    params.setValue(
      "min_match_score",
      0.8,
      "Minimum score for a match to be considered valid in `matchSpectrum()`."
    );
    params.setMinFloat("min_match_score", 0.0);
    params.setMaxFloat("min_match_score", 1.0);

    params.setValue("min_fragment_mz", 0.0, "Minimal m/z of a fragment ion choosen as a transition");
    params.setValue("max_fragment_mz", 2000.0, "Maximal m/z of a fragment ion choosen as a transition");
    params.setValue("relative_allowable_product_mass", 10.0, "Threshold m/z of a product relatively to the precurosor m/z (can be negative)");

    params.setValue("deisotoping:use_deisotoper", "false", "Use Deisotoper (if no fragment annotation is used)");
    params.setValue("deisotoping:fragment_tolerance", 1.0, "Tolerance used to match isotopic peaks");
    params.setValue("deisotoping:fragment_unit", "ppm", "Unit of the fragment tolerance");
    params.setValidStrings("deisotoping:fragment_unit", {"ppm","Da"});
    params.setValue("deisotoping:min_charge", 1, "The minimum charge considered");
    params.setMinInt("deisotoping:min_charge", 1);
    params.setValue("deisotoping:max_charge", 1, "The maximum charge considered");
    params.setMinInt("deisotoping:max_charge", 1);
    params.setValue("deisotoping:min_isopeaks", 2, "The minimum number of isotopic peaks (at least 2) required for an isotopic cluster");
    params.setMinInt("deisotoping:min_isopeaks", 2);
    params.setValue("deisotoping:max_isopeaks", 3, "The maximum number of isotopic peaks (at least 2) considered for an isotopic cluster");
    params.setMinInt("deisotoping:max_isopeaks", 3);
    params.setValue("deisotoping:keep_only_deisotoped", "false", "Only monoisotopic peaks of fragments with isotopic pattern are retained");
    params.setValue("deisotoping:annotate_charge", "false", "Annotate the charge to the peaks");

    params.setValue("max_precursor_mass_threashold", 10.0, "Tolerance used to set intensity to zero for peaks with mz higher than precursor mz");
  }

  void TargetedSpectraExtractor::annotateSpectra(
      const std::vector<MSSpectrum>& spectra,
      const FeatureMap& ms1_features,
      FeatureMap& ms2_features,
      std::vector<MSSpectrum>& annotated_spectra) const
  {
    annotated_spectra.clear();
    for (const auto& spectrum : spectra)
    {
      if (spectrum.getMSLevel() == 1)
      {
        continue; // we want to annotate MS2 spectra only
      }

      const double spectrum_rt = spectrum.getRT();
      const std::vector<Precursor>& precursors = spectrum.getPrecursors();
      if (precursors.empty())
      {
        OPENMS_LOG_WARN << "annotateSpectra(): No precursor MZ found. Setting spectrum_mz to 0." << std::endl;
      }
      const double spectrum_mz = precursors.empty() ? 0.0 : precursors.front().getMZ();

      // Lambda to check the mz/rt thresholds
      auto checkRtAndMzTol = [](const double& spectrum_mz, const double& spectrum_rt, 
        const double& target_mz, const double& target_rt, const double& mz_window, const double& rt_window) 
      {
        const double rt_left_lim = spectrum_rt - rt_window / 2.0;
        const double rt_right_lim = spectrum_rt + rt_window / 2.0;
        const double mz_left_lim = spectrum_mz - mz_window / 2.0;
        const double mz_right_lim = spectrum_mz + mz_window / 2.0;
        if (spectrum_mz != 0.0)
        {
          return target_rt >= rt_left_lim && target_rt <= rt_right_lim && target_mz >= mz_left_lim && target_mz <= mz_right_lim;
        }
        else
        {
          return target_rt >= rt_left_lim && target_rt <= rt_right_lim;
        }
      };

      // Lambda to create the feature
      auto construct_feature = [checkRtAndMzTol, spectrum_rt, spectrum_mz, &ms2_features, &annotated_spectra, &spectrum](
        const OpenMS::Feature& feature, const double& mz_tol, const double& rt_win)
      {
        const auto& peptide_ref_s = feature.getMetaValue("PeptideRef");
        const auto& native_id_s = feature.getMetaValue("native_id");
        // check for null annotations resulting from unnanotated features
        if (peptide_ref_s != "null")
        {
          const double target_mz = feature.getMZ();
          const double target_rt = feature.getRT();
          if (checkRtAndMzTol(spectrum_mz, spectrum_rt, target_mz, target_rt, mz_tol, rt_win))
          {
            OPENMS_LOG_DEBUG << "annotateSpectra(): " << peptide_ref_s << "]";
            OPENMS_LOG_DEBUG << " (target_rt: " << target_rt << ") (target_mz: " << target_mz << ")" << std::endl;
            MSSpectrum annotated_spectrum = spectrum;
            annotated_spectrum.setName(peptide_ref_s);
            annotated_spectra.push_back(std::move(annotated_spectrum));
            // fill the ms2 features map
            Feature ms2_feature;
            ms2_feature.setUniqueId();
            ms2_feature.setRT(spectrum_rt);
            ms2_feature.setMZ(spectrum_mz);
            ms2_feature.setIntensity(feature.getIntensity());
            ms2_feature.setMetaValue("native_id", native_id_s);
            ms2_feature.setMetaValue("PeptideRef", peptide_ref_s);
            ms2_features.push_back(std::move(ms2_feature));
          }
        }
      };

      for (const auto& feature : ms1_features)
      {
        if (!feature.getSubordinates().empty())
        {
          // iterate through the subordinate level
          for (const auto& subordinate : feature.getSubordinates())
          {
            construct_feature(subordinate, mz_tolerance_, rt_window_);
          }
        }
        else
        {
          construct_feature(feature, mz_tolerance_, rt_window_);
        }
      }
    }
  }

  void TargetedSpectraExtractor::searchSpectrum(OpenMS::FeatureMap& feat_map, 
    OpenMS::FeatureMap& feat_map_output, bool add_unidentified_features) const
  {
    OpenMS::AccurateMassSearchEngine ams;
    OpenMS::MzTab output;
    ams.setParameters(param_.copy("AccurateMassSearchEngine:", true));
    ams.init();
    ams.run(feat_map, output);
    // Remake the feature map replacing the peptide hits as features/sub-features
    feat_map_output.clear();
    for (const OpenMS::Feature& feature : feat_map)
    {
      const auto& peptide_identifications = feature.getPeptideIdentifications();
      if (peptide_identifications.size())
      {
        for (const auto& ident : peptide_identifications)
        {
          for (const auto& hit : ident.getHits())
          {
            OpenMS::Feature f;
            OpenMS::Feature s = feature;
            f.setUniqueId();
            s.setUniqueId();
            if (hit.getMetaValue("identifier").toStringList().at(0) != "null")
            {
              f.setMetaValue("PeptideRef", hit.getMetaValue("identifier").toStringList().at(0));
              s.setMetaValue("PeptideRef", hit.getMetaValue("identifier").toStringList().at(0));
              std::string native_id = hit.getMetaValue("chemical_formula").toString() + ";" + hit.getMetaValue("modifications").toString();
              s.setMetaValue("native_id", native_id);
              s.setMetaValue("identifier", hit.getMetaValue("identifier"));
              s.setMetaValue("description", hit.getMetaValue("description"));
              s.setMetaValue("modifications", hit.getMetaValue("modifications"));
              std::string adducts;
              try
              {
                // Extract adduct: the first letter stands for the actual metabolite and then everything after are the abducts up until the ";"
                // For example, M-H;1- will give -H
                std::string str = hit.getMetaValue("modifications").toString();
                std::string delimiter = ";";
                adducts = str.substr(1, str.find(delimiter) - 1);
              }
              catch (const std::exception& e)
              {
                OPENMS_LOG_ERROR << e.what();
              }
              s.setMetaValue("adducts", adducts);
              OpenMS::EmpiricalFormula chemform(hit.getMetaValue("chemical_formula").toString());
              double adduct_mass = s.getMZ() * std::abs(hit.getCharge()) + static_cast<double>(hit.getMetaValue("mz_error_Da")) - chemform.getMonoWeight();
              s.setMetaValue("dc_charge_adduct_mass", adduct_mass);
              s.setMetaValue("chemical_formula", hit.getMetaValue("chemical_formula"));
              s.setMetaValue("mz_error_ppm", hit.getMetaValue("mz_error_ppm"));
              s.setMetaValue("mz_error_Da", hit.getMetaValue("mz_error_Da"));
              s.setCharge(hit.getCharge());
              f.setSubordinates({s});
              feat_map_output.push_back(f);
            }
            else if (add_unidentified_features)
            {
              // "PeptideRef" metavalue should have been set during peak picking, but if not...
              std::ostringstream mass_of_the_peak;
              mass_of_the_peak << s.getMZ();

              // Fill in accurateMassSearch metavalues
              DataValue identifiers(std::vector<std::string>({mass_of_the_peak.str()}));
              s.setMetaValue("identifier", identifiers);
              s.setMetaValue("description", "");
              s.setMetaValue("modifications", "");
              s.setMetaValue("adducts", "");
              s.setMetaValue("dc_charge_adduct_mass", 0);
              s.setMetaValue("chemical_formula", "");
              s.setMetaValue("mz_error_ppm", 0);
              s.setMetaValue("mz_error_Da", 0);
              // s.setCharge(hit.getCharge()); // The polarity should have been set during peak picking
              f.setSubordinates({s});
              feat_map_output.push_back(f);
            }
          }
        }
      }
    }
  }

  void TargetedSpectraExtractor::annotateSpectra(
    const std::vector<MSSpectrum>& spectra,
    const TargetedExperiment& targeted_exp,
    std::vector<MSSpectrum>& annotated_spectra,
    FeatureMap& features,
    const bool compute_features
  ) const
  {
    annotated_spectra.clear();
    features.clear(true);
    const std::vector<ReactionMonitoringTransition>& transitions = targeted_exp.getTransitions();
    for (Size i = 0; i < spectra.size(); ++i)
    {
      const MSSpectrum& spectrum = spectra[i];
      const double spectrum_rt = spectrum.getRT();
      const double rt_left_lim = spectrum_rt - rt_window_ / 2.0;
      const double rt_right_lim = spectrum_rt + rt_window_ / 2.0;
      const std::vector<Precursor>& precursors = spectrum.getPrecursors();
      if (precursors.empty())
      {
        OPENMS_LOG_WARN << "annotateSpectra(): No precursor MZ found. Setting spectrum_mz to 0." << std::endl;
      }
      const double spectrum_mz = precursors.empty() ? 0.0 : precursors.front().getMZ();
      const double mz_tolerance = mz_unit_is_Da_ ? mz_tolerance_ : mz_tolerance_ / 1e6;

      // When spectrum_mz is 0, the mz check on transitions is inhibited
      const double mz_left_lim = spectrum_mz ? spectrum_mz - mz_tolerance : std::numeric_limits<double>::min();
      const double mz_right_lim = spectrum_mz ? spectrum_mz + mz_tolerance : std::numeric_limits<double>::max();

      OPENMS_LOG_DEBUG << "annotateSpectra(): [" << i << "] (RT: " << spectrum_rt << ") (MZ: " << spectrum_mz << ")" << std::endl;

      for (Size j = 0; j < transitions.size(); ++j)
      {
        const TargetedExperimentHelper::Peptide& peptide = targeted_exp.getPeptideByRef(transitions[j].getPeptideRef());
        double target_rt = peptide.getRetentionTime();
        if (peptide.getRetentionTimeUnit() == TargetedExperimentHelper::RetentionTime::RTUnit::MINUTE)
        {
          target_rt *= 60.0;
        }
        const double target_mz = transitions[j].getPrecursorMZ();
        if (target_rt >= rt_left_lim && target_rt <= rt_right_lim &&
            target_mz >= mz_left_lim && target_mz <= mz_right_lim)
        {
          OPENMS_LOG_DEBUG << "annotateSpectra(): [" << j << "][" << transitions[j].getPeptideRef() << "]";
          OPENMS_LOG_DEBUG << " (target_rt: " << target_rt << ") (target_mz: " << target_mz << ")" << std::endl << std::endl;
          MSSpectrum annotated_spectrum = spectrum;
          annotated_spectrum.setName(transitions[j].getPeptideRef());
          annotated_spectra.push_back(annotated_spectrum);
          if (compute_features)
          {
            Feature feature;
            feature.setRT(spectrum_rt);
            feature.setMZ(spectrum_mz);
            feature.setMetaValue("transition_name", transitions[j].getPeptideRef());
            features.push_back(feature);
          }
        }
      }
    }
    OPENMS_LOG_DEBUG << "annotateSpectra(): (input size: " << spectra.size() << ") (annotated spectra: " << annotated_spectra.size() << ")\n" << std::endl;
  }

  void TargetedSpectraExtractor::annotateSpectra(
    const std::vector<MSSpectrum>& spectra,
    const TargetedExperiment& targeted_exp,
    std::vector<MSSpectrum>& annotated_spectra
  ) const
  {
    FeatureMap features;
    const bool compute_features { false };
    annotateSpectra(spectra, targeted_exp, annotated_spectra, features, compute_features);
  }

  void TargetedSpectraExtractor::pickSpectrum(const MSSpectrum& spectrum, MSSpectrum& picked_spectrum) const
  {
    if (!spectrum.isSorted())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "Spectrum must be sorted by position");
    }

    // Smooth the spectrum
    MSSpectrum smoothed_spectrum = spectrum;
    if (use_gauss_)
    {
      GaussFilter gauss;
      Param filter_parameters = gauss.getParameters();
      filter_parameters.update(param_.copy("GaussFilter:", true));
      gauss.setParameters(filter_parameters);
      gauss.filter(smoothed_spectrum);
    }
    else
    {
      SavitzkyGolayFilter sgolay;
      Param filter_parameters = sgolay.getParameters();
      filter_parameters.update(param_.copy("SavitzkyGolayFilter:", true));
      sgolay.setParameters(filter_parameters);
      sgolay.filter(smoothed_spectrum);
    }

    // Find initial seeds (peak picking)
    Param pepi_param = PeakPickerHiRes().getDefaults();
    pepi_param.update(param_.copy("PeakPickerHiRes:", true));
    // disable spacing constraints, since we're dealing with spectrum
    pepi_param.setValue("spacing_difference", 0.0);
    pepi_param.setValue("spacing_difference_gap", 0.0);
    pepi_param.setValue("report_FWHM", "true");
    pepi_param.setValue("report_FWHM_unit", "absolute");
    picked_spectrum.clear(true);
    PeakPickerHiRes pp;
    pp.setParameters(pepi_param);
    pp.pick(smoothed_spectrum, picked_spectrum);

    std::vector<Int> peaks_pos_to_erase;
    const double fwhm_threshold = mz_unit_is_Da_ ? fwhm_threshold_ : fwhm_threshold_ / 1e6;
    for (Int i = picked_spectrum.size() - 1; i >= 0; --i)
    {
      if (picked_spectrum[i].getIntensity() < peak_height_min_ ||
          picked_spectrum[i].getIntensity() > peak_height_max_ ||
          picked_spectrum.getFloatDataArrays()[0][i] < fwhm_threshold)
      {
        peaks_pos_to_erase.push_back(i);
      }
    }

    if (peaks_pos_to_erase.size() != picked_spectrum.size()) // if not all peaks are to be removed
    {
      for (Int i : peaks_pos_to_erase) // then keep only the valid peaks (and fwhm)
      {
        picked_spectrum.erase(picked_spectrum.begin() + i);
        picked_spectrum.getFloatDataArrays()[0].erase(picked_spectrum.getFloatDataArrays()[0].begin() + i);
      }
    }
    else // otherwise output an empty picked_spectrum
    {
      picked_spectrum.clear(true);
    }

    OPENMS_LOG_DEBUG << "pickSpectrum(): " << spectrum.getName() << " (input size: " <<
      spectrum.size() << ") (picked: " << picked_spectrum.size() << ")\n" << std::endl;
  }

  void TargetedSpectraExtractor::scoreSpectra(
    const std::vector<MSSpectrum>& annotated_spectra,
    const std::vector<MSSpectrum>& picked_spectra,
    FeatureMap& features,
    std::vector<MSSpectrum>& scored_spectra,
    const bool compute_features
  ) const
  {
    scored_spectra.clear();
    scored_spectra.resize(annotated_spectra.size());
    if (compute_features && scored_spectra.size() != features.size())
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    for (Size i = 0; i < annotated_spectra.size(); ++i)
    {
      double total_tic { 0 };
      for (Size j = 0; j < annotated_spectra[i].size(); ++j)
      {
        total_tic += annotated_spectra[i][j].getIntensity();
      }

      double avgFWHM { 0 };
      if (!picked_spectra[i].getFloatDataArrays().empty())
      {
        for (Size j = 0; j < picked_spectra[i].getFloatDataArrays()[0].size(); ++j)
        {
          avgFWHM += picked_spectra[i].getFloatDataArrays()[0][j];
        }
        avgFWHM /= picked_spectra[i].getFloatDataArrays()[0].size();
      }
      SignalToNoiseEstimatorMedian<MSSpectrum> sne;
      Param p;
      p.setValue("win_len", 40.0);
      p.setValue("noise_for_empty_window", 2.0);
      p.setValue("min_required_elements", 10);
      sne.setParameters(p);
      sne.init(annotated_spectra[i]);
      double avgSNR { 0 };
      for (Size j = 0; j < annotated_spectra[i].size(); ++j)
      {
        avgSNR += sne.getSignalToNoise(j);
      }
      avgSNR /= annotated_spectra[i].size();

      const double log10_total_tic = log10(total_tic);
      const double inverse_avgFWHM = 1.0 / avgFWHM;
      const double score = log10_total_tic * tic_weight_ + inverse_avgFWHM * fwhm_weight_ + avgSNR * snr_weight_;

      scored_spectra[i] = annotated_spectra[i];
      scored_spectra[i].getFloatDataArrays().resize(5);
      scored_spectra[i].getFloatDataArrays()[1].setName("score");
      scored_spectra[i].getFloatDataArrays()[1].push_back(score);
      scored_spectra[i].getFloatDataArrays()[2].setName("log10_total_tic");
      scored_spectra[i].getFloatDataArrays()[2].push_back(log10_total_tic);
      scored_spectra[i].getFloatDataArrays()[3].setName("inverse_avgFWHM");
      scored_spectra[i].getFloatDataArrays()[3].push_back(inverse_avgFWHM);
      scored_spectra[i].getFloatDataArrays()[4].setName("avgSNR");
      scored_spectra[i].getFloatDataArrays()[4].push_back(avgSNR);

      if (compute_features)
      {
        // The intensity of a feature is (proportional to) its total ion count
        // http://www.openms.de/documentation/classOpenMS_1_1Feature.html
        features[i].setIntensity(score);
        features[i].setMetaValue("log10_total_tic", log10_total_tic);
        features[i].setMetaValue("inverse_avgFWHM", inverse_avgFWHM);
        features[i].setMetaValue("avgFWHM", avgFWHM);
        features[i].setMetaValue("avgSNR", avgSNR);
        std::vector<Feature> subordinates(picked_spectra[i].size());
        for (Size j = 0; j < picked_spectra[i].size(); ++j)
        {
          subordinates[j].setMZ(picked_spectra[i][j].getMZ());
          subordinates[j].setIntensity(picked_spectra[i][j].getIntensity());
          subordinates[j].setMetaValue("FWHM", picked_spectra[i].getFloatDataArrays()[0][j]);
        }
        features[i].setSubordinates(subordinates);
      }
    }
  }

  void TargetedSpectraExtractor::scoreSpectra(
    const std::vector<MSSpectrum>& annotated_spectra,
    const std::vector<MSSpectrum>& picked_spectra,
    std::vector<MSSpectrum>& scored_spectra
  ) const
  {
    FeatureMap features;
    const bool compute_features { false };
    scoreSpectra(annotated_spectra, picked_spectra, features, scored_spectra, compute_features);
  }

  void TargetedSpectraExtractor::selectSpectra(
    const std::vector<MSSpectrum>& scored_spectra,
    const FeatureMap& features,
    std::vector<MSSpectrum>& selected_spectra,
    FeatureMap& selected_features,
    const bool compute_features
  ) const
  {
    if (compute_features && scored_spectra.size() != features.size())
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    std::map<std::string,UInt> transition_best_spec;
    for (UInt i = 0; i < scored_spectra.size(); ++i)
    {
      if (scored_spectra[i].getFloatDataArrays()[1][0] < min_select_score_)
      {
        continue;
      }
      const std::string& transition_name = scored_spectra[i].getName();
      std::map<std::string,UInt>::const_iterator it = transition_best_spec.find(transition_name);
      if (it == transition_best_spec.cend())
      {
        transition_best_spec.emplace(transition_name, i);
      }
      else if (scored_spectra[it->second].getFloatDataArrays()[1][0] <
               scored_spectra[i].getFloatDataArrays()[1][0])
      {
        transition_best_spec.erase(transition_name);
        transition_best_spec.emplace(transition_name, i);
      }
    }

    selected_spectra.clear();
    selected_features.clear(true);

    for (const auto& m : transition_best_spec)
    {
      selected_spectra.push_back(scored_spectra[m.second]);
      if (compute_features) selected_features.push_back(features[m.second]);
    }
  }

  void TargetedSpectraExtractor::selectSpectra(
    const std::vector<MSSpectrum>& scored_spectra,
    std::vector<MSSpectrum>& selected_spectra
  ) const
  {
    FeatureMap dummy_features;
    FeatureMap dummy_selected_features;
    const bool compute_features { false };
    selectSpectra(scored_spectra, dummy_features, selected_spectra, dummy_selected_features, compute_features);
  }

  void TargetedSpectraExtractor::extractSpectra(
    const MSExperiment& experiment,
    const TargetedExperiment& targeted_exp,
    std::vector<MSSpectrum>& extracted_spectra,
    FeatureMap& extracted_features,
    const bool compute_features
  ) const
  {
    // get the spectra from the experiment
    const std::vector<MSSpectrum>& spectra = experiment.getSpectra();

    // annotate spectra
    std::vector<MSSpectrum> annotated;
    FeatureMap features;
    annotateSpectra(spectra, targeted_exp, annotated, features, compute_features);

    // pick peaks from annotate spectra
    std::vector<MSSpectrum> picked(annotated.size());
    for (Size i = 0; i < annotated.size(); ++i)
    {
      pickSpectrum(annotated[i], picked[i]);
    }

    // remove empty picked<> spectra, and accordingly update annotated<> and features
    for (Int i = annotated.size() - 1; i >= 0; --i)
    {
      if (picked[i].empty())
      {
        annotated.erase(annotated.begin() + i);
        picked.erase(picked.begin() + i);
        if (compute_features) features.erase(features.begin() + i);
      }
    }

    // score spectra
    std::vector<MSSpectrum> scored;
    scoreSpectra(annotated, picked, features, scored, compute_features);

    // select the best spectrum for each group of spectra having the same name
    selectSpectra(scored, features, extracted_spectra, extracted_features, compute_features);
  }

  void TargetedSpectraExtractor::extractSpectra(
    const MSExperiment& experiment,
    const TargetedExperiment& targeted_exp,
    std::vector<MSSpectrum>& extracted_spectra
  ) const
  {
    FeatureMap extracted_features;
    const bool compute_features { false };
    extractSpectra(experiment, targeted_exp, extracted_spectra, extracted_features, compute_features);
  }

  void TargetedSpectraExtractor::extractSpectra(
    const MSExperiment& experiment,
    const FeatureMap& ms1_features,
    std::vector<MSSpectrum>& extracted_spectra
  ) const
  {
    FeatureMap extracted_features;
    extractSpectra(experiment, ms1_features, extracted_spectra, extracted_features, false);
  }

  void TargetedSpectraExtractor::extractSpectra(
    const MSExperiment& experiment,
    const FeatureMap& ms1_features,
    std::vector<MSSpectrum>& extracted_spectra,
    FeatureMap& extracted_features
  ) const
  {
    extractSpectra(experiment, ms1_features, extracted_spectra, extracted_features, true);
  }

  void TargetedSpectraExtractor::extractSpectra(
    const MSExperiment& experiment,
    const FeatureMap& ms1_features,
    std::vector<MSSpectrum>& extracted_spectra,
    FeatureMap& extracted_features,
    const bool compute_features
  ) const
  {
    // annotate spectra
    std::vector<OpenMS::MSSpectrum> annotated_spectra;
    OpenMS::FeatureMap ms2_features;
    annotateSpectra(experiment.getSpectra(), ms1_features, ms2_features, annotated_spectra);

    // pickSpectra
    std::vector<MSSpectrum> picked_spectra(annotated_spectra.size());
    for (Size i = 0; i < annotated_spectra.size(); ++i)
    {
      pickSpectrum(annotated_spectra[i], picked_spectra[i]);
    }

    // score and select
    std::vector<OpenMS::MSSpectrum> scored_spectra;
    scoreSpectra(annotated_spectra, picked_spectra, scored_spectra);

    // select the best spectrum for each group of spectra having the same name
    // NOTE: It maybe needed to take the top N instead of the top 1 spectra in the future
    selectSpectra(scored_spectra, ms2_features, extracted_spectra, extracted_features, compute_features);
  }

  void TargetedSpectraExtractor::matchSpectrum(
    const MSSpectrum& input_spectrum,
    const Comparator& cmp,
    std::vector<Match>& matches
  ) const
  {
    // TODO: remove times debug info
    // std::clock_t start;
    // start = std::clock();
    matches.clear();
    std::vector<std::pair<Size,double>> scores;

    cmp.generateScores(input_spectrum, scores, min_match_score_);

    // Sort the vector of scores
    std::sort(scores.begin(), scores.end(),
      [](const std::pair<Size,double>& a, const std::pair<Size,double>& b)
      {
        return a.second > b.second;
      });

    // Set the number of best matches to return
    const Size n = std::min(top_matches_to_report_, scores.size());

    // Construct a vector of n `Match`es
    for (Size i = 0; i < n; ++i)
    {
      const Size spec_idx { scores[i].first };
      const double spec_score { scores[i].second };
      matches.emplace_back(cmp.getLibrary()[spec_idx], spec_score);
    }
  }

  void TargetedSpectraExtractor::targetedMatching(
    const std::vector<MSSpectrum>& spectra,
    const Comparator& cmp,
    FeatureMap& features
  )
  {
    if (spectra.size() != features.size())
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    std::vector<Size> no_matches_idx; // to keep track of those features without a match
    const Size tmp = top_matches_to_report_;
    top_matches_to_report_ = 1;

    for (Size i = 0; i < spectra.size(); ++i)
    {
      std::vector<Match> matches;
      matchSpectrum(spectra[i], cmp, matches);
      if (!matches.empty())
      {
        features[i].setMetaValue("spectral_library_name", matches[0].spectrum.getName());
        features[i].setMetaValue("spectral_library_score", matches[0].score);
        const String& comments = matches[0].spectrum.metaValueExists("Comments") ?
          matches[0].spectrum.getMetaValue("Comments") : "";
        features[i].setMetaValue("spectral_library_comments", comments);
      }
      else
      {
        no_matches_idx.push_back(i);
        features[i].setMetaValue("spectral_library_name", "");
        features[i].setMetaValue("spectral_library_score", 0.0);
        features[i].setMetaValue("spectral_library_comments", "");
      }
    }

    top_matches_to_report_ = tmp;

    if (!no_matches_idx.empty())
    {
      String warn_msg = "No match was found for " + std::to_string(no_matches_idx.size()) + " `Feature`s. Indices: ";
      for (const Size idx : no_matches_idx)
      {
        warn_msg += std::to_string(idx) + " ";
      }
      OPENMS_LOG_WARN << std:: endl << warn_msg << std::endl;
    }
  }

  void TargetedSpectraExtractor::untargetedMatching(
    const std::vector<MSSpectrum>& spectra,
    const Comparator& cmp,
    FeatureMap& features
  )
  {
    features.clear(true);

    std::vector<MSSpectrum> picked(spectra.size());
    for (Size i = 0; i < spectra.size(); ++i)
    {
      pickSpectrum(spectra[i], picked[i]);
    }

    // remove empty picked<> spectra
    for (Int i = spectra.size() - 1; i >= 0; --i)
    {
      if (picked[i].empty())
      {
        picked.erase(picked.begin() + i);
      }
    }

    for (const MSSpectrum& spectrum : picked)
    {
      const std::vector<Precursor>& precursors = spectrum.getPrecursors();
      if (precursors.empty())
      {
        OPENMS_LOG_WARN << "untargetedMatching(): No precursor MZ found. Setting spectrum_mz to 0." << std::endl;
      }
      const double spectrum_mz = precursors.empty() ? 0.0 : precursors.front().getMZ();
      Feature feature;
      feature.setRT(spectrum.getRT());
      feature.setMZ(spectrum_mz);
      features.push_back(feature);
    }

    targetedMatching(picked, cmp, features);
  }

  void TargetedSpectraExtractor::constructTransitionsList(const OpenMS::FeatureMap& ms1_features, const OpenMS::FeatureMap& ms2_features, TargetedExperiment& t_exp) const
  { 
    // Create a map based on PeptideRef between MS1 and MS2 features
    std::map<std::string, std::vector<const Feature*>> ms1_to_ms2;
    for (const auto& feature : ms2_features)
    {
      for (const auto& subordinate : feature.getSubordinates())
      {
        ms1_to_ms2[subordinate.getMetaValue("PeptideRef")].push_back(&subordinate);
      }
    }

    // Create individual maps for the transitions, peptides and proteins to account for the same IDs but different RTs
    std::map<std::string, ReactionMonitoringTransition> rmt_map;
    std::map<std::string, OpenMS::TargetedExperiment::Peptide> peptides_map;
    std::map<std::string, OpenMS::TargetedExperiment::Protein> proteins_map;
    for (const auto& ms1_feature : ms1_features)
    {
      std::string peptide_ref = ms1_feature.getMetaValue("PeptideRef");
      OpenMS::TargetedExperiment::Protein protein;
      protein.id = peptide_ref;
      protein.addMetaValues(ms1_feature);
      proteins_map.emplace(peptide_ref, protein); // OK to reject duplicate keys

      OpenMS::ReactionMonitoringTransition::RetentionTime rt_f;
      rt_f.setRT(ms1_feature.getRT());

      OpenMS::TargetedExperiment::Peptide peptide;
      peptide.id = peptide_ref;
      peptide.setChargeState(ms1_feature.getCharge());
      peptide.addMetaValues(ms1_feature);
      peptide.protein_refs.emplace_back(peptide_ref);
      peptide.rts.push_back(rt_f);
      auto found_peptide = peptides_map.emplace(peptide_ref, peptide);
      if (!found_peptide.second)
      {
        peptides_map.at(peptide_ref).rts.push_back(rt_f);
      }
      
      for (const auto& ms2_feature : ms1_to_ms2[peptide_ref])
      {
        auto current_mz = ms2_feature->getMZ();
        if ((current_mz > min_fragment_mz_ && current_mz < max_fragment_mz_) &&
            (current_mz < ms1_feature.getMZ() + relative_allowable_product_mass_))
        {
          OpenMS::ReactionMonitoringTransition rmt;
          rmt.setLibraryIntensity(ms2_feature->getIntensity());
          rmt.setName(ms2_feature->getMetaValue("native_id"));

          OpenMS::ReactionMonitoringTransition::RetentionTime rt_s;
          rt_s.setRT(ms2_feature->getRT());
          rmt.setRetentionTime(rt_s);

          std::ostringstream os;
          os << peptide_ref << "_" << ms2_feature->getMetaValue("native_id") << "_" << ms2_feature->getRT();
          rmt.setNativeID(os.str());
          rmt.setPeptideRef(peptide_ref);
          rmt.setPrecursorMZ(ms1_feature.getMZ());
          rmt.setProductMZ(ms2_feature->getMZ());
          rmt.addMetaValues(*ms2_feature);
          rmt_map.emplace(os.str(), rmt); // OK to reject duplicate keys
        }
      }
    }

    // Reconstruct the final vectors from the maps
    std::vector<TargetedExperiment::Peptide> peptides;
    for (const auto& p : peptides_map) {
      peptides.push_back(p.second);
    }
    std::vector<TargetedExperiment::Protein> proteins;
    for (const auto& p : proteins_map)
    {
      proteins.push_back(p.second);
    }
    std::vector<ReactionMonitoringTransition> rmt_vec;
    for (const auto& p : rmt_map)
    {
      rmt_vec.push_back(p.second);
    }
    t_exp.setProteins(proteins);
    t_exp.setPeptides(peptides);
    t_exp.setTransitions(rmt_vec);

    // validate
    OpenMS::TransitionTSVFile tsv_file;
    tsv_file.validateTargetedExperiment(t_exp);
  }

  void TargetedSpectraExtractor::mergeFeatures(const OpenMS::FeatureMap& fmap_input, OpenMS::FeatureMap& fmap_output) const
  {
    try
    {
      // Pass 1: organize into a map by combining features and subordinates with the same `identifier`
      std::map<OpenMS::String, std::vector<OpenMS::Feature>> fmapmap;
      organizeMapWithSameIdentifier(fmap_input, fmapmap);

      // Pass 2: compute the consensus manually
      for (const auto& f_map : fmapmap)
      {
        // compute the total intensity for weighting
        double total_intensity = 0;
        for (const auto& f : f_map.second)
        {
          if (f.metaValueExists("peak_apex_int"))
            total_intensity += (double) f.getMetaValue("peak_apex_int");
          else
            total_intensity += f.getIntensity();
        }

        // compute the weighted averages
        double rt = 0.0, m = 0.0, intensity = 0.0, peak_apex_int = 0.0;
        double weighting_factor = 1.0 / f_map.second.size();// will be updated
        for (const auto& f : f_map.second)
        {
          // compute the weighting factor
          if (f.metaValueExists("peak_apex_int"))
            weighting_factor = (double) f.getMetaValue("peak_apex_int") / total_intensity;
          else
            weighting_factor = f.getIntensity() / total_intensity;

          // compute the weighted averages
          rt += f.getRT() * weighting_factor;
          m += f.getMZ() * weighting_factor;
          intensity += f.getIntensity();
          if (f.metaValueExists("peak_apex_int"))
            peak_apex_int += (double) f.getMetaValue("peak_apex_int");
        }

        // make the feature map and assign subordinates
        OpenMS::Feature f;
        f.setUniqueId();

        // parse the identifier
        std::string id_f;
        try
        {
          id_f = f_map.first.prefix('_');
        }
        catch (const std::exception& e)
        {
          OPENMS_LOG_ERROR << e.what();
        }
        f.setMetaValue("PeptideRef", id_f);
        f.setMZ(m);
        f.setRT(rt);
        f.setMetaValue("scan_polarity", f_map.second.front().getMetaValue("scan_polarity"));
        f.setIntensity(intensity);
        f.setMetaValue("peak_apex_int", peak_apex_int);
        f.setSubordinates(f_map.second);
        fmap_output.push_back(f);
      }
    }
    catch (const std::exception& e)
    {
      OPENMS_LOG_ERROR << e.what();
    }
  }

  void TargetedSpectraExtractor::storeSpectraMSP(const String& filename, MSExperiment& experiment) const
  {
    if (deisotoping_use_deisotoper_)
    {
      deisotopeMS2Spectra_(experiment);
    }

    removeMS2SpectraPeaks_(experiment);

    // Store
    FileHandler().storeExperiment(filename, experiment, {FileTypes::MSP});
  }
  
  void TargetedSpectraExtractor::deisotopeMS2Spectra_(MSExperiment& experiment) const
  {
    for (auto& peakmap_it : experiment.getSpectra())
    {
      MSSpectrum& spectrum = peakmap_it;
      if (spectrum.getMSLevel() == 1)
      {
        continue;
      }
      bool fragment_unit_ppm = deisotoping_fragment_unit_ == "ppm" ? true : false;
      bool make_single_charged = false;
      Deisotoper::deisotopeAndSingleCharge(spectrum,
                                           deisotoping_fragment_tolerance_,
                                           fragment_unit_ppm,
                                           deisotoping_min_charge_,
                                           deisotoping_max_charge_,
                                           deisotoping_keep_only_deisotoped_,
                                           deisotoping_min_isopeaks_,
                                           deisotoping_max_isopeaks_,
                                           make_single_charged,
                                           deisotoping_annotate_charge_);
    }  
  }

  void TargetedSpectraExtractor::removeMS2SpectraPeaks_(MSExperiment& experiment) const
  {
    // remove peaks form MS2 which are at a higher mz than the precursor + 10 ppm
    for (auto& peakmap_it : experiment.getSpectra())
    {
      MSSpectrum& spectrum = peakmap_it;
      if (spectrum.getMSLevel() == 1)
      {
        continue;
      }
      // if peak mz higher than precursor mz set intensity to zero
      double prec_mz = spectrum.getPrecursors()[0].getMZ();
      double mass_diff = Math::ppmToMass(max_precursor_mass_threashold_, prec_mz);
      for (auto& spec : spectrum)
      {
        if (spec.getMZ() > prec_mz + mass_diff)
        {
          spec.setIntensity(0);
        }
      }
      spectrum.erase(remove_if(spectrum.begin(),
                               spectrum.end(),
                               InIntensityRange<PeakMap::PeakType>(1,
                                                                   std::numeric_limits<PeakMap::PeakType::IntensityType>::max(),
                                                                   true)),
                     spectrum.end());
    }
  }

  void TargetedSpectraExtractor::organizeMapWithSameIdentifier(const OpenMS::FeatureMap& fmap_input, std::map<OpenMS::String, std::vector<OpenMS::Feature>>& fmapmap) const
  {
    auto construct_feature = [&fmapmap](const OpenMS::Feature& feature)
    {
      if (feature.metaValueExists("PeptideRef") && feature.metaValueExists("identifier"))
      {
        std::string id = std::string(feature.getMetaValue("PeptideRef")) + std::string("_") + std::string(feature.getMetaValue("identifier").toStringList().at(0));
        std::string id_f = id + std::string("_") + std::to_string(feature.getRT());
        auto found_f = fmapmap.emplace(id_f, std::vector<OpenMS::Feature>({feature}));
        if (!found_f.second)
        {
          fmapmap.at(id_f).push_back(feature);
        }
      }
    };
    for (const OpenMS::Feature& f : fmap_input)
    {
      construct_feature(f);
      for (const OpenMS::Feature& s : f.getSubordinates())
      {
        construct_feature(s);
      }
    }
  }

  void TargetedSpectraExtractor::BinnedSpectrumComparator::init(const std::vector<MSSpectrum>& library, const std::map<String,DataValue>& options)
  {
    if (options.count("bin_size"))
    {
      bin_size_ = options.at("bin_size");
    }
    if (options.count("peak_spread"))
    {
      peak_spread_ = options.at("peak_spread");
    }
    if (options.count("bin_offset"))
    {
      bin_offset_ = options.at("bin_offset");
    }
    library_.clear();
    bs_library_.clear();
    for (const MSSpectrum& s : library)
    {
      library_.push_back(s);
      bs_library_.emplace_back(s, bin_size_, false, peak_spread_, bin_offset_);
    }
    OPENMS_LOG_INFO << "The library contains " << bs_library_.size() << " spectra." << std::endl;
  }

}// namespace OpenMS
