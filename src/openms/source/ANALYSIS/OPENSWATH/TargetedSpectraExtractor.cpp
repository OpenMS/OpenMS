// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/TargetedSpectraExtractor.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>
#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>
#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>

namespace OpenMS
{
  TargetedSpectraExtractor::TargetedSpectraExtractor() :
    DefaultParamHandler("TargetedSpectraExtractor")
  {
    getDefaultParameters(defaults_);

    subsections_.push_back("SavitzkyGolayFilter");
    defaults_.setValue("SavitzkyGolayFilter:frame_length", 15);
    defaults_.setValue("SavitzkyGolayFilter:polynomial_order", 3);

    subsections_.push_back("GaussFilter");
    defaults_.setValue("GaussFilter:gaussian_width", 0.2);

    subsections_.push_back("PeakPickerHiRes");
    defaults_.setValue("PeakPickerHiRes:signal_to_noise", 1.0);

    defaultsToParam_(); // write defaults into Param object param_
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
    params.setValidStrings("mz_unit_is_Da", ListUtils::create<String>("false,true"));

    params.setValue("use_gauss", "true", "Use Gaussian filter for smoothing (alternative is Savitzky-Golay filter)");
    params.setValidStrings("use_gauss", ListUtils::create<String>("false,true"));

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
      for (Size j = 0; j < picked_spectra[i].getFloatDataArrays()[0].size(); ++j)
      {
        avgFWHM += picked_spectra[i].getFloatDataArrays()[0][j];
      }
      avgFWHM /= picked_spectra[i].getFloatDataArrays()[0].size();

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

    for (auto it = transition_best_spec.cbegin(); it != transition_best_spec.cend(); ++it)
    {
      selected_spectra.push_back(scored_spectra[it->second]);
      if (compute_features) selected_features.push_back(features[it->second]);
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

  void TargetedSpectraExtractor::matchSpectrum(
    const MSSpectrum& input_spectrum,
    const Comparator& cmp,
    std::vector<Match>& matches
  )
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

    // std::cout << "MATCH TIME: " << ((std::clock() - start) / (double)CLOCKS_PER_SEC) << std::endl;
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
      if (matches.size())
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

    if (no_matches_idx.size())
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
}
