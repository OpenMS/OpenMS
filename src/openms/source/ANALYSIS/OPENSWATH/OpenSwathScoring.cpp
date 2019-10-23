// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathScoring.h>

#include <OpenMS/KERNEL/ComparatorUtils.h>
#include <OpenMS/CONCEPT/Macros.h>

// scoring
#include <OpenMS/OPENSWATHALGO/ALGO/Scoring.h>
#include <OpenMS/OPENSWATHALGO/ALGO/MRMScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SONARScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/IonMobilityScoring.h>

// auxiliary
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SpectrumAddition.h>

// basic file operations

namespace OpenMS
{

  void sortSpectrumByMZ(OpenSwath::Spectrum& spec)
  {
    //sort index list
    std::vector<std::pair<double, Size> > sorted_indices;
    sorted_indices.reserve(spec.getMZArray()->data.size());
    auto mz_it = spec.getMZArray()->data.begin();
    for (Size i = 0; i < spec.getMZArray()->data.size(); ++i)
    {
      sorted_indices.emplace_back(*mz_it, i);
      ++mz_it;
    }
    std::stable_sort(sorted_indices.begin(), sorted_indices.end(), PairComparatorFirstElement<std::pair<double, Size> >());

    // extract list of indices
    std::vector<Size> select_indices;
    select_indices.reserve(sorted_indices.size());
    for (const auto& sidx : sorted_indices)
    {
      select_indices.push_back(sidx.second);
    }

    for (auto& da : spec.getDataArrays() )
    {
      if (da->data.empty()) continue;
      OpenSwath::BinaryDataArrayPtr tmp(new OpenSwath::BinaryDataArray);
      tmp->description = da->description;
      tmp->data.reserve(select_indices.size());
      for (Size i = 0; i < select_indices.size(); ++i)
      {
        tmp->data.push_back( da->data[ select_indices[i] ] );
      }
      da = tmp;
    }

    OPENMS_POSTCONDITION( std::adjacent_find(spec.getMZArray()->data.begin(),
           spec.getMZArray()->data.end(), std::greater<double>()) == spec.getMZArray()->data.end(),
           "Postcondition violated: m/z vector needs to be sorted!" )
  }
}

namespace OpenMS
{

  /// Constructor
  OpenSwathScoring::OpenSwathScoring() :
    rt_normalization_factor_(1.0),
    spacing_for_spectra_resampling_(0.005),
    add_up_spectra_(1),
    spectra_addition_method_("simple"),
    im_drift_extra_pcnt_(0.0)
  {
  }

  /// Destructor
  OpenSwathScoring::~OpenSwathScoring()
  {
  }

  void OpenSwathScoring::initialize(double rt_normalization_factor,
                                    int add_up_spectra,
                                    double spacing_for_spectra_resampling,
                                    const double drift_extra,
                                    const OpenSwath_Scores_Usage & su,
                                    const std::string& spectrum_addition_method)
  {
    this->rt_normalization_factor_ = rt_normalization_factor;
    this->add_up_spectra_ = add_up_spectra;
    this->spectra_addition_method_ = spectrum_addition_method;
    this->im_drift_extra_pcnt_ = drift_extra;
    this->spacing_for_spectra_resampling_ = spacing_for_spectra_resampling;
    this->su_ = su;
  }

  void OpenSwathScoring::calculateDIAScores(OpenSwath::IMRMFeature* imrmfeature,
                                            const std::vector<TransitionType>& transitions,
                                            const std::vector<OpenSwath::SwathMap>& swath_maps,
                                            OpenSwath::SpectrumAccessPtr ms1_map,
                                            OpenMS::DIAScoring& diascoring,
                                            const CompoundType& compound,
                                            OpenSwath_Scores& scores,
                                            std::vector<double>& masserror_ppm,
                                            const double drift_lower,
                                            const double drift_upper,
                                            const double drift_target)
  {
    OPENMS_PRECONDITION(imrmfeature != nullptr, "Feature to be scored cannot be null");
    OPENMS_PRECONDITION(transitions.size() > 0, "There needs to be at least one transition.");
    OPENMS_PRECONDITION(swath_maps.size() > 0, "There needs to be at least one swath map.");

    // Identify corresponding SONAR maps (if more than one map is used)
    std::vector<OpenSwath::SwathMap> used_swath_maps;
    if (swath_maps.size() > 1 || transitions.empty())
    {
      double precursor_mz = transitions[0].getPrecursorMZ();
      for (size_t i = 0; i < swath_maps.size(); ++i)
      {
        if (swath_maps[i].ms1) {continue;} // skip MS1
        if (precursor_mz > swath_maps[i].lower && precursor_mz < swath_maps[i].upper)
        {
          used_swath_maps.push_back(swath_maps[i]);
        }
      }
    }
    else
    {
      used_swath_maps = swath_maps;
    }

    std::vector<double> normalized_library_intensity;
    getNormalized_library_intensities_(transitions, normalized_library_intensity);

    // find spectrum that is closest to the apex of the peak using binary search
    OpenSwath::SpectrumPtr spectrum = fetchSpectrumSwath(used_swath_maps, imrmfeature->getRT(), add_up_spectra_, drift_lower, drift_upper);

    // calculate drift extraction width for current spectrum (with some extra for cross-correlation)
    double drift_width = fabs(drift_upper - drift_lower);
    double drift_lower_used = drift_lower - drift_width * im_drift_extra_pcnt_;
    double drift_upper_used = drift_upper + drift_width * im_drift_extra_pcnt_;

    // score drift time dimension
    if (drift_upper > 0 && su_.use_im_scores)
    {
      double dia_extract_window_ = (double)diascoring.getParameters().getValue("dia_extraction_window");
      bool dia_extraction_ppm_ = diascoring.getParameters().getValue("dia_extraction_unit") == "ppm";
      auto drift_spectrum = fetchSpectrumSwath(used_swath_maps, imrmfeature->getRT(), add_up_spectra_, drift_lower_used, drift_upper_used);
      IonMobilityScoring::driftScoring(drift_spectrum, transitions, scores,
                                       drift_lower, drift_upper, drift_target,
                                       dia_extract_window_, dia_extraction_ppm_,
                                       false, im_drift_extra_pcnt_);
    }

    // Mass deviation score
    diascoring.dia_massdiff_score(transitions, spectrum, normalized_library_intensity, scores.massdev_score, scores.weighted_massdev_score, masserror_ppm);

    // DIA dotproduct and manhattan score based on library intensity
    diascoring.score_with_isotopes(spectrum, transitions, scores.dotprod_score_dia, scores.manhatt_score_dia);

    // Isotope correlation / overlap score: Is this peak part of an
    // isotopic pattern or is it the monoisotopic peak in an isotopic
    // pattern?
    // Currently this is computed for an averagine model of a peptide so its
    // not optimal for metabolites - but better than nothing, given that for
    // most fragments we dont really know their composition
    diascoring.dia_isotope_scores(transitions, spectrum, imrmfeature, scores.isotope_correlation, scores.isotope_overlap);

    // Peptide-specific scores
    if (compound.isPeptide())
    {
      // Presence of b/y series score
      OpenMS::AASequence aas;
      int by_charge_state = 1; // for which charge states should we check b/y series
      OpenSwathDataAccessHelper::convertPeptideToAASequence(compound, aas);
      diascoring.dia_by_ion_score(spectrum, aas, by_charge_state, scores.bseries_score, scores.yseries_score);
    }

    if (ms1_map && ms1_map->getNrSpectra() > 0) 
    {
      double precursor_mz = transitions[0].precursor_mz;
      double rt = imrmfeature->getRT();

      calculatePrecursorDIAScores(ms1_map, diascoring, precursor_mz, rt, compound, scores, drift_lower, drift_upper);

      if (drift_upper > 0 && su_.use_im_scores)
      {
        double dia_extract_window_ = (double)diascoring.getParameters().getValue("dia_extraction_window");
        bool dia_extraction_ppm_ = diascoring.getParameters().getValue("dia_extraction_unit") == "ppm";
        IonMobilityScoring::driftScoringMS1( fetchSpectrumSwath(ms1_map, imrmfeature->getRT(), add_up_spectra_, drift_lower_used, drift_upper_used),
            transitions, scores, drift_lower, drift_upper, drift_target, dia_extract_window_, dia_extraction_ppm_, false, im_drift_extra_pcnt_);

        IonMobilityScoring::driftScoringMS1Contrast(
            fetchSpectrumSwath(used_swath_maps, imrmfeature->getRT(), add_up_spectra_, drift_lower_used, drift_upper_used),
            fetchSpectrumSwath(ms1_map, imrmfeature->getRT(), add_up_spectra_, drift_lower, drift_upper),
            transitions, scores, drift_lower, drift_upper, dia_extract_window_, dia_extraction_ppm_, im_drift_extra_pcnt_);
      }
    }

  }

  void OpenSwathScoring::calculatePrecursorDIAScores(OpenSwath::SpectrumAccessPtr ms1_map, 
                                   OpenMS::DIAScoring & diascoring, 
                                   double precursor_mz, 
                                   double rt, 
                                   const CompoundType& compound, 
                                   OpenSwath_Scores & scores,
                                   double drift_lower, double drift_upper)
  {
    // Compute precursor-level scores:
    // - compute mass difference in ppm
    // - compute isotopic pattern score
    if (ms1_map && ms1_map->getNrSpectra() > 0)
    {
      OpenSwath::SpectrumPtr ms1_spectrum = fetchSpectrumSwath(ms1_map, rt, add_up_spectra_, drift_lower, drift_upper);
      diascoring.dia_ms1_massdiff_score(precursor_mz, ms1_spectrum, scores.ms1_ppm_score);

      // derive precursor charge state (get from data if possible)
      int precursor_charge = 1;
      if (compound.getChargeState() != 0) 
      {
        precursor_charge = compound.getChargeState();
      }

      if (compound.isPeptide())
      {
        diascoring.dia_ms1_isotope_scores(precursor_mz, ms1_spectrum,
                                          precursor_charge, scores.ms1_isotope_correlation,
                                          scores.ms1_isotope_overlap);
      }
      else
      {
        diascoring.dia_ms1_isotope_scores(precursor_mz, ms1_spectrum,
                                          precursor_charge, scores.ms1_isotope_correlation,
                                          scores.ms1_isotope_overlap, compound.sum_formula);
      }
    }
  }

  void OpenSwathScoring::calculateDIAIdScores(OpenSwath::IMRMFeature* imrmfeature,
                                              const TransitionType & transition,
                                              const std::vector<OpenSwath::SwathMap> swath_maps,
                                              OpenMS::DIAScoring & diascoring,
                                              OpenSwath_Scores & scores,
                                              double drift_lower, double drift_upper)
  {
    OPENMS_PRECONDITION(imrmfeature != nullptr, "Feature to be scored cannot be null");
    OPENMS_PRECONDITION(swath_maps.size() > 0, "There needs to be at least one swath map.");

    // Identify corresponding SONAR maps (if more than one map is used)
    std::vector<OpenSwath::SwathMap> used_swath_maps;
    if (swath_maps.size() > 1)
    {
      double precursor_mz = transition.getPrecursorMZ();
      for (size_t i = 0; i < swath_maps.size(); ++i)
      {
        if (swath_maps[i].ms1) {continue;} // skip MS1
        if (precursor_mz > swath_maps[i].lower && precursor_mz < swath_maps[i].upper)
        {
          used_swath_maps.push_back(swath_maps[i]);
        }
      }
    }
    else
    {
      used_swath_maps = swath_maps;
    }

    // find spectrum that is closest to the apex of the peak using binary search
    OpenSwath::SpectrumPtr spectrum = fetchSpectrumSwath(used_swath_maps, imrmfeature->getRT(), add_up_spectra_, drift_lower, drift_upper);

    // If no charge is given, we assume it to be 1
    int putative_product_charge = 1;
    if (transition.getProductChargeState() > 0)
    {
      putative_product_charge = transition.getProductChargeState();
    }

    // Isotope correlation / overlap score: Is this peak part of an
    // isotopic pattern or is it the monoisotopic peak in an isotopic
    // pattern?
    diascoring.dia_ms1_isotope_scores(transition.getProductMZ(), spectrum, putative_product_charge, scores.isotope_correlation, scores.isotope_overlap);
    // Mass deviation score
    diascoring.dia_ms1_massdiff_score(transition.getProductMZ(), spectrum, scores.massdev_score);
  }

  void OpenSwathScoring::calculateChromatographicScores(
        OpenSwath::IMRMFeature* imrmfeature,
        const std::vector<std::string>& native_ids,
        const std::vector<std::string>& precursor_ids,
        const std::vector<double>& normalized_library_intensity,
        std::vector<OpenSwath::ISignalToNoisePtr>& signal_noise_estimators,
        OpenSwath_Scores & scores)
  {
    OPENMS_PRECONDITION(imrmfeature != nullptr, "Feature to be scored cannot be null");
    OpenSwath::MRMScoring mrmscore_;
    if (su_.use_coelution_score_ || su_.use_shape_score_ || (imrmfeature->getPrecursorIDs().size() > 0 && su_.use_ms1_correlation))
      mrmscore_.initializeXCorrMatrix(imrmfeature, native_ids);

    // XCorr score (coelution)
    if (su_.use_coelution_score_)
    {
      scores.xcorr_coelution_score = mrmscore_.calcXcorrCoelutionScore();
      scores.weighted_coelution_score = mrmscore_.calcXcorrCoelutionWeightedScore(normalized_library_intensity);
    }

    // XCorr score (shape)
    // mean over the intensities at the max of the crosscorrelation
    // FEATURE : weigh by the intensity as done by mQuest
    // FEATURE : normalize with the intensity at the peak group apex?
    if (su_.use_shape_score_)
    {
      scores.xcorr_shape_score = mrmscore_.calcXcorrShapeScore();
      scores.weighted_xcorr_shape = mrmscore_.calcXcorrShapeWeightedScore(normalized_library_intensity);
    }

    // check that the MS1 feature is present and that the MS1 correlation should be calculated
    if (imrmfeature->getPrecursorIDs().size() > 0 && su_.use_ms1_correlation)
    {
      // we need at least two precursor isotopes
      if (precursor_ids.size() > 1)
      {
        mrmscore_.initializeXCorrPrecursorMatrix(imrmfeature, precursor_ids);
        scores.ms1_xcorr_coelution_score = mrmscore_.calcXcorrPrecursorCoelutionScore();
        scores.ms1_xcorr_shape_score = mrmscore_.calcXcorrPrecursorShapeScore();
      }
      mrmscore_.initializeXCorrPrecursorContrastMatrix(imrmfeature, precursor_ids, native_ids); // perform cross-correlation on monoisotopic precursor
      scores.ms1_xcorr_coelution_contrast_score = mrmscore_.calcXcorrPrecursorContrastCoelutionScore();
      scores.ms1_xcorr_shape_contrast_score = mrmscore_.calcXcorrPrecursorContrastShapeScore();

      mrmscore_.initializeXCorrPrecursorCombinedMatrix(imrmfeature, precursor_ids, native_ids); // perform cross-correlation on monoisotopic precursor
      scores.ms1_xcorr_coelution_combined_score = mrmscore_.calcXcorrPrecursorCombinedCoelutionScore();
      scores.ms1_xcorr_shape_combined_score = mrmscore_.calcXcorrPrecursorCombinedShapeScore();
    }

    if (su_.use_nr_peaks_score_)
    {
      scores.nr_peaks = boost::numeric_cast<int>(imrmfeature->size());
    }

    // Signal to noise scoring
    if (su_.use_sn_score_)
    {
      scores.sn_ratio = mrmscore_.calcSNScore(imrmfeature, signal_noise_estimators);
      // everything below S/N 1 can be set to zero (and the log safely applied)
      if (scores.sn_ratio < 1)
      { 
        scores.log_sn_score = 0;
      }
      else
      {
        scores.log_sn_score = std::log(scores.sn_ratio);
      }
    }

    // Mutual information scoring
    if (su_.use_mi_score_)
    {
      mrmscore_.initializeMIMatrix(imrmfeature, native_ids);
      scores.mi_score = mrmscore_.calcMIScore();
      scores.weighted_mi_score = mrmscore_.calcMIWeightedScore(normalized_library_intensity);
    }

    // check that the MS1 feature is present and that the MS1 MI should be calculated
    if (imrmfeature->getPrecursorIDs().size() > 0 && su_.use_ms1_mi)
    {
      // we need at least two precursor isotopes
      if (precursor_ids.size() > 1)
      {
        mrmscore_.initializeMIPrecursorMatrix(imrmfeature, precursor_ids);
        scores.ms1_mi_score = mrmscore_.calcMIPrecursorScore();
      }
      mrmscore_.initializeMIPrecursorContrastMatrix(imrmfeature, precursor_ids, native_ids);
      scores.ms1_mi_contrast_score = mrmscore_.calcMIPrecursorContrastScore();

      mrmscore_.initializeMIPrecursorCombinedMatrix(imrmfeature, precursor_ids, native_ids);
      scores.ms1_mi_combined_score = mrmscore_.calcMIPrecursorCombinedScore();
    }
  }

  void OpenSwathScoring::calculateChromatographicIdScores(
        OpenSwath::IMRMFeature* imrmfeature,
        const std::vector<std::string>& native_ids_identification,
        const std::vector<std::string>& native_ids_detection,
        std::vector<OpenSwath::ISignalToNoisePtr>& signal_noise_estimators,
        OpenSwath_Ind_Scores & idscores)
  {
    OPENMS_PRECONDITION(imrmfeature != nullptr, "Feature to be scored cannot be null");
    OpenSwath::MRMScoring mrmscore_;
    mrmscore_.initializeXCorrContrastMatrix(imrmfeature, native_ids_identification, native_ids_detection);

    if (su_.use_coelution_score_)
    {
      idscores.ind_xcorr_coelution_score = mrmscore_.calcSeparateXcorrContrastCoelutionScore();
    }

    if (su_.use_shape_score_)
    {
      idscores.ind_xcorr_shape_score = mrmscore_.calcSeparateXcorrContrastShapeScore();
    }

    // Signal to noise scoring
    if (su_.use_sn_score_)
    {
      idscores.ind_log_sn_score = mrmscore_.calcSeparateSNScore(imrmfeature, signal_noise_estimators);
    }

    // Mutual information scoring
    if (su_.use_mi_score_)
    {
      mrmscore_.initializeMIContrastMatrix(imrmfeature, native_ids_identification, native_ids_detection);
      idscores.ind_mi_score = mrmscore_.calcSeparateMIContrastScore();
    }
  }

  void OpenSwathScoring::calculateLibraryScores(
        OpenSwath::IMRMFeature* imrmfeature,
        const std::vector<TransitionType> & transitions,
        const CompoundType& pep,
        const double normalized_feature_rt,
        OpenSwath_Scores & scores)
  {
    OPENMS_PRECONDITION(imrmfeature != nullptr, "Feature to be scored cannot be null");
    std::vector<double> normalized_library_intensity;
    getNormalized_library_intensities_(transitions, normalized_library_intensity);

    std::vector<std::string> native_ids;
    OpenSwath::MRMScoring mrmscore_;
    for (Size i = 0; i < transitions.size(); i++) {native_ids.push_back(transitions[i].getNativeID());}

    if (su_.use_library_score_)
    {
      mrmscore_.calcLibraryScore(imrmfeature, transitions,
          scores.library_corr, scores.library_norm_manhattan, scores.library_manhattan,
          scores.library_dotprod, scores.library_sangle, scores.library_rootmeansquare);
    }

    // Retention time score
    if (su_.use_rt_score_)
    {
      // rt score is delta iRT
      double normalized_experimental_rt = normalized_feature_rt;
      double rt_score = mrmscore_.calcRTScore(pep, normalized_experimental_rt);

      scores.normalized_experimental_rt = normalized_experimental_rt;
      scores.raw_rt_score = rt_score;
      scores.norm_rt_score = rt_score / rt_normalization_factor_;
    }
  }

  void OpenSwathScoring::getNormalized_library_intensities_(const std::vector<TransitionType> & transitions,
                                                            std::vector<double>& normalized_library_intensity)
  {
    normalized_library_intensity.clear();
    for (Size i = 0; i < transitions.size(); i++)
    {
      normalized_library_intensity.push_back(transitions[i].getLibraryIntensity());
    }
    for (Size i = 0; i < normalized_library_intensity.size(); i++)
    {
      // the library intensity should never be below zero
      if (normalized_library_intensity[i] < 0.0) { normalized_library_intensity[i] = 0.0; }
    }
    OpenSwath::Scoring::normalize_sum(&normalized_library_intensity[0], boost::numeric_cast<int>(normalized_library_intensity.size()));
  }

  OpenSwath::SpectrumPtr OpenSwathScoring::fetchSpectrumSwath(OpenSwath::SpectrumAccessPtr swath_map,
                                                              double RT, int nr_spectra_to_add, const double drift_lower, const double drift_upper)
  {
    return getAddedSpectra_(swath_map, RT, nr_spectra_to_add, drift_lower, drift_upper);
  }

  OpenSwath::SpectrumPtr OpenSwathScoring::fetchSpectrumSwath(std::vector<OpenSwath::SwathMap> swath_maps,
                                                              double RT, int nr_spectra_to_add, const double drift_lower, const double drift_upper)
  {
    if (swath_maps.size() == 1)
    {
      return getAddedSpectra_(swath_maps[0].sptr, RT, nr_spectra_to_add, drift_lower, drift_upper);
    }
    else
    {
      // multiple SWATH maps for a single precursor -> this is SONAR data
      std::vector<OpenSwath::SpectrumPtr> all_spectra;
      for (size_t i = 0; i < swath_maps.size(); ++i)
      {
        OpenSwath::SpectrumPtr spec = getAddedSpectra_(swath_maps[i].sptr, RT, nr_spectra_to_add, drift_lower, drift_upper);
        all_spectra.push_back(spec);
      }
      OpenSwath::SpectrumPtr spectrum_ = SpectrumAddition::addUpSpectra(all_spectra, spacing_for_spectra_resampling_, true);
      return spectrum_;
    }
  }

  OpenSwath::SpectrumPtr filterByDrift(const OpenSwath::SpectrumPtr input, const double drift_lower, const double drift_upper)
  {
    OPENMS_PRECONDITION(drift_upper > 0, "Cannot filter by drift time if upper value is less or equal to zero");
    //OPENMS_PRECONDITION(input->getDriftTimeArray() != nullptr, "Cannot filter by drift time if no drift time is available.");

    if (input->getDriftTimeArray() == nullptr)
    {
      std::cerr << "Warning: Cannot filter by drift time if no drift time is available.\n";
      return input;
    }
      
    OpenSwath::SpectrumPtr output(new OpenSwath::Spectrum);

    OpenSwath::BinaryDataArrayPtr mz_arr = input->getMZArray();
    OpenSwath::BinaryDataArrayPtr int_arr = input->getIntensityArray();
    OpenSwath::BinaryDataArrayPtr im_arr = input->getDriftTimeArray();

    std::vector<double>::const_iterator mz_it = mz_arr->data.begin();
    std::vector<double>::const_iterator int_it = int_arr->data.begin();
    std::vector<double>::const_iterator im_it = im_arr->data.begin();
    std::vector<double>::const_iterator mz_end = mz_arr->data.end();

    OpenSwath::BinaryDataArrayPtr mz_arr_out(new OpenSwath::BinaryDataArray);
    OpenSwath::BinaryDataArrayPtr intens_arr_out(new OpenSwath::BinaryDataArray);
    OpenSwath::BinaryDataArrayPtr im_arr_out(new OpenSwath::BinaryDataArray);
    im_arr_out->description = im_arr->description;

    size_t n = mz_arr->data.size();
    im_arr_out->data.reserve(n);
    while (mz_it != mz_end)
    {
      if (*im_it > drift_lower && *im_it < drift_upper)
      {
        mz_arr_out->data.push_back( *mz_it );
        intens_arr_out->data.push_back( *int_it );
        im_arr_out->data.push_back( *im_it );
      }
      ++mz_it;
      ++int_it;
      ++im_it;
    }
    output->setMZArray(mz_arr_out);
    output->setIntensityArray(intens_arr_out);
    output->getDataArrays().push_back(im_arr_out);
    return output;
  }


  OpenSwath::SpectrumPtr OpenSwathScoring::getAddedSpectra_(OpenSwath::SpectrumAccessPtr swath_map,
                                                            double RT, int nr_spectra_to_add, const double drift_lower, const double drift_upper)
  {
    std::vector<std::size_t> indices = swath_map->getSpectraByRT(RT, 0.0);
    OpenSwath::SpectrumPtr added_spec(new OpenSwath::Spectrum);
    added_spec->getDataArrays().push_back( OpenSwath::BinaryDataArrayPtr(new OpenSwath::BinaryDataArray) );
    added_spec->getDataArrays().back()->description = "Ion Mobility";

    if (indices.empty() )
    {
      return added_spec;
    }
    int closest_idx = boost::numeric_cast<int>(indices[0]);
    if (indices[0] != 0 &&
        std::fabs(swath_map->getSpectrumMetaById(boost::numeric_cast<int>(indices[0]) - 1).RT - RT) <
        std::fabs(swath_map->getSpectrumMetaById(boost::numeric_cast<int>(indices[0])).RT - RT))
    {
      closest_idx--;
    }

    if (nr_spectra_to_add == 1)
    {
      added_spec = swath_map->getSpectrumById(closest_idx);
      if (drift_upper > 0) 
      {
        added_spec = filterByDrift(added_spec, drift_lower, drift_upper);
      }
    }
    else
    {
      std::vector<OpenSwath::SpectrumPtr> all_spectra;
      // always add the spectrum 0, then add those right and left
      all_spectra.push_back(swath_map->getSpectrumById(closest_idx));
      for (int i = 1; i <= nr_spectra_to_add / 2; i++) // cast to int is intended!
      {
        if (closest_idx - i >= 0)
        {
          all_spectra.push_back(swath_map->getSpectrumById(closest_idx - i));
        }
        if (closest_idx + i < (int)swath_map->getNrSpectra())
        {
          all_spectra.push_back(swath_map->getSpectrumById(closest_idx + i));
        }
      }

      // Filter all spectra by drift time before further processing
      if (drift_upper > 0) 
      {
        for (auto& s: all_spectra) s = filterByDrift(s, drift_lower, drift_upper);
      }

      // add up all spectra
      if (spectra_addition_method_ == "simple")
      {
        // Ensure that we have the same number of data arrays as in the input spectrum
        if (!all_spectra.empty() && all_spectra[0]->getDataArrays().size() > 2)
        {
          for (Size k = 2; k < all_spectra[0]->getDataArrays().size(); k++)
          {
            OpenSwath::BinaryDataArrayPtr tmp (new OpenSwath::BinaryDataArray());
            tmp->description = all_spectra[0]->getDataArrays()[k]->description;
            added_spec->getDataArrays().push_back(tmp);
          }
        }

        // Simply add up data and sort in the end
        for (const auto& s : all_spectra)
        {
          for (Size k = 0; k < s->getDataArrays().size(); k++)
          {
            auto& v1 = added_spec->getDataArrays()[k]->data;
            auto& v2 = s->getDataArrays()[k]->data;

            v1.reserve( v1.size() + v2.size() ); 
            v1.insert( v1.end(), v2.begin(), v2.end() );
          }
        }
        sortSpectrumByMZ(*added_spec);
      }
      else
      {
        added_spec = SpectrumAddition::addUpSpectra(all_spectra, spacing_for_spectra_resampling_, true);
      }
    }

    OPENMS_POSTCONDITION( std::adjacent_find(added_spec->getMZArray()->data.begin(),
           added_spec->getMZArray()->data.end(), std::greater<double>()) == added_spec->getMZArray()->data.end(),
           "Postcondition violated: m/z vector needs to be sorted!" )

    return added_spec;
  }

}

