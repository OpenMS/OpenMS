// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathScoring.h>

#include <OpenMS/CONCEPT/Macros.h>

// scoring
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathScores.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DIAScoring.h>
#include <OpenMS/OPENSWATHALGO/ALGO/Scoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/IonMobilityScoring.h>

// auxiliary
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SpectrumAddition.h>

#include <utility>
namespace OpenMS
{

  /// Constructor
  OpenSwathScoring::OpenSwathScoring() :
    rt_normalization_factor_(1.0),
    spacing_for_spectra_resampling_(0.005),
    add_up_spectra_(1),
    spectra_addition_method_(SpectrumAdditionMethod::ADDITION),
    im_drift_extra_pcnt_(0.0)
  {
  }

  /// Destructor
  OpenSwathScoring::~OpenSwathScoring() = default;

  void OpenSwathScoring::initialize(double rt_normalization_factor,
                                    int add_up_spectra,
                                    double spacing_for_spectra_resampling,
                                    const double drift_extra,
                                    const OpenSwath_Scores_Usage & su,
                                    const std::string& spectrum_addition_method,
                                    bool use_ms1_ion_mobility)
  {
    this->rt_normalization_factor_ = rt_normalization_factor;
    this->add_up_spectra_ = add_up_spectra;
    if (spectrum_addition_method == "simple")
    {
      this->spectra_addition_method_ = SpectrumAdditionMethod::ADDITION;
    }
    else if (spectrum_addition_method == "resample")
    {
      this->spectra_addition_method_ = SpectrumAdditionMethod::RESAMPLE;
    }
    else
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "spectrum_addition_method must be simple or resample", spectrum_addition_method);
    }

    this->im_drift_extra_pcnt_ = drift_extra;
    this->spacing_for_spectra_resampling_ = spacing_for_spectra_resampling;
    this->su_ = su;
    this->use_ms1_ion_mobility_ = use_ms1_ion_mobility;
  }

  void OpenSwathScoring::calculateDIAScores(OpenSwath::IMRMFeature* imrmfeature,
                                            const std::vector<TransitionType>& transitions,
                                            const std::vector<OpenSwath::SwathMap>& swath_maps,
                                            const OpenSwath::SpectrumAccessPtr& ms1_map,
                                            const OpenMS::DIAScoring& diascoring,
                                            const CompoundType& compound,
                                            OpenSwath_Scores& scores,
                                            std::vector<double>& masserror_ppm,
                                            const double drift_target,// TODO is this needed
                                            const RangeMobility& im_range)
  {
    OPENMS_PRECONDITION(imrmfeature != nullptr, "Feature to be scored cannot be null");
    OPENMS_PRECONDITION(transitions.size() > 0, "There needs to be at least one transition.");
    OPENMS_PRECONDITION(swath_maps.size() > 0, "There needs to be at least one swath map.");

    std::vector<double> normalized_library_intensity;
    getNormalized_library_intensities_(transitions, normalized_library_intensity);

    // find spectrum that is closest to the apex of the peak using binary search
    std::vector<OpenSwath::SpectrumPtr> spectra = fetchSpectrumSwath(swath_maps, imrmfeature->getRT(), add_up_spectra_, im_range);

    // set the DIA parameters
    // TODO Cache these parameters
    double dia_extract_window_ = (double)diascoring.getParameters().getValue("dia_extraction_window");
    bool dia_extraction_ppm_ = diascoring.getParameters().getValue("dia_extraction_unit") == "ppm";

    // score drift time dimension
    if ( su_.use_im_scores)
    {
      IonMobilityScoring::driftScoring(spectra, transitions, scores,
                                       drift_target, im_range,
                                       dia_extract_window_, dia_extraction_ppm_,
                                       false, im_drift_extra_pcnt_);
    }


    // Mass deviation score
    diascoring.dia_massdiff_score(transitions, spectra, normalized_library_intensity, im_range, scores.massdev_score, scores.weighted_massdev_score, masserror_ppm);

    //TODO this score and the next, both rely on the CoarseIsotope of the PeptideAveragine. Maybe we could
    // DIA dotproduct and manhattan score based on library intensity and sum formula if present
    if (su_.use_ms2_isotope_scores)
    {
      diascoring.score_with_isotopes(spectra, transitions, im_range, scores.dotprod_score_dia, scores.manhatt_score_dia);

      // Isotope correlation / overlap score: Is this peak part of an
      // isotopic pattern or is it the monoisotopic peak in an isotopic
      // pattern?
      // Currently this is computed for an averagine model of a peptide so its
      // not optimal for metabolites - but better than nothing, given that for
      // most fragments we don't really know their composition
      diascoring
          .dia_isotope_scores(transitions, spectra, imrmfeature, im_range, scores.isotope_correlation, scores.isotope_overlap);
    }

    // Peptide-specific scores (only useful, when product transitions are REAL fragments, e.g. not in FFID)
    // and only if sequence is known (non-empty)
    if (compound.isPeptide() && !compound.sequence.empty() && su_.use_ionseries_scores)
    {
      // Presence of b/y series score
      OpenMS::AASequence aas;
      int by_charge_state = 1; // for which charge states should we check b/y series
      OpenSwathDataAccessHelper::convertPeptideToAASequence(compound, aas);
      diascoring.dia_by_ion_score(spectra, aas, by_charge_state, im_range, scores.bseries_score, scores.yseries_score);
    }


    RangeMobility im_range_ms1;
    if (use_ms1_ion_mobility_)
    {
      im_range_ms1 = im_range;
    }
    else // do not extract across IM in MS1
    {
      im_range_ms1 = RangeMobility();
    }

    if (ms1_map && ms1_map->getNrSpectra() > 0)
    {
      double precursor_mz = transitions[0].precursor_mz;
      double rt = imrmfeature->getRT();

      calculatePrecursorDIAScores(ms1_map, diascoring, precursor_mz, rt, compound, im_range_ms1, scores);
    }


    if ( (ms1_map && ms1_map->getNrSpectra() > 0) && ( su_.use_im_scores) )  // IM MS1 scores
    {
        double dia_extract_window_ = (double)diascoring.getParameters().getValue("dia_extraction_window");
        bool dia_extraction_ppm_ = diascoring.getParameters().getValue("dia_extraction_unit") == "ppm";
        double rt = imrmfeature->getRT();

        std::vector<OpenSwath::SpectrumPtr> ms1_spectrum = fetchSpectrumSwath(ms1_map, rt, add_up_spectra_, im_range_ms1);
        IonMobilityScoring::driftScoringMS1(ms1_spectrum,
            transitions, scores, drift_target, im_range_ms1, dia_extract_window_, dia_extraction_ppm_, false, im_drift_extra_pcnt_);

        IonMobilityScoring::driftScoringMS1Contrast(spectra, ms1_spectrum,
            transitions, scores, im_range_ms1, dia_extract_window_, dia_extraction_ppm_, im_drift_extra_pcnt_);
    }
  }

  void OpenSwathScoring::calculatePrecursorDIAScores(const OpenSwath::SpectrumAccessPtr& ms1_map,
                                   const OpenMS::DIAScoring & diascoring,
                                   double precursor_mz,
                                   double rt,
                                   const CompoundType& compound,
                                   RangeMobility im_range,
                                   OpenSwath_Scores & scores)
  {
    // change im_range based on ms1 settings
    if (!use_ms1_ion_mobility_)
    {
      im_range.clear();
    }

    // Compute precursor-level scores:
    // - compute mass difference in ppm
    // - compute isotopic pattern score
    if (ms1_map && ms1_map->getNrSpectra() > 0)
    {
      std::vector<OpenSwath::SpectrumPtr> ms1_spectrum = fetchSpectrumSwath(ms1_map, rt, add_up_spectra_, im_range);
      diascoring.dia_ms1_massdiff_score(precursor_mz, ms1_spectrum, im_range, scores.ms1_ppm_score);

      // derive precursor charge state (get from data if possible)
      int precursor_charge = 1;
      if (compound.getChargeState() != 0)
      {
        precursor_charge = compound.getChargeState();
      }

      if (compound.isPeptide())
      {
        if (!compound.sequence.empty())
        {
          diascoring.dia_ms1_isotope_scores(precursor_mz, ms1_spectrum, im_range, scores.ms1_isotope_correlation,
                                            scores.ms1_isotope_overlap,
                                            AASequence::fromString(compound.sequence).getFormula(Residue::Full, precursor_charge));
        }
        else
        {
          diascoring.dia_ms1_isotope_scores_averagine(precursor_mz, ms1_spectrum, precursor_charge, im_range,
                                                      scores.ms1_isotope_correlation,
                                                      scores.ms1_isotope_overlap);
        }
      }
      else
      {
        if (!compound.sequence.empty())
        {
          EmpiricalFormula empf{compound.sequence};
          //Note: this only sets the charge to be extracted again in the following function.
          // It is not really used in EmpiricalFormula. Also the m/z of the formula is not used since
          // it is shadowed by the exact precursor_mz.
          //TODO check if charges are the same (in case the charge was actually present in the sum_formula?)
          empf.setCharge(precursor_charge);
          diascoring.dia_ms1_isotope_scores(precursor_mz, ms1_spectrum, im_range, scores.ms1_isotope_correlation,
                                            scores.ms1_isotope_overlap,
                                            empf);
        }
        else
        {
          diascoring.dia_ms1_isotope_scores_averagine(precursor_mz, ms1_spectrum, precursor_charge, im_range,
                                                      scores.ms1_isotope_correlation,
                                                      scores.ms1_isotope_overlap);
        }
      }
    }
  }

  void OpenSwathScoring::calculateDIAIdScores(OpenSwath::IMRMFeature* imrmfeature,
                                              const TransitionType & transition,
                                              const std::vector<OpenSwath::SwathMap>& swath_maps,
                                              RangeMobility& im_range,
                                              const OpenMS::DIAScoring & diascoring,
                                              OpenSwath_Scores & scores)
  {
    OPENMS_PRECONDITION(imrmfeature != nullptr, "Feature to be scored cannot be null");
    OPENMS_PRECONDITION(swath_maps.size() > 0, "There needs to be at least one swath map.");

    if (!use_ms1_ion_mobility_)
    {
      im_range.clear();
    }

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
    std::vector<OpenSwath::SpectrumPtr> spectrum = fetchSpectrumSwath(used_swath_maps, imrmfeature->getRT(), add_up_spectra_, im_range);

    // If no charge is given, we assume it to be 1
    int putative_product_charge = 1;
    if (transition.getProductChargeState() != 0)
    {
      putative_product_charge = transition.getProductChargeState();
    }

    // Isotope correlation / overlap score: Is this peak part of an
    // isotopic pattern or is it the monoisotopic peak in an isotopic
    // pattern?
    diascoring.dia_ms1_isotope_scores_averagine(transition.getProductMZ(),
                                                spectrum,
                                                putative_product_charge,
                                                im_range,
                                                scores.isotope_correlation,
                                                scores.isotope_overlap);
    // Mass deviation score
    diascoring.dia_ms1_massdiff_score(transition.getProductMZ(), spectrum, im_range, scores.massdev_score);
  }

  void OpenSwathScoring::calculateChromatographicScores(
        OpenSwath::IMRMFeature* imrmfeature,
        const std::vector<std::string>& native_ids,
        const std::vector<std::string>& precursor_ids,
        const std::vector<double>& normalized_library_intensity,
        std::vector<OpenSwath::ISignalToNoisePtr>& signal_noise_estimators,
        OpenSwath_Scores & scores) const
  {
    OPENMS_PRECONDITION(imrmfeature != nullptr, "Feature to be scored cannot be null");
    OpenSwath::MRMScoring mrmscore_;
    if (su_.use_coelution_score_ || su_.use_shape_score_ || (!imrmfeature->getPrecursorIDs().empty() && su_.use_ms1_correlation))
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
    if (!imrmfeature->getPrecursorIDs().empty() && su_.use_ms1_correlation)
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
    if (!imrmfeature->getPrecursorIDs().empty() && su_.use_ms1_mi)
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
        OpenSwath_Ind_Scores & idscores) const
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
    native_ids.reserve(transitions.size());
    for (const auto& trans : transitions)
    {
      native_ids.push_back(trans.getNativeID());
    }

    if (su_.use_library_score_)
    {
      OpenSwath::MRMScoring::calcLibraryScore(imrmfeature, transitions,
          scores.library_corr, scores.library_norm_manhattan, scores.library_manhattan,
          scores.library_dotprod, scores.library_sangle, scores.library_rootmeansquare);
    }

    // Retention time score
    if (su_.use_rt_score_)
    {
      // rt score is delta iRT
      double normalized_experimental_rt = normalized_feature_rt;
      double rt_score = OpenSwath::MRMScoring::calcRTScore(pep, normalized_experimental_rt);

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

  SpectrumSequence OpenSwathScoring::fetchSpectrumSwath(OpenSwath::SpectrumAccessPtr swathmap, double RT, int nr_spectra_to_add, const RangeMobility& im_range)
  {

    SpectrumSequence all_spectra = swathmap->getMultipleSpectra(RT, nr_spectra_to_add);
    if (spectra_addition_method_ == SpectrumAdditionMethod::ADDITION)
    {
      return all_spectra; // return vector, addition is done later
    }
    else // (spectra_addition_method_ == SpectrumAdditionMethod::RESAMPLE)
    {
      std::vector<OpenSwath::SpectrumPtr> spectrum_out;
      //added_spec = SpectrumAddition::addUpSpectra(all_spectra, spacing_for_spectra_resampling_, true);
      spectrum_out.push_back(SpectrumAddition::addUpSpectra(all_spectra, im_range, spacing_for_spectra_resampling_, true));
      return spectrum_out;
    }
  }

  SpectrumSequence OpenSwathScoring::fetchSpectrumSwath(std::vector<OpenSwath::SwathMap> swath_maps, double RT, int nr_spectra_to_add, const RangeMobility& im_range)
  {
    OPENMS_PRECONDITION(nr_spectra_to_add >= 1, "nr_spectra_to_add must be at least 1.")
    OPENMS_PRECONDITION(!swath_maps.empty(), "swath_maps vector cannot be empty")

    // This is not SONAR data
    if (swath_maps.size() == 1)
    {
      return fetchSpectrumSwath(swath_maps[0].sptr, RT, nr_spectra_to_add, im_range);
    }
    else
    {
      // data is not IM enhanced
      if (!im_range.isEmpty())
      {
        // multiple SWATH maps for a single precursor -> this is SONAR data, in all cases only return a single spectrum
        SpectrumSequence all_spectra;

        if (spectra_addition_method_ == SpectrumAdditionMethod::ADDITION)
        {
          for (size_t i = 0; i < swath_maps.size(); ++i)
          {
            SpectrumSequence spectrumSequence = swath_maps[i].sptr->getMultipleSpectra(RT, nr_spectra_to_add, im_range.getMin(), im_range.getMax());
          }
        }
        else // (spectra_addition_method_ == SpectrumAdditionMethod::RESAMPLE)
        {
          for (size_t i = 0; i < swath_maps.size(); ++i)
          {
            SpectrumSequence spectrumSequence = swath_maps[i].sptr->getMultipleSpectra(RT, nr_spectra_to_add, im_range.getMin(), im_range.getMax());
            all_spectra.push_back(SpectrumAddition::addUpSpectra(spectrumSequence, spacing_for_spectra_resampling_, true));
          }
        }
        return { SpectrumAddition::addUpSpectra(all_spectra, spacing_for_spectra_resampling_, true) };
      }
      else // im_range.isEmpty()
      {
        // multiple SWATH maps for a single precursor -> this is SONAR data, in all cases only return a single spectrum
        SpectrumSequence all_spectra;

        if (spectra_addition_method_ == SpectrumAdditionMethod::ADDITION)
        {
          for (size_t i = 0; i < swath_maps.size(); ++i)
          {
            SpectrumSequence spectrumSequence = swath_maps[i].sptr->getMultipleSpectra(RT, nr_spectra_to_add);
            all_spectra.push_back(SpectrumAddition::concatenateSpectra(spectrumSequence));
          }
        }
        else // (spectra_addition_method_ == SpectrumAdditionMethod::RESAMPLE)
        {
          for (size_t i = 0; i < swath_maps.size(); ++i)
          {
            SpectrumSequence spectrumSequence = swath_maps[i].sptr->getMultipleSpectra(RT, nr_spectra_to_add);
            all_spectra.push_back(SpectrumAddition::addUpSpectra(spectrumSequence, spacing_for_spectra_resampling_, true));
          }
        }
        return { SpectrumAddition::addUpSpectra(all_spectra, spacing_for_spectra_resampling_, true) };
      }
    }
  }
}
