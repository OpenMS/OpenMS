// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/IonMobilityScoring.h>

#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/CONCEPT/LogStream.h>

// scoring
#include <OpenMS/OPENSWATHALGO/ALGO/Scoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SONARScoring.h>

// auxiliary
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SpectrumAddition.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DIAHelper.h>

// #define DEBUG_IMSCORING

namespace OpenMS
{
  std::vector<double> IonMobilityScoring::computeGrid_(const std::vector< IonMobilogram >& mobilograms, double eps)
  {
    // Extract all ion mobility values across all transitions and produce a
    // grid of all permitted ion mobility values
    std::vector<double> im_grid;
    std::vector< double > mobilityValues;
    for (const auto & im_profile : mobilograms)
    {
      mobilityValues.reserve(mobilityValues.size() + im_profile.size());
      for (const auto & k : im_profile) mobilityValues.push_back(k.im);
    }

    // sort all extracted values
    std::sort(mobilityValues.begin(), mobilityValues.end());

    // Reduce mobility values to grid (consider equal if closer than eps)
    //
    // In some cases there are not enough datapoints available (one of the
    // transitions has no datapoints)
    if (!mobilityValues.empty())
    {
      im_grid.push_back( mobilityValues[0] );
      for (Size k = 1; k < mobilityValues.size(); k++)
      {
        double diff = fabs(mobilityValues[k] - mobilityValues[k-1]);
        if (diff > eps)
        {
          im_grid.push_back( mobilityValues[k] );
        }
      }
    }
    return im_grid;
  }

  /*
   @brief Extracts ion mobility values projected onto a grid

   For a given ion mobility profile and a grid, compute an ion mobilogram
   across the grid for each ion mobility data point. Returns two data arrays
   for the ion mobilogram: intensity (y) and ion mobility (x). Zero values are
   inserted if no data point was found for a given grid value.

   @param profile The ion mobility data
   @param im_grid The grid to be used
   @param al_int_values The intensity vector (y)
   @param al_im_values The ion mobility vector (x)
   @param eps Epsilon used for computing the ion mobility grid
   @param max_peak_idx The grid position of the maximum

  */
  void IonMobilityScoring::alignToGrid_(const IonMobilogram& profile,
               const std::vector<double>& im_grid,
               std::vector< double >& al_int_values,
               std::vector< double >& al_im_values,
               double eps,
               Size & max_peak_idx)
  {
    auto pr_it = profile.begin();
    max_peak_idx = 0;
    double max_int = 0;
    for (Size k = 0; k < im_grid.size(); k++)
    {
      // In each iteration, the IM value of pr_it should be equal to or
      // larger than the master container. If it is equal, we add the current
      // data point, if it is larger we add zero and advance the counter k.
      if (pr_it != profile.end() && fabs(pr_it->im - im_grid[k] ) < eps*10)
      {
        al_int_values.push_back(pr_it->intensity);
        al_im_values.push_back(pr_it->im);
        ++pr_it;
      }
      else
      {
        al_int_values.push_back(0.0);
        al_im_values.push_back( im_grid[k] );
      }
      // OPENMS_LOG_DEBUG << "grid position " << im_grid[k] << " profile position " << pr_it->first << std::endl;

      // check that we did not advance past
      if (pr_it != profile.end() && (im_grid[k] - pr_it->im) > eps*10)
      {
        std::cout << " This should never happen, pr_it has advanced past the master container: " << im_grid[k]  << "  / " <<  pr_it->im  << std::endl;
        throw Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
      }

      // collect maxima
      if (pr_it != profile.end() && pr_it->intensity > max_int)
      {
        max_int = pr_it->intensity;
        max_peak_idx = k;
      }
    }
  }

  // compute ion mobilogram as well as im weighted average. This is based off of integrateWindows() in DIAHelper.cpp
  void IonMobilityScoring::computeIonMobilogram(const SpectrumSequence& spectra,
                              const RangeMZ& mz_range,
                              const RangeMobility& im_range,
                              double & im,
                              double & intensity,
                              IonMobilogram& res,
                              double eps)
  {

    // rounding multiplier for the ion mobility value
    // TODO: how to improve this -- will work up to 42949.67296
    double IM_IDX_MULT = 1/eps;

    // We need to store all values that map to the same ion mobility in the
    // same spot in the ion mobilogram (they are not sorted by ion mobility in
    // the input data), therefore create a map to map to bins.
    std::map< int, double> im_chrom;

    for (const auto& spectrum : spectra)
    {
      OPENMS_PRECONDITION(spectrum->getDriftTimeArray() != nullptr, "Cannot filter by drift time if no drift time is available.");
      OPENMS_PRECONDITION(spectrum->getMZArray()->data.size() == spectrum->getIntensityArray()->data.size(), "MZ and Intensity array need to have the same length.");
      OPENMS_PRECONDITION(spectrum->getMZArray()->data.size() == spectrum->getDriftTimeArray()->data.size(), "MZ and Drift Time array need to have the same length.");

      auto mz_arr_end = spectrum->getMZArray()->data.end();
      auto int_it = spectrum->getIntensityArray()->data.begin();
      auto im_it = spectrum->getDriftTimeArray()->data.begin();

      // this assumes that the spectra are sorted!
      auto mz_it = std::lower_bound(spectrum->getMZArray()->data.begin(), mz_arr_end, mz_range.getMin());
      // auto mz_it_end = std::lower_bound(mz_it, mz_arr_end, mz_end);

      // also advance intensity and ion mobility iterator now
      auto iterator_pos = std::distance(spectrum->getMZArray()->data.begin(), mz_it);
      std::advance(int_it, iterator_pos);
      std::advance(im_it, iterator_pos);

      // Start iteration from mz start, end iteration when mz value is larger than mz_end, only store only storing ion mobility values that are in the range
      double mz_end = mz_range.getMax();
      while ( ( *mz_it < mz_end ) && (mz_it < mz_arr_end) )
      {
        if (im_range.contains(*im_it))
        {
          intensity += (*int_it);
          im += (*int_it) * (*im_it);
          im_chrom[ int((*im_it)*IM_IDX_MULT) ] += *int_it;
        }
        ++mz_it;
        ++int_it;
        ++im_it;
      }
    }

    // compute the weighted average ion mobility
    if (intensity > 0.)
    {
      im /= intensity;
    }
    else
    {
      im = -1;
      intensity = 0;
    }

    res.reserve(res.size() + im_chrom.size());
    for (const auto& k : im_chrom)
    {
      res.emplace_back(k.first / IM_IDX_MULT, k.second );
    }
  }

  /// Constructor
  IonMobilityScoring::IonMobilityScoring() = default;

  /// Destructor
  IonMobilityScoring::~IonMobilityScoring() = default;

  void IonMobilityScoring::driftScoringMS1Contrast(const SpectrumSequence& spectra, const SpectrumSequence& ms1spectrum,
                                                   const std::vector<TransitionType> & transitions,
                                                   OpenSwath_Scores & scores,
                                                   RangeMobility im_range,
                                                   const double dia_extract_window_,
                                                   const bool dia_extraction_ppm_,
                                                   const double drift_extra)
  {
    OPENMS_PRECONDITION(!spectra.empty(), "Spectra cannot be empty")
    OPENMS_PRECONDITION(!ms1spectrum.empty(), "MS1 spectrum cannot be empty")
    OPENMS_PRECONDITION(!transitions.empty(), "Need at least one transition");

    //TODO not sure what error format is best
    for (const auto& s:spectra)
    {
      if (s->getDriftTimeArray() == nullptr)
      {
        OPENMS_LOG_DEBUG << " ERROR: Drift time is missing in ion mobility spectrum!" << std::endl;
        return;
      }
    }

    for (const auto& s:ms1spectrum)
    {
      if (s->getDriftTimeArray() == nullptr)
      {
        OPENMS_LOG_DEBUG << " ERROR: Drift time is missing in MS1 ion mobility spectrum!" << std::endl;
        return;
      }
    }

    double eps = 1e-5; // eps for two grid cells to be considered equal

    // extend IM range by drift_extra
    im_range.scaleBy(drift_extra * 2. + 1); // multiple by 2 because want drift extra to be extended by that amount on either side

    std::vector< IonMobilogram > mobilograms;

    // Step 1: MS2 extraction
    for (std::size_t k = 0; k < transitions.size(); k++)
    {
      double im(0), intensity(0);
      IonMobilogram res;
      const TransitionType transition = transitions[k];
      // Calculate the difference of the theoretical ion mobility and the actually measured ion mobility
      RangeMZ mz_range = DIAHelpers::createMZRangePPM(transition.getProductMZ(), dia_extract_window_, dia_extraction_ppm_);

      computeIonMobilogram(spectra, mz_range, im_range, im, intensity, res, eps);
      mobilograms.push_back( std::move(res) );
    }

    // Step 2: MS1 extraction
    double im(0), intensity(0);
    IonMobilogram ms1_profile;
    RangeMZ mz_range = DIAHelpers::createMZRangePPM(transitions[0].getPrecursorMZ(), dia_extract_window_, dia_extraction_ppm_);

    computeIonMobilogram(ms1spectrum, mz_range, im_range, im, intensity, ms1_profile, eps); // TODO: aggregate over isotopes
    mobilograms.push_back(ms1_profile);

    std::vector<double> im_grid = computeGrid_(mobilograms, eps); // ensure grid is based on all profiles!
    mobilograms.pop_back();

    // Step 3: Align the IonMobilogram vectors to the grid
    std::vector< std::vector< double > > aligned_mobilograms;
    for (const auto & mobilogram : mobilograms)
    {
      std::vector< double > arrInt, arrIM;
      Size max_peak_idx = 0;
      alignToGrid_(mobilogram, im_grid, arrInt, arrIM, eps, max_peak_idx);
      aligned_mobilograms.push_back(arrInt);
    }

    std::vector< double > ms1_int_values, ms1_im_values;
    Size max_peak_idx = 0;
    alignToGrid_(ms1_profile, im_grid, ms1_int_values, ms1_im_values, eps, max_peak_idx);

    // Step 4: MS1 contrast scores
    {
      OpenSwath::MRMScoring mrmscore_;
      mrmscore_.initializeXCorrPrecursorContrastMatrix({ms1_int_values}, aligned_mobilograms);
      OPENMS_LOG_DEBUG << "all-all: Contrast Scores : coelution precursor : " << mrmscore_.calcXcorrPrecursorContrastCoelutionScore() << " / shape  precursor " <<
        mrmscore_.calcXcorrPrecursorContrastShapeScore() << std::endl;
      scores.im_ms1_contrast_coelution = mrmscore_.calcXcorrPrecursorContrastCoelutionScore();
      scores.im_ms1_contrast_shape = mrmscore_.calcXcorrPrecursorContrastShapeScore();
    }

    // Step 5: contrast precursor vs summed fragment ions
    std::vector<double> fragment_values;
    fragment_values.resize(ms1_int_values.size(), 0);
    for (Size k = 0; k < fragment_values.size(); k++)
    {
      for (Size i = 0; i < aligned_mobilograms.size(); i++)
      {
        fragment_values[k] += aligned_mobilograms[i][k];
      }
    }

    OpenSwath::MRMScoring mrmscore_;
    // horribly broken: provides vector of length 1, but expects at least length 2 in calcXcorrPrecursorContrastCoelutionScore()
    mrmscore_.initializeXCorrPrecursorContrastMatrix({ms1_int_values}, {fragment_values});
    OPENMS_LOG_DEBUG << "Contrast Scores : coelution precursor : " << mrmscore_.calcXcorrPrecursorContrastSumFragCoelutionScore() << " / shape  precursor " <<
       mrmscore_.calcXcorrPrecursorContrastSumFragShapeScore() << std::endl;

    // in order to prevent assertion error call calcXcorrPrecursorContrastSumFragCoelutionScore, same as calcXcorrPrecursorContrastCoelutionScore() however different assertion
    scores.im_ms1_sum_contrast_coelution = mrmscore_.calcXcorrPrecursorContrastSumFragCoelutionScore();

    // in order to prevent assertion error call calcXcorrPrecursorContrastSumFragShapeScore(), same as calcXcorrPrecursorContrastShapeScore() however different assertion.
    scores.im_ms1_sum_contrast_shape = mrmscore_.calcXcorrPrecursorContrastSumFragShapeScore();
  }

  void IonMobilityScoring::driftScoringMS1(const SpectrumSequence & spectra,
                                           const std::vector<TransitionType> & transitions,
                                           OpenSwath_Scores & scores,
                                           const double drift_target,
                                           RangeMobility im_range,
                                           const double dia_extract_window_,
                                           const bool dia_extraction_ppm_,
                                           const bool /* use_spline */,
                                           const double drift_extra)
  {
    OPENMS_PRECONDITION(!spectra.empty(), "Spectra cannot be empty")
    OPENMS_PRECONDITION(!transitions.empty(), "Need at least one transition");

    for (auto s:spectra){
      if (s->getDriftTimeArray() == nullptr)
      {
        OPENMS_LOG_DEBUG << " ERROR: Drift time is missing in ion mobility spectrum!" << std::endl;
        return;
      }
    }

    im_range.scaleBy(drift_extra * 2. + 1); // multiple by 2 because want drift extra to be extended by that amount on either side

    double im(0), intensity(0), mz(0);
    RangeMZ mz_range = DIAHelpers::createMZRangePPM(transitions[0].getPrecursorMZ(), dia_extract_window_, dia_extraction_ppm_);

    DIAHelpers::integrateWindow(spectra, mz, im, intensity, mz_range, im_range);

    // Record the measured ion mobility
    scores.im_ms1_drift = im;

    // Calculate the difference of the theoretical ion mobility and the actually measured ion mobility
    scores.im_ms1_delta_score = fabs(drift_target - im);
    scores.im_ms1_delta = drift_target - im;
  }

  void IonMobilityScoring::driftScoring(const SpectrumSequence& spectra,
                                        const std::vector<TransitionType> & transitions,
                                        OpenSwath_Scores & scores,
                                        const double drift_target,
                                        RangeMobility im_range,
                                        const double dia_extract_window_,
                                        const bool dia_extraction_ppm_,
                                        const bool /* use_spline */,
                                        const double drift_extra)
  {
    OPENMS_PRECONDITION(!spectra.empty(), "Spectra cannot be empty");
    for (auto s:spectra)
    {
      if (s->getDriftTimeArray() == nullptr)
      {
        OPENMS_LOG_DEBUG << " ERROR: Drift time is missing in ion mobility spectrum!" << std::endl;
        return;
      }
    }

    double eps = 1e-5; // eps for two grid cells to be considered equal

    im_range.scaleBy(drift_extra * 2. + 1); // multiple by 2 because want drift extra to be extended by that amount on either side

    double delta_drift = 0;
    double delta_drift_abs = 0;
    // IonMobilogram: a data structure that holds points <im_value, intensity>
    std::vector< IonMobilogram > mobilograms;
    double computed_im = 0;
    double computed_im_weighted = 0;
    double sum_intensity = 0;
    int tr_used = 0;

    // Step 1: MS2 extraction
    for (std::size_t k = 0; k < transitions.size(); k++)
    {
      const TransitionType transition = transitions[k];
      IonMobilogram res;
      double im(0), intensity(0);

      // Calculate the difference of the theoretical ion mobility and the actually measured ion mobility
      RangeMZ mz_range = DIAHelpers::createMZRangePPM(transition.getProductMZ(), dia_extract_window_, dia_extraction_ppm_);

      //double left(transition.getProductMZ()), right(transition.getProductMZ());
      //DIAHelpers::adjustExtractionWindow(right, left, dia_extract_window_, dia_extraction_ppm_);
      computeIonMobilogram(spectra, mz_range, im_range, im, intensity, res, eps);
      mobilograms.push_back(res);

      // TODO what do to about those that have no signal ?
      if (intensity <= 0.0) {continue;} // note: im is -1 then

      tr_used++;

      delta_drift_abs += fabs(drift_target - im);
      delta_drift += drift_target - im;
      OPENMS_LOG_DEBUG << "  -- have delta drift time " << fabs(drift_target -im ) << " with im " << im << std::endl;
      computed_im += im;
      computed_im_weighted += im * intensity;
      sum_intensity += intensity;
      // delta_drift_weighted += delta_drift * normalized_library_intensity[k];
      // weights += normalized_library_intensity[k];
    }

    if (tr_used != 0)
    {
      delta_drift /= tr_used;
      delta_drift_abs /= tr_used;
      computed_im /= tr_used;
      computed_im_weighted /= sum_intensity;
    }
    else
    {
      delta_drift = -1;
      delta_drift_abs = -1;
      computed_im = -1;
      computed_im_weighted = -1;
    }

    OPENMS_LOG_DEBUG << " Scoring delta drift time " << delta_drift << std::endl;
    OPENMS_LOG_DEBUG << " Scoring weighted delta drift time " << computed_im_weighted << " -> get difference " << std::fabs(computed_im_weighted - drift_target)<< std::endl;
    scores.im_delta_score = delta_drift_abs;
    scores.im_delta = delta_drift;
    scores.im_drift = computed_im;
    scores.im_drift_weighted = computed_im_weighted;

    // Step 2: Align the IonMobilogram vectors to the grid
    std::vector<double> im_grid = computeGrid_(mobilograms, eps);
    std::vector< std::vector< double > > aligned_mobilograms;
    for (const auto & mobilogram : mobilograms)
    {
      std::vector< double > arr_int, arr_IM;
      Size max_peak_idx = 0;
      alignToGrid_(mobilogram, im_grid, arr_int, arr_IM, eps, max_peak_idx);
      if (!arr_int.empty()) aligned_mobilograms.push_back(arr_int);
    }

    // Step 3: Compute cross-correlation scores based on ion mobilograms
    if (aligned_mobilograms.size() < 2)
    {
      scores.im_xcorr_coelution_score = 0;
      scores.im_xcorr_shape_score = std::numeric_limits<double>::quiet_NaN();
      return;
    }


    OpenSwath::MRMScoring mrmscore_;
    mrmscore_.initializeXCorrMatrix(aligned_mobilograms);

    double xcorr_coelution_score = mrmscore_.calcXcorrCoelutionScore();
    double xcorr_shape_score = mrmscore_.calcXcorrShapeScore(); // can be nan!

    scores.im_xcorr_coelution_score = xcorr_coelution_score;
    scores.im_xcorr_shape_score = xcorr_shape_score;
  }
}
