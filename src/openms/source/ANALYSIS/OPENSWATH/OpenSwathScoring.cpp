// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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

// scoring
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/Scoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/MRMScoring.h>

// auxiliary
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SpectrumAddition.h>

namespace OpenMS
{

  /// Constructor
  OpenSwathScoring::OpenSwathScoring() :
    rt_normalization_factor_(1.0),
    add_up_spectra_(1),
    spacing_for_spectra_resampling_(0.005)
  {
  }

  /// Destructor
  OpenSwathScoring::~OpenSwathScoring()
  {
  }

  void OpenSwathScoring::initialize(double rt_normalization_factor,
    int add_up_spectra, double spacing_for_spectra_resampling,
    OpenSwath_Scores_Usage & su)
  {
    this->rt_normalization_factor_ = rt_normalization_factor;
    this->add_up_spectra_ = add_up_spectra;
    this->spacing_for_spectra_resampling_ = spacing_for_spectra_resampling;
    this->su_ = su;
  }

  void OpenSwathScoring::calculateDIAScores(OpenSwath::IMRMFeature* imrmfeature, const std::vector<TransitionType> & transitions,
      OpenSwath::SpectrumAccessPtr swath_map, OpenSwath::SpectrumAccessPtr ms1_map, OpenMS::DIAScoring & diascoring,
      const PeptideType& pep, OpenSwath_Scores & scores)
  {
    OPENMS_PRECONDITION(transitions.size() > 0, "There needs to be at least one transition.");

    std::vector<double> normalized_library_intensity;
    getNormalized_library_intensities_(transitions, normalized_library_intensity);

    // parameters
    int by_charge_state = 1; // for which charge states should we check b/y series
    double precursor_mz = transitions[0].precursor_mz;

    // find spectrum that is closest to the apex of the peak using binary search
    OpenSwath::SpectrumPtr spectrum_ = getAddedSpectra_(swath_map, imrmfeature->getRT(), add_up_spectra_);
    OpenSwath::SpectrumPtr* spectrum = &spectrum_;

    // Isotope correlation / overlap score: Is this peak part of an
    // isotopic pattern or is it the monoisotopic peak in an isotopic
    // pattern?
    diascoring.dia_isotope_scores(transitions, (*spectrum), imrmfeature, scores.isotope_correlation, scores.isotope_overlap);
    // Mass deviation score
    diascoring.dia_massdiff_score(transitions, (*spectrum), normalized_library_intensity,
        scores.massdev_score, scores.weighted_massdev_score);

    // Presence of b/y series score
    OpenMS::AASequence aas;
    OpenSwathDataAccessHelper::convertPeptideToAASequence(pep, aas);
    diascoring.dia_by_ion_score((*spectrum), aas, by_charge_state, scores.bseries_score, scores.yseries_score);

    // FEATURE we should not punish so much when one transition is missing!
    scores.massdev_score = scores.massdev_score / transitions.size();

    // DIA dotproduct and manhattan score
    diascoring.score_with_isotopes((*spectrum), transitions, scores.dotprod_score_dia, scores.manhatt_score_dia);

    // MS1 ppm score : check that the map is not NULL and contains spectra
    if (ms1_map && ms1_map->getNrSpectra() > 0) 
    {
      OpenSwath::SpectrumPtr ms1_spectrum = getAddedSpectra_(ms1_map, imrmfeature->getRT(), add_up_spectra_);
      diascoring.dia_ms1_massdiff_score(precursor_mz, ms1_spectrum, scores.ms1_ppm_score);
      diascoring.dia_ms1_isotope_scores(precursor_mz, ms1_spectrum, pep.getChargeState(), scores.ms1_isotope_correlation, scores.ms1_isotope_overlap);
    }
  }

  void OpenSwathScoring::calculateDIAIdScores(OpenSwath::IMRMFeature* imrmfeature, const TransitionType & transition,
      OpenSwath::SpectrumAccessPtr swath_map, OpenMS::DIAScoring & diascoring,
      OpenSwath_Scores & scores)
  {
    // find spectrum that is closest to the apex of the peak using binary search
    OpenSwath::SpectrumPtr spectrum_ = getAddedSpectra_(swath_map, imrmfeature->getRT(), add_up_spectra_);
    OpenSwath::SpectrumPtr* spectrum = &spectrum_;

    // If no charge is given, we assume it to be 1
    int putative_product_charge = 1;
    if (transition.getProductChargeState() > 0)
    {
      putative_product_charge = transition.getProductChargeState();
    }

    // Isotope correlation / overlap score: Is this peak part of an
    // isotopic pattern or is it the monoisotopic peak in an isotopic
    // pattern?
    diascoring.dia_ms1_isotope_scores(transition.getProductMZ(), (*spectrum), putative_product_charge, scores.isotope_correlation, scores.isotope_overlap);
    // Mass deviation score
    diascoring.dia_ms1_massdiff_score(transition.getProductMZ(), (*spectrum), scores.massdev_score);
  }

  void OpenSwathScoring::calculateChromatographicScores(
        OpenSwath::IMRMFeature* imrmfeature,
        const std::vector<std::string>& native_ids,
        const std::vector<double>& normalized_library_intensity,
        std::vector<OpenSwath::ISignalToNoisePtr>& signal_noise_estimators,
        OpenSwath_Scores & scores)
  {
    OpenSwath::MRMScoring mrmscore_;
    mrmscore_.initializeXCorrMatrix(imrmfeature, native_ids);

    // XCorr score (coelution)
    if (su_.use_coelution_score_)
    {
      scores.xcorr_coelution_score = mrmscore_.calcXcorrCoelutionScore();
      scores.weighted_coelution_score = mrmscore_.calcXcorrCoelutionScore_weighted(normalized_library_intensity);
    }

    // XCorr score (shape)
    // mean over the intensities at the max of the crosscorrelation
    // FEATURE : weigh by the intensity as done by mQuest
    // FEATURE : normalize with the intensity at the peak group apex?
    if (su_.use_shape_score_)
    {
      scores.xcorr_shape_score = mrmscore_.calcXcorrShape_score();
      scores.weighted_xcorr_shape = mrmscore_.calcXcorrShape_score_weighted(normalized_library_intensity);
    }

    // check that the MS1 feature is present and that the MS1 correlation should be calculated
    if (imrmfeature->getPrecursorIDs().size() > 0 && su_.use_ms1_correlation)
    {
      mrmscore_.initializeMS1XCorr(imrmfeature, native_ids, "Precursor_i0"); // perform cross-correlation on monoisotopic precursor
      scores.xcorr_ms1_coelution_score = mrmscore_.calcMS1XcorrCoelutionScore();
      scores.xcorr_ms1_shape_score = mrmscore_.calcMS1XcorrShape_score();
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
      if (scores.sn_ratio < 1) { scores.log_sn_score = 0; }
      else { scores.log_sn_score = std::log(scores.sn_ratio); }
    }
  }

  void OpenSwathScoring::calculateChromatographicIdScores(
        OpenSwath::IMRMFeature* imrmfeature,
        const std::vector<std::string>& native_ids_identification,
        const std::vector<std::string>& native_ids_detection,
        std::vector<OpenSwath::ISignalToNoisePtr>& signal_noise_estimators,
        OpenSwath_Scores & idscores)
  { 
    OpenSwath::MRMScoring mrmscore_;
    mrmscore_.initializeXCorrIdMatrix(imrmfeature, native_ids_identification, native_ids_detection);

    if (su_.use_coelution_score_)
    {
      idscores.ind_xcorr_coelution_score = mrmscore_.calcIndXcorrIdCoelutionScore();
    }

    if (su_.use_shape_score_)
    {
      idscores.ind_xcorr_shape_score = mrmscore_.calcIndXcorrIdShape_score();
    }

    // Signal to noise scoring
    if (su_.use_sn_score_)
    {
      idscores.ind_log_sn_score = mrmscore_.calcIndSNScore(imrmfeature, signal_noise_estimators);
    }
  }

  void OpenSwathScoring::calculateLibraryScores(
        OpenSwath::IMRMFeature* imrmfeature,
        const std::vector<TransitionType> & transitions,
        const PeptideType& pep,
        const double normalized_feature_rt,
        OpenSwath_Scores & scores)
  {
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

  OpenSwath::SpectrumPtr OpenSwathScoring::getAddedSpectra_(OpenSwath::SpectrumAccessPtr swath_map, 
      double RT, int nr_spectra_to_add)
  {
    std::vector<std::size_t> indices = swath_map->getSpectraByRT(RT, 0.0);
    if (indices.empty() ) 
    {
      OpenSwath::SpectrumPtr sptr(new OpenSwath::Spectrum);
      return sptr;
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
      OpenSwath::SpectrumPtr spectrum_ = swath_map->getSpectrumById(closest_idx);
      return spectrum_;
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
      OpenSwath::SpectrumPtr spectrum_ = SpectrumAddition::addUpSpectra(all_spectra, spacing_for_spectra_resampling_, true);
      return spectrum_;
    }
  }

}

