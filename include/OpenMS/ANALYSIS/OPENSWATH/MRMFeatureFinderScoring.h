// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

#ifndef OPENMS_ANALYSIS_OPENSWATH_MRMFEATUREFINDERSCORING_H
#define OPENMS_ANALYSIS_OPENSWATH_MRMFEATUREFINDERSCORING_H

#define run_identifier "unique_run_identifier"
#define USE_SP_INTERFACE

// move to TOPPTool
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/FORMAT/TraMLFile.h>

#include <OpenMS/KERNEL/MRMTransitionGroup.h>
#include <OpenMS/KERNEL/MRMFeature.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>

// peak picking & noise estimation
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h>

// data access
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/DataStructures.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/MRMFeatureAccessOpenMS.h>

// scoring
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/Scoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/MRMScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DIAScoring.h>

// auxiliary
#include <OpenMS/ANALYSIS/OPENSWATH/SpectrumAddition.h>

#ifdef _OPENMP
#include <omp.h>
#endif

bool SortDoubleDoublePairFirst(const std::pair<double, double>& left, const std::pair<double, double>& right);

namespace OpenMS
{

  /**
  @brief A structure to store which scores should be used by the Algorithm
  */
  struct OpenSwath_Scores_Usage
  {
    // Which scores to use
    bool use_coelution_score_;
    bool use_shape_score_;
    bool use_rt_score_;
    bool use_library_score_;
    bool use_elution_model_score_;
    bool use_intensity_score_;
    bool use_total_xic_score_;
    bool use_nr_peaks_score_;
    bool use_sn_score_;
  };

  /**
  @brief A structure to hold the different scores computed by the FeatureFinder
  */
  struct OpenSwath_Scores
  {
    double elution_model_fit_score;
    double library_corr;
    double library_rmsd;
    double norm_rt_score;
    double isotope_correlation;
    double isotope_overlap;
    double massdev_score;
    double xcorr_coelution_score;
    double xcorr_shape_score;
    double yseries_score;
    double bseries_score;
    double log_sn_score;

    double weighted_coelution_score;
    double weighted_xcorr_shape;
    double weighted_massdev_score;

    double library_manhattan;
    double library_dotprod;
    double intensity;
    double total_xic;
    double nr_peaks;
    double sn_ratio;

    double rt_difference;
    double normalized_experimental_rt;
    double raw_rt_score;

    double dotprod_score_dia;
    double manhatt_score_dia;

    OpenSwath_Scores() :
      elution_model_fit_score(0),
      library_corr(0),
      library_rmsd(0),
      norm_rt_score(0),
      isotope_correlation(0),
      isotope_overlap(0),
      massdev_score(0),
      xcorr_coelution_score(0),
      xcorr_shape_score(0),
      yseries_score(0),
      bseries_score(0),
      log_sn_score(0),
      weighted_coelution_score(0),
      weighted_xcorr_shape(0),
      weighted_massdev_score(0),
      library_manhattan(0),
      library_dotprod(0),
      intensity(0),
      total_xic(0),
      nr_peaks(0),
      sn_ratio(0),
      dotprod_score_dia(0),
      manhatt_score_dia(0)
    {
    }


    double get_quick_lda_score(double library_corr, double library_rmsd, double norm_rt_score, double xcorr_coelution_score,
                               double xcorr_shape_score, double log_sn_score)
    {
      // some scores based on manual evaluation of 80 chromatograms
      // quick LDA average model on 100 2xCrossvalidated runs (0.85 TPR/0.17 FDR)
      // true: mean 4.2 with sd 1.055
      // false: mean -0.07506772  with sd 1.055
      // below -0.5 removes around 30% of the peaks
      // below 0    removes around 50% of the peaks
      // below 0.5  removes around 70% of the peaks
      // below 1.0  removes around 85% of the peaks
      // below 1.5  removes around 93% of the peaks
      // below 2.0  removes around 97% of the peaks
      //
      // NOTE this score means "better" if it is more negative!
      double lda_quick_score =
        library_corr                    * -0.5319046 +
        library_rmsd                    *  2.1643962 +
        norm_rt_score                   *  8.0353047 +
        xcorr_coelution_score           *  0.1458914 +
        xcorr_shape_score               * -1.6901925 +
        log_sn_score                    * -0.8002824;
      return lda_quick_score;
    }

    double calculate_lda_prescore(OpenSwath_Scores scores)
    {

      // LDA average model on 100 2xCrossvalidated runs (0.91 TPR/0.20 FDR)
      /*
      double xx_old_lda_prescore =
      intensity_score       * -2.296679          +
      library_corr          * -0.1223876         +
      library_rmsd          *  2.013638          +
      nr_peaks_score        *  0.01683357        +
      rt_score              *  0.00143999        +
      sn_score              * -0.1619762         +
      total_xic_score       *  0.00000003697898  +
      xcorr_coelution_score *  0.05909583        +
      xcorr_shape_score     * -0.4699841;

      // NOTE this score means "better" if it is more negative!
      */

      return scores.library_corr                     * -0.34664267 +
             scores.library_rmsd                     *  2.98700722 +
             scores.norm_rt_score                    *  7.05496384 +
             scores.xcorr_coelution_score            *  0.09445371 +
             scores.xcorr_shape_score                * -5.71823862 +
             scores.log_sn_score                     * -0.72989582 +
             scores.elution_model_fit_score          *  1.88443209;
    }

    double calculate_swath_lda_prescore(OpenSwath_Scores scores)
    {

      // Swath - LDA average model on 100 2xCrossvalidated runs (0.76 TPR/0.20 FDR) [without elution model]
      /*
      double xx_old_swath_prescore =
      intensity_score              * -3.148838e+00  +
      library_corr                 * -7.562403e-02  +
      library_rmsd                 *  1.786286e+00  +
      nr_peaks_score               * -7.674263e-03  +
      rt_score                     *  1.748377e-03  +
      sn_score                     * -1.372636e-01  +
      total_xic_score              *  7.278437e-08  +
      xcorr_coelution_score        *  1.181813e-01  +
      weighted_coelution_score     * -7.661783e-02  +
      xcorr_shape_score            * -6.903933e-02  +
      weighted_xcorr_shape         * -4.234820e-01  +
      bseries_score                * -2.022380e-02  +
      massdev_score                *  2.844948e-02  +
      massdev_score_weighted       *  1.133209e-02  +
      yseries_score                * -9.510874e-02  +
      isotope_corr                 * -1.619902e+00  +
      isotope_overlap              *  2.890688e-01  ;

      // NOTE this score means "better" if it is more negative!
      */

      return scores.library_corr              * -0.19011762 +
             scores.library_rmsd              *  2.47298914 +
             scores.norm_rt_score             *  5.63906731 +
             scores.isotope_correlation       * -0.62640133 +
             scores.isotope_overlap           *  0.36006925 +
             scores.massdev_score             *  0.08814003 +
             scores.xcorr_coelution_score     *  0.13978311 +
             scores.xcorr_shape_score         * -1.16475032 +
             scores.yseries_score             * -0.19267813 +
             scores.log_sn_score              * -0.61712054;

/*


Gold standard, best sample
 main_var_xx_swath_prelim_score  0.291440015642621
 var_bseries_score 0.0496492555026149
 var_dotprod_score -0.522561744728316
 var_elution_model_fit_score -1.99429446109581
 var_intensity_score 1.70915451039584
 var_isotope_correlation_score 0.966260829910062
 var_isotope_overlap_score -14.216079147368
 var_library_corr  0.061432632721274
 var_library_dotprod -3.79958938222036
 var_library_manhattan -1.36520528433508
 var_library_rmsd  -6.44998534845163
 var_log_sn_score  -0.0389995774588385
 var_manhatt_score -0.0944805864772705
 var_massdev_score 0.0144460056621709
 var_massdev_score_weighted  -0.0494772144218002
 var_norm_rt_score -9.04596725429934
 var_xcorr_coelution -0.141763244951207
 var_xcorr_coelution_weighted  0.00261409408565438
 var_xcorr_shape 4.89741810577371
 var_xcorr_shape_weighted  0.342723332762697
 var_yseries_score -0.188316503432445


Strep  Strep0_Repl2_R02/runlogs_mprophet.tar.gz 
main_var_xx_swath_prelim_score  0.231523019269729
var_bseries_score   -0.0488528503276347
var_elution_model_fit_score -0.47977060647858
var_intensity_score -0.80664074459128
var_isotope_correlation_score   2.34488326031997
var_isotope_overlap_score   -2.14735763746488
var_library_corr    -0.395167010986141
var_library_rmsd    -13.1295053007338
var_log_sn_score    0.265784828465348
var_massdev_score   0.0150193500103614
var_massdev_score_weighted  -0.109859906028132
var_norm_rt_score   -25.7107556062008
var_xcorr_coelution 0.244590396074410
var_xcorr_coelution_weighted    -0.918578472543494
var_xcorr_shape 2.18720521365230
var_xcorr_shape_weighted    -0.815295893352108
var_yseries_score   -0.0620070175846356

Strep10_Repl2_R02/runlogs_mprophet.tar.gz 
main_var_xx_swath_prelim_score  0.293470108599468
var_bseries_score   -0.0129641361717189
var_elution_model_fit_score -0.44993587229358
var_intensity_score -0.828540564651968
var_isotope_correlation_score   2.76284687671386
var_isotope_overlap_score   -2.26460097307479
var_library_corr    -0.445369627383142
var_library_rmsd    -13.2905041886848
var_log_sn_score    0.224626177093898
var_massdev_score   0.0185003919755981
var_massdev_score_weighted  -0.0899477179756381
var_norm_rt_score   -24.4807649346717
var_xcorr_coelution 0.218195211767293
var_xcorr_coelution_weighted    -0.91949559943762
var_xcorr_shape 1.77358514815991
var_xcorr_shape_weighted    -0.616535104461374
var_yseries_score   -0.0652111196389966




// FINAL AQUA gold standard classifier
human
main_var_xx_swath_prelim_score  0.4384384475524
var_bseries_score   0.00227405501436837
var_elution_model_fit_score -2.06412570248571
var_intensity_score -1.26021147555789
var_isotope_correlation_score   1.21887083303546
var_isotope_overlap_score   -1.60051046353231
var_library_corr    -0.33958843974352
var_library_rmsd    -5.20235596662978
var_log_sn_score    0.24021015633787
var_massdev_score   0.0399855393620327
var_massdev_score_weighted  -0.0907785715261295
var_norm_rt_score   -16.2155920223681
var_xcorr_coelution 0.0805852135076143
var_xcorr_coelution_weighted    -0.387927719728573
var_xcorr_shape 1.885899937033
var_xcorr_shape_weighted    2.45579580649067
var_yseries_score   0.138306574987678

yeast
main_var_xx_swath_prelim_score  0.369009421609329
var_bseries_score   0.0157508674154482
var_elution_model_fit_score -1.67348268698707
var_intensity_score -1.11972743418717
var_isotope_correlation_score   1.68717154416093
var_isotope_overlap_score   -1.38410070381813
var_library_corr    -0.454409692201745
var_library_rmsd    -6.08160902837145
var_log_sn_score    0.157259477914274
var_massdev_score   0.0543919580711367
var_massdev_score_weighted  -0.137296627160332
var_norm_rt_score   -28.4381743938298
var_xcorr_coelution 0.0256469469673884
var_xcorr_coelution_weighted    -0.362865323100099
var_xcorr_shape 1.88863198062243
var_xcorr_shape_weighted    1.3518953353109
var_yseries_score   0.115472572686466

water
main_var_xx_swath_prelim_score  0.174880281226536
var_bseries_score   -0.0606466737704899
var_elution_model_fit_score -0.123252502705892
var_intensity_score 1.91714146537607
var_isotope_correlation_score   0.914387652486204
var_isotope_overlap_score   -1.46521560409083
var_library_corr    -0.485498555013885
var_library_rmsd    -8.3847526088391
var_log_sn_score    0.00644514889704832
var_massdev_score   0.0177435175558717
var_massdev_score_weighted  -0.0899451169038299
var_norm_rt_score   -15.1458716759687
var_xcorr_coelution -0.370050235089866
var_xcorr_coelution_weighted    0.21512520647974
var_xcorr_shape 0.563413547839886
var_xcorr_shape_weighted    -0.270773625703933
var_yseries_score   -0.0327896378737766



*/
    }

  };

  /**
  @brief A class that calls the scoring routines
  */
  class SwathScorer 
  {
    typedef OpenSwath::LightPeptide PeptideType;
    typedef OpenSwath::LightTransition TransitionType;

    DoubleReal rt_normalization_factor_;
    int add_up_spectra_;
    DoubleReal spacing_for_spectra_resampling_;
    OpenSwath_Scores_Usage su_;

  public:

    /**
    @brief Initialize the scoring object
    */
    void initialize( DoubleReal rt_normalization_factor_,
      int add_up_spectra_, DoubleReal spacing_for_spectra_resampling_,
      OpenSwath_Scores_Usage & su_)
    {
      this->rt_normalization_factor_ = rt_normalization_factor_;
      this->add_up_spectra_ = add_up_spectra_;
      this->spacing_for_spectra_resampling_ = spacing_for_spectra_resampling_;
      this->su_ = su_;
    }

    /** @brief Score a single chromatographic feature.
     *
     * The scores are returned in the OpenSwath_Scores object. Only those
     * scores specified in the OpenSwath_Scores_Usage object are computed.
     *
    */
    void scoreMRMFeature(OpenSwath::ITransitionGroup* itransition_group,
          const PeptideType& pep,
          std::vector<OpenSwath::ISignalToNoisePtr>& signal_noise_estimators,
          OpenSwath::IMRMFeature* imrmfeature,
          const std::vector<TransitionType> & transitions,
          TransformationDescription & trafo,
          OpenSwath::SpectrumAccessPtr swath_map, 
          OpenMS::DIAScoring& diascoring_t,
          OpenSwath_Scores & scores)
    {
        int group_size = boost::numeric_cast<int>(transitions.size());

        // calculate the normalized library intensity (expected value of the intensities)
        std::vector<double> normalized_library_intensity;
        itransition_group->getLibraryIntensities(normalized_library_intensity);
        // std::vector<double> normalized_library_intensity(normalized_library_intensity_.size() );
        OpenSwath::Scoring::normalize_sum(&normalized_library_intensity[0], boost::numeric_cast<int>(normalized_library_intensity.size()));

        // calcxcorr -> for each lag do the correlation, normally use lag 0
        // xcorr_matrix  => correlate chromatogram i with chromatogram j
        bool normalize = true;
        OpenSwath::MRMScoring mrmscore_;
        mrmscore_.initializeXCorrMatrix(imrmfeature, itransition_group, normalize);

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

        // FEATURE : how should we best calculate correlation between library and experiment?
        // FEATURE : spectral angle
        if (su_.use_library_score_)
        {
          mrmscore_.calcLibraryScore(imrmfeature, transitions, scores.library_corr, scores.library_rmsd, scores.library_manhattan, scores.library_dotprod);
        }

        // Retention time score
        if (su_.use_rt_score_)
        {
          // rt score is delta iRT
          double normalized_experimental_rt = trafo.apply(imrmfeature->getRT());
          double rt_score = mrmscore_.calcRTScore(pep, normalized_experimental_rt);

          scores.normalized_experimental_rt = normalized_experimental_rt;
          scores.raw_rt_score = rt_score;
          scores.norm_rt_score = rt_score / rt_normalization_factor_;
        }


        if (su_.use_nr_peaks_score_)
        {
          scores.nr_peaks = group_size;
        }

        if (su_.use_sn_score_)
        {
          scores.sn_ratio = mrmscore_.calcSNScore(imrmfeature, signal_noise_estimators);
          if (scores.sn_ratio < 1) // fix to make sure, that log(sn_score = 0) = -inf does not occur
          {
            scores.log_sn_score = 0;
          }
          else
          {
            scores.log_sn_score = std::log(scores.sn_ratio);
          }
        }

        double quick_lda_dismiss = 0;
        double lda_quick_score = -scores.get_quick_lda_score(scores.library_corr,
            scores.library_rmsd, scores.norm_rt_score, scores.xcorr_coelution_score, scores.xcorr_shape_score, scores.log_sn_score);

        if (lda_quick_score < quick_lda_dismiss)
        {
          // continue;
        }

        if (swath_map->getNrSpectra() > 0)
        {
          calculateSwathScores_(imrmfeature, swath_map, normalized_library_intensity, scores, transitions, diascoring_t, pep);
        }
    }

    /** @brief Score a single chromatographic feature using DIA / SWATH scores.
     *
     * The scores are returned in the OpenSwath_Scores object. 
     *
    */
    void calculateSwathScores_(OpenSwath::IMRMFeature* imrmfeature, OpenSwath::SpectrumAccessPtr swath_map,
        std::vector<double>& normalized_library_intensity, OpenSwath_Scores & scores, const std::vector<TransitionType> & transitions, 
        OpenMS::DIAScoring & diascoring, const PeptideType& pep)
    {
      // parameters
      int by_charge_state = 1; // for which charge states should we check b/y series

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

#ifdef DEBUG_MRMPEAKPICKER
      cout << "added corr isotope_score " << scores.isotope_correlation << endl;
      cout << "added overlap isotope_score " << scores.isotope_overlap << endl;
      cout << "added score massdev_score " << scores.massdev_score << endl;
      cout << "added score weighted massdev_score " << scores.weighted_massdev_score << endl;
      cout << "added score bseries_score " << scores.bseries_score << endl;
      cout << "added score yseries_score " << scores.yseries_score << endl;
#endif
    }

    /// Returns the addition of "nr_spectra_to_add" spectra around the given RT
    OpenSwath::SpectrumPtr getAddedSpectra_(OpenSwath::SpectrumAccessPtr swath_map, double RT, int nr_spectra_to_add)
    {
      std::vector<std::size_t> indices = swath_map->getSpectraByRT(RT, 0.0);
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
          all_spectra.push_back(swath_map->getSpectrumById(closest_idx - i));
          all_spectra.push_back(swath_map->getSpectrumById(closest_idx + i));
        }
        OpenSwath::SpectrumPtr spectrum_ = SpectrumAddition::addUpSpectra(all_spectra, spacing_for_spectra_resampling_, true);
        return spectrum_;
      }
    }

  };

  /**
  @brief The MRMFeatureFinder finds and scores peaks of transitions that coelute.

  It does so using an internal peakpicker for each chromatogram and then
  creating consensus / meta-peaks (MRMFeatures) that contain the information of
  all corresponding chromatograms at the peak-position. It then goes on to
  score those MRMFeatures using different criteria described in the
  MRMScoring class.

  @htmlinclude OpenMS_MRMFeatureFinderScoring.parameters

  */
  class OPENMS_DLLAPI MRMFeatureFinderScoring :
    public DefaultParamHandler,
    public ProgressLogger
  {

public:

    ///Type definitions
    //@{

    // All the filters expect MSSpectrum<PeakT>, thus we give it an "MSSpectrum"
    // but filled with Chromatogram Peaks.
    
    // this is the type in which we store the chromatograms for this analysis
    typedef MSSpectrum<ChromatogramPeak> RichPeakChromatogram;
    typedef OpenSwath::LightTransition TransitionType;
    typedef OpenSwath::LightTargetedExperiment TargetedExpType;
    typedef OpenSwath::LightPeptide PeptideType;
    typedef OpenSwath::LightProtein ProteinType;
    typedef OpenSwath::LightModification ModificationType;
    // a transition group holds the MSSpectra with the Chromatogram peaks from above
    typedef MRMTransitionGroup<MSSpectrum <ChromatogramPeak>, TransitionType> MRMTransitionGroupType; 
    typedef std::map<String, MRMTransitionGroupType> TransitionGroupMapType;
    //@}

    /// Constructor
    MRMFeatureFinderScoring();

    /// Destructor
    ~MRMFeatureFinderScoring();

    /** @brief Pick features in one experiment containing chromatogram
     *
     * Function for for wrapping in Python, only uses OpenMS datastructures 
     * and does not return the map.
     *
    */
    void pickExperiment(MSExperiment<Peak1D> & chromatograms, FeatureMap<Feature>& output, TargetedExperiment& transition_exp_,
                        TransformationDescription trafo, MSExperiment<Peak1D>& swath_map)
    {
      OpenSwath::LightTargetedExperiment transition_exp;
      OpenSwathDataAccessHelper::convertTargetedExp(transition_exp_, transition_exp);
      TransitionGroupMapType transition_group_map;

      OpenSwath::SpectrumAccessPtr chromatogram_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(chromatograms);
      OpenSwath::SpectrumAccessPtr empty_swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);

      pickExperiment(chromatogram_ptr, output, transition_exp, trafo, empty_swath_ptr, transition_group_map);
    }

    /** @brief Pick features in one experiment containing chromatogram
    */
    void pickExperiment(OpenSwath::SpectrumAccessPtr input, FeatureMap<Feature>& output, OpenSwath::LightTargetedExperiment& transition_exp,
                        TransformationDescription trafo, OpenSwath::SpectrumAccessPtr swath_map, TransitionGroupMapType& transition_group_map)
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
      for (TransitionGroupMapType::iterator trgroup_it = transition_group_map.begin(); trgroup_it != transition_group_map.end(); trgroup_it++)
      {
        if (trgroup_it->second.getChromatograms().size() > 0) {counter++; }
      }
      std::cout << "Will analyse " << counter << " peptides with a total of " << transition_exp.getTransitions().size() << " transitions " << std::endl;

      //
      // Step 3
      //
      // Go through all transition groups: first create consensus features, then score them
      Size progress = 0;
      startProgress(0, transition_group_map.size(), "picking peaks");
      for (TransitionGroupMapType::iterator trgroup_it = transition_group_map.begin(); trgroup_it != transition_group_map.end(); trgroup_it++)
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
        scorePeakgroups(trgroup_it->second, trafo, swath_map, output);

      }
      endProgress();

      //output.sortByPosition(); // if the exact same order is needed
      return;
    }

    /** @brief Prepares the internal mappings of peptides and proteins.
    */
    void prepareProteinPeptideMaps_(OpenSwath::LightTargetedExperiment& transition_exp)
    {
      for (Size i = 0; i < transition_exp.getPeptides().size(); i++)
      {
        PeptideRefMap_[transition_exp.getPeptides()[i].id] = &transition_exp.getPeptides()[i];
      }

      for (Size i = 0; i < transition_exp.getProteins().size(); i++)
      {
        ProteinRefMap_[transition_exp.getProteins()[i].id] = &transition_exp.getProteins()[i];
      }
    }

    /** @brief Map the chromatograms to the transitions.
     *
     * Map an input experiment (mzML) and transition list (TraML) onto each other
     * when they share identifiers, e.g. if the transition id is the same as the
     * chromatogram native id.
    */
    void mapExperimentToTransitionList(OpenSwath::SpectrumAccessPtr input, OpenSwath::LightTargetedExperiment& transition_exp,
                                       TransitionGroupMapType& transition_group_map, TransformationDescription trafo, double rt_extraction_window);

    /** @brief Set the flag for strict mapping
    */
    void setStrictFlag(bool f)
    {
      strict_ = f;
    }

    /** @brief Score all peak groups of a transition group
     *
     * Iterate through all features found along the chromatograms of the
     * transition group and score each one individually.
     *
    */
    void scorePeakgroups(MRMTransitionGroupType& transition_group, TransformationDescription & trafo,
                         OpenSwath::SpectrumAccessPtr swath_map, FeatureMap<Feature>& output)
    {
      typedef typename MRMTransitionGroupType::PeakType PeakT;
      std::vector<OpenSwath::ISignalToNoisePtr> signal_noise_estimators;
      std::vector<MRMFeature> feature_list;

      DoubleReal sn_win_len_ = (DoubleReal)param_.getValue("TransitionGroupPicker:sn_win_len");
      DoubleReal sn_bin_count_ = (DoubleReal)param_.getValue("TransitionGroupPicker:sn_bin_count");
      for (Size k = 0; k < transition_group.getChromatograms().size(); k++)
      {
        OpenSwath::ISignalToNoisePtr snptr(new OpenMS::SignalToNoiseOpenMS< PeakT >(transition_group.getChromatograms()[k], sn_win_len_, sn_bin_count_));
        signal_noise_estimators.push_back(snptr);
      }

      const PeptideType* pep = PeptideRefMap_[transition_group.getTransitionGroupID()];
      const ProteinType* prot = ProteinRefMap_[pep->protein_ref];

      // get the expected rt value for this peptide
      double expected_rt = pep->rt;
      TransformationDescription newtr = trafo;
      newtr.invert();
      expected_rt = newtr.apply(expected_rt);

      SwathScorer scorer;
      scorer.initialize(rt_normalization_factor_, add_up_spectra_, spacing_for_spectra_resampling_, su_);

      // Go through all peak groups (found MRM features) and score them
      for (std::vector<MRMFeature>::iterator mrmfeature = transition_group.getFeaturesMuteable().begin();
           mrmfeature != transition_group.getFeaturesMuteable().end(); mrmfeature++)
      {
        OpenSwath::IMRMFeature* imrmfeature;
        imrmfeature = new MRMFeatureOpenMS(*mrmfeature);
        OpenSwath::ITransitionGroup* itransition_group;
        itransition_group = new TransitionGroupOpenMS<MRMTransitionGroupType::SpectrumT,
                          MRMTransitionGroupType::TransitionT>(transition_group);

        LOG_DEBUG << "000000000000000000000000000000000000000000000000000000000000000000000000000 " << std::endl;
        LOG_DEBUG << "scoring feature " << (*mrmfeature) << " == " << mrmfeature->getMetaValue("PeptideRef") <<
        "[ expected RT " << PeptideRefMap_[mrmfeature->getMetaValue("PeptideRef")]->rt << " / " << expected_rt << " ]" <<
        " with " << transition_group.size()  << " nr transitions and nr chromats " << transition_group.getChromatograms().size() << std::endl;

        int group_size = boost::numeric_cast<int>(transition_group.size());
        if (group_size == 0)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
            "Error: Transition group " + transition_group.getTransitionGroupID() + " has no chromatograms.");
        }
        if (group_size < 2)
        {
          std::cerr << "Error: transition group " << transition_group.getTransitionGroupID() 
            << " has less than 2 chromatograms. It has " << group_size << std::endl;
          continue;
        }

        ///////////////////////////////////
        // call the scoring
        ///////////////////////////////////

        OpenSwath_Scores scores;
        scorer.scoreMRMFeature(itransition_group, *pep,
          signal_noise_estimators, imrmfeature, transition_group.getTransitions(),
          trafo, swath_map, diascoring_, scores);

        if (su_.use_coelution_score_) { 
          mrmfeature->addScore("var_xcorr_coelution", scores.xcorr_coelution_score);
          mrmfeature->addScore("var_xcorr_coelution_weighted ", scores.weighted_coelution_score); }
        if (su_.use_shape_score_) { 
          mrmfeature->addScore("var_xcorr_shape", scores.xcorr_shape_score);
          mrmfeature->addScore("var_xcorr_shape_weighted", scores.weighted_xcorr_shape); }
        if (su_.use_library_score_) { 
          mrmfeature->addScore("var_library_corr", scores.library_corr);
          mrmfeature->addScore("var_library_rmsd", scores.library_rmsd);
          mrmfeature->addScore("var_library_manhattan", scores.library_manhattan);
          mrmfeature->addScore("var_library_dotprod", scores.library_dotprod); }
        if (su_.use_rt_score_) { 
          mrmfeature->addScore("delta_rt", mrmfeature->getRT() - expected_rt);
          mrmfeature->addScore("assay_rt", expected_rt);
          mrmfeature->addScore("norm_RT", scores.normalized_experimental_rt);
          mrmfeature->addScore("rt_score", scores.raw_rt_score);
          mrmfeature->addScore("var_norm_rt_score", scores.norm_rt_score); }
        // TODO do we really want these intensity scores ?
        if (su_.use_intensity_score_) { mrmfeature->addScore("var_intensity_score", mrmfeature->getIntensity() / (double)mrmfeature->getMetaValue("total_xic")); }
        if (su_.use_total_xic_score_) { mrmfeature->addScore("total_xic", (double)mrmfeature->getMetaValue("total_xic")); }
        if (su_.use_nr_peaks_score_) { mrmfeature->addScore("nr_peaks", scores.nr_peaks); }
        if (su_.use_sn_score_) { mrmfeature->addScore("sn_ratio", scores.sn_ratio); mrmfeature->addScore("var_log_sn_score", scores.log_sn_score); }
        // TODO get it working with imrmfeature
        if (su_.use_elution_model_score_) { 
          scores.elution_model_fit_score = emgscoring_.calcElutionFitScore((*mrmfeature), transition_group);
          mrmfeature->addScore("var_elution_model_fit_score", scores.elution_model_fit_score); }

        double xx_lda_prescore = -scores.calculate_lda_prescore(scores);
        bool swath_present = (swath_map->getNrSpectra() > 0);
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
        if (swath_present)
        {
          mrmfeature->addScore("var_isotope_correlation_score", scores.isotope_correlation);
          mrmfeature->addScore("var_isotope_overlap_score", scores.isotope_overlap);
          mrmfeature->addScore("var_massdev_score", scores.massdev_score);
          mrmfeature->addScore("var_massdev_score_weighted", scores.weighted_massdev_score);
          mrmfeature->addScore("var_bseries_score", scores.bseries_score);
          mrmfeature->addScore("var_yseries_score", scores.yseries_score);
          mrmfeature->addScore("var_dotprod_score", scores.dotprod_score_dia);
          mrmfeature->addScore("var_manhatt_score", scores.manhatt_score_dia);

          double xx_swath_prescore = -scores.calculate_swath_lda_prescore(scores);
          mrmfeature->addScore("main_var_xx_swath_prelim_score", xx_swath_prescore);
          mrmfeature->setOverallQuality(xx_swath_prescore);
        }

        ///////////////////////////////////////////////////////////////////////////
        // add the peptide hit information to the feature
        ///////////////////////////////////////////////////////////////////////////
        PeptideIdentification pep_id_ = PeptideIdentification();
        PeptideHit pep_hit_ = PeptideHit();

        if (pep->getChargeState() != -1)
        {
          pep_hit_.setCharge(pep->getChargeState());
        }
        pep_hit_.setScore(xx_lda_prescore);
        if (swath_present)
        {
          pep_hit_.setScore(mrmfeature->getScore("xx_swath_prelim_score"));
        }
        pep_hit_.setSequence((String)pep->sequence);
        pep_hit_.addProteinAccession(prot->id);
        pep_id_.insertHit(pep_hit_);
        pep_id_.setIdentifier(run_identifier);

        mrmfeature->getPeptideIdentifications().push_back(pep_id_);
        mrmfeature->ensureUniqueId();
        mrmfeature->setMetaValue("PrecursorMZ", transition_group.getTransitions()[0].getPrecursorMZ());
        mrmfeature->setSubordinates(mrmfeature->getFeatures()); // add all the subfeatures as subordinates
        double total_intensity = 0, total_peak_apices = 0;
        for (std::vector<Feature>::iterator sub_it = mrmfeature->getSubordinates().begin(); sub_it != mrmfeature->getSubordinates().end(); sub_it++)
        {
          if (!write_convex_hull_) {sub_it->getConvexHulls().clear(); }
          sub_it->ensureUniqueId();
          if (sub_it->getMZ() > quantification_cutoff_)
          {
            total_intensity += sub_it->getIntensity();
            total_peak_apices += (DoubleReal)sub_it->getMetaValue("peak_apex_int");
          }
        }
        // overwrite the reported intensities with those above the m/z cutoff
        mrmfeature->setIntensity(total_intensity);
        mrmfeature->setMetaValue("peak_apices_sum", total_peak_apices);
        feature_list.push_back((*mrmfeature));

        delete imrmfeature;
        delete itransition_group;
      }

      // Order by quality
      std::sort(feature_list.begin(), feature_list.end(), OpenMS::Feature::OverallQualityLess());
      std::reverse(feature_list.begin(), feature_list.end());

      for (Size i = 0; i < feature_list.size(); i++)
      {
        if (stop_report_after_feature_ >= 0 && i >= (Size)stop_report_after_feature_) {break; }
        output.push_back(feature_list[i]);
      }
    }

private:

    /// Synchronize members with param class
    void updateMembers_();

    // Variables
    DoubleReal rt_extraction_window_;
    DoubleReal quantification_cutoff_;


    int stop_report_after_feature_;

    // bool do_local_fdr_;
    bool write_convex_hull_;
    bool strict_;

    // TODO
    DoubleReal rt_normalization_factor_;
    int add_up_spectra_;
    DoubleReal spacing_for_spectra_resampling_;
    // TODO

    OpenSwath_Scores_Usage su_;

    std::map<OpenMS::String, const PeptideType*> PeptideRefMap_;
    std::map<OpenMS::String, const ProteinType*> ProteinRefMap_;

    OpenMS::DIAScoring diascoring_;
    OpenMS::EmgScoring emgscoring_;
  };
}

#undef run_identifier
#endif
