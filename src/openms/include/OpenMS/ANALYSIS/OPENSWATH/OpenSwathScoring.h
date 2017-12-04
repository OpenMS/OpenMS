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

#ifndef OPENMS_ANALYSIS_OPENSWATH_OPENSWATHSCORING_H
#define OPENMS_ANALYSIS_OPENSWATH_OPENSWATHSCORING_H

// data access
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/DataStructures.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ITransition.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/SwathMap.h>

// scoring
#include <OpenMS/ANALYSIS/OPENSWATH/DIAScoring.h>

#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

namespace OpenMS
{

  /** @brief A structure to store which scores should be used by the Algorithm
   *
   * This can be used to turn on/off individual scores.
  */
  struct OPENMS_DLLAPI OpenSwath_Scores_Usage
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
    bool use_dia_scores_;
    bool use_sonar_scores;
    bool use_ms1_correlation;
    bool use_ms1_fullscan;
    bool use_uis_scores;
    
    OpenSwath_Scores_Usage() :
      use_coelution_score_(true),
      use_shape_score_(true),
      use_rt_score_(true),
      use_library_score_(true),
      use_elution_model_score_(true),
      use_intensity_score_(true),
      use_total_xic_score_(true),
      use_nr_peaks_score_(true),
      use_sn_score_(true),
      use_dia_scores_(true),
      use_sonar_scores(true),
      use_ms1_correlation(true),
      use_ms1_fullscan(true),
      use_uis_scores(true)
    {}

  };

  /** @brief A structure to hold the different scores computed by OpenSWATH
   *
   * This struct is used to store the individual OpenSWATH (sub-)scores. It
   * also allows to compute some preliminary quality score for a feature by
   * using a predefined combination of the individual scores determined using
   * LDA.
   *
  */
  struct OPENMS_DLLAPI OpenSwath_Scores
  {
    double elution_model_fit_score;
    double library_corr;
    double library_norm_manhattan;
    double library_rootmeansquare;
    double library_sangle;
    double norm_rt_score;
    double isotope_correlation;
    std::string ind_isotope_correlation;
    double isotope_overlap;
    std::string ind_isotope_overlap;
    double massdev_score;
    std::string ind_massdev_score;
    double xcorr_coelution_score;
    std::string ind_xcorr_coelution_score;
    double xcorr_shape_score;
    std::string ind_xcorr_shape_score;
    double yseries_score;
    double bseries_score;
    double log_sn_score;
    std::string ind_log_sn_score;
    int ind_num_transitions;
    std::string ind_transition_names;
    std::string ind_area_intensity;
    std::string ind_apex_intensity;
    std::string ind_log_intensity;

    double weighted_coelution_score;
    double weighted_xcorr_shape;
    double weighted_massdev_score;
   
    double xcorr_ms1_coelution_score;
    double xcorr_ms1_shape_score;
    double ms1_ppm_score;
    double ms1_isotope_correlation;
    double ms1_isotope_overlap;

    double sonar_sn;
    double sonar_diff;
    double sonar_trend;
    double sonar_rsq;
    double sonar_shape;
    double sonar_lag;

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
      library_norm_manhattan(0),
      library_rootmeansquare(0),
      library_sangle(0),
      norm_rt_score(0),
      isotope_correlation(0),
      ind_isotope_correlation(""),
      isotope_overlap(0),
      ind_isotope_overlap(""),
      massdev_score(0),
      ind_massdev_score(""),
      xcorr_coelution_score(0),
      ind_xcorr_coelution_score(""),
      xcorr_shape_score(0),
      ind_xcorr_shape_score(""),
      yseries_score(0),
      bseries_score(0),
      log_sn_score(0),
      ind_log_sn_score(""),
      ind_num_transitions(0),
      ind_transition_names(""),
      ind_log_intensity(""),
      weighted_coelution_score(0),
      weighted_xcorr_shape(0),
      weighted_massdev_score(0),
      xcorr_ms1_coelution_score(0),
      xcorr_ms1_shape_score(0),
      ms1_ppm_score(0),
      ms1_isotope_correlation(0),
      ms1_isotope_overlap(0),
      sonar_sn(0),
      sonar_diff(0),
      sonar_trend(0),
      sonar_rsq(0),
      sonar_shape(0),
      sonar_lag(0),
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


    double get_quick_lda_score(double library_corr_, double library_norm_manhattan_, double norm_rt_score_, double xcorr_coelution_score_,
                               double xcorr_shape_score_, double log_sn_score_)
    {
      // some scores based on manual evaluation of 80 chromatograms
      // quick LDA average model on 100 2 x Crossvalidated runs (0.85 TPR/0.17 FDR)
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
        library_corr_                    * -0.5319046 +
        library_norm_manhattan_          *  2.1643962 +
        norm_rt_score_                   *  8.0353047 +
        xcorr_coelution_score_           *  0.1458914 +
        xcorr_shape_score_               * -1.6901925 +
        log_sn_score_                    * -0.8002824;
      return lda_quick_score;
    }

    double calculate_lda_prescore(OpenSwath_Scores scores)
    {

      // LDA average model on 100 2 x Crossvalidated runs (0.91 TPR/0.20 FDR)
      /*
      double xx_old_lda_prescore =
      intensity_score       * -2.296679          +
      library_corr          * -0.1223876         +
      library_norm_manhattan*  2.013638          +
      nr_peaks_score        *  0.01683357        +
      rt_score              *  0.00143999        +
      sn_score              * -0.1619762         +
      total_xic_score       *  0.00000003697898  +
      xcorr_coelution_score *  0.05909583        +
      xcorr_shape_score     * -0.4699841;

      // NOTE this score means "better" if it is more negative!
      */

      return scores.library_corr                     * -0.34664267 +
             scores.library_norm_manhattan           *  2.98700722 +
             scores.norm_rt_score                    *  7.05496384 +
             scores.xcorr_coelution_score            *  0.09445371 +
             scores.xcorr_shape_score                * -5.71823862 +
             scores.log_sn_score                     * -0.72989582 +
             scores.elution_model_fit_score          *  1.88443209;
    }

    double calculate_swath_lda_prescore(OpenSwath_Scores scores)
    {

      // Swath - LDA average model on 100 2 x Crossvalidated runs (0.76 TPR/0.20 FDR) [without elution model]
      /*
      double xx_old_swath_prescore =
      intensity_score              * -3.148838e+00  +
      library_corr                 * -7.562403e-02  +
      library_norm_manhattan       *  1.786286e+00  +
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
             scores.library_norm_manhattan    *  2.47298914 +
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
 var_library_norm_manhattan  -6.44998534845163
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
var_library_norm_manhattan    -13.1295053007338
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
var_library_norm_manhattan    -13.2905041886848
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
var_library_norm_manhattan    -5.20235596662978
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
var_library_norm_manhattan    -6.08160902837145
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
var_library_norm_manhattan    -8.3847526088391
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

  /** @brief A class that calls the scoring routines
   *
   * Use this class to invoke the individual OpenSWATH scoring routines.
   * 
  */
  class OPENMS_DLLAPI OpenSwathScoring 
  {
    typedef OpenSwath::LightCompound CompoundType;
    typedef OpenSwath::LightTransition TransitionType;

    double rt_normalization_factor_;
    int add_up_spectra_;
    double spacing_for_spectra_resampling_;
    OpenSwath_Scores_Usage su_;

  public:

    /// Constructor
    OpenSwathScoring();

    /// Destructor
    ~OpenSwathScoring();

    /** @brief Initialize the scoring object
     *
     * Sets the parameters for the scoring.
     *
     * @param rt_normalization_factor Specifies the range of the normalized retention time space
     * @param add_up_spectra How many spectra to add up (default 1)
     * @param spacing_for_spectra_resampling Spacing factor for spectra addition
     * @param su Which scores to actually compute
     *
    */
    void initialize(double rt_normalization_factor,
      int add_up_spectra, double spacing_for_spectra_resampling,
      OpenSwath_Scores_Usage & su);

    /** @brief Score a single peakgroup in a chromatogram using only chromatographic properties.
     *
     * This function only uses the chromatographic properties (coelution,
     * signal to noise, etc.) of a peakgroup in a chromatogram to compute
     * scores. If more information is available, also consider using the
     * library based scoring and the full-spectrum based scoring.
     *
     * The scores are returned in the OpenSwath_Scores object. Only those
     * scores specified in the OpenSwath_Scores_Usage object are computed.
     *
     * @param imrmfeature The feature to be scored
     * @param native_ids The list of native ids (giving a canonical ordering of the transitions)
     * @param normalized_library_intensity The weights to be used for each transition (e.g. normalized library intensities)
     * @param signal_noise_estimators The signal-to-noise estimators for each transition
     * @param scores The object to store the result
     *
    */
    void calculateChromatographicScores(
          OpenSwath::IMRMFeature* imrmfeature,
          const std::vector<std::string>& native_ids,
          const std::string& precursor_chrom_id,
          const std::vector<double>& normalized_library_intensity,
          std::vector<OpenSwath::ISignalToNoisePtr>& signal_noise_estimators,
          OpenSwath_Scores & scores);

    /** @brief Score identification transitions against detection transitions of a single peakgroup 
     * in a chromatogram using only chromatographic properties.
     *
     * This function only uses the chromatographic properties (coelution,
     * signal to noise, etc.) of a peakgroup in a chromatogram to compute
     * scores. The scores are computed by scoring identification against detection
     * transitions.
     *
     * The scores are returned in the OpenSwath_Scores object. Only those
     * scores specified in the OpenSwath_Scores_Usage object are computed.
     *
     * @param imrmfeature The feature to be scored
     * @param native_ids_identification The list of identification native ids (giving a canonical ordering of the transitions)
     * @param native_ids_detection The list of detection native ids (giving a canonical ordering of the transitions)
     * @param signal_noise_estimators The signal-to-noise estimators for each transition
     * @param scores The object to store the result
     *
    */
    void calculateChromatographicIdScores(
          OpenSwath::IMRMFeature* imrmfeature,
          const std::vector<std::string>& native_ids_identification,
          const std::vector<std::string>& native_ids_detection,
          std::vector<OpenSwath::ISignalToNoisePtr>& signal_noise_estimators,
          OpenSwath_Scores & scores);

    /** @brief Score a single chromatographic feature against a spectral library
     *
     * The spectral library is provided in a set of transition objects and a
     * peptide object. Both contain information about the expected elution time
     * on the chromatography and the relative intensity of the transitions.
     *
     * The scores are returned in the OpenSwath_Scores object. 
     *
     * @param imrmfeature The feature to be scored
     * @param transitions The library transition to score the feature against
     * @param pep The peptide corresponding to the library transitions
     * @param normalized_feature_rt The retention time of the feature in normalized space
     * @param scores The object to store the result
     *
    */
    void calculateLibraryScores(
          OpenSwath::IMRMFeature* imrmfeature,
          const std::vector<TransitionType> & transitions,
          const CompoundType& compound,
          const double normalized_feature_rt,
          OpenSwath_Scores & scores);

    /** @brief Score a single chromatographic feature using DIA / SWATH scores.
     *
     * The scores are returned in the OpenSwath_Scores object. 
     *
     * @param imrmfeature The feature to be scored
     * @param transitions The library transition to score the feature against
     * @param swath_maps The SWATH-MS (DIA) maps from which to retrieve full MS/MS spectra at the chromatographic peak apices
     * @param ms1_map The corresponding MS1 (precursor ion map) from which the precursor spectra can be retrieved (optional, may be NULL)
     * @param diascoring DIA Scoring object to use for scoring
     * @param pep The peptide corresponding to the library transitions
     * @param scores The object to store the result
     *
    */
    void calculateDIAScores(OpenSwath::IMRMFeature* imrmfeature, 
        const std::vector<TransitionType> & transitions,
        std::vector<OpenSwath::SwathMap> swath_maps,
        OpenSwath::SpectrumAccessPtr ms1_map,
        OpenMS::DIAScoring & diascoring,
        const CompoundType& compound,
        OpenSwath_Scores & scores);

    /** @brief Score a single chromatographic feature using the precursor map.
     *
     * The scores are returned in the OpenSwath_Scores object. 
     *
     * @param ms1_map The MS1 (precursor ion map) from which the precursor spectra can be retrieved
     * @param diascoring DIA Scoring object to use for scoring
     * @param precursor_mz The m/z ratio of the precursor
     * @param rt The compound retention time
     * @param scores The object to store the result
     *
    */
    void calculatePrecursorDIAScores(OpenSwath::SpectrumAccessPtr ms1_map, 
                                     OpenMS::DIAScoring & diascoring, 
                                     double precursor_mz, 
                                     double rt, 
                                     const CompoundType& compound, 
                                     OpenSwath_Scores & scores);

    /** @brief Score a single chromatographic feature using DIA / SWATH scores.
     *
     * The scores are returned in the OpenSwath_Scores object. 
     *
     * @param imrmfeature The feature to be scored
     * @param transitions The library transition to score the feature against
     * @param swath_maps The SWATH-MS (DIA) maps from which to retrieve full MS/MS spectra at the chromatographic peak apices
     * @param diascoring DIA Scoring object to use for scoring
     * @param scores The object to store the result
     *
    */
    void calculateDIAIdScores(OpenSwath::IMRMFeature* imrmfeature,
        const TransitionType & transition,
        std::vector<OpenSwath::SwathMap> swath_maps,
        OpenMS::DIAScoring & diascoring,
        OpenSwath_Scores & scores);

    /** @brief Computing the normalized library intensities from the transition objects
     *
     * The intensities are normalized such that the sum to one.
     *
     * @param transitions The library transition to score the feature against
     * @param normalized_library_intensity The resulting normalized library intensities
     *
    */
    void getNormalized_library_intensities_(const std::vector<TransitionType> & transitions,
        std::vector<double>& normalized_library_intensity);

    /** @brief Returns an averaged spectrum
     *
     * This function will sum up (add) the intensities of multiple spectra
     * around the given retention time and return an "averaged" spectrum which
     * may contain less noise.
     *
     * @param swath_map The map containing the spectra
     * @param RT The target retention time
     * @param nr_spectra_to_add How many spectra to add up
     *
    */
    OpenSwath::SpectrumPtr getAddedSpectra_(OpenSwath::SpectrumAccessPtr swath_map, 
        double RT, int nr_spectra_to_add);

    /** @brief Returns an averaged spectrum
     *
     * This function will sum up (add) the intensities of multiple spectra from
     * multiple swath maps (assuming these are SONAR maps of shifted precursor
     * isolation windows) around the given retention time and return an
     * "averaged" spectrum which may contain less noise.
     *
     * @param swath_maps The maps containing the spectra
     * @param RT The target retention time
     * @param nr_spectra_to_add How many spectra to add up
     *
    */
    OpenSwath::SpectrumPtr getAddedSpectra_(std::vector<OpenSwath::SwathMap> swath_maps,
                                            double RT, int nr_spectra_to_add);

  };
}

#endif // OPENMS_ANALYSIS_OPENSWATH_OPENSWATHSCORING_H
