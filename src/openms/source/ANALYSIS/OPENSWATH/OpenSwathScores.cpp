// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathScores.h>

namespace OpenMS
{

    double OpenSwath_Scores::get_quick_lda_score(double library_corr_,
                                                 double library_norm_manhattan_,
                                                 double norm_rt_score_,
                                                 double xcorr_coelution_score_,
                                                 double xcorr_shape_score_,
                                                 double log_sn_score_) const
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

    double OpenSwath_Scores::calculate_lda_prescore(const OpenSwath_Scores& scores) const
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

    double OpenSwath_Scores::calculate_lda_single_transition(const OpenSwath_Scores& scores) const
    {
      // Manually derived scoring model for single transition peakgroups
      return scores.norm_rt_score                    *  7.05496384 +
             scores.log_sn_score                     * -0.72989582 +
             scores.elution_model_fit_score          *  -1.08443209;
    }

    double OpenSwath_Scores::calculate_swath_lda_prescore(const OpenSwath_Scores& scores) const
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
}

