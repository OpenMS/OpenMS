// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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

#ifndef OPENMS_ANALYSIS_OPENSWATH_MRMPEAKPICKER_H
#define OPENMS_ANALYSIS_OPENSWATH_MRMPEAKPICKER_H

#define run_identifier "unique_run_identifier"
#define USE_SP_INTERFACE

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/FORMAT/TraMLFile.h>

#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>
#include <OpenMS/KERNEL/MRMTransitionGroup.h>
#include <OpenMS/KERNEL/MRMFeature.h>

// peak picking
#include <OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h>

// data access
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/DataStructures.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/MRMFeatureAccessOpenMS.h>

// scoring
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/Scoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/MRMScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DIAScoring.h>

#ifdef _OPENMP
#include <omp.h>
#endif

bool SortDoubleDoublePairFirst(const std::pair<double, double> & left, const std::pair<double, double> & right);

namespace OpenMS
{

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
    double log_sn_score;

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
    }

  };

  /**
  @brief The MRMFeatureFinder finds and scores peaks of transitions that coelute.

  It does so using an internal peakpicker for each chromatogram and then
  creating consensus / meta-peaks (MRMFeatures) that contain the information of
  all corresponding chromatograms at the peak-position. It then goes on to
  score those MRMFeatures using different criteria described in the
  MRMScoring class.

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
    typedef MSSpectrum<ChromatogramPeak> RichPeakChromatogram; // this is the type in which we store the chromatograms for this analysis
    typedef OpenSwath::LightTransition TransitionType;
    typedef OpenSwath::LightTargetedExperiment TargetedExpType;
    typedef OpenSwath::LightPeptide PeptideType;
    typedef OpenSwath::LightProtein ProteinType;
    typedef OpenSwath::LightModification ModificationType;
    typedef MRMTransitionGroup<MSSpectrum, ChromatogramPeak, TransitionType> MRMTransitionGroupType; // a transition group holds the MSSpectra with the Chromatogram peaks from above
    typedef std::map<String, MRMTransitionGroupType> TransitionGroupMapType;
    //@}

    /// Constructor
    MRMFeatureFinderScoring();

    /// Destructor
    ~MRMFeatureFinderScoring();

    // pick features in one experiment containing chromatograms
    void pickExperiment(OpenSwath::SpectrumAccessPtr input, FeatureMap<Feature> & output, OpenSwath::LightTargetedExperiment & transition_exp,
                        TransformationDescription & trafo, OpenSwath::SpectrumAccessPtr swath_map, TransitionGroupMapType & transition_group_map)
    {
      handle_params();

      //
      // Step 1
      //
      // Store the peptide retention times in an intermediate map
      PeptideRTMap.clear();
      for (Size i = 0; i < transition_exp.getPeptides().size(); i++)
      {
        PeptideType pep = transition_exp.getPeptides()[i];
        PeptideRTMap[pep.id] = pep.rt;
        PeptideRefMap[pep.id] = &transition_exp.getPeptides()[i];
      }

      // Store the proteins from the input in the output feature map
      std::vector<ProteinHit> protein_hits;
      for (Size i = 0; i < transition_exp.getProteins().size(); i++)
      {
        const ProteinType & prot = transition_exp.getProteins()[i];
        ProteinRefMap[transition_exp.getProteins()[i].id] = &transition_exp.getProteins()[i];
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

      // Local FDR is a comparison whether there is any other assay that could
      // explain the current peak group. Thus we compare this MRMFeature to all
      // the other assays we have in the database and calculate the score.
      if (do_local_fdr_)
      {
        String tr_library = param_.getValue("local_fdr_lib");
        //loadLibraryForLocalFDR(tr_library);
      }

      //
      // Step 3
      //
      // Go through all transition groups: first create consensus features, then score them
      Size progress = 0;
      startProgress(0, transition_group_map.size(), "picking peaks");
      for (TransitionGroupMapType::iterator trgroup_it = transition_group_map.begin(); trgroup_it != transition_group_map.end(); trgroup_it++)
      {

        setProgress(++progress);
        MRMTransitionGroupType & transition_group = trgroup_it->second;
        if (transition_group.getChromatograms().size() == 0 || transition_group.getTransitions().size() == 0)
        {
          continue;
        }

        MRMTransitionGroupPicker trgroup_picker;
        Param trgroup_picker_parameters = trgroup_picker.getParameters();
        trgroup_picker_parameters.setValue("sgolay_frame_length", sgolay_frame_length_);
        trgroup_picker_parameters.setValue("sgolay_polynomial_order", sgolay_polynomial_order_);
        trgroup_picker_parameters.setValue("gauss_width", gauss_width_);
        trgroup_picker_parameters.setValue("peak_width", peak_width_);
        trgroup_picker_parameters.setValue("signal_to_noise", signal_to_noise_);
        trgroup_picker_parameters.setValue("sn_win_len", sn_win_len_);
        trgroup_picker_parameters.setValue("sn_bin_count", sn_bin_count_);
        trgroup_picker_parameters.setValue("stop_after_feature", stop_after_feature_);
        trgroup_picker_parameters.setValue("stop_report_after_feature", stop_report_after_feature_);
        trgroup_picker_parameters.setValue("stop_after_intensity_ratio", stop_after_intensity_ratio_);
        trgroup_picker.setParameters(trgroup_picker_parameters);
        trgroup_picker.handle_params();

        trgroup_picker.pickTransitionGroup(transition_group);

        scorePeakgroups(trgroup_it->second, trafo, swath_map, output);

      }
      endProgress();

      //output.sortByPosition(); // if the exact same order is needed
      return;
    }

    // Map an input experiment (mzML) and transition list (TraML) onto each other
    // when they share identifiers, e.g. if the transition id is the same as the
    // chromatogram native id.
    void mapExperimentToTransitionList(OpenSwath::SpectrumAccessPtr input, OpenSwath::LightTargetedExperiment & transition_exp,
                                       TransitionGroupMapType & transition_group_map, TransformationDescription trafo, double rt_extraction_window);

    void setStrictFlag(bool f)
    {
      strict = f;
    }

private:

    /// Score all peak groups
    template <template <typename> class SpectrumT, typename PeakT, typename TransitionT>
    void scorePeakgroups(MRMTransitionGroup<SpectrumT, PeakT, TransitionT> & transition_group, TransformationDescription & trafo,
                         OpenSwath::SpectrumAccessPtr  swath_map, FeatureMap<Feature> & output)
    {
      //std::vector<SignalToNoiseEstimatorMedian<RichPeakChromatogram> > signal_noise_estimators;
      std::vector<OpenSwath::ISignalToNoisePtr> signal_noise_estimators;
      std::vector<MRMFeature> feature_list;

      for (Size k = 0; k < transition_group.getChromatograms().size(); k++)
      {
        OpenSwath::ISignalToNoisePtr snptr(new OpenMS::SignalToNoiseOpenMS<PeakT>(transition_group.getChromatograms()[k], sn_win_len_, sn_bin_count_));
        signal_noise_estimators.push_back(snptr);
      }

      // get the expected rt value for this peptide
      double expected_rt = PeptideRTMap[transition_group.getTransitionGroupID()];
      TransformationDescription newtr = trafo;
      newtr.invert();
      expected_rt = newtr.apply(expected_rt);

      // Go through all peak groups (found MRM features) and score them
      for (std::vector<MRMFeature>::iterator mrmfeature = transition_group.getFeaturesMuteable().begin();
           mrmfeature != transition_group.getFeaturesMuteable().end(); mrmfeature++)
      {

        OpenSwath::IMRMFeature * imrmfeature;
        imrmfeature = new MRMFeatureOpenMS(*mrmfeature);

        OpenSwath::ITransitionGroup * itransition_group;
        itransition_group = new TransitionGroupOpenMS<SpectrumT, PeakT, TransitionT>(transition_group);

#ifdef DEBUG_MRMPEAKPICKER
        std::cout << "000000000000000000000000000000000000000000000000000000000000000000000000000 " << std::endl;
        std::cout << "scoring feature " << (*mrmfeature) << " == " << mrmfeature->getMetaValue("PeptideRef") <<
        "[ expected RT " << PeptideRTMap[mrmfeature->getMetaValue("PeptideRef")] << " / " << expected_rt << " ]" <<
        " with " << transition_group.size()  << " nr transitions and nr chromats " << transition_group.getChromatograms().size() << std::endl;
#endif

        int group_size = transition_group.size();
        if (group_size == 0)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                           "Error: Transition group " + transition_group.getTransitionGroupID() + " has no chromatograms.");
        }
        if (group_size < 2)
        {
          std::cerr << "Error: transition group " << transition_group.getTransitionGroupID() << " has less than 2 chromatograms. It has " << group_size << std::endl;
          continue;
          //throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Error: transition group " + transition_group.getTransitionGroupID() + " has less than 2 chromatograms.");
        }

        // calculate the normalized library intensity (expected value of the intensities)
        std::vector<double> normalized_library_intensity;
        transition_group.getLibraryIntensity(normalized_library_intensity);
        Scoring::normalize_sum(&normalized_library_intensity[0], normalized_library_intensity.size());

        // calcxcorr -> for each lag do the correlation, normally use lag 0
        // xcorr_matrix  => correlate chromatogram i with chromatogram j
        bool normalize = true;
        mrmscore.initializeXCorrMatrix(imrmfeature, itransition_group, normalize);

        // XCorr score (coelution)
        double xcorr_coelution_score = 0;
        if (use_coelution_score_)
        {
          xcorr_coelution_score = mrmscore.calcXcorrCoelutionScore();
          mrmfeature->addScore("var_xcorr_coelution", xcorr_coelution_score);
        }

        double weighted_coelution_score = 0;
        if (use_coelution_score_)
        {
          weighted_coelution_score = mrmscore.calcXcorrCoelutionScore_weighted(normalized_library_intensity);
          mrmfeature->addScore("var_xcorr_coelution_weighted ", weighted_coelution_score);
        }

        // XCorr score (shape)
        // mean over the intensities at the max of the crosscorrelation
        // FEATURE : weigh by the intensity as done by mQuest
        // FEATURE : normalize with the intensity at the peak group apex?
        double xcorr_shape_score = 0;
        if (use_shape_score_)
        {
          xcorr_shape_score = mrmscore.calcXcorrShape_score();
          mrmfeature->addScore("var_xcorr_shape", xcorr_shape_score);
        }

        double weighted_xcorr_shape = 0;
        if (use_shape_score_)
        {
          weighted_xcorr_shape = mrmscore.calcXcorrShape_score_weighted(normalized_library_intensity);
          mrmfeature->addScore("var_xcorr_shape_weighted", weighted_xcorr_shape);
        }

        // FEATURE : how should we best calculate correlation between library and experiment?
        // FEATURE : spectral angle
        double library_corr = 0, library_rmsd = 0;
        double library_manhattan, library_dotprod;

        if (use_library_score_)
        {
          mrmscore.calcLibraryScore(imrmfeature, transition_group.getTransitions(), library_corr, library_rmsd, library_manhattan, library_dotprod);
          mrmfeature->addScore("var_library_corr", library_corr);
          mrmfeature->addScore("var_library_rmsd", library_rmsd);
          mrmfeature->addScore("var_library_manhattan", library_manhattan); // new score
          mrmfeature->addScore("var_library_dotprod", library_dotprod); // new score
        }

        // Retention time score
        double rt_score = 0, norm_rt_score = 0;
        if (use_rt_score_)
        {
          // get the id, then get the expected and the experimental retention time
          String native_id = transition_group.getChromatograms()[0].getNativeID();
          TransitionType tr = transition_group.getTransition(native_id);
          const PeptideType * pep = PeptideRefMap[tr.getPeptideRef()];
          double experimental_rt = mrmfeature->getFeature(native_id).getRT();
          double normalized_experimental_rt = trafo.apply(experimental_rt);
          // rt score is delta iRT
          rt_score = mrmscore.calcRTScore(*pep, normalized_experimental_rt);
          norm_rt_score = rt_score / rt_normalization_factor_;
          mrmfeature->addScore("delta_rt", mrmfeature->getRT() - expected_rt);
          mrmfeature->addScore("assay_rt", expected_rt);
          mrmfeature->addScore("norm_RT", normalized_experimental_rt);
          mrmfeature->addScore("rt_score", rt_score);
          mrmfeature->addScore("var_norm_rt_score", norm_rt_score);
        }

        // Intensity score
        double intensity_score = 0;
        if (use_intensity_score_)
        {
          intensity_score = mrmfeature->getIntensity() / (double)mrmfeature->getMetaValue("total_xic");
          mrmfeature->addScore("var_intensity_score", intensity_score);
        }

        double total_xic_score = 0;
        if (use_total_xic_score_)
        {
          total_xic_score = (double)mrmfeature->getMetaValue("total_xic");
          mrmfeature->addScore("total_xic", total_xic_score);
        }

        double nr_peaks_score = 0;
        if (use_nr_peaks_score_)
        {
          nr_peaks_score = group_size;
          mrmfeature->addScore("nr_peaks", nr_peaks_score);
        }

        double sn_score = 0, log_sn_score = 0;
        if (use_sn_score_)
        {
          sn_score = mrmscore.calcSNScore(imrmfeature, signal_noise_estimators);
          if (sn_score < 1) // fix to make sure, that log(sn_score = 0) = -inf does not occur
          {
            log_sn_score = 0;
          }
          else
          {
            log_sn_score = std::log(sn_score);
          }
          mrmfeature->addScore("sn_ratio", sn_score);
          mrmfeature->addScore("var_log_sn_score", log_sn_score);
        }

        OpenSwath_Scores scores;
        double quick_lda_dismiss = 0;
        double lda_quick_score = -scores.get_quick_lda_score(library_corr, library_rmsd, norm_rt_score, xcorr_coelution_score, xcorr_shape_score, log_sn_score);

        if (lda_quick_score < quick_lda_dismiss)
        {
          // continue;
        }

        double elution_model_fit_score = 0;
        if (use_elution_model_score_)
        {
          elution_model_fit_score = emgscoring.calcElutionFitScore((*mrmfeature), transition_group);
          mrmfeature->addScore("var_elution_model_fit_score", elution_model_fit_score);
        }

        double xx_lda_prescore;
        scores.library_corr              = library_corr;
        scores.library_rmsd              = library_rmsd;
        scores.norm_rt_score             = norm_rt_score;
        scores.elution_model_fit_score   = elution_model_fit_score;
        scores.log_sn_score              = log_sn_score;
        scores.xcorr_coelution_score     = xcorr_coelution_score;
        scores.xcorr_shape_score         = xcorr_shape_score;
        xx_lda_prescore = -scores.calculate_lda_prescore(scores);

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

        if (swath_present)
        {
          calculate_swath_scores(transition_group, *mrmfeature, swath_map, normalized_library_intensity, scores);
        }

#if 0
        if (do_local_fdr_)
        {
          calculate_local_fdr_scores(transition_group, *mrmfeature, trafo);
        }
#endif

        ///////////////////////////////////////////////////////////////////////////
        // add the peptide hit information to the feature
        ///////////////////////////////////////////////////////////////////////////

        const PeptideType * pep = PeptideRefMap[transition_group.getTransitions()[0].getPeptideRef()];
        const ProteinType * prot = ProteinRefMap[pep->protein_ref];

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
        if (stop_report_after_feature_ >= 0 && i >= (Size)stop_report_after_feature_ ) {break;} 
        output.push_back(feature_list[i]);
      }
    }

    template <template <typename> class SpectrumT, typename PeakT, typename TransitionT>
    void calculate_swath_scores(MRMTransitionGroup<SpectrumT, PeakT, TransitionT> & transition_group, MRMFeature & mrmfeature_,
                                OpenSwath::SpectrumAccessPtr swath_map, std::vector<double> & normalized_library_intensity, OpenSwath_Scores scores)
    {
      MRMFeature * mrmfeature = &mrmfeature_;

      // parameters
      int nr_isotopes = 4; // how many isotopes to consider
      int nr_charges = 4; // how many charges to look at for lower m/z
      int by_charge_state = 1; // for which charge states should we check b/y series

      int group_size = transition_group.size();
      diascoring.set_dia_parameters(dia_extract_window_, dia_centroided_, dia_byseries_intensity_min_, dia_byseries_ppm_diff_, nr_isotopes, nr_charges);

      // find spectrum that is closest to the apex of the peak using binary search
      // TODO (hroest): add up multiple spectra (as optional parameter)
      std::vector<std::size_t> indices = swath_map->getSpectraByRT(mrmfeature->getRT(), 0.0);
      int closest_idx = indices[0];
      OpenSwath::SpectrumPtr spectrum_;
      if (indices[0] != 0 &&
          std::fabs(swath_map->getSpectrumMetaById(indices[0] - 1).RT - mrmfeature->getRT()) <
          std::fabs(swath_map->getSpectrumMetaById(indices[0]).RT - mrmfeature->getRT()))
      { closest_idx--; }
      spectrum_ = swath_map->getSpectrumById(closest_idx);
      OpenSwath::SpectrumPtr * spectrum = &spectrum_;
#if 0
      MSSpectrum<SwathPeakType> * spectrum = &(*swath_map.RTBegin(mrmfeature->getRT()));
      if (spectrum != &(*swath_map.begin()) &&
          std::fabs((spectrum - 1)->getRT() - mrmfeature->getRT()) <
          std::fabs((spectrum)->getRT() - mrmfeature->getRT()))
      { spectrum--; }
#endif

      // Isotope correlation / overlap score: Is this peak part of an
      // isotopic pattern or is it the monoisotopic peak in an isotopic
      // pattern?
      OpenSwath::IMRMFeature * imrmfeature = new MRMFeatureOpenMS(*mrmfeature);
      double isotope_corr = 0, isotope_overlap = 0;
      diascoring.dia_isotope_scores(transition_group.getTransitions(),
        (*spectrum), imrmfeature, isotope_corr, isotope_overlap);
      // Mass deviation score
      double ppm_score = 0, ppm_score_weighted = 0;
      diascoring.dia_massdiff_score(transition_group.getTransitions(),
        (*spectrum), normalized_library_intensity, ppm_score, ppm_score_weighted);

      // Presence of b/y series score
      double bseries_score = 0, yseries_score = 0;
      const PeptideType * pep = PeptideRefMap[transition_group.getTransitions()[0].getPeptideRef()];
      OpenMS::AASequence aas = (String)pep->sequence;
      for (std::vector<ModificationType>::const_iterator it = pep->modifications.begin(); it != pep->modifications.end(); ++it)
      {
        aas.setModification(it->location, "UniMod:" + it->unimod_id);
      }
      diascoring.dia_by_ion_score((*spectrum), aas, by_charge_state, bseries_score, yseries_score);
      mrmfeature->addScore("var_isotope_correlation_score", isotope_corr);
      mrmfeature->addScore("var_isotope_overlap_score", isotope_overlap);
#ifdef DEBUG_MRMPEAKPICKER
      cout << "added corr isotope_score " << isotope_corr << endl;
      cout << "added overlap isotope_score " << isotope_overlap << endl;
#endif

      // FEATURE we should not punish so much when one transition is missing!
      double massdev_score = ppm_score / group_size;
      double massdev_score_weighted = ppm_score_weighted;
      mrmfeature->addScore("var_massdev_score", massdev_score);
      mrmfeature->addScore("var_massdev_score_weighted", massdev_score_weighted);
#ifdef DEBUG_MRMPEAKPICKER
      cout << "added score massdev_score " << massdev_score << endl;
      cout << "added score weighted massdev_score " << massdev_score_weighted << endl;
#endif

      mrmfeature->addScore("var_bseries_score", bseries_score);
      mrmfeature->addScore("var_yseries_score", yseries_score);
#ifdef DEBUG_MRMPEAKPICKER
      cout << "added score bseries_score " << bseries_score << endl;
      cout << "added score yseries_score " << yseries_score << endl;
#endif

      double dotprod_score_dia;
      double manhatt_score_dia;

      diascoring.score_with_isotopes((*spectrum), transition_group.getTransitions(), dotprod_score_dia, manhatt_score_dia);

      mrmfeature->addScore("var_dotprod_score", dotprod_score_dia);
      mrmfeature->addScore("var_manhatt_score", manhatt_score_dia);

      scores.yseries_score             = yseries_score;
      scores.isotope_correlation       = isotope_corr;
      scores.isotope_overlap           = isotope_overlap;
      scores.massdev_score             = massdev_score;
      double xx_swath_prescore = -scores.calculate_swath_lda_prescore(scores);
      mrmfeature->addScore("main_var_xx_swath_prelim_score", xx_swath_prescore);
      mrmfeature->setOverallQuality(xx_swath_prescore);
#ifdef DEBUG_MRMPEAKPICKER
      cout << "added xx_swath_prescore (everything above 2 is good) " << xx_swath_prescore << endl;
#endif
      delete imrmfeature;
    }

#if 0
    template <template <typename> class SpectrumT, typename PeakT, typename TransitionT>
    void calculate_local_fdr_scores(MRMTransitionGroup<SpectrumT, PeakT, TransitionT> & transition_group, MRMFeature & mrmfeature_, TransformationDescription & trafo)
    {
      int group_size = transition_group.size();
      MRMFeature * mrmfeature = &mrmfeature_;

      double local_fdr_score = 0;

      std::vector<std::pair<double, double> > library_intensity;
      std::vector<double> lib_int;
      std::vector<double> exp_int;
      double Tr_predicted, rmsd, delta, assay_score, normalized_experimental_rt;

      // compute the measured experimental intensities
      {
        std::vector<std::pair<double, double> > experimental_intensity;
        String native_id;
        for (int k = 0; k < group_size; k++)
        {
          native_id = transition_group.getChromatograms()[k].getNativeID();
          Feature f = mrmfeature->getFeature(native_id);
          experimental_intensity.push_back(std::make_pair(f.getMZ(), f.getIntensity()));
        }
        std::sort(experimental_intensity.begin(), experimental_intensity.end(), SortDoubleDoublePairFirst); // sort by m/z
        for (int k = 0; k < group_size; k++)
        {
          exp_int.push_back(experimental_intensity[k].second);
        }
      }

      // compute the score for the real transition
      {
        double experimental_rt, inten;
        experimental_rt = mrmfeature->getRT();
        normalized_experimental_rt = trafo.apply(experimental_rt);

        Tr_predicted = PeptideRTMap[transition_group.getTransitions()[0].getPeptideRef()];
        for (int k = 0; k < group_size; k++)
        {
          inten = transition_group.getTransitions()[k].getLibraryIntensity();
          library_intensity.push_back(std::make_pair(transition_group.getTransitions()[k].getProductMZ(), inten));
        }
        std::sort(library_intensity.begin(), library_intensity.end(), SortDoubleDoublePairFirst); // sort by m/z
        for (Size k = 0; k < transition_group.size(); k++)
        {
          lib_int.push_back(library_intensity[k].second);
        }

        rmsd = Scoring::RMSD(&exp_int[0], &lib_int[0], group_size);
        delta = std::fabs(Tr_predicted - normalized_experimental_rt);
        delta *= delta;
        assay_score = 4.699867 - 0.0001334994 * delta - 17.16685 * rmsd;
      }

      // now loop over all assays that we have and compare them to our
      // *real* assay and count how many of them score better than
      // the real one.
      for (TransitionGroupMapType::const_iterator all_itr = lib_transition_group_map.begin(); all_itr != lib_transition_group_map.end(); all_itr++)
      {

        library_intensity.clear();
        lib_int.clear();
        Tr_predicted = 0;

        // skip the real assay (assuming they have the same transition group id)
        if (all_itr->second.getTransitionGroupID() == transition_group.getTransitionGroupID())
        {
          continue;
        }

        // get the peptide retention time
        Tr_predicted = PeptideRTMap[transition_group.getTransitions()[0].getPeptideRef()];

        // TODO (weissrh) some mapping magic on how to map this assay on the one we have
        for (Size k = 0; k < all_itr->second.getTransitions().size(); k++)
        {
          double inten = all_itr->second.getTransitions()[k].getLibraryIntensity();
          library_intensity.push_back(std::make_pair(all_itr->second.getTransitions()[k].getProductMZ(), inten));
        }
        std::sort(library_intensity.begin(), library_intensity.end(), SortDoubleDoublePairFirst); // sort by m/z
        for (Size k = 0; k < all_itr->second.getTransitions().size(); k++)
        {
          lib_int.push_back(library_intensity[k].second);
        }

        int n = std::min(group_size, (int)all_itr->second.getTransitions().size());
        rmsd = Scoring::RMSD(&exp_int[0], &lib_int[0], n);

        // double correlation = gsl_stats_correlation(normx, 1, normy, 1, group_size);

        // apply GLM score and add punishment if another random assay scores
        // better than the target assay.
        //
        //     my $glm = 4.699867-0.0001334994*$delta-17.16685*$rmsd;
        //     $glm = -200 if $glm < -200; # things fall apart when the value is
        //     really low
        //     $self->{_model_prob} = 1/(1+1/exp($glm));
        //
        //
        //     where $delta is (deltatime)^2:
        //     $delta =
        //     ($self->get_ntime()-$A->get_ntime())*($self->get_ntime()-$A->get_ntime()) if $delta == -1;
        //
        //     and $rmsd is
        //
        //     $residuals +=
        //     ($TRANSPEAK->get_rel_area()-$TRANSITION->get_rel_area())*($TRANSPEAK->get_rel_area()-$TRANSITION->get_rel_area());
        //
        delta = std::fabs(Tr_predicted - normalized_experimental_rt);
        delta *= delta;
        // double this_score = 3*delta + 3*rmsd ; //+ 3*correlation;
        double this_score = 4.699867 - 0.0001334994 * delta - 17.16685 * rmsd;

        // TODO (weisserh) do something smarter than counting how many assays have a better score!
        if (this_score > assay_score)
        {
          //cout << " glm score " << this_score << " ( 3*" << delta<< " + 3*" <<  rmsd<< " ) " << lib_int.size( ) << " " << all_itr->second.size()  << " " << exp_int.size() <<  " vs assay score " << assay_score << endl;
          local_fdr_score += 1;
        }
      }

      mrmfeature->addScore("var_local_fdr_score", local_fdr_score);
#ifdef DEBUG_MRMPEAKPICKER
      cout << "added score local_fdr_score " << local_fdr_score << endl;
#endif

    }
#endif

    void handle_params();

    // Variables

    UInt sgolay_frame_length_;
    UInt sgolay_polynomial_order_;
    DoubleReal gauss_width_;

    DoubleReal peak_width_;
    DoubleReal signal_to_noise_;

    DoubleReal sn_win_len_;
    UInt sn_bin_count_;

    int stop_after_feature_;
    int stop_report_after_feature_;
    DoubleReal stop_after_intensity_ratio_;

    int emgfitter_maxiterations_;

    DoubleReal rt_extraction_window_;
    DoubleReal quantification_cutoff_;

    // for the SWATH analysis
    DoubleReal dia_byseries_ppm_diff_;
    DoubleReal dia_byseries_intensity_min_;
    DoubleReal dia_extract_window_;
    bool dia_centroided_;

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

    bool do_local_fdr_;
    bool write_convex_hull_;

    bool strict;

    DoubleReal rt_normalization_factor_;

    std::map<OpenMS::String, double> PeptideRTMap;

    std::map<OpenMS::String, const PeptideType *> PeptideRefMap;
    std::map<OpenMS::String, const ProteinType *> ProteinRefMap;

    TransitionGroupMapType lib_transition_group_map;

    MRMScoring mrmscore;
    OpenSwath::DIAScoring diascoring;
    EmgScoring emgscoring;

  };
}

#undef run_identifier
#endif
