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

#ifndef OPENMS_ANALYSIS_OPENSWATH_OPENSWATHALGO_MRMFEATURESCORING_H
#define OPENMS_ANALYSIS_OPENSWATH_OPENSWATHALGO_MRMFEATURESCORING_H

#include <string>
#include <boost/math/special_functions/fpclassify.hpp> // for isnan

#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h"
#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ITransition.h"
#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h"
#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/corr.h"
#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/Scoring.h"
#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DIAHelpers.h"

namespace OpenMS
{

  /**
      @brief This class implements different scores for peaks found in SRM/MRM.

      It uses scores based on different parameters of the peaks from the
      individual transitions and stores them individually. The idea and the
      scores are based on the following paper:
        Reiter L, Rinner O, Picotti P, Hüttenhain R, Beck M, Brusniak MY,
        Hengartner MO, Aebersold R.  mProphet: automated data processing and
        statistical validation for large-scale SRM experiments.  Nat Methods.
        2011 May;8(5):430-5. Epub 2011 Mar 20.

      The currently implemented scores include:
      - xcorr_coelution: Crosscorrelation of the different transitions
      - xcorr_shape: Crosscorrelation shape score (whether the maximal
                     crosscorrelation coincides with the maximal intensity)
      - library_rmsd: normalized RMSD of the measured intensities to the expected intensities
      - library_correlation: correlation of the measured intensities to the expected intensities
      - rt_score: deviation from the expected retention time
      - elution_fit_score: how well the elution profile fits a theoretical elution profile

  */
  class MRMFeatureScoring
  {

public:

    ///Type definitions
    //@{
    /// Cross Correlation array
    typedef std::map<int, double> XCorrArrayType;
    /// Cross Correlation matrix
    typedef std::vector<std::vector<XCorrArrayType> > XCorrMatrixType;

    typedef std::string String;

    typedef OpenSwath::SpectrumPtr SpectrumType;
    typedef OpenSwath::LightTransition TransitionType;
    typedef OpenSwath::LightPeptide PeptideType;
    typedef OpenSwath::LightProtein ProteinType;

    typedef boost::shared_ptr<OpenSwath::IFeature> FeatureType;
    //@}

    /** @name Accessors */
    //@{
    /// non-muteable access to the Cross-correlation matrix
    const XCorrMatrixType & getXCorrMatrix() const;
    //@}

    /** @name Scores */
    //@{
    /// Initialize the scoring object and building the cross-correlation matrix
    void initializeXCorrMatrix(OpenSwath::IMRMFeature * mrmfeature, OpenSwath::ITransitionGroup * transition_group, bool normalize)
    {
      std::vector<double> intensityi, intensityj;
      xcorr_matrix_.resize(transition_group->size());
      for (std::size_t i = 0; i < transition_group->size(); i++)
      {
        String native_id = transition_group->getNativeIDs()[i];
        FeatureType fi = mrmfeature->getFeature(native_id);
        xcorr_matrix_[i].resize(transition_group->size());
        intensityi.clear();
        fi->getIntensity(intensityi);
        for (std::size_t j = i; j < transition_group->size(); j++)
        {
          String native_id2 = transition_group->getNativeIDs()[j];
          FeatureType fj = mrmfeature->getFeature(native_id2);
          intensityj.clear();
          fj->getIntensity(intensityj);
          //std::cout << " = Computing crosscorrelation " << i << " / "  << j << " or " << native_id << " vs " << native_id2 << std::endl;
          if (normalize)
          {
            xcorr_matrix_[i][j] = Scoring::normalizedCalcxcorr(intensityi, intensityj, intensityi.size(), 1);
          }
          else
          {
            throw "not implemented";
          }

        }
      }
    }

    /// calculate the cross-correlation score
    double calcXcorrCoelutionScore();

    /// calculate the cross-correlation shape score
    double calcXcorrShape_score();

    /// calculate the weighted cross-correlation shape score
    double calcXcorrShape_score_weighted(const std::vector<double> & normalized_library_intensity);

    /// calculate the weighted cross-correlation score
    double calcXcorrCoelutionScore_weighted(const std::vector<double> & normalized_library_intensity);

    /// calculate the library correlation score (correlation and rmsd)
    static void calcLibraryScore(OpenSwath::IMRMFeature * mrmfeature, const std::vector<TransitionType> & transitions,
        double & correlation, double & rmsd, double & manhattan, double & dotprod)
    {
      std::vector<double> library_intensity;
      std::vector<double> experimental_intensity;
      String native_id;

      for (std::size_t k = 0; k < transitions.size(); k++)
      {
        // TODO get apex or integrated intensity?
        native_id = transitions[k].getNativeID();
        double intensity = transitions[k].getLibraryIntensity();
        // the library intensity should never be below zero
        if (intensity < 0.0) 
        {
          intensity = 0.0;
        }
        experimental_intensity.push_back(mrmfeature->getFeature(native_id)->getIntensity());
        library_intensity.push_back(intensity);
      }

      //OPENMS_PRECONDITION(library_intensity.size() == experimental_intensity.size(), "Both vectors need to have the same size");

#ifdef MRMSCORING_TESTING
      for (std::size_t k = 0; k < transitions.size(); k++)
      {
        native_id = transitions[k].getNativeID();
        std::cout << native_id << " Lib vs exp " << library_intensity[k] << " " << experimental_intensity[k] << std::endl;
      }
#endif

      manhattan = OpenSwath::manhattanScoring(experimental_intensity,library_intensity);
      dotprod = OpenSwath::dotprodScoring(experimental_intensity,library_intensity);


      Scoring::normalize_sum(&experimental_intensity[0], transitions.size());
      Scoring::normalize_sum(&library_intensity[0], transitions.size());

      rmsd = Scoring::RMSD(&experimental_intensity[0], &library_intensity[0], transitions.size());
      correlation = OpenSwath::cor_pearson(experimental_intensity.begin(), experimental_intensity.end(),library_intensity.begin());


      //double c0, c1, cov00, cov01, cov11, sumsq;
      //int ret = gsl_fit_linear(normx, 1, normy, 1, transitions.size(), &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
      if (boost::math::isnan(correlation)) correlation = -1.0;
    }

    /// calculate the retention time correlation score
    static double calcRTScore(const PeptideType & peptide, double normalized_experimental_rt)
    {
      double experimental_rt, expected_rt;
      expected_rt = peptide.rt;

      if (expected_rt <= -1000)
      {
        return 0;
      }

      // use the transformed experimental retention time and then take the difference.
      double rt_score = std::fabs(normalized_experimental_rt - expected_rt);
      return rt_score;
    }

    /// calculate the Signal to Noise ratio
    //  using a vector of SignalToNoiseEstimatorMedian that were calculated for
    //  each chromatogram of the transition_group.
    static double calcSNScore(OpenSwath::IMRMFeature * mrmfeature, std::vector<OpenSwath::ISignalToNoisePtr> & signal_noise_estimators)
    {
      double sn_score = 0;

      for (std::size_t k = 0; k < signal_noise_estimators.size(); k++)
      {
      	sn_score += signal_noise_estimators[k]->getValueAtRT(mrmfeature->getRT());
      }
      return sn_score / signal_noise_estimators.size();
    }

    //@}

    /** @name Helper functions */
    //@{
/*
    double RMSD(double x[], double y[], int n);

    XCorrArrayType calcxcorr(std::vector<double> & data1, std::vector<double> & data2, bool normalize);

    /// Calculate crosscorrelation on std::vector data (which is first normalized)
    MRMFeatureScoring::XCorrArrayType normalizedCalcxcorr(std::vector<double> & data1, std::vector<double> & data2, int maxdelay, int lag);

    /// Calculate crosscorrelation on std::vector data without normalization
    MRMFeatureScoring::XCorrArrayType calcxcorr_new(std::vector<double> & data1, std::vector<double> & data2, int maxdelay, int lag);

    MRMFeatureScoring::XCorrArrayType calcxcorr_lag1(std::vector<double> & data1, std::vector<double> & data2, int maxdelay);

    /// Find best peak in an cross-correlation (highest apex)
    XCorrArrayType::iterator xcorrArrayGetMaxPeak(XCorrArrayType array);

    /// Standardize a vector (subtract mean, divide by standard deviation)
    void standardize_data(std::vector<double> & data);

    /// divide each element of x by the sum of the vector
    void normalize_sum(double x[], int n);
*/

     //@}


private:



    /// Calculate crosscorrelation on Convex Hulls

    /// Fxn from mQuest to calculate similarity between library intensity and experimental ones
    //double deltaRatioSum(double x[], double y[], int n);

    /// Fxn from FeatureFinderAlgorithmMRM

    /** @name Members */
    //@{

    /// An Emg fitter for the elution profile score
    //TODO remove EmgFitter and related methods to own class.

    //EmgFitter1D fitter_emg1D;

    /// the precomputed cross correlation matrix
    XCorrMatrixType xcorr_matrix_;



    //@}

  };
}

#endif
