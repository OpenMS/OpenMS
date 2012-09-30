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

//#define MRMSCORING_TESTING
#include <algorithm>
#include <algorithm>
#include <iterator>
#include <iostream>

#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/MRMScoring.h"
#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/Scoring.h"
#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/meanAndSd.h"

#ifdef OPENMS_ASSERTIONS
#define OPENMS_PRECONDITION(condition, message)\
	if (!(condition))\
    { throw std::runtime_error(message); }
#else
#define OPENMS_PRECONDITION(condition, message)
#endif

namespace OpenMS
{

  const MRMScoring::XCorrMatrixType & MRMScoring::getXCorrMatrix() const
  {
    return xcorr_matrix_;
  }

  void MRMScoring::initializeXCorrMatrix(OpenSwath::IMRMFeature * mrmfeature, OpenSwath::ITransitionGroup * transition_group, bool normalize)
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

  // see /IMSB/users/reiterl/bin/code/biognosys/trunk/libs/mrm_libs/MRM_pgroup.pm
  // _calc_xcorr_coelution_score
  //
  //   for each i,j get xcorr_matrix array => find max of the crosscorrelation
  //   store the delta to the retention time
  // return $deltascore_mean + $deltascore_stdev
  double MRMScoring::calcXcorrCoelutionScore()
  {
    OPENMS_PRECONDITION(xcorr_matrix_.size() > 1, "Expect cross-correlation matrix of at least 2x2");

    std::vector<int> deltas;
    for (std::size_t i = 0; i < xcorr_matrix_.size(); i++)
    {
      for (std::size_t  j = i; j < xcorr_matrix_.size(); j++)
      {
        // first is the X value (RT), should be an int
        deltas.push_back(std::fabs(Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][j])->first));
#ifdef MRMSCORING_TESTING
        std::cout << "&&_xcoel append " << std::fabs(Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][j])->first) << std::endl;
#endif
      }
    }

    OpenSwath::mean_and_stddev msc;
    msc = std::for_each(deltas.begin(), deltas.end(), msc);
    double deltas_mean = msc.mean();
    double deltas_stdv = msc.sample_stddev();
    //double deltas_mean = gsl_stats_int_mean(&deltas[0], 1, deltas.size());
    //double deltas_stdv = gsl_stats_int_sd(&deltas[0], 1, deltas.size());

    double xcorr_coelution_score = deltas_mean + deltas_stdv;
    return xcorr_coelution_score;
  }

  double MRMScoring::calcXcorrCoelutionScore_weighted(
    const std::vector<double> & normalized_library_intensity)
  {
    OPENMS_PRECONDITION(xcorr_matrix_.size() > 1, "Expect cross-correlation matrix of at least 2x2");

#ifdef MRMSCORING_TESTING
    double weights = 0;
#endif
    std::vector<double> deltas;
    for (std::size_t i = 0; i < xcorr_matrix_.size(); i++)
    {
      deltas.push_back(
        std::fabs(Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][i])->first)
        * normalized_library_intensity[i]
        * normalized_library_intensity[i]);
#ifdef MRMSCORING_TESTING
      std::cout << "_xcoel_weighted " << i << " " << i << " " << Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][i])->first << " weight " <<
      normalized_library_intensity[i] * normalized_library_intensity[i] << std::endl;
      weights += normalized_library_intensity[i] * normalized_library_intensity[i];
#endif
      for (std::size_t j = i + 1; j < xcorr_matrix_.size(); j++)
      {
        // first is the X value (RT), should be an int
        deltas.push_back(
          std::fabs(Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][j])->first)
          * normalized_library_intensity[i]
          * normalized_library_intensity[j] * 2);
#ifdef MRMSCORING_TESTING
        std::cout << "_xcoel_weighted " << i << " " << j << " " << Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][j])->first << " weight " <<
        normalized_library_intensity[i] * normalized_library_intensity[j] * 2 << std::endl;
        weights += normalized_library_intensity[i] * normalized_library_intensity[j];
#endif

      }
    }

#ifdef MRMSCORING_TESTING
    std::cout << " all weights sum " << weights << std::endl;
#endif

    /*
     double deltas_mean = gsl_stats_int_mean(&deltas[0],1,deltas.size() );
     double deltas_stdv = gsl_stats_int_sd(  &deltas[0],1,deltas.size() );

     double xcorr_coelution_score = deltas_mean + deltas_stdv;
     return xcorr_coelution_score;
     */
    return std::accumulate(deltas.begin(), deltas.end(), 0.0);
  }

  // see /IMSB/users/reiterl/bin/code/biognosys/trunk/libs/mrm_libs/MRM_pgroup.pm
  // _calc_xcorr_shape_score
  //
  //   for each i,j get xcorr_matrix array => find max of the crosscorrelation
  //   calculate whether the maximal crosscorrelation coincides with the maximal intensity
  ///
  double MRMScoring::calcXcorrShape_score()
  {
    OPENMS_PRECONDITION(xcorr_matrix_.size() > 1, "Expect cross-correlation matrix of at least 2x2");

    std::vector<double> intensities;
    for (std::size_t i = 0; i < xcorr_matrix_.size(); i++)
    {
      for (std::size_t j = i; j < xcorr_matrix_.size(); j++)
      {
        // second is the Y value (intensity)
        intensities.push_back(Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][j])->second);
      }
    }
    OpenSwath::mean_and_stddev msc;
    msc = std::for_each(intensities.begin(), intensities.end(), msc);
    return msc.mean();
  }

  double MRMScoring::calcXcorrShape_score_weighted(
    const std::vector<double> & normalized_library_intensity)
  {
    OPENMS_PRECONDITION(xcorr_matrix_.size() > 1, "Expect cross-correlation matrix of at least 2x2");

    // TODO (hroest) : check implementation
    //         see _calc_weighted_xcorr_shape_score in MRM_pgroup.pm
    //         -- they only multiply up the intensity once
    std::vector<double> intensities;
    for (std::size_t i = 0; i < xcorr_matrix_.size(); i++)
    {
      intensities.push_back(
        Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][i])->second
        * normalized_library_intensity[i]
        * normalized_library_intensity[i]);
#ifdef MRMSCORING_TESTING
      std::cout << "_xcorr_weighted " << i << " " << i << " " << Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][i])->second << " weight " <<
      normalized_library_intensity[i] * normalized_library_intensity[i] << std::endl;
#endif
      for (std::size_t j = i + 1; j < xcorr_matrix_.size(); j++)
      {
        intensities.push_back(
          Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][j])->second
          * normalized_library_intensity[i]
          * normalized_library_intensity[j] * 2);
#ifdef MRMSCORING_TESTING
        std::cout << "_xcorr_weighted " << i << " " << j << " " << Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][j])->second << " weight " <<
        normalized_library_intensity[i] * normalized_library_intensity[j] * 2 << std::endl;
#endif
      }
    }
    return std::accumulate(intensities.begin(), intensities.end(), 0.0);
  }

  void MRMScoring::calcLibraryScore(OpenSwath::IMRMFeature * mrmfeature, const std::vector<TransitionType> & transitions,
                                           double & correlation, double & rmsd, double & manhattan, double & dotprod)
  {
    std::vector<double> library_intensity;
    std::vector<double> experimental_intensity;
    String native_id;

    for (std::size_t k = 0; k < transitions.size(); k++)
    {
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

    OPENMS_PRECONDITION(library_intensity.size() == experimental_intensity.size(), "Both vectors need to have the same size");

#ifdef MRMSCORING_TESTING
    for (std::size_t k = 0; k < transitions.size(); k++)
    {
      native_id = transitions[k].getNativeID();
      std::cout << native_id << " Lib vs exp " << library_intensity[k] << " " << experimental_intensity[k] << std::endl;
    }
#endif

    manhattan = OpenSwath::manhattanScoring(experimental_intensity, library_intensity);
    dotprod = OpenSwath::dotprodScoring(experimental_intensity, library_intensity);


    Scoring::normalize_sum(&experimental_intensity[0], transitions.size());
    Scoring::normalize_sum(&library_intensity[0], transitions.size());

    rmsd = Scoring::RMSD(&experimental_intensity[0], &library_intensity[0], transitions.size());
    correlation = OpenSwath::cor_pearson(experimental_intensity.begin(), experimental_intensity.end(), library_intensity.begin());


    //double c0, c1, cov00, cov01, cov11, sumsq;
    //int ret = gsl_fit_linear(normx, 1, normy, 1, transitions.size(), &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
    if (boost::math::isnan(correlation))
      correlation = -1.0;
  }

  double MRMScoring::calcRTScore(const PeptideType & peptide, double normalized_experimental_rt)
  {
    double expected_rt;
    expected_rt = peptide.rt;

    if (expected_rt <= -1000)
    {
      return 0;
    }

    // use the transformed experimental retention time and then take the difference.
    double rt_score = std::fabs(normalized_experimental_rt - expected_rt);
    return rt_score;
  }

  double MRMScoring::calcSNScore(OpenSwath::IMRMFeature * mrmfeature, std::vector<OpenSwath::ISignalToNoisePtr> & signal_noise_estimators)
  {
    double sn_score = 0;

    for (std::size_t k = 0; k < signal_noise_estimators.size(); k++)
    {
      sn_score += signal_noise_estimators[k]->getValueAtRT(mrmfeature->getRT());
    }
    return sn_score / signal_noise_estimators.size();
  }

}
