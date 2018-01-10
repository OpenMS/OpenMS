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

#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/MRMScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/StatsHelpers.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/Macros.h>
//#define MRMSCORING_TESTING
#include <algorithm>
#include <iostream>
#include <iterator>


namespace OpenSwath
{

  const MRMScoring::XCorrMatrixType& MRMScoring::getXCorrMatrix() const
  {
    return xcorr_matrix_;
  }

  void MRMScoring::initializeXCorrMatrix(OpenSwath::IMRMFeature* mrmfeature, std::vector<String> native_ids)
  {
    std::vector<double> intensityi, intensityj;
    xcorr_matrix_.resize(native_ids.size());
    for (std::size_t i = 0; i < native_ids.size(); i++)
    {
      String native_id = native_ids[i];
      FeatureType fi = mrmfeature->getFeature(native_id);
      xcorr_matrix_[i].resize(native_ids.size());
      intensityi.clear();
      fi->getIntensity(intensityi);
      for (std::size_t j = i; j < native_ids.size(); j++)
      {
        String native_id2 = native_ids[j];
        FeatureType fj = mrmfeature->getFeature(native_id2);
        intensityj.clear();
        fj->getIntensity(intensityj);
        // compute normalized cross correlation
        xcorr_matrix_[i][j] = Scoring::normalizedCrossCorrelation(intensityi, intensityj, boost::numeric_cast<int>(intensityi.size()), 1);
      }
    }
  }

  void MRMScoring::initializeMS1XCorr(OpenSwath::IMRMFeature* mrmfeature, std::vector<String> native_ids, std::string precursor_id)
  {
    std::vector<double> intensityi, intensity_ms1;
    mrmfeature->getPrecursorFeature(precursor_id)->getIntensity(intensity_ms1);
    ms1_xcorr_vector_.resize(native_ids.size());
    for (std::size_t i = 0; i < native_ids.size(); i++)
    {
      String native_id = native_ids[i];
      FeatureType fi = mrmfeature->getFeature(native_id);
      intensityi.clear();
      fi->getIntensity(intensityi);
      ms1_xcorr_vector_[i] = Scoring::normalizedCrossCorrelation(
        intensityi, intensity_ms1, boost::numeric_cast<int>(intensityi.size()), 1);
    }
  }

  void MRMScoring::initializeXCorrIdMatrix(OpenSwath::IMRMFeature* mrmfeature, std::vector<String> native_ids_identification, std::vector<String> native_ids_detection)
  { 
    std::vector<double> intensityi, intensityj;
    xcorr_matrix_.resize(native_ids_identification.size());
    for (std::size_t i = 0; i < native_ids_identification.size(); i++)
    { 
      String native_id = native_ids_identification[i];
      FeatureType fi = mrmfeature->getFeature(native_id);
      xcorr_matrix_[i].resize(native_ids_detection.size());
      intensityi.clear();
      fi->getIntensity(intensityi);
      for (std::size_t j = 0; j < native_ids_detection.size(); j++)
      {
        String native_id2 = native_ids_detection[j];
        FeatureType fj = mrmfeature->getFeature(native_id2);
        intensityj.clear();
        fj->getIntensity(intensityj);
        // compute normalized cross correlation
        xcorr_matrix_[i][j] = Scoring::normalizedCrossCorrelation(intensityi, intensityj, boost::numeric_cast<int>(intensityi.size()), 1);
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
    OPENSWATH_PRECONDITION(xcorr_matrix_.size() > 1, "Expect cross-correlation matrix of at least 2x2");

    std::vector<int> deltas;
    for (std::size_t i = 0; i < xcorr_matrix_.size(); i++)
    {
      for (std::size_t  j = i; j < xcorr_matrix_.size(); j++)
      {
        // first is the X value (RT), should be an int
        deltas.push_back(std::abs(Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][j])->first));
#ifdef MRMSCORING_TESTING
        std::cout << "&&_xcoel append " << std::abs(Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][j])->first) << std::endl;
#endif
      }
    }

    OpenSwath::mean_and_stddev msc;
    msc = std::for_each(deltas.begin(), deltas.end(), msc);
    double deltas_mean = msc.mean();
    double deltas_stdv = msc.sample_stddev();

    double xcorr_coelution_score = deltas_mean + deltas_stdv;
    return xcorr_coelution_score;
  }

  std::string MRMScoring::calcIndXcorrIdCoelutionScore()
  {
    OPENSWATH_PRECONDITION(xcorr_matrix_.size() > 0 && xcorr_matrix_[0].size() > 1, "Expect cross-correlation matrix of at least 2x1");

    std::vector<double> deltas;
    for (std::size_t i = 0; i < xcorr_matrix_.size(); i++)
    {
      double deltas_id = 0;
      for (std::size_t  j = 0; j < xcorr_matrix_[0].size(); j++)
      {
        // first is the X value (RT), should be an int
        deltas_id += std::abs(Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][j])->first);
#ifdef MRMSCORING_TESTING
        std::cout << "&&_xcoel append " << std::abs(Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][j])->first) << std::endl;
#endif
      }
      deltas.push_back(deltas_id / xcorr_matrix_[0].size());
    }

    std::stringstream ss;
    for (size_t i = 0; i < deltas.size(); i++)
    {
      if (i != 0)
        ss << ";";
      ss << deltas[i];
    }

    return ss.str();
  }

  double MRMScoring::calcXcorrCoelutionScore_weighted(
    const std::vector<double>& normalized_library_intensity)
  {
    OPENSWATH_PRECONDITION(xcorr_matrix_.size() > 1, "Expect cross-correlation matrix of at least 2x2");

#ifdef MRMSCORING_TESTING
    double weights = 0;
#endif
    std::vector<double> deltas;
    for (std::size_t i = 0; i < xcorr_matrix_.size(); i++)
    {
      deltas.push_back(
        std::abs(Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][i])->first)
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
          std::abs(Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][j])->first)
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
    OPENSWATH_PRECONDITION(xcorr_matrix_.size() > 1, "Expect cross-correlation matrix of at least 2x2");

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

  std::string MRMScoring::calcIndXcorrIdShape_score()
  {
    OPENSWATH_PRECONDITION(xcorr_matrix_.size() > 0 && xcorr_matrix_[0].size() > 1, "Expect cross-correlation matrix of at least 2x1");

    std::vector<double> intensities;
    for (std::size_t i = 0; i < xcorr_matrix_.size(); i++)
    {
      double intensities_id = 0;
      for (std::size_t j = 0; j < xcorr_matrix_[0].size(); j++)
      {
        // second is the Y value (intensity)
        intensities_id += Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][j])->second;
      }
      intensities.push_back(intensities_id / xcorr_matrix_[0].size());
    }

    std::stringstream ss;
    for (size_t i = 0; i <intensities.size(); i++)
    {
      if (i != 0)
        ss << ";";
      ss << intensities[i];
    }

    return ss.str();
  }

  double MRMScoring::calcXcorrShape_score_weighted(
    const std::vector<double>& normalized_library_intensity)
  {
    OPENSWATH_PRECONDITION(xcorr_matrix_.size() > 1, "Expect cross-correlation matrix of at least 2x2");

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

  double MRMScoring::calcMS1XcorrCoelutionScore()
  {
    OPENSWATH_PRECONDITION(ms1_xcorr_vector_.size() > 1, "Expect cross-correlation vector of a size of least 2");

    std::vector<int> deltas;
    for (std::size_t i = 0; i < ms1_xcorr_vector_.size(); i++)
    {
      // first is the X value (RT), should be an int
      deltas.push_back(std::abs(Scoring::xcorrArrayGetMaxPeak(ms1_xcorr_vector_[i])->first));
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

  double MRMScoring::calcMS1XcorrShape_score()
  {
    OPENSWATH_PRECONDITION(ms1_xcorr_vector_.size() > 1, "Expect cross-correlation vector of a size of least 2");

    std::vector<double> intensities;
    for (std::size_t i = 0; i < ms1_xcorr_vector_.size(); i++)
    {
      // second is the Y value (intensity)
      intensities.push_back(Scoring::xcorrArrayGetMaxPeak(ms1_xcorr_vector_[i])->second);
    }
    OpenSwath::mean_and_stddev msc;
    msc = std::for_each(intensities.begin(), intensities.end(), msc);
    return msc.mean();
  }

  void MRMScoring::calcLibraryScore(OpenSwath::IMRMFeature* mrmfeature, const std::vector<TransitionType>& transitions,
                                    double& correlation, double& norm_manhattan, double& manhattan, double& dotprod, double& spectral_angle, double& rmsd)
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
      experimental_intensity.push_back(static_cast<double>(mrmfeature->getFeature(native_id)->getIntensity()));
      library_intensity.push_back(intensity);
    }

    OPENSWATH_PRECONDITION(library_intensity.size() == experimental_intensity.size(), "Both vectors need to have the same size");

#ifdef MRMSCORING_TESTING
    for (std::size_t k = 0; k < transitions.size(); k++)
    {
      native_id = transitions[k].getNativeID();
      std::cout << native_id << " Lib vs exp " << library_intensity[k] << " " << experimental_intensity[k] << std::endl;
    }
#endif

    manhattan = OpenSwath::manhattanScoring(experimental_intensity, library_intensity);
    dotprod = OpenSwath::dotprodScoring(experimental_intensity, library_intensity);

    spectral_angle = Scoring::SpectralAngle(&experimental_intensity[0], &library_intensity[0], boost::numeric_cast<unsigned int>(transitions.size()));

    Scoring::normalize_sum(&experimental_intensity[0], boost::numeric_cast<unsigned int>(transitions.size()));
    Scoring::normalize_sum(&library_intensity[0], boost::numeric_cast<unsigned int>(transitions.size()));

    norm_manhattan = Scoring::NormalizedManhattanDist(&experimental_intensity[0], &library_intensity[0], boost::numeric_cast<unsigned int>(transitions.size()));
    rmsd = Scoring::RootMeanSquareDeviation(&experimental_intensity[0], &library_intensity[0], boost::numeric_cast<unsigned int>(transitions.size()));
    correlation = OpenSwath::cor_pearson(experimental_intensity.begin(), experimental_intensity.end(), library_intensity.begin());

    if (boost::math::isnan(correlation))
    {
      correlation = -1.0;
    }
  }

  double MRMScoring::calcRTScore(const PeptideType& peptide, double normalized_experimental_rt)
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

  double MRMScoring::calcSNScore(OpenSwath::IMRMFeature* mrmfeature, std::vector<OpenSwath::ISignalToNoisePtr>& signal_noise_estimators)
  {
    OPENSWATH_PRECONDITION(signal_noise_estimators.size() > 0, "Input S/N estimators needs to be larger than 0");

    double sn_score = 0;
    if (signal_noise_estimators.size() == 0)
    {
      return 0;
    }

    for (std::size_t k = 0; k < signal_noise_estimators.size(); k++)
    {
      sn_score += signal_noise_estimators[k]->getValueAtRT(mrmfeature->getRT());
    }
    return sn_score / signal_noise_estimators.size();
  }

  std::string MRMScoring::calcIndSNScore(OpenSwath::IMRMFeature* mrmfeature, std::vector<OpenSwath::ISignalToNoisePtr>& signal_noise_estimators)
  {
    OPENSWATH_PRECONDITION(signal_noise_estimators.size() > 0, "Input S/N estimators needs to be larger than 0");

    std::vector<double> sn_score;
    if (signal_noise_estimators.size() == 0)
    {
      return std::string();
    }

    for (std::size_t k = 0; k < signal_noise_estimators.size(); k++)
    {
      sn_score.push_back(signal_noise_estimators[k]->getValueAtRT(mrmfeature->getRT()));
    }

    std::stringstream ss;
    for (size_t i = 0; i <sn_score.size(); i++)
    {
      if (i != 0)
        ss << ";";

      if (sn_score[i] < 1) // everything below S/N 1 can be set to zero (and the log safely applied)
      {
        ss << 0;
      }
      else
      {
        ss << std::log(sn_score[i]);
      }
    }

    return ss.str();
  }
}
