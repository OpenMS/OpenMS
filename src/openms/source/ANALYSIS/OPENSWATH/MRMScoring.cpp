// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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

#include <OpenMS/ANALYSIS/OPENSWATH/MRMScoring.h>
#include <OpenMS/OPENSWATHALGO/ALGO/StatsHelpers.h>
#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/OPENSWATHALGO/Macros.h>
//#define MRMSCORING_TESTING
#include <algorithm>
#include <iostream>
#include <iterator>

#include <boost/math/special_functions/fpclassify.hpp> // for isnan
#include <boost/numeric/conversion/cast.hpp>

namespace OpenSwath
{

    const MRMScoring::XCorrMatrixType& MRMScoring::getXCorrMatrix() const
    {
      return xcorr_matrix_;
    }

    void MRMScoring::initializeXCorrMatrix(const std::vector< std::vector< double > >& data)
    {
      xcorr_matrix_.resize(data.size(), data.size());
      xcorr_matrix_max_peak_.resize(data.size(), data.size());
      xcorr_matrix_max_peak_sec_.resize(data.size(), data.size());

      for (std::size_t i = 0; i < data.size(); i++)
      {
        std::vector< double > tmp1(data[i]);
        for (std::size_t j = i; j < data.size(); j++)
        {
          // compute normalized cross correlation
          std::vector< double > tmp2(data[j]);
          xcorr_matrix_.setValue(i, j, Scoring::normalizedCrossCorrelation(tmp1, tmp2, boost::numeric_cast<int>(data[i].size()), 1));
          auto x = Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_.getValue(i, j));
          xcorr_matrix_max_peak_.setValue(i, j, std::abs(x->first));
          xcorr_matrix_max_peak_sec_.setValue(i, j, x->second);
        }
      }
    }

    const MRMScoring::XCorrMatrixType& MRMScoring::getXCorrContrastMatrix() const
    {
      return xcorr_contrast_matrix_;
    }

    const MRMScoring::XCorrMatrixType& MRMScoring::getXCorrPrecursorContrastMatrix() const
    {
      return xcorr_precursor_contrast_matrix_;
    }

    const MRMScoring::XCorrMatrixType& MRMScoring::getXCorrPrecursorCombinedMatrix() const
    {
      return xcorr_precursor_combined_matrix_;
    }

    void MRMScoring::initializeXCorrMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<String>& native_ids)
    {
      std::vector<double> intensityi, intensityj;
      xcorr_matrix_.resize(native_ids.size(), native_ids.size());
      xcorr_matrix_max_peak_.resize(native_ids.size(), native_ids.size());
      xcorr_matrix_max_peak_sec_.resize(native_ids.size(), native_ids.size());
      for (std::size_t i = 0; i < native_ids.size(); i++)
      {
        FeatureType fi = mrmfeature->getFeature(native_ids[i]);
        intensityi.clear();
        fi->getIntensity(intensityi);
        for (std::size_t j = i; j < native_ids.size(); j++)
        {
          FeatureType fj = mrmfeature->getFeature(native_ids[j]);
          intensityj.clear();
          fj->getIntensity(intensityj);
          // compute normalized cross correlation
          xcorr_matrix_.setValue(i, j, Scoring::normalizedCrossCorrelation(intensityi, intensityj, boost::numeric_cast<int>(intensityi.size()), 1));
          auto x = Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_.getValue(i, j));
          xcorr_matrix_max_peak_.setValue(i, j, std::abs(x->first));
          xcorr_matrix_max_peak_sec_.setValue(i, j, x->second);
        }
      }
    }

    void MRMScoring::initializeXCorrContrastMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<String>& native_ids_set1, const std::vector<String>& native_ids_set2)
    {
      std::vector<double> intensityi, intensityj;
      xcorr_contrast_matrix_.resize(native_ids_set1.size(), native_ids_set2.size());
      xcorr_contrast_matrix_max_peak_.resize(native_ids_set1.size(), native_ids_set2.size());
      xcorr_contrast_matrix_max_peak_sec_.resize(native_ids_set1.size(), native_ids_set2.size());
      for (std::size_t i = 0; i < native_ids_set1.size(); i++)
      {
        FeatureType fi = mrmfeature->getFeature(native_ids_set1[i]);
        intensityi.clear();
        fi->getIntensity(intensityi);
        for (std::size_t j = 0; j < native_ids_set2.size(); j++)
        {
          FeatureType fj = mrmfeature->getFeature(native_ids_set2[j]);
          intensityj.clear();
          fj->getIntensity(intensityj);
          // compute normalized cross correlation
          xcorr_contrast_matrix_.setValue(i, j, Scoring::normalizedCrossCorrelation(intensityi, intensityj, boost::numeric_cast<int>(intensityi.size()), 1));
          auto x = Scoring::xcorrArrayGetMaxPeak(xcorr_contrast_matrix_.getValue(i, j));
          xcorr_contrast_matrix_max_peak_.setValue(i, j, std::abs(x->first));
          xcorr_contrast_matrix_max_peak_sec_.setValue(i, j, x->second);
        }
      }
    }

    void MRMScoring::initializeXCorrPrecursorMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<String>& precursor_ids)
    {
      std::vector<double> intensityi, intensityj;
      xcorr_precursor_matrix_.resize(precursor_ids.size(), precursor_ids.size());
      xcorr_precursor_matrix_max_peak_.resize(precursor_ids.size(), precursor_ids.size());
      xcorr_precursor_matrix_max_peak_sec_.resize(precursor_ids.size(), precursor_ids.size());
      for (std::size_t i = 0; i < precursor_ids.size(); i++)
      {
        FeatureType fi = mrmfeature->getPrecursorFeature(precursor_ids[i]);
        intensityi.clear();
        fi->getIntensity(intensityi);
        for (std::size_t j = i; j < precursor_ids.size(); j++)
        {
          FeatureType fj = mrmfeature->getPrecursorFeature(precursor_ids[j]);
          intensityj.clear();
          fj->getIntensity(intensityj);
          // compute normalized cross correlation
          xcorr_precursor_matrix_.setValue(i, j, Scoring::normalizedCrossCorrelation(intensityi, intensityj, boost::numeric_cast<int>(intensityi.size()), 1));
          auto x = Scoring::xcorrArrayGetMaxPeak(xcorr_precursor_matrix_.getValue(i, j));
          xcorr_precursor_matrix_max_peak_.setValue(i, j, std::abs(x->first));
          xcorr_precursor_matrix_max_peak_sec_.setValue(i, j, x->second);
        }
      }
    }

    void MRMScoring::initializeXCorrPrecursorContrastMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<String>& precursor_ids, const std::vector<String>& native_ids)
    {
      std::vector<double> intensityi, intensityj;
      xcorr_precursor_contrast_matrix_.resize(precursor_ids.size(), native_ids.size());
      xcorr_precursor_contrast_matrix_max_peak_.resize(precursor_ids.size(), native_ids.size());
      xcorr_precursor_contrast_matrix_max_peak_sec_.resize(precursor_ids.size(), native_ids.size());
      for (std::size_t i = 0; i < precursor_ids.size(); i++)
      {
        FeatureType fi = mrmfeature->getPrecursorFeature(precursor_ids[i]);
        intensityi.clear();
        fi->getIntensity(intensityi);
        for (std::size_t j = 0; j < native_ids.size(); j++)
        {
          FeatureType fj = mrmfeature->getFeature(native_ids[j]);
          intensityj.clear();
          fj->getIntensity(intensityj);
          // compute normalized cross correlation
          xcorr_precursor_contrast_matrix_.setValue(i, j, Scoring::normalizedCrossCorrelation(intensityi, intensityj, boost::numeric_cast<int>(intensityi.size()), 1));
          auto x = Scoring::xcorrArrayGetMaxPeak(xcorr_precursor_contrast_matrix_.getValue(i, j));
          xcorr_precursor_contrast_matrix_max_peak_.setValue(i, j, std::abs(x->first));
          xcorr_precursor_contrast_matrix_max_peak_sec_.setValue(i, j, x->second);
        }
      }
    }

    void MRMScoring::initializeXCorrPrecursorContrastMatrix(const std::vector< std::vector< double > >& data_precursor, const std::vector< std::vector< double > >& data_fragments)
    {
      xcorr_precursor_contrast_matrix_.resize(data_precursor.size(), data_fragments.size());
      xcorr_precursor_contrast_matrix_max_peak_.resize(data_precursor.size(), data_fragments.size());
      xcorr_precursor_contrast_matrix_max_peak_sec_.resize(data_precursor.size(), data_fragments.size());
      for (std::size_t i = 0; i < data_precursor.size(); i++)
      {
        std::vector< double > tmp1(data_precursor[i]);
        for (std::size_t j = 0; j < data_fragments.size(); j++)
        {
          // compute normalized cross correlation
          std::vector< double > tmp2(data_fragments[j]);
          xcorr_precursor_contrast_matrix_.setValue(i, j, Scoring::normalizedCrossCorrelation(tmp1, tmp2, boost::numeric_cast<int>(tmp1.size()), 1));
          auto x = Scoring::xcorrArrayGetMaxPeak(xcorr_precursor_contrast_matrix_.getValue(i, j));
          xcorr_precursor_contrast_matrix_max_peak_.setValue(i, j, std::abs(x->first));
          xcorr_precursor_contrast_matrix_max_peak_sec_.setValue(i, j, x->second);
#ifdef MRMSCORING_TESTING
          std::cout << " fill xcorr_precursor_contrast_matrix_ "<< tmp1.size() << " / " << tmp2.size() << " : " << xcorr_precursor_contrast_matrix_[i][j].data.size() << std::endl;
#endif
        }
      }
    }

    void MRMScoring::initializeXCorrPrecursorCombinedMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<String>& precursor_ids, const std::vector<String>& native_ids)
    {
      std::vector<double> intensityi, intensityj;
      std::vector<FeatureType> features;

      for (std::size_t i = 0; i < precursor_ids.size(); i++)
      {
        FeatureType fi = mrmfeature->getPrecursorFeature(precursor_ids[i]);
        features.push_back(fi);
      }
      for (std::size_t j = 0; j < native_ids.size(); j++)
      {
        FeatureType fj = mrmfeature->getFeature(native_ids[j]);
        features.push_back(fj);
      }

      xcorr_precursor_combined_matrix_.resize(features.size(), features.size());
      xcorr_precursor_combined_matrix_max_peak_.resize(features.size(), features.size());
      xcorr_precursor_combined_matrix_max_peak_sec_.resize(features.size(), features.size());
      for (std::size_t i = 0; i < features.size(); i++)
      {
        FeatureType fi = features[i];
        intensityi.clear();
        fi->getIntensity(intensityi);
        for (std::size_t j = 0; j < features.size(); j++)
        {
          FeatureType fj = features[j];
          intensityj.clear();
          fj->getIntensity(intensityj);
          // compute normalized cross correlation
          xcorr_precursor_combined_matrix_.setValue(i, j, Scoring::normalizedCrossCorrelation(intensityi, intensityj, boost::numeric_cast<int>(intensityi.size()), 1));
          auto x = Scoring::xcorrArrayGetMaxPeak(xcorr_precursor_combined_matrix_.getValue(i, j));
          xcorr_precursor_combined_matrix_max_peak_.setValue(i, j, std::abs(x->first));
          xcorr_precursor_combined_matrix_max_peak_sec_.setValue(i, j, x->second);
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
      OPENSWATH_PRECONDITION(xcorr_matrix_max_peak_.rows() > 1, "Expect cross-correlation matrix of at least 2x2");

      std::vector<int> deltas;
      for (std::size_t i = 0; i < xcorr_matrix_max_peak_.rows(); i++)
      {
        for (std::size_t  j = i; j < xcorr_matrix_max_peak_.rows(); j++)
        {
          // first is the X value (RT), should be an int
          //deltas.push_back(std::abs(Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_.getValue(i, j))->first));
          deltas.push_back(xcorr_matrix_max_peak_.getValue(i,j));
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

    double MRMScoring::calcXcorrCoelutionWeightedScore(
            const std::vector<double>& normalized_library_intensity)
    {
      OPENSWATH_PRECONDITION(xcorr_matrix_max_peak_.rows() > 1, "Expect cross-correlation matrix of at least 2x2");

#ifdef MRMSCORING_TESTING
      double weights = 0;
#endif
      double deltas{0};
      for (std::size_t i = 0; i < xcorr_matrix_max_peak_.rows(); i++)
      {
        deltas += (xcorr_matrix_max_peak_.getValue(i, i)//std::abs(Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_.getValue(i, i))->first)
                   * normalized_library_intensity[i]
                   * normalized_library_intensity[i]);
#ifdef MRMSCORING_TESTING
        std::cout << "_xcoel_weighted " << i << " " << i << " " << Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][i])->first << " weight " <<
        normalized_library_intensity[i] * normalized_library_intensity[i] << std::endl;
      weights += normalized_library_intensity[i] * normalized_library_intensity[i];
#endif
        for (std::size_t j = i + 1; j < xcorr_matrix_max_peak_.rows(); j++)
        {
          // first is the X value (RT), should be an int
          deltas += (xcorr_matrix_max_peak_.getValue(i, j)//std::abs(Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_.getValue(i, j))->first)
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

      return deltas;
    }

    double MRMScoring::calcXcorrContrastCoelutionScore()
    {
      OPENSWATH_PRECONDITION(xcorr_contrast_matrix_max_peak_.rows() > 0 && xcorr_contrast_matrix_max_peak_.cols() > 1, "Expect cross-correlation matrix of at least 1x2");

      std::vector<int> deltas;
      for (auto e : xcorr_contrast_matrix_max_peak_)
      {
        // first is the X value (RT), should be an int
        deltas.push_back(e);          //(std::abs(Scoring::xcorrArrayGetMaxPeak(e)->first));
#ifdef MRMSCORING_TESTING
        std::cout << "&&_xcoel append " << std::abs(Scoring::xcorrArrayGetMaxPeak(xcorr_contrast_matrix_[i][j])->first) << std::endl;
#endif
      }

      OpenSwath::mean_and_stddev msc;
      msc = std::for_each(deltas.begin(), deltas.end(), msc);
      double deltas_mean = msc.mean();
      double deltas_stdv = msc.sample_stddev();

      double xcorr_coelution_score = deltas_mean + deltas_stdv;
      return xcorr_coelution_score;
    }

    std::vector<double> MRMScoring::calcSeparateXcorrContrastCoelutionScore()
    {
      OPENSWATH_PRECONDITION(xcorr_contrast_matrix_max_peak_.rows() > 0 && xcorr_contrast_matrix_max_peak_.cols() > 1, "Expect cross-correlation matrix of at least 1x2");

      std::vector<double > deltas;
      for (std::size_t i = 0; i < xcorr_contrast_matrix_max_peak_.rows(); i++)
      {
        double deltas_id = 0;
        for (std::size_t  j = 0; j < xcorr_contrast_matrix_max_peak_.cols(); j++)
        {
          // first is the X value (RT), should be an int
          deltas_id += xcorr_contrast_matrix_max_peak_.getValue(i, j);
#ifdef MRMSCORING_TESTING
          std::cout << "&&_xcoel append " << xcorr_contrast_matrix_max_peak_getValue(i, j) << std::endl;
#endif
        }
        deltas.push_back(deltas_id / xcorr_contrast_matrix_max_peak_.cols());
      }

      return deltas;
    }

    double MRMScoring::calcXcorrPrecursorCoelutionScore()
    {
      OPENSWATH_PRECONDITION(xcorr_precursor_matrix_max_peak_.rows() > 1, "Expect cross-correlation matrix of at least 2x2");

      std::vector<int> deltas;
      for (std::size_t i = 0; i < xcorr_precursor_matrix_max_peak_.rows(); i++)
      {
        for (std::size_t  j = i; j < xcorr_precursor_matrix_max_peak_.rows(); j++)
        {
          // first is the X value (RT), should be an int
          deltas.push_back(xcorr_precursor_matrix_max_peak_.getValue(i, j));
#ifdef MRMSCORING_TESTING
          std::cout << "&&_xcoel append " << std::abs(Scoring::xcorrArrayGetMaxPeak(xcorr_precursor_matrix_[i][j])->first) << std::endl;
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

    double MRMScoring::calcXcorrPrecursorContrastCoelutionScore()
    {
      OPENSWATH_PRECONDITION(xcorr_precursor_contrast_matrix_max_peak_.rows() > 0 && xcorr_precursor_contrast_matrix_max_peak_.cols() > 1, "Expect cross-correlation matrix of at least 1x2");

      std::vector<int> deltas;
      for (auto e : xcorr_precursor_contrast_matrix_max_peak_)
      {
        // first is the X value (RT), should be an int
        deltas.push_back(e);
#ifdef MRMSCORING_TESTING
        std::cout << "&&_xcoel append " << std::abs(Scoring::xcorrArrayGetMaxPeak(xcorr_precursor_contrast_matrix_[i][j])->first) << std::endl;
#endif
      }

      OpenSwath::mean_and_stddev msc;
      msc = std::for_each(deltas.begin(), deltas.end(), msc);
      double deltas_mean = msc.mean();
      double deltas_stdv = msc.sample_stddev();

      double xcorr_coelution_score = deltas_mean + deltas_stdv;
      return xcorr_coelution_score;
    }

    double MRMScoring::calcXcorrPrecursorCombinedCoelutionScore()
    {
      OPENSWATH_PRECONDITION(xcorr_precursor_combined_matrix_max_peak_.rows() > 1, "Expect cross-correlation matrix of at least 2x2");

      std::vector<int> deltas;
      for (std::size_t i = 0; i < xcorr_precursor_combined_matrix_max_peak_.rows(); i++)
      {
        for (std::size_t  j = i; j < xcorr_precursor_combined_matrix_max_peak_.rows(); j++)
        {
          // first is the X value (RT), should be an int
          deltas.push_back(xcorr_precursor_combined_matrix_max_peak_.getValue(i, j));
#ifdef MRMSCORING_TESTING
          std::cout << "&&_xcoel append " << std::abs(Scoring::xcorrArrayGetMaxPeak(xcorr_precursor_combined_matrix_[i][j])->first) << std::endl;
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

    // see /IMSB/users/reiterl/bin/code/biognosys/trunk/libs/mrm_libs/MRM_pgroup.pm
    // _calc_xcorr_shape_score
    //
    //   for each i,j get xcorr_matrix array => find max of the crosscorrelation
    //   calculate whether the maximal crosscorrelation coincides with the maximal intensity
    ///
    double MRMScoring::calcXcorrShapeScore()
    {
      OPENSWATH_PRECONDITION(xcorr_matrix_max_peak_sec_.rows() > 1, "Expect cross-correlation matrix of at least 2x2");

      size_t element_number{0};
      double intensities{0};
      for (std::size_t i = 0; i < xcorr_matrix_max_peak_sec_.rows(); i++)
      {
        for (std::size_t j = i; j < xcorr_matrix_max_peak_sec_.rows(); j++)
        {
          // second is the Y value (intensity)
          intensities += xcorr_matrix_max_peak_sec_.getValue(i, j);
          element_number++;
        }
      }
      return intensities / element_number;
    }

    double MRMScoring::calcXcorrShapeWeightedScore(
            const std::vector<double>& normalized_library_intensity)
    {
      OPENSWATH_PRECONDITION(xcorr_matrix_max_peak_sec_.rows() > 1, "Expect cross-correlation matrix of at least 2x2");

      // TODO (hroest) : check implementation
      //         see _calc_weighted_xcorr_shape_score in MRM_pgroup.pm
      //         -- they only multiply up the intensity once
      double intensities{0};
      for (std::size_t i = 0; i < xcorr_matrix_max_peak_sec_.rows(); i++)
      {
        intensities += (xcorr_matrix_max_peak_sec_.getValue(i, i)
                        * normalized_library_intensity[i]
                        * normalized_library_intensity[i]);
#ifdef MRMSCORING_TESTING
        std::cout << "_xcorr_weighted " << i << " " << i << " " << Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][i])->second << " weight " <<
        normalized_library_intensity[i] * normalized_library_intensity[i] << std::endl;
#endif
        for (std::size_t j = i + 1; j < xcorr_matrix_max_peak_sec_.rows(); j++)
        {
          intensities += (xcorr_matrix_max_peak_sec_.getValue(i, j)
                          * normalized_library_intensity[i]
                          * normalized_library_intensity[j] * 2);
#ifdef MRMSCORING_TESTING
          std::cout << "_xcorr_weighted " << i << " " << j << " " << Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][j])->second << " weight " <<
          normalized_library_intensity[i] * normalized_library_intensity[j] * 2 << std::endl;
#endif
        }
      }
      return intensities;
    }

    double MRMScoring::calcXcorrContrastShapeScore()
    {
      OPENSWATH_PRECONDITION(xcorr_contrast_matrix_max_peak_sec_.rows() > 0 && xcorr_contrast_matrix_max_peak_sec_.cols() > 1, "Expect cross-correlation matrix of at least 1x2");

      double intensities{0};
      for(auto e : xcorr_contrast_matrix_max_peak_sec_)
      {
        intensities += e;
      }
      return intensities;
    }

    std::vector<double> MRMScoring::calcSeparateXcorrContrastShapeScore()
    {
      OPENSWATH_PRECONDITION(xcorr_contrast_matrix_max_peak_sec_.rows() > 0 && xcorr_contrast_matrix_max_peak_sec_.cols() > 1, "Expect cross-correlation matrix of at least 1x2");

      std::vector<double> intensities;
      for (std::size_t i = 0; i < xcorr_contrast_matrix_max_peak_sec_.rows(); i++)
      {
        double intensities_id = 0;
        for (std::size_t j = 0; j < xcorr_contrast_matrix_max_peak_sec_.cols(); j++)
        {
          // second is the Y value (intensity)
          intensities_id += xcorr_contrast_matrix_max_peak_sec_.getValue(i,j);
        }
        intensities.push_back(intensities_id / xcorr_contrast_matrix_max_peak_sec_.cols());
      }

      return intensities;
    }

    double MRMScoring::calcXcorrPrecursorShapeScore()
    {
      OPENSWATH_PRECONDITION(xcorr_precursor_matrix_max_peak_sec_.rows() > 1, "Expect cross-correlation matrix of at least 2x2");

      double intensities{0};
      for(auto e : xcorr_precursor_matrix_max_peak_sec_)
      {
        intensities += e;
      }
      //xcorr_precursor_matrix_ is a triangle matrix
      size_t element_number = xcorr_precursor_matrix_max_peak_sec_.rows()*xcorr_precursor_matrix_.rows()/2 + (xcorr_precursor_matrix_max_peak_sec_.rows()+1)/2;
      return intensities / element_number;
    }

    double MRMScoring::calcXcorrPrecursorContrastShapeScore()
    {
      OPENSWATH_PRECONDITION(xcorr_precursor_contrast_matrix_max_peak_sec_.rows() > 0 && xcorr_precursor_contrast_matrix_max_peak_sec_.cols() > 1, "Expect cross-correlation matrix of at least 1x2");


      double intensities{0};
      for(auto e : xcorr_precursor_contrast_matrix_max_peak_sec_)
      {
        intensities += e;
      }
      return intensities / xcorr_precursor_contrast_matrix_max_peak_sec_.size();
    }

    double MRMScoring::calcXcorrPrecursorCombinedShapeScore()
    {
      OPENSWATH_PRECONDITION(xcorr_precursor_combined_matrix_max_peak_sec_.rows() > 1, "Expect cross-correlation matrix of at least 2x2");


      double intensities{0};
      for(size_t i = 0; i < xcorr_precursor_combined_matrix_max_peak_sec_.rows(); i++)
      {
        for(size_t j = i; j < xcorr_precursor_combined_matrix_max_peak_sec_.cols(); j++)
        {
          intensities += xcorr_precursor_combined_matrix_max_peak_sec_.getValue(i, j);
        }
      }
      //xcorr_precursor-combined_matrix_ is a triangle matrix
      size_t element_number = xcorr_precursor_combined_matrix_max_peak_sec_.rows()*xcorr_precursor_combined_matrix_max_peak_sec_.rows()/2 + (xcorr_precursor_combined_matrix_max_peak_sec_.rows()+1)/2;
      return intensities / element_number;
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

      if (boost::math::isnan(spectral_angle))
      {
        spectral_angle = 0.0;
      }

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

    std::vector<double> MRMScoring::calcSeparateSNScore(OpenSwath::IMRMFeature* mrmfeature, std::vector<OpenSwath::ISignalToNoisePtr>& signal_noise_estimators)
    {
      OPENSWATH_PRECONDITION(signal_noise_estimators.size() > 0, "Input S/N estimators needs to be larger than 0");

      std::vector<double> sn_scores;
      if (signal_noise_estimators.size() == 0)
      {
        return {};
      }

      for (std::size_t k = 0; k < signal_noise_estimators.size(); k++)
      {
        if (signal_noise_estimators[k]->getValueAtRT(mrmfeature->getRT()) < 1)
          // everything below S/N 1 can be set to zero (and the log safely applied)
        {
          sn_scores.push_back(0);
        }
        else
        {
          sn_scores.push_back(std::log(signal_noise_estimators[k]->getValueAtRT(mrmfeature->getRT())));
        }
      }

      return sn_scores;
    }

    const OpenMS::Matrix<double> & MRMScoring::getMIMatrix() const
    {
      return mi_matrix_;
    }

    const OpenMS::Matrix<double> & MRMScoring::getMIContrastMatrix() const
    {
      return mi_contrast_matrix_;
    }

    const OpenMS::Matrix<double> & MRMScoring::getMIPrecursorContrastMatrix() const
    {
      return mi_precursor_contrast_matrix_;
    }

    const OpenMS::Matrix<double> & MRMScoring::getMIPrecursorCombinedMatrix() const
    {
      return mi_precursor_combined_matrix_;
    }

    void MRMScoring::initializeMIMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<String>& native_ids)
    {
      std::vector<double> intensityi, intensityj;
      mi_matrix_.resize(native_ids.size(),native_ids.size());
      std::vector<unsigned int> rank_vec1, rank_vec2;
      for (std::size_t i = 0; i < native_ids.size(); i++)
      {
        FeatureType fi = mrmfeature->getFeature(native_ids[i]);

        intensityi.clear();
        fi->getIntensity(intensityi);
        Scoring::computeRank(intensityi, rank_vec1);
        for (std::size_t j = i; j < native_ids.size(); j++)
        {
          FeatureType fj = mrmfeature->getFeature(native_ids[j]);
          intensityj.clear();
          fj->getIntensity(intensityj);
          Scoring::computeRank(intensityj, rank_vec2);
          // compute ranked mutual information
          mi_matrix_.setValue(i,j,Scoring::rankedMutualInformation(rank_vec1, rank_vec2));
        }
      }
    }

    void MRMScoring::initializeMIContrastMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<String>& native_ids_set1, const std::vector<String>& native_ids_set2)
    {
      std::vector<double> intensityi, intensityj;
      mi_contrast_matrix_.resize(native_ids_set1.size(), native_ids_set2.size());
      std::vector<unsigned int> rank_vec1, rank_vec2;
      for (std::size_t i = 0; i < native_ids_set1.size(); i++)
      {
        FeatureType fi = mrmfeature->getFeature(native_ids_set1[i]);
        //mi_contrast_matrix_[i].resize(native_ids_set2.size());
        intensityi.clear();
        fi->getIntensity(intensityi);
        Scoring::computeRank(intensityi, rank_vec1);
        for (std::size_t j = 0; j < native_ids_set2.size(); j++)
        {
          FeatureType fj = mrmfeature->getFeature(native_ids_set2[j]);
          intensityj.clear();
          fj->getIntensity(intensityj);
          Scoring::computeRank(intensityj, rank_vec2);
          // compute ranked mutual information
          mi_contrast_matrix_.setValue(i, j, Scoring::rankedMutualInformation(rank_vec1, rank_vec2));
        }
      }
    }

    void MRMScoring::initializeMIPrecursorMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<String>& precursor_ids)
    {
      std::vector<double> intensityi, intensityj;
      mi_precursor_matrix_.resize(precursor_ids.size(),precursor_ids.size());
      std::vector<unsigned int> rank_vec1, rank_vec2;
      for (std::size_t i = 0; i < precursor_ids.size(); i++)
      {
        FeatureType fi = mrmfeature->getPrecursorFeature(precursor_ids[i]);
        intensityi.clear();
        fi->getIntensity(intensityi);
        Scoring::computeRank(intensityi, rank_vec1);
        for (std::size_t j = i; j < precursor_ids.size(); j++)
        {
          FeatureType fj = mrmfeature->getPrecursorFeature(precursor_ids[j]);
          intensityj.clear();
          fj->getIntensity(intensityj);
          Scoring::computeRank(intensityj, rank_vec2);
          // compute ranked mutual information
          mi_precursor_matrix_.setValue(i, j, Scoring::rankedMutualInformation(rank_vec1, rank_vec2));
        }
      }
    }

    void MRMScoring::initializeMIPrecursorContrastMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<String>& precursor_ids, const std::vector<String>& native_ids)
    {
      std::vector<double> intensityi, intensityj;
      mi_precursor_contrast_matrix_.resize(precursor_ids.size(), native_ids.size());
      std::vector<unsigned int> rank_vec1, rank_vec2;
      for (std::size_t i = 0; i < precursor_ids.size(); i++)
      {
        FeatureType fi = mrmfeature->getPrecursorFeature(precursor_ids[i]);
        //mi_precursor_contrast_matrix_[i].resize(native_ids.size());
        intensityi.clear();
        fi->getIntensity(intensityi);
        Scoring::computeRank(intensityi, rank_vec1);
        for (std::size_t j = 0; j < native_ids.size(); j++)
        {
          FeatureType fj = mrmfeature->getFeature(native_ids[j]);
          intensityj.clear();
          fj->getIntensity(intensityj);
          Scoring::computeRank(intensityj, rank_vec2);
          // compute ranked mutual information
          mi_precursor_contrast_matrix_.setValue(i, j, Scoring::rankedMutualInformation(rank_vec1, rank_vec2));
        }
      }
    }

    void MRMScoring::initializeMIPrecursorCombinedMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<String>& precursor_ids, const std::vector<String>& native_ids)
    {
      std::vector<double> intensityi, intensityj;
      std::vector<FeatureType> features;

      for (std::size_t i = 0; i < precursor_ids.size(); i++)
      {
        FeatureType fi = mrmfeature->getPrecursorFeature(precursor_ids[i]);
        features.push_back(fi);
      }
      for (std::size_t j = 0; j < native_ids.size(); j++)
      {
        FeatureType fj = mrmfeature->getFeature(native_ids[j]);
        features.push_back(fj);
      }
      std::vector<unsigned int> rank_vec1, rank_vec2;
      mi_precursor_combined_matrix_.resize(features.size(), features.size());
      for (std::size_t i = 0; i < features.size(); i++)
      {
        FeatureType fi = features[i];
        intensityi.clear();
        fi->getIntensity(intensityi);
        Scoring::computeRank(intensityi, rank_vec1);
        for (std::size_t j = 0; j < features.size(); j++)
        {
          FeatureType fj = features[j];
          intensityj.clear();
          fj->getIntensity(intensityj);
          Scoring::computeRank(intensityj, rank_vec2);
          // compute ranked mutual information
          mi_precursor_combined_matrix_.setValue(i ,j, Scoring::rankedMutualInformation(rank_vec1, rank_vec2));
        }
      }
    }

    double MRMScoring::calcMIScore()
    {
      OPENSWATH_PRECONDITION(mi_matrix_.rows() > 1, "Expect mutual information matrix of at least 2x2");

      double mi_scores{0};
      for(auto e : mi_matrix_)
      {
        mi_scores += e;
      }
      //mi_matrix_ is a triangular matrix
      size_t element_number = mi_matrix_.rows() * mi_matrix_.rows() / 2 + (mi_matrix_.rows() + 1) / 2;
      return mi_scores / element_number;
    }

    double MRMScoring::calcMIWeightedScore(
            const std::vector<double>& normalized_library_intensity)
    {
      OPENSWATH_PRECONDITION(mi_matrix_.rows() > 1, "Expect mutual information matrix of at least 2x2");

      double mi_scores{0};
      for (std::size_t i = 0; i < mi_matrix_.rows(); i++)
      {
        mi_scores += mi_matrix_.getValue(i,i)
                     * normalized_library_intensity[i]
                     * normalized_library_intensity[i];
#ifdef MRMSCORING_TESTING
        std::cout << "_mi_weighted " << i << " " << i << " " << mi_matrix_[i][i] << " weight " <<
        normalized_library_intensity[i] * normalized_library_intensity[i] << std::endl;
#endif
        for (std::size_t j = i + 1; j < mi_matrix_.rows(); j++)
        {
          mi_scores += mi_matrix_.getValue(i,j)
                       * normalized_library_intensity[i]
                       * normalized_library_intensity[j] * 2;
#ifdef MRMSCORING_TESTING
          std::cout << "_mi_weighted " << i << " " << j << " " << mi_matrix_[i][j] << " weight " <<
          normalized_library_intensity[i] * normalized_library_intensity[j] * 2 << std::endl;
#endif
        }
      }
      return mi_scores;
    }

    double MRMScoring::calcMIPrecursorScore()
    {
      OPENSWATH_PRECONDITION(mi_precursor_matrix_.rows() > 1, "Expect mutual information matrix of at least 2x2");

      double mi_scores{0};
      for(auto e : mi_precursor_matrix_)
      {
        mi_scores += e;
      }
      //mi_precursor_matrix_ is a triangular matrix
      size_t element_number = mi_precursor_matrix_.rows()*mi_precursor_matrix_.rows()/2 + (mi_precursor_matrix_.rows()+1)/2;
      return mi_scores / element_number;
    }

    double MRMScoring::calcMIPrecursorContrastScore()
    {
      OPENSWATH_PRECONDITION(mi_precursor_contrast_matrix_.rows() > 0 && mi_precursor_contrast_matrix_.cols() > 1, "Expect mutual information matrix of at least 1x2");

      double mi_scores{0};
      for(auto e : mi_precursor_contrast_matrix_)
      {
        mi_scores += e;
      }
      return mi_scores / mi_precursor_contrast_matrix_.size();
    }

    double MRMScoring::calcMIPrecursorCombinedScore()
    {
      OPENSWATH_PRECONDITION(mi_precursor_combined_matrix_.rows() > 1, "Expect mutual information matrix of at least 2x2");

      double mi_scores{0};

      for(auto e: mi_precursor_combined_matrix_)
      {
        mi_scores += e;
      }
      return mi_scores / mi_precursor_combined_matrix_.size();
    }

    std::vector<double> MRMScoring::calcSeparateMIContrastScore()
    {
      OPENSWATH_PRECONDITION(mi_contrast_matrix_.rows() > 0 && mi_contrast_matrix_.cols() > 1, "Expect mutual information matrix of at least 1x2");

      std::vector<double> mi_scores;
      mi_scores.resize(mi_contrast_matrix_.rows());
      for (std::size_t i = 0; i < mi_contrast_matrix_.rows(); i++)
      {
        double mi_scores_id = 0;
        for (std::size_t j = 0; j < mi_contrast_matrix_.cols(); j++)
        {
          mi_scores_id += mi_contrast_matrix_.getValue(i,j);
        }
        mi_scores[i] = mi_scores_id / mi_contrast_matrix_.cols();
      }
      return mi_scores;
    }
}