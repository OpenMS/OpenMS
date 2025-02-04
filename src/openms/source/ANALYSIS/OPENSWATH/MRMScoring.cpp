// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
#include <cmath> // for isnan

namespace OpenSwath
{

    const MRMScoring::XCorrMatrixType& MRMScoring::getXCorrMatrix() const
    {
      return xcorr_matrix_;
    }

    void MRMScoring::initializeXCorrMatrix(const std::vector< std::vector< double > >& data)
    {
      xcorr_matrix_.getEigenMatrix().resize(data.size(), data.size());
      xcorr_matrix_max_peak_.getEigenMatrix().resize(data.size(), data.size());
      xcorr_matrix_max_peak_sec_.getEigenMatrix().resize(data.size(), data.size());

      std::vector< std::vector< double > > tmp_data = data;
      for (std::size_t i = 0; i < tmp_data.size(); i++)
      {
        Scoring::standardize_data(tmp_data[i]);
      }

      for (std::size_t i = 0; i < data.size(); i++)
      {
        for (std::size_t j = i; j < data.size(); j++)
        {
          // compute normalized cross correlation
          xcorr_matrix_(i, j) = Scoring::normalizedCrossCorrelationPost(tmp_data[i], tmp_data[j], static_cast<int>(data[i].size()), 1);
          auto x = Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_(i, j));
          xcorr_matrix_max_peak_(i, j) = std::abs(x->first);
          xcorr_matrix_max_peak_sec_(i, j) = x->second;
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

    void fillIntensityFromFeature(OpenSwath::IMRMFeature* mrmfeature, const std::vector<std::string>& ids, std::vector<std::vector<double>>& intensity)
    {
      intensity.resize(ids.size());
      for (std::size_t i = 0; i < intensity.size(); i++)
      {
        MRMScoring::FeatureType fi = mrmfeature->getFeature(ids[i]);
        fi->getIntensity(intensity[i]);
      }
    }

    void fillIntensityFromPrecursorFeature(OpenSwath::IMRMFeature* mrmfeature, const std::vector<std::string>& ids, std::vector<std::vector<double>>& intensity)
    {
      intensity.resize(ids.size());
      for (std::size_t i = 0; i < intensity.size(); i++)
      {
        MRMScoring::FeatureType fi = mrmfeature->getPrecursorFeature(ids[i]);
        fi->getIntensity(intensity[i]);
      }
    }

    void MRMScoring::initializeXCorrMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<std::string>& native_ids)
    {
      std::vector<std::vector<double>> intensity;
      fillIntensityFromFeature(mrmfeature, native_ids, intensity);
      for (std::size_t i = 0; i < intensity.size(); i++)
      {
        Scoring::standardize_data(intensity[i]);
      }

      xcorr_matrix_.getEigenMatrix().resize(native_ids.size(), native_ids.size());
      xcorr_matrix_max_peak_.getEigenMatrix().resize(native_ids.size(), native_ids.size());
      xcorr_matrix_max_peak_sec_.getEigenMatrix().resize(native_ids.size(), native_ids.size());

      for (std::size_t i = 0; i < native_ids.size(); i++)
      {
        for (std::size_t j = i; j < native_ids.size(); j++)
        {
          // compute normalized cross correlation
          xcorr_matrix_(i, j) = Scoring::normalizedCrossCorrelationPost(intensity[i], intensity[j], static_cast<int>(intensity[i].size()), 1);
          auto x = Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_(i, j));
          xcorr_matrix_max_peak_(i, j) = std::abs(x->first);
          xcorr_matrix_max_peak_sec_(i, j) = x->second;
        }
      }
    }

    void MRMScoring::initializeXCorrContrastMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<std::string>& native_ids_set1, const std::vector<std::string>& native_ids_set2)
    {
      std::vector<std::vector<double>> intensityi, intensityj;
      fillIntensityFromFeature(mrmfeature, native_ids_set1, intensityi);
      for (std::size_t i = 0; i < intensityi.size(); i++)
      {
        Scoring::standardize_data(intensityi[i]);
      }
      fillIntensityFromFeature(mrmfeature, native_ids_set2, intensityj);
      for (std::size_t i = 0; i < intensityj.size(); i++)
      {
        Scoring::standardize_data(intensityj[i]);
      }

      xcorr_contrast_matrix_.getEigenMatrix().resize(native_ids_set1.size(), native_ids_set2.size());
      xcorr_contrast_matrix_max_peak_sec_.getEigenMatrix().resize(native_ids_set1.size(), native_ids_set2.size());
            
      for (std::size_t i = 0; i < native_ids_set1.size(); i++)
      {
        for (std::size_t j = 0; j < native_ids_set2.size(); j++)
        {
          // compute normalized cross correlation
          xcorr_contrast_matrix_(i, j) = Scoring::normalizedCrossCorrelationPost(intensityi[i], intensityj[j], static_cast<int>(intensityi[i].size()), 1);
          auto x = Scoring::xcorrArrayGetMaxPeak(xcorr_contrast_matrix_(i, j));
          xcorr_contrast_matrix_max_peak_sec_(i, j) = x->second;
        }
      }
    }

    void MRMScoring::initializeXCorrPrecursorMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<std::string>& precursor_ids)
    {
      std::vector<std::vector<double>> intensity;
      fillIntensityFromPrecursorFeature(mrmfeature, precursor_ids, intensity);
      for (std::size_t i = 0; i < intensity.size(); i++)
      {
        Scoring::standardize_data(intensity[i]);
      }

      xcorr_precursor_matrix_.getEigenMatrix().resize(precursor_ids.size(), precursor_ids.size());
      for (std::size_t i = 0; i < precursor_ids.size(); i++)
      {
        for (std::size_t j = i; j < precursor_ids.size(); j++)
        {
          // compute normalized cross correlation
          xcorr_precursor_matrix_(i, j) = Scoring::normalizedCrossCorrelationPost(intensity[i], intensity[j], static_cast<int>(intensity[i].size()), 1);
        }
      }
    }

    void MRMScoring::initializeXCorrPrecursorContrastMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<std::string>& precursor_ids, const std::vector<std::string>& native_ids)
    {
      std::vector<std::vector<double>> intensityi, intensityj;
      fillIntensityFromPrecursorFeature(mrmfeature, precursor_ids, intensityi);
      for (std::size_t i = 0; i < intensityi.size(); i++)
      {
        Scoring::standardize_data(intensityi[i]);
      }
      fillIntensityFromFeature(mrmfeature, native_ids, intensityj);
      for (std::size_t i = 0; i < intensityj.size(); i++)
      {
        Scoring::standardize_data(intensityj[i]);
      }

      xcorr_precursor_contrast_matrix_.getEigenMatrix().resize(precursor_ids.size(), native_ids.size());
      for (std::size_t i = 0; i < precursor_ids.size(); i++)
      {
        for (std::size_t j = 0; j < native_ids.size(); j++)
        {
          // compute normalized cross correlation
          xcorr_precursor_contrast_matrix_(i, j) = Scoring::normalizedCrossCorrelationPost(intensityi[i], intensityj[j], static_cast<int>(intensityi[i].size()), 1);
        }
      }
    }

    void MRMScoring::initializeXCorrPrecursorContrastMatrix(const std::vector< std::vector< double > >& data_precursor, const std::vector< std::vector< double > >& data_fragments)
    {
      xcorr_precursor_contrast_matrix_.getEigenMatrix().resize(data_precursor.size(), data_fragments.size());
      std::vector< std::vector< double > > tmp_data_precursor = data_precursor;
      std::vector< std::vector< double > > tmp_data_fragments = data_fragments;
      for (std::size_t i = 0; i < tmp_data_precursor.size(); i++)
      {
        Scoring::standardize_data(tmp_data_precursor[i]);
      }
      for (std::size_t i = 0; i < tmp_data_fragments.size(); i++)
      {
        Scoring::standardize_data(tmp_data_fragments[i]);
      }

      for (std::size_t i = 0; i < data_precursor.size(); i++)
      {
        for (std::size_t j = 0; j < data_fragments.size(); j++)
        {
          // compute normalized cross correlation
          xcorr_precursor_contrast_matrix_(i, j) = Scoring::normalizedCrossCorrelationPost(tmp_data_precursor[i], tmp_data_fragments[j], static_cast<int>(tmp_data_precursor[i].size()), 1);
#ifdef MRMSCORING_TESTING
          std::cout << " fill xcorr_precursor_contrast_matrix_ "<< tmp_data_precursor[i].size() << " / " << tmp_data_fragments[j].size() << " : " << xcorr_precursor_contrast_matrix_[i][j].data.size() << std::endl;
#endif
        }
      }
    }

    void MRMScoring::initializeXCorrPrecursorCombinedMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<std::string>& precursor_ids, const std::vector<std::string>& native_ids)
    {
      std::vector<std::vector<double>> intensityi, intensityj;
      fillIntensityFromPrecursorFeature(mrmfeature, precursor_ids, intensityi);
      fillIntensityFromFeature(mrmfeature, native_ids, intensityj);
      std::vector<std::vector<double>> combined_intensity;
      for (std::size_t i = 0; i < intensityi.size(); i++)
      {
        combined_intensity.push_back(intensityi[i]);
      }
      for (std::size_t j = 0; j < intensityj.size(); j++)
      {
        combined_intensity.push_back(intensityj[j]);
      }
      for (std::size_t i = 0; i < combined_intensity.size(); i++)
      {
        Scoring::standardize_data(combined_intensity[i]);
      }

      xcorr_precursor_combined_matrix_.getEigenMatrix().resize(combined_intensity.size(), combined_intensity.size());
      for (std::size_t i = 0; i < combined_intensity.size(); i++)
      {
        for (std::size_t j = i; j < combined_intensity.size(); j++)
        {
          // compute normalized cross correlation
          xcorr_precursor_combined_matrix_(i, j) = Scoring::normalizedCrossCorrelationPost(combined_intensity[i], combined_intensity[j], static_cast<int>(combined_intensity[i].size()), 1);
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
    
      OpenSwath::mean_and_stddev msc;
      for (long int i = 0; i < xcorr_matrix_max_peak_.rows(); i++)
      {
        for (long int  j = i; j < xcorr_matrix_max_peak_.rows(); j++)
        {
          // first is the X value (RT), should be an int
          //deltas.push_back(std::abs(Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_.getValue(i, j))->first));
          msc(xcorr_matrix_max_peak_(i,j));
#ifdef MRMSCORING_TESTING
          std::cout << "&&_xcoel append " << std::abs(Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][j])->first) << std::endl;
#endif
        }
      }

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
      for (long int i = 0; i < xcorr_matrix_max_peak_.rows(); i++)
      {
        deltas += (xcorr_matrix_max_peak_(i, i)//std::abs(Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_.getValue(i, i))->first)
                   * normalized_library_intensity[i]
                   * normalized_library_intensity[i]);
#ifdef MRMSCORING_TESTING
        std::cout << "_xcoel_weighted " << i << " " << i << " " << Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][i])->first << " weight " <<
        normalized_library_intensity[i] * normalized_library_intensity[i] << std::endl;
      weights += normalized_library_intensity[i] * normalized_library_intensity[i];
#endif
        for (long int j = i + 1; j < xcorr_matrix_max_peak_.rows(); j++)
        {
          // first is the X value (RT), should be an int
          deltas += (xcorr_matrix_max_peak_(i, j)//std::abs(Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_.getValue(i, j))->first)
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

    std::vector<double> MRMScoring::calcSeparateXcorrContrastCoelutionScore()
    {
      OPENSWATH_PRECONDITION(xcorr_contrast_matrix_.rows() > 0 && xcorr_contrast_matrix_.cols() > 1, "Expect cross-correlation matrix of at least 1x2");

      std::vector<double > deltas;
      for (long int i = 0; i < xcorr_contrast_matrix_.rows(); i++)
      {
        double deltas_id = 0;
        for (long int  j = 0; j < xcorr_contrast_matrix_.cols(); j++)
        {
          // first is the X value (RT), should be an int
          auto x = Scoring::xcorrArrayGetMaxPeak(xcorr_contrast_matrix_(i, j));
          deltas_id += std::abs(x->first);
#ifdef MRMSCORING_TESTING
          std::cout << "&&_xcoel append " << xcorr_contrast_matrix_max_peak_getValue(i, j) << std::endl;
#endif
        }
        deltas.push_back(deltas_id / xcorr_contrast_matrix_.cols());
      }

      return deltas;
    }

    double MRMScoring::calcXcorrPrecursorCoelutionScore()
    {
      OPENSWATH_PRECONDITION(xcorr_precursor_matrix_.rows() > 1, "Expect cross-correlation matrix of at least 2x2");

      OpenSwath::mean_and_stddev msc;
      for (long int i = 0; i < xcorr_precursor_matrix_.rows(); i++)
      {
        for (long int  j = i; j < xcorr_precursor_matrix_.rows(); j++)
        {
          // first is the X value (RT), should be an int
          auto x = Scoring::xcorrArrayGetMaxPeak(xcorr_precursor_matrix_(i, j));
          msc(std::abs(x->first));
#ifdef MRMSCORING_TESTING
          std::cout << "&&_xcoel append " << std::abs(Scoring::xcorrArrayGetMaxPeak(xcorr_precursor_matrix_[i][j])->first) << std::endl;
#endif
        }
      }

      double deltas_mean = msc.mean();
      double deltas_stdv = msc.sample_stddev();

      double xcorr_coelution_score = deltas_mean + deltas_stdv;
      return xcorr_coelution_score;
    }

    double MRMScoring::calcXcorrPrecursorContrastCoelutionScore()
    {
      OPENSWATH_PRECONDITION(xcorr_precursor_contrast_matrix_.rows() > 0 && xcorr_precursor_contrast_matrix_.cols() > 1, "Expect cross-correlation matrix of at least 1x2");

      OpenSwath::mean_and_stddev msc;
      size_t n_entries = xcorr_precursor_contrast_matrix_.getEigenMatrix().size();
      const auto& em = xcorr_precursor_contrast_matrix_.getEigenMatrix();

      for (size_t i = 0; i < n_entries; i++)
      {
        // first is the X value (RT), should be an int
        auto e = *(em.data() + i);
        msc(std::abs(Scoring::xcorrArrayGetMaxPeak(e)->first));
#ifdef MRMSCORING_TESTING
        std::cout << "&&_xcoel append " << std::abs(Scoring::xcorrArrayGetMaxPeak(xcorr_precursor_contrast_matrix_[i][j])->first) << std::endl;
#endif
      }

      double deltas_mean = msc.mean();
      double deltas_stdv = msc.sample_stddev();

      double xcorr_coelution_score = deltas_mean + deltas_stdv;
      return xcorr_coelution_score;
    }

    double MRMScoring::calcXcorrPrecursorContrastSumFragCoelutionScore()
    {
      OPENSWATH_PRECONDITION(xcorr_precursor_contrast_matrix_.rows() > 0 && xcorr_precursor_contrast_matrix_.cols() > 0, "Expect cross-correlation matrix of at least 1x1");

      OpenSwath::mean_and_stddev msc;
      size_t n_entries = xcorr_precursor_contrast_matrix_.getEigenMatrix().size();
      auto& em = xcorr_precursor_contrast_matrix_.getEigenMatrix();
      for (size_t i = 0; i < n_entries; i++)
      {
        // first is the X value (RT), should be an int
        auto e = *(em.data() + i);
        msc(std::abs(Scoring::xcorrArrayGetMaxPeak(e)->first));

#ifdef MRMSCORING_TESTING
        std::cout << "&&_xcoel append " << std::abs(Scoring::xcorrArrayGetMaxPeak(xcorr_precursor_contrast_matrix_[i][j])->first) << std::endl;
#endif
      }

      double deltas_mean = msc.mean();
      double deltas_stdv = msc.sample_stddev();

      double xcorr_coelution_score = deltas_mean + deltas_stdv;
      return xcorr_coelution_score;
    }

    double MRMScoring::calcXcorrPrecursorCombinedCoelutionScore()
    {
      OPENSWATH_PRECONDITION(xcorr_precursor_combined_matrix_.rows() > 1, "Expect cross-correlation matrix of at least 2x2");

      OpenSwath::mean_and_stddev msc;
      for (long int i = 0; i < xcorr_precursor_combined_matrix_.rows(); i++)
      {
        for (long int  j = i; j < xcorr_precursor_combined_matrix_.rows(); j++)
        {
          // first is the X value (RT), should be an int
          auto x = Scoring::xcorrArrayGetMaxPeak(xcorr_precursor_combined_matrix_(i, j));
          msc(std::abs(x->first));
#ifdef MRMSCORING_TESTING
          std::cout << "&&_xcoel append " << std::abs(Scoring::xcorrArrayGetMaxPeak(xcorr_precursor_combined_matrix_[i][j])->first) << std::endl;
#endif
        }
      }
      
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
      for (long int i = 0; i < xcorr_matrix_max_peak_sec_.rows(); i++)
      {
        for (long int j = i; j < xcorr_matrix_max_peak_sec_.rows(); j++)
        {
          // second is the Y value (intensity)
          intensities += xcorr_matrix_max_peak_sec_(i, j);
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
      for (long int i = 0; i < xcorr_matrix_max_peak_sec_.rows(); i++)
      {
        intensities += (xcorr_matrix_max_peak_sec_(i, i)
                        * normalized_library_intensity[i]
                        * normalized_library_intensity[i]);
#ifdef MRMSCORING_TESTING
        std::cout << "_xcorr_weighted " << i << " " << i << " " << Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][i])->second << " weight " <<
        normalized_library_intensity[i] * normalized_library_intensity[i] << std::endl;
#endif
        for (long int j = i + 1; j < xcorr_matrix_max_peak_sec_.rows(); j++)
        {
          intensities += (xcorr_matrix_max_peak_sec_(i, j)
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
      const auto& em = xcorr_contrast_matrix_max_peak_sec_.getEigenMatrix();
      return em.sum();
    }

    std::vector<double> MRMScoring::calcSeparateXcorrContrastShapeScore()
    {
      OPENSWATH_PRECONDITION(xcorr_contrast_matrix_max_peak_sec_.rows() > 0 && xcorr_contrast_matrix_max_peak_sec_.cols() > 1, "Expect cross-correlation matrix of at least 1x2");

      std::vector<double> intensities;
      for (long int i = 0; i < xcorr_contrast_matrix_max_peak_sec_.rows(); i++)
      {
        double intensities_id = 0;
        for (long int j = 0; j < xcorr_contrast_matrix_max_peak_sec_.cols(); j++)
        {
          // second is the Y value (intensity)
          intensities_id += xcorr_contrast_matrix_max_peak_sec_(i,j);
        }
        intensities.push_back(intensities_id / xcorr_contrast_matrix_max_peak_sec_.cols());
      }

      return intensities;
    }

    double MRMScoring::calcXcorrPrecursorShapeScore()
    {
      OPENSWATH_PRECONDITION(xcorr_precursor_matrix_.rows() > 1, "Expect cross-correlation matrix of at least 2x2");

      double intensities{0};
      for(long int i = 0; i < xcorr_precursor_matrix_.rows(); i++)
      {
        for(long int j = i; j < xcorr_precursor_matrix_.cols(); j++)
        {
          auto x = Scoring::xcorrArrayGetMaxPeak(xcorr_precursor_matrix_(i, j));
          intensities += x->second;
        }
      }
      //xcorr_precursor_matrix_ is a triangle matrix
      size_t element_number = xcorr_precursor_matrix_.rows()*xcorr_precursor_matrix_.rows()/2 + (xcorr_precursor_matrix_.rows()+1)/2;
      return intensities / element_number;
    }

    double MRMScoring::calcXcorrPrecursorContrastSumFragShapeScore()
    {
      OPENSWATH_PRECONDITION(xcorr_precursor_contrast_matrix_.rows() > 0 && xcorr_precursor_contrast_matrix_.cols() > 0, "Expect cross-correlation matrix of at least 1x1");

      double intensities{0};
      const auto& em = xcorr_precursor_contrast_matrix_.getEigenMatrix();
      size_t n_elements = em.size();
      for (size_t i = 0; i != n_elements; ++i)
      {
        const auto& e = *(em.data() + i); 
        intensities += Scoring::xcorrArrayGetMaxPeak(e)->second;;
      }

      return intensities / (double)n_elements;
    }

    double MRMScoring::calcXcorrPrecursorContrastShapeScore()
    {
      OPENSWATH_PRECONDITION(xcorr_precursor_contrast_matrix_.rows() > 0 && xcorr_precursor_contrast_matrix_.cols() > 1, "Expect cross-correlation matrix of at least 1x2");

      double intensities{0};

      const auto& em = xcorr_precursor_contrast_matrix_.getEigenMatrix();
      size_t n_elements = em.size();
      for (size_t i = 0; i != n_elements; ++i)
      {
        const auto& e = *(em.data() + i); 
        intensities += Scoring::xcorrArrayGetMaxPeak(e)->second;
      }
      return intensities / (double)n_elements;      
    }

    double MRMScoring::calcXcorrPrecursorCombinedShapeScore()
    {
      OPENSWATH_PRECONDITION(xcorr_precursor_combined_matrix_.rows() > 1, "Expect cross-correlation matrix of at least 2x2");

      double intensities{0};
      for(long int i = 0; i < xcorr_precursor_combined_matrix_.rows(); i++)
      {
        for(long int j = i; j < xcorr_precursor_combined_matrix_.cols(); j++)
        {
          auto x = Scoring::xcorrArrayGetMaxPeak(xcorr_precursor_combined_matrix_(i, j));
          intensities += x->second;
        }
      }
      //xcorr_precursor-combined_matrix_ is a triangle matrix
      size_t element_number = xcorr_precursor_combined_matrix_.rows()*xcorr_precursor_combined_matrix_.rows()/2 + (xcorr_precursor_combined_matrix_.rows()+1)/2;
      return intensities / element_number;
    }

    void MRMScoring::calcLibraryScore(OpenSwath::IMRMFeature* mrmfeature, const std::vector<TransitionType>& transitions,
                                      double& correlation, double& norm_manhattan, double& manhattan, double& dotprod, double& spectral_angle, double& rmsd)
    {
      std::vector<double> library_intensity;
      std::vector<double> experimental_intensity;
      std::string native_id;

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

      spectral_angle = Scoring::SpectralAngle(&experimental_intensity[0], &library_intensity[0], static_cast<unsigned int>(transitions.size()));

      if (std::isnan(spectral_angle))
      {
        spectral_angle = 0.0;
      }

      Scoring::normalize_sum(&experimental_intensity[0], static_cast<unsigned int>(transitions.size()));
      Scoring::normalize_sum(&library_intensity[0], static_cast<unsigned int>(transitions.size()));

      norm_manhattan = Scoring::NormalizedManhattanDist(&experimental_intensity[0], &library_intensity[0], static_cast<unsigned int>(transitions.size()));
      rmsd = Scoring::RootMeanSquareDeviation(&experimental_intensity[0], &library_intensity[0], static_cast<unsigned int>(transitions.size()));
      correlation = OpenSwath::cor_pearson(experimental_intensity.begin(), experimental_intensity.end(), library_intensity.begin());

      if (std::isnan(correlation))
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
      if (signal_noise_estimators.empty())
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
      if (signal_noise_estimators.empty())
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

    void MRMScoring::initializeMIMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<std::string>& native_ids)
    {
      std::vector<std::vector<double>> intensity;
      std::vector<std::vector<unsigned int>> rank_vec{};
      fillIntensityFromFeature(mrmfeature, native_ids, intensity);
      std::vector<unsigned int> max_rank_vec = Scoring::computeRankVector(intensity, rank_vec);

      mi_matrix_.getEigenMatrix().resize(native_ids.size(), native_ids.size());
      mi_matrix_.getEigenMatrix().setZero(); 
      for (std::size_t i = 0; i < native_ids.size(); i++)
      {
        for (std::size_t j = i; j < native_ids.size(); j++)
        {
          // compute ranked mutual information
          mi_matrix_(i, j) = Scoring::rankedMutualInformation(rank_vec[i], rank_vec[j], max_rank_vec[i], max_rank_vec[j]);
        }
      }
    }

    void MRMScoring::initializeMIContrastMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<std::string>& native_ids_set1, const std::vector<std::string>& native_ids_set2)
    { 
      std::vector<std::vector<double>> intensityi, intensityj;
      std::vector<std::vector<unsigned int>> rank_vec1{}, rank_vec2{};
      fillIntensityFromFeature(mrmfeature, native_ids_set1, intensityi);
      fillIntensityFromFeature(mrmfeature, native_ids_set2, intensityj);
      std::vector<unsigned int> max_rank_vec1 = Scoring::computeRankVector(intensityi, rank_vec1);
      std::vector<unsigned int> max_rank_vec2 = Scoring::computeRankVector(intensityj, rank_vec2);

      mi_contrast_matrix_.getEigenMatrix().resize(native_ids_set1.size(), native_ids_set2.size());
      for (std::size_t i = 0; i < native_ids_set1.size(); i++)
      {
        for (std::size_t j = 0; j < native_ids_set2.size(); j++)
        {
          // compute ranked mutual information
          mi_contrast_matrix_(i, j) = Scoring::rankedMutualInformation(rank_vec1[i], rank_vec2[j], max_rank_vec1[i], max_rank_vec2[j]);
        }
      }
    }

    void MRMScoring::initializeMIPrecursorMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<std::string>& precursor_ids)
    {
      std::vector<std::vector<double>> intensity;
      std::vector<std::vector<unsigned int>> rank_vec;
      fillIntensityFromPrecursorFeature(mrmfeature, precursor_ids, intensity);
      std::vector<unsigned int> max_rank_vec = Scoring::computeRankVector(intensity, rank_vec);

      mi_precursor_matrix_.getEigenMatrix().resize(precursor_ids.size(), precursor_ids.size());
      mi_precursor_matrix_.getEigenMatrix().setZero();

      for (std::size_t i = 0; i < precursor_ids.size(); i++)
      {
        for (std::size_t j = i; j < precursor_ids.size(); j++)
        {
          // compute ranked mutual information
          mi_precursor_matrix_(i, j) = Scoring::rankedMutualInformation(rank_vec[i], rank_vec[j], max_rank_vec[i], max_rank_vec[j]);
        }
      }
    }

    void MRMScoring::initializeMIPrecursorContrastMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<std::string>& precursor_ids, const std::vector<std::string>& native_ids)
    {
      std::vector<std::vector<double>> intensityi, intensityj;
      std::vector<std::vector<unsigned int>> rank_vec1{}, rank_vec2{};
      fillIntensityFromPrecursorFeature(mrmfeature, precursor_ids, intensityi);
      fillIntensityFromFeature(mrmfeature, native_ids, intensityj);
      std::vector<unsigned int> max_rank_vec1 = Scoring::computeRankVector(intensityi, rank_vec1);
      std::vector<unsigned int> max_rank_vec2 = Scoring::computeRankVector(intensityj, rank_vec2);

      mi_precursor_contrast_matrix_.getEigenMatrix().resize(precursor_ids.size(), native_ids.size());
      for (std::size_t i = 0; i < precursor_ids.size(); i++)
      {
        for (std::size_t j = 0; j < native_ids.size(); j++)
        {
          // compute ranked mutual information
          mi_precursor_contrast_matrix_(i, j) = Scoring::rankedMutualInformation(rank_vec1[i], rank_vec2[j], max_rank_vec1[i], max_rank_vec2[j]);
        }
      }
    }

    void MRMScoring::initializeMIPrecursorCombinedMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<std::string>& precursor_ids, const std::vector<std::string>& native_ids)
    {
      std::vector<std::vector<unsigned int>> rank_vec{};
      std::vector<std::vector<double>> intensity;
      fillIntensityFromPrecursorFeature(mrmfeature, precursor_ids, intensity);
      std::vector<unsigned int> max_rank_vec = Scoring::computeRankVector(intensity, rank_vec);
      intensity.clear();
      fillIntensityFromFeature(mrmfeature, native_ids, intensity);
      std::vector<unsigned int> max_rank_vec_tmp = Scoring::computeRankVector(intensity, rank_vec);
      max_rank_vec.reserve(max_rank_vec.size() + native_ids.size());
      max_rank_vec.insert(max_rank_vec.end(), max_rank_vec_tmp.begin(), max_rank_vec_tmp.end());
      
      mi_precursor_combined_matrix_.getEigenMatrix().resize(rank_vec.size(), rank_vec.size());
      for (std::size_t i = 0; i < rank_vec.size(); i++)
      { 
        for (std::size_t j = i; j < rank_vec.size(); j++)
        {
          // compute ranked mutual information
          double curr_mutual_score = Scoring::rankedMutualInformation(rank_vec[i], rank_vec[j], max_rank_vec[i], max_rank_vec[j]);
          mi_precursor_combined_matrix_(i, j) = curr_mutual_score;
          if (i != j) mi_precursor_combined_matrix_(j, i) = curr_mutual_score;
        }
      }
    }
    
    double MRMScoring::calcMIScore()
    {
      OPENSWATH_PRECONDITION(mi_matrix_.rows() > 1, "Expect mutual information matrix of at least 2x2");
      const auto& em = mi_matrix_.getEigenMatrix();
      double mi_scores = em.sum();
      //mi_matrix_ is a triangular matrix
      size_t element_number = mi_matrix_.rows() * mi_matrix_.rows() / 2 + (mi_matrix_.rows() + 1) / 2;
      return mi_scores / element_number;
    }

    double MRMScoring::calcMIWeightedScore(
            const std::vector<double>& normalized_library_intensity)
    {
      OPENSWATH_PRECONDITION(mi_matrix_.rows() > 1, "Expect mutual information matrix of at least 2x2");

      double mi_scores{0};
      for (long int i = 0; i < mi_matrix_.rows(); i++)
      {
        mi_scores += mi_matrix_(i, i)
                     * normalized_library_intensity[i]
                     * normalized_library_intensity[i];
#ifdef MRMSCORING_TESTING
        std::cout << "_mi_weighted " << i << " " << i << " " << mi_matrix_[i][i] << " weight " <<
        normalized_library_intensity[i] * normalized_library_intensity[i] << std::endl;
#endif
        for (long int j = i + 1; j < mi_matrix_.rows(); j++)
        {
          mi_scores += mi_matrix_(i, j)
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

      const auto& em = mi_precursor_matrix_.getEigenMatrix();
      double mi_scores = em.sum();
      //mi_precursor_matrix_ is a triangular matrix
      size_t element_number = mi_precursor_matrix_.rows()*mi_precursor_matrix_.rows()/2 + (mi_precursor_matrix_.rows()+1)/2;
      return mi_scores / (double)element_number;
    }

    double MRMScoring::calcMIPrecursorContrastScore()
    {
      OPENSWATH_PRECONDITION(mi_precursor_contrast_matrix_.rows() > 0 && mi_precursor_contrast_matrix_.cols() > 1, "Expect mutual information matrix of at least 1x2");

      const auto& em = mi_precursor_contrast_matrix_.getEigenMatrix();
      size_t n_entries = em.size();
      double mi_scores = em.sum();

      return mi_scores / (double)n_entries;
    }

    double MRMScoring::calcMIPrecursorCombinedScore()
    {
      OPENSWATH_PRECONDITION(mi_precursor_combined_matrix_.rows() > 1, "Expect mutual information matrix of at least 2x2");

      const auto& em = mi_precursor_combined_matrix_.getEigenMatrix();
      size_t n_entries = em.size();

      double mi_scores = em.sum();

      return mi_scores / (double)n_entries;
    }

    std::vector<double> MRMScoring::calcSeparateMIContrastScore()
    {
      OPENSWATH_PRECONDITION(mi_contrast_matrix_.rows() > 0 && mi_contrast_matrix_.cols() > 1, "Expect mutual information matrix of at least 1x2");

      std::vector<double> mi_scores;
      mi_scores.resize(mi_contrast_matrix_.rows());      
      for (long int i = 0; i < mi_contrast_matrix_.rows(); i++)
      {
        double mi_scores_id = 0;
        for (long int j = 0; j < mi_contrast_matrix_.cols(); j++)
        {
          mi_scores_id += mi_contrast_matrix_(i, j);
        }
        mi_scores[i] = mi_scores_id / mi_contrast_matrix_.cols();
      }
      return mi_scores;
    }
}
