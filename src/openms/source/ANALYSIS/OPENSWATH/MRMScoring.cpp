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
#include <OpenMS/DATASTRUCTURES/MatrixEigen.h>
#include <OpenMS/OPENSWATHALGO/Macros.h>
#include <OpenMS/MATH/StatisticFunctions.h>
//#define MRMSCORING_TESTING
#include <algorithm>
#include <iostream>
#include <iterator>
#include <cmath> // for isnan

namespace OpenSwath
{
    // Maximum lag (in data points) for cross-correlation computation.
    // Lags beyond this range are physically meaningless for coelution scoring
    // and computing them wastes O(N) work per lag per pair.
    static constexpr int XCORR_MAX_DELAY = 10;

    /// Compute w^T * M_sym * w for an upper-triangular Matrix<double>.
    /// Uses selfadjointView to treat the upper triangle as symmetric without copying.
    static double weightedTriangularSum(const OpenMS::Matrix<double>& mat, const std::vector<double>& w)
    {
      OPENSWATH_PRECONDITION(mat.rows() == mat.cols() && static_cast<size_t>(mat.rows()) == w.size(), "Matrix must be square and match weight vector size");
      auto em = OpenMS::eigenView(mat);
      Eigen::Map<const Eigen::VectorXd> ew(w.data(), w.size());
      return ew.dot(em.selfadjointView<Eigen::Upper>() * ew);
    }

    /// Overload for Matrix<int> (used by xcorr_matrix_max_peak_).
    static double weightedTriangularSum(const OpenMS::Matrix<int>& mat, const std::vector<double>& w)
    {
      OPENSWATH_PRECONDITION(mat.rows() == mat.cols() && static_cast<size_t>(mat.rows()) == w.size(), "Matrix must be square and match weight vector size");
      auto em = OpenMS::eigenView(mat);
      Eigen::MatrixXd sym = em.cast<double>();
      Eigen::Map<const Eigen::VectorXd> ew(w.data(), w.size());
      return ew.dot(sym.selfadjointView<Eigen::Upper>() * ew);
    }

    const MRMScoring::XCorrMatrixType& MRMScoring::getXCorrMatrix() const
    {
      return xcorr_matrix_;
    }

    void MRMScoring::initializeXCorrMatrix(const std::vector< std::vector< double > >& data)
    {
      xcorr_matrix_.resize(data.size(), data.size());
      xcorr_matrix_max_peak_.resize(data.size(), data.size());
      xcorr_matrix_max_peak_sec_.resize(data.size(), data.size());

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
          xcorr_matrix_(i, j) = Scoring::normalizedCrossCorrelationPost(tmp_data[i], tmp_data[j], std::min(XCORR_MAX_DELAY, static_cast<int>(data[i].size())), 1);
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

      xcorr_matrix_.resize(native_ids.size(), native_ids.size());
      xcorr_matrix_max_peak_.resize(native_ids.size(), native_ids.size());
      xcorr_matrix_max_peak_sec_.resize(native_ids.size(), native_ids.size());

      for (std::size_t i = 0; i < native_ids.size(); i++)
      {
        for (std::size_t j = i; j < native_ids.size(); j++)
        {
          // compute normalized cross correlation
          xcorr_matrix_(i, j) = Scoring::normalizedCrossCorrelationPost(intensity[i], intensity[j], std::min(XCORR_MAX_DELAY, static_cast<int>(intensity[i].size())), 1);
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

      xcorr_contrast_matrix_.resize(native_ids_set1.size(), native_ids_set2.size());
      xcorr_contrast_matrix_max_peak_sec_.resize(native_ids_set1.size(), native_ids_set2.size());
            
      for (std::size_t i = 0; i < native_ids_set1.size(); i++)
      {
        for (std::size_t j = 0; j < native_ids_set2.size(); j++)
        {
          // compute normalized cross correlation
          xcorr_contrast_matrix_(i, j) = Scoring::normalizedCrossCorrelationPost(intensityi[i], intensityj[j], std::min(XCORR_MAX_DELAY, static_cast<int>(intensityi[i].size())), 1);
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

      xcorr_precursor_matrix_.resize(precursor_ids.size(), precursor_ids.size());
      for (std::size_t i = 0; i < precursor_ids.size(); i++)
      {
        for (std::size_t j = i; j < precursor_ids.size(); j++)
        {
          // compute normalized cross correlation
          xcorr_precursor_matrix_(i, j) = Scoring::normalizedCrossCorrelationPost(intensity[i], intensity[j], std::min(XCORR_MAX_DELAY, static_cast<int>(intensity[i].size())), 1);
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

      xcorr_precursor_contrast_matrix_.resize(precursor_ids.size(), native_ids.size());
      for (std::size_t i = 0; i < precursor_ids.size(); i++)
      {
        for (std::size_t j = 0; j < native_ids.size(); j++)
        {
          // compute normalized cross correlation
          xcorr_precursor_contrast_matrix_(i, j) = Scoring::normalizedCrossCorrelationPost(intensityi[i], intensityj[j], std::min(XCORR_MAX_DELAY, static_cast<int>(intensityi[i].size())), 1);
        }
      }
    }

    void MRMScoring::initializeXCorrPrecursorContrastMatrix(const std::vector< std::vector< double > >& data_precursor, const std::vector< std::vector< double > >& data_fragments)
    {
      xcorr_precursor_contrast_matrix_.resize(data_precursor.size(), data_fragments.size());
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
          xcorr_precursor_contrast_matrix_(i, j) = Scoring::normalizedCrossCorrelationPost(tmp_data_precursor[i], tmp_data_fragments[j], std::min(XCORR_MAX_DELAY, static_cast<int>(tmp_data_precursor[i].size())), 1);
#ifdef MRMSCORING_TESTING
          std::cout << " fill xcorr_precursor_contrast_matrix_ "<< tmp_data_precursor[i].size() << " / " << tmp_data_fragments[j].size() << " : " << xcorr_precursor_contrast_matrix_[i][j].data.size() << '\n';
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

      xcorr_precursor_combined_matrix_.resize(combined_intensity.size(), combined_intensity.size());
      for (std::size_t i = 0; i < combined_intensity.size(); i++)
      {
        for (std::size_t j = i; j < combined_intensity.size(); j++)
        {
          // compute normalized cross correlation
          xcorr_precursor_combined_matrix_(i, j) = Scoring::normalizedCrossCorrelationPost(combined_intensity[i], combined_intensity[j], std::min(XCORR_MAX_DELAY, static_cast<int>(combined_intensity[i].size())), 1);
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
          std::cout << "&&_xcoel append " << std::abs(Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][j])->first) << '\n';
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
      return weightedTriangularSum(xcorr_matrix_max_peak_, normalized_library_intensity);
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
          std::cout << "&&_xcoel append " << xcorr_contrast_matrix_max_peak_getValue(i, j) << '\n';
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
          std::cout << "&&_xcoel append " << std::abs(Scoring::xcorrArrayGetMaxPeak(xcorr_precursor_matrix_[i][j])->first) << '\n';
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
      size_t n_entries = xcorr_precursor_contrast_matrix_.size();
      auto em = OpenMS::eigenView(xcorr_precursor_contrast_matrix_);

      for (size_t i = 0; i < n_entries; i++)
      {
        // first is the X value (RT), should be an int
        auto e = *(em.data() + i);
        msc(std::abs(Scoring::xcorrArrayGetMaxPeak(e)->first));
#ifdef MRMSCORING_TESTING
        std::cout << "&&_xcoel append " << std::abs(Scoring::xcorrArrayGetMaxPeak(xcorr_precursor_contrast_matrix_[i][j])->first) << '\n';
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
      size_t n_entries = xcorr_precursor_contrast_matrix_.size();
      auto em = OpenMS::eigenView(xcorr_precursor_contrast_matrix_);
      for (size_t i = 0; i < n_entries; i++)
      {
        // first is the X value (RT), should be an int
        auto e = *(em.data() + i);
        msc(std::abs(Scoring::xcorrArrayGetMaxPeak(e)->first));

#ifdef MRMSCORING_TESTING
        std::cout << "&&_xcoel append " << std::abs(Scoring::xcorrArrayGetMaxPeak(xcorr_precursor_contrast_matrix_[i][j])->first) << '\n';
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
          std::cout << "&&_xcoel append " << std::abs(Scoring::xcorrArrayGetMaxPeak(xcorr_precursor_combined_matrix_[i][j])->first) << '\n';
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

      // xcorr_matrix_max_peak_sec_ is upper-triangular (lower triangle is zero from resize())
      auto em = OpenMS::eigenView(xcorr_matrix_max_peak_sec_);
      double intensities = em.sum();
      size_t n = xcorr_matrix_max_peak_sec_.rows();
      size_t element_number = n * n / 2 + (n + 1) / 2;
      return intensities / element_number;
    }

    double MRMScoring::calcXcorrShapeWeightedScore(
            const std::vector<double>& normalized_library_intensity)
    {
      OPENSWATH_PRECONDITION(xcorr_matrix_max_peak_sec_.rows() > 1, "Expect cross-correlation matrix of at least 2x2");
      return weightedTriangularSum(xcorr_matrix_max_peak_sec_, normalized_library_intensity);
    }

    double MRMScoring::calcXcorrContrastShapeScore()
    {
      OPENSWATH_PRECONDITION(xcorr_contrast_matrix_max_peak_sec_.rows() > 0 && xcorr_contrast_matrix_max_peak_sec_.cols() > 1, "Expect cross-correlation matrix of at least 1x2");
      auto em = OpenMS::eigenView(xcorr_contrast_matrix_max_peak_sec_);
      return em.sum();
    }

    std::vector<double> MRMScoring::calcSeparateXcorrContrastShapeScore()
    {
      OPENSWATH_PRECONDITION(xcorr_contrast_matrix_max_peak_sec_.rows() > 0 && xcorr_contrast_matrix_max_peak_sec_.cols() > 1, "Expect cross-correlation matrix of at least 1x2");

      auto em = OpenMS::eigenView(xcorr_contrast_matrix_max_peak_sec_);
      Eigen::VectorXd row_means = em.rowwise().sum() / em.cols();
      return std::vector<double>(row_means.data(), row_means.data() + row_means.size());
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
      auto em = OpenMS::eigenView(xcorr_precursor_contrast_matrix_);
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

      auto em = OpenMS::eigenView(xcorr_precursor_contrast_matrix_);
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
      std::cout << native_id << " Lib vs exp " << library_intensity[k] << " " << experimental_intensity[k] << '\n';
    }
#endif

      manhattan = OpenSwath::manhattanScoring(experimental_intensity, library_intensity);
      dotprod = OpenSwath::dotprodScoring(experimental_intensity, library_intensity);

      spectral_angle = Scoring::SpectralAngle(experimental_intensity, library_intensity);

      if (std::isnan(spectral_angle))
      {
        spectral_angle = 0.0;
      }

      Scoring::normalize_sum(experimental_intensity);
      Scoring::normalize_sum(library_intensity);

      norm_manhattan = Scoring::NormalizedManhattanDist(experimental_intensity, library_intensity);
      rmsd = Scoring::RootMeanSquareDeviation(experimental_intensity, library_intensity);
      correlation = OpenMS::Math::pearsonCorrelationCoefficient(experimental_intensity.begin(), experimental_intensity.end(), library_intensity.begin(), library_intensity.end());

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

      mi_matrix_.resize(native_ids.size(), native_ids.size());
      mi_matrix_.fill(0.0); 
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

      mi_contrast_matrix_.resize(native_ids_set1.size(), native_ids_set2.size());
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

      mi_precursor_matrix_.resize(precursor_ids.size(), precursor_ids.size());
      mi_precursor_matrix_.fill(0.0);

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

      mi_precursor_contrast_matrix_.resize(precursor_ids.size(), native_ids.size());
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
      
      mi_precursor_combined_matrix_.resize(rank_vec.size(), rank_vec.size());
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
      auto em = OpenMS::eigenView(mi_matrix_);
      double mi_scores = em.sum();
      //mi_matrix_ is a triangular matrix
      size_t element_number = mi_matrix_.rows() * mi_matrix_.rows() / 2 + (mi_matrix_.rows() + 1) / 2;
      return mi_scores / element_number;
    }

    double MRMScoring::calcMIWeightedScore(
            const std::vector<double>& normalized_library_intensity)
    {
      OPENSWATH_PRECONDITION(mi_matrix_.rows() > 1, "Expect mutual information matrix of at least 2x2");
      return weightedTriangularSum(mi_matrix_, normalized_library_intensity);
    }

    double MRMScoring::calcMIPrecursorScore()
    {
      OPENSWATH_PRECONDITION(mi_precursor_matrix_.rows() > 1, "Expect mutual information matrix of at least 2x2");

      auto em = OpenMS::eigenView(mi_precursor_matrix_);
      double mi_scores = em.sum();
      //mi_precursor_matrix_ is a triangular matrix
      size_t element_number = mi_precursor_matrix_.rows()*mi_precursor_matrix_.rows()/2 + (mi_precursor_matrix_.rows()+1)/2;
      return mi_scores / (double)element_number;
    }

    double MRMScoring::calcMIPrecursorContrastScore()
    {
      OPENSWATH_PRECONDITION(mi_precursor_contrast_matrix_.rows() > 0 && mi_precursor_contrast_matrix_.cols() > 1, "Expect mutual information matrix of at least 1x2");

      auto em = OpenMS::eigenView(mi_precursor_contrast_matrix_);
      size_t n_entries = em.size();
      double mi_scores = em.sum();

      return mi_scores / (double)n_entries;
    }

    double MRMScoring::calcMIPrecursorCombinedScore()
    {
      OPENSWATH_PRECONDITION(mi_precursor_combined_matrix_.rows() > 1, "Expect mutual information matrix of at least 2x2");

      auto em = OpenMS::eigenView(mi_precursor_combined_matrix_);
      size_t n_entries = em.size();

      double mi_scores = em.sum();

      return mi_scores / (double)n_entries;
    }

    std::vector<double> MRMScoring::calcSeparateMIContrastScore()
    {
      OPENSWATH_PRECONDITION(mi_contrast_matrix_.rows() > 0 && mi_contrast_matrix_.cols() > 1, "Expect mutual information matrix of at least 1x2");

      auto em = OpenMS::eigenView(mi_contrast_matrix_);
      Eigen::VectorXd row_means = em.rowwise().sum() / em.cols();
      return std::vector<double>(row_means.data(), row_means.data() + row_means.size());
    }
}
