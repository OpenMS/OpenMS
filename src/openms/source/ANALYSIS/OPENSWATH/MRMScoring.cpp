// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/MRMScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/MRMFeatureAccessOpenMS.h>
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
  namespace
  {
    struct MRMScoringPool
    {
      std::vector<std::vector<double>> vv;
      // small 1-D temporaries reused across calls to avoid frequent heap allocs
      std::vector<double> tmp_double_a;
      std::vector<double> tmp_double_b;
      void ensure_size(std::size_t outer)
      {
        if (vv.size() < outer) vv.resize(outer);
        for (std::size_t i = 0; i < outer; ++i) vv[i].clear();
      }
      void ensure_1d_capacity(std::size_t n)
      {
        tmp_double_a.clear();
        tmp_double_b.clear();
        tmp_double_a.reserve(n);
        tmp_double_b.reserve(n);
      }
    };
    thread_local MRMScoringPool mrm_poolA;
    thread_local MRMScoringPool mrm_poolB;
  }
    // Maximum lag (in data points) for cross-correlation computation.
    // Lags beyond this range are physically meaningless for coelution scoring
    // and computing them wastes O(N) work per lag per pair.
    static constexpr int XCORR_MAX_DELAY = 10;

    /// Compute w^T * M_sym * w for an upper-triangular matrix.
    template <typename MatrixType>
    static double weightedTriangularSum(const MatrixType& mat, const std::vector<double>& w)
    {
      OPENSWATH_PRECONDITION(mat.rows() == mat.cols() && static_cast<std::size_t>(mat.rows()) == w.size(), "Matrix must be square and match weight vector size");
      double weighted_sum = 0.0;
      for (std::size_t i = 0; i < w.size(); ++i)
      {
        const double wi = w[i];
        weighted_sum += wi * wi * static_cast<double>(mat(i, i));
        for (std::size_t j = i + 1; j < w.size(); ++j)
        {
          weighted_sum += 2.0 * wi * w[j] * static_cast<double>(mat(i, j));
        }
      }
      return weighted_sum;
    }

    static double meanStdevUpperTriangle(const OpenMS::Matrix<int>& mat)
    {
      OpenSwath::mean_and_stddev msc;
      for (std::size_t i = 0; i < mat.rows(); ++i)
      {
        for (std::size_t j = i; j < mat.cols(); ++j)
        {
          msc(mat(i, j));
        }
      }
      return msc.mean() + msc.sample_stddev();
    }

    static double meanStdevAll(const OpenMS::Matrix<int>& mat)
    {
      OpenSwath::mean_and_stddev msc;
      const std::size_t n_entries = mat.size();
      const int* data = mat.data();
      for (std::size_t i = 0; i < n_entries; ++i)
      {
        msc(data[i]);
      }
      return msc.mean() + msc.sample_stddev();
    }

    static double meanUpperTriangle(const OpenMS::Matrix<double>& mat)
    {
      double sum = 0.0;
      std::size_t n_entries = 0;
      for (std::size_t i = 0; i < mat.rows(); ++i)
      {
        for (std::size_t j = i; j < mat.cols(); ++j)
        {
          sum += mat(i, j);
          ++n_entries;
        }
      }
      return n_entries > 0 ? sum / static_cast<double>(n_entries) : 0.0;
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

      // reuse thread-local pool to hold standardized copies (preserve capacity)
      mrm_poolA.ensure_size(data.size());
      std::vector<std::vector<double>>& tmp_data = mrm_poolA.vv;
      for (std::size_t i = 0; i < data.size(); i++)
      {
        tmp_data[i] = data[i];
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
      if (const auto* openms_feature = dynamic_cast<const OpenMS::MRMFeatureOpenMS*>(mrmfeature))
      {
        openms_feature->getFeatureIntensities(ids, intensity);
        return;
      }
      intensity.resize(ids.size());
      for (std::size_t i = 0; i < intensity.size(); i++)
      {
        MRMScoring::FeatureType fi = mrmfeature->getFeature(ids[i]);
        fi->getIntensity(intensity[i]);
      }
    }

    void fillIntensityFromPrecursorFeature(OpenSwath::IMRMFeature* mrmfeature, const std::vector<std::string>& ids, std::vector<std::vector<double>>& intensity)
    {
      if (const auto* openms_feature = dynamic_cast<const OpenMS::MRMFeatureOpenMS*>(mrmfeature))
      {
        openms_feature->getPrecursorFeatureIntensities(ids, intensity);
        return;
      }
      intensity.resize(ids.size());
      for (std::size_t i = 0; i < intensity.size(); i++)
      {
        MRMScoring::FeatureType fi = mrmfeature->getPrecursorFeature(ids[i]);
        fi->getIntensity(intensity[i]);
      }
    }

    void MRMScoring::initializeXCorrMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<std::string>& native_ids)
    {
      mrm_poolA.ensure_size(native_ids.size());
      std::vector<std::vector<double>>& intensity = mrm_poolA.vv;
      fillIntensityFromFeature(mrmfeature, native_ids, intensity);
      for (std::size_t i = 0; i < native_ids.size(); i++)
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
      mrm_poolA.ensure_size(native_ids_set1.size());
      std::vector<std::vector<double>>& intensityi = mrm_poolA.vv;
      mrm_poolB.ensure_size(native_ids_set2.size());
      std::vector<std::vector<double>>& intensityj = mrm_poolB.vv;
      fillIntensityFromFeature(mrmfeature, native_ids_set1, intensityi);
      for (std::size_t i = 0; i < native_ids_set1.size(); i++)
      {
        Scoring::standardize_data(intensityi[i]);
      }
      // ensure pool for intensityj
      mrm_poolB.ensure_size(native_ids_set2.size());
      fillIntensityFromFeature(mrmfeature, native_ids_set2, intensityj);
      for (std::size_t i = 0; i < native_ids_set2.size(); i++)
      {
        Scoring::standardize_data(intensityj[i]);
      }

      xcorr_contrast_matrix_.resize(native_ids_set1.size(), native_ids_set2.size());
      xcorr_contrast_matrix_max_peak_.resize(native_ids_set1.size(), native_ids_set2.size());
      xcorr_contrast_matrix_max_peak_sec_.resize(native_ids_set1.size(), native_ids_set2.size());
            
      for (std::size_t i = 0; i < native_ids_set1.size(); i++)
      {
        for (std::size_t j = 0; j < native_ids_set2.size(); j++)
        {
          // compute normalized cross correlation
          xcorr_contrast_matrix_(i, j) = Scoring::normalizedCrossCorrelationPost(intensityi[i], intensityj[j], std::min(XCORR_MAX_DELAY, static_cast<int>(intensityi[i].size())), 1);
          auto x = Scoring::xcorrArrayGetMaxPeak(xcorr_contrast_matrix_(i, j));
          xcorr_contrast_matrix_max_peak_(i, j) = std::abs(x->first);
          xcorr_contrast_matrix_max_peak_sec_(i, j) = x->second;
        }
      }
    }

    void MRMScoring::initializeXCorrPrecursorMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<std::string>& precursor_ids)
    {
      mrm_poolA.ensure_size(precursor_ids.size());
      std::vector<std::vector<double>>& intensity = mrm_poolA.vv;
      fillIntensityFromPrecursorFeature(mrmfeature, precursor_ids, intensity);
      for (std::size_t i = 0; i < precursor_ids.size(); i++)
      {
        Scoring::standardize_data(intensity[i]);
      }

      xcorr_precursor_matrix_.resize(precursor_ids.size(), precursor_ids.size());
      xcorr_precursor_matrix_max_peak_.resize(precursor_ids.size(), precursor_ids.size());
      xcorr_precursor_matrix_max_peak_sec_.resize(precursor_ids.size(), precursor_ids.size());
      for (std::size_t i = 0; i < precursor_ids.size(); i++)
      {
        for (std::size_t j = i; j < precursor_ids.size(); j++)
        {
          // compute normalized cross correlation
          xcorr_precursor_matrix_(i, j) = Scoring::normalizedCrossCorrelationPost(intensity[i], intensity[j], std::min(XCORR_MAX_DELAY, static_cast<int>(intensity[i].size())), 1);
          auto x = Scoring::xcorrArrayGetMaxPeak(xcorr_precursor_matrix_(i, j));
          xcorr_precursor_matrix_max_peak_(i, j) = std::abs(x->first);
          xcorr_precursor_matrix_max_peak_sec_(i, j) = x->second;
        }
      }
    }

    void MRMScoring::initializeXCorrPrecursorContrastMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<std::string>& precursor_ids, const std::vector<std::string>& native_ids)
    {
      mrm_poolA.ensure_size(precursor_ids.size());
      std::vector<std::vector<double>>& intensityi = mrm_poolA.vv;
      mrm_poolB.ensure_size(native_ids.size());
      std::vector<std::vector<double>>& intensityj = mrm_poolB.vv; // reuse pool
      fillIntensityFromPrecursorFeature(mrmfeature, precursor_ids, intensityi);
      for (std::size_t i = 0; i < precursor_ids.size(); i++)
      {
        Scoring::standardize_data(intensityi[i]);
      }
      mrm_poolB.ensure_size(native_ids.size());
      fillIntensityFromFeature(mrmfeature, native_ids, intensityj);
      for (std::size_t i = 0; i < native_ids.size(); i++)
      {
        Scoring::standardize_data(intensityj[i]);
      }

      xcorr_precursor_contrast_matrix_.resize(precursor_ids.size(), native_ids.size());
      xcorr_precursor_contrast_matrix_max_peak_.resize(precursor_ids.size(), native_ids.size());
      xcorr_precursor_contrast_matrix_max_peak_sec_.resize(precursor_ids.size(), native_ids.size());
      for (std::size_t i = 0; i < precursor_ids.size(); i++)
      {
        for (std::size_t j = 0; j < native_ids.size(); j++)
        {
          // compute normalized cross correlation
          xcorr_precursor_contrast_matrix_(i, j) = Scoring::normalizedCrossCorrelationPost(intensityi[i], intensityj[j], std::min(XCORR_MAX_DELAY, static_cast<int>(intensityi[i].size())), 1);
          auto x = Scoring::xcorrArrayGetMaxPeak(xcorr_precursor_contrast_matrix_(i, j));
          xcorr_precursor_contrast_matrix_max_peak_(i, j) = std::abs(x->first);
          xcorr_precursor_contrast_matrix_max_peak_sec_(i, j) = x->second;
        }
      }
    }

    void MRMScoring::initializeXCorrPrecursorContrastMatrix(const std::vector< std::vector< double > >& data_precursor, const std::vector< std::vector< double > >& data_fragments)
    {
      xcorr_precursor_contrast_matrix_.resize(data_precursor.size(), data_fragments.size());
      xcorr_precursor_contrast_matrix_max_peak_.resize(data_precursor.size(), data_fragments.size());
      xcorr_precursor_contrast_matrix_max_peak_sec_.resize(data_precursor.size(), data_fragments.size());
      // reuse thread-local pools to hold temp copies (preserve capacity)
      mrm_poolA.ensure_size(data_precursor.size());
      std::vector<std::vector<double>>& tmp_data_precursor = mrm_poolA.vv;
      for (std::size_t i = 0; i < data_precursor.size(); ++i)
      {
        tmp_data_precursor[i] = data_precursor[i];
        Scoring::standardize_data(tmp_data_precursor[i]);
      }
      mrm_poolB.ensure_size(data_fragments.size());
      std::vector<std::vector<double>>& tmp_data_fragments = mrm_poolB.vv;
      for (std::size_t i = 0; i < data_fragments.size(); ++i)
      {
        tmp_data_fragments[i] = data_fragments[i];
        Scoring::standardize_data(tmp_data_fragments[i]);
      }

      for (std::size_t i = 0; i < data_precursor.size(); i++)
      {
        for (std::size_t j = 0; j < data_fragments.size(); j++)
        {
          // compute normalized cross correlation
          xcorr_precursor_contrast_matrix_(i, j) = Scoring::normalizedCrossCorrelationPost(tmp_data_precursor[i], tmp_data_fragments[j], std::min(XCORR_MAX_DELAY, static_cast<int>(tmp_data_precursor[i].size())), 1);
          auto x = Scoring::xcorrArrayGetMaxPeak(xcorr_precursor_contrast_matrix_(i, j));
          xcorr_precursor_contrast_matrix_max_peak_(i, j) = std::abs(x->first);
          xcorr_precursor_contrast_matrix_max_peak_sec_(i, j) = x->second;
#ifdef MRMSCORING_TESTING
          std::cout << " fill xcorr_precursor_contrast_matrix_ "<< tmp_data_precursor[i].size() << " / " << tmp_data_fragments[j].size() << " : " << xcorr_precursor_contrast_matrix_[i][j].data.size() << '\n';
#endif
        }
      }
    }

    void MRMScoring::initializeXCorrPrecursorCombinedMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<std::string>& precursor_ids, const std::vector<std::string>& native_ids)
    {
      // use separate pools for precursor and fragment intensity buffers
      mrm_poolA.ensure_size(precursor_ids.size());
      std::vector<std::vector<double>>& intensityi = mrm_poolA.vv;
      fillIntensityFromPrecursorFeature(mrmfeature, precursor_ids, intensityi);
      mrm_poolB.ensure_size(native_ids.size());
      std::vector<std::vector<double>>& intensityj = mrm_poolB.vv;
      fillIntensityFromFeature(mrmfeature, native_ids, intensityj);
      // avoid allocating a combined vector by indexing into the two pooled buffers
      const std::size_t ni = intensityi.size();
      const std::size_t nj = intensityj.size();
      const std::size_t combined_size = ni + nj;

      // standardize in-place for both pools
      for (std::size_t i = 0; i < ni; ++i) Scoring::standardize_data(intensityi[i]);
      for (std::size_t j = 0; j < nj; ++j) Scoring::standardize_data(intensityj[j]);

      xcorr_precursor_combined_matrix_.resize(combined_size, combined_size);
      xcorr_precursor_combined_matrix_max_peak_.resize(combined_size, combined_size);
      xcorr_precursor_combined_matrix_max_peak_sec_.resize(combined_size, combined_size);
      for (std::size_t ii = 0; ii < combined_size; ++ii)
      {
        std::vector<double>* pvi = (ii < ni) ? &intensityi[ii] : &intensityj[ii - ni];
        std::vector<double>& vi = *pvi;
        for (std::size_t jj = ii; jj < combined_size; ++jj)
        {
          std::vector<double>* pvj = (jj < ni) ? &intensityi[jj] : &intensityj[jj - ni];
          std::vector<double>& vj = *pvj;
          xcorr_precursor_combined_matrix_(ii, jj) = Scoring::normalizedCrossCorrelationPost(vi, vj, std::min(XCORR_MAX_DELAY, static_cast<int>(vi.size())), 1);
          auto x = Scoring::xcorrArrayGetMaxPeak(xcorr_precursor_combined_matrix_(ii, jj));
          xcorr_precursor_combined_matrix_max_peak_(ii, jj) = std::abs(x->first);
          xcorr_precursor_combined_matrix_max_peak_sec_(ii, jj) = x->second;
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
      return meanStdevUpperTriangle(xcorr_matrix_max_peak_);
    }

    double MRMScoring::calcXcorrCoelutionWeightedScore(
            const std::vector<double>& normalized_library_intensity)
    {
      OPENSWATH_PRECONDITION(xcorr_matrix_max_peak_.rows() > 1, "Expect cross-correlation matrix of at least 2x2");
      return weightedTriangularSum(xcorr_matrix_max_peak_, normalized_library_intensity);
    }

    std::vector<double> MRMScoring::calcSeparateXcorrContrastCoelutionScore()
    {
      OPENSWATH_PRECONDITION(xcorr_contrast_matrix_max_peak_.rows() > 0 && xcorr_contrast_matrix_max_peak_.cols() > 1, "Expect cross-correlation matrix of at least 1x2");

      std::vector<double > deltas;
      for (std::size_t i = 0; i < xcorr_contrast_matrix_max_peak_.rows(); ++i)
      {
        double deltas_id = 0;
        for (std::size_t j = 0; j < xcorr_contrast_matrix_max_peak_.cols(); ++j)
        {
          deltas_id += xcorr_contrast_matrix_max_peak_(i, j);
#ifdef MRMSCORING_TESTING
          std::cout << "&&_xcoel append " << xcorr_contrast_matrix_max_peak_getValue(i, j) << '\n';
#endif
        }
        deltas.push_back(deltas_id / xcorr_contrast_matrix_max_peak_.cols());
      }

      return deltas;
    }

    double MRMScoring::calcXcorrPrecursorCoelutionScore()
    {
      OPENSWATH_PRECONDITION(xcorr_precursor_matrix_max_peak_.rows() > 1, "Expect cross-correlation matrix of at least 2x2");
      return meanStdevUpperTriangle(xcorr_precursor_matrix_max_peak_);
    }

    double MRMScoring::calcXcorrPrecursorContrastCoelutionScore()
    {
      OPENSWATH_PRECONDITION(xcorr_precursor_contrast_matrix_max_peak_.rows() > 0 && xcorr_precursor_contrast_matrix_max_peak_.cols() > 1, "Expect cross-correlation matrix of at least 1x2");
      return meanStdevAll(xcorr_precursor_contrast_matrix_max_peak_);
    }

    double MRMScoring::calcXcorrPrecursorContrastSumFragCoelutionScore()
    {
      OPENSWATH_PRECONDITION(xcorr_precursor_contrast_matrix_max_peak_.rows() > 0 && xcorr_precursor_contrast_matrix_max_peak_.cols() > 0, "Expect cross-correlation matrix of at least 1x1");
      return meanStdevAll(xcorr_precursor_contrast_matrix_max_peak_);
    }

    double MRMScoring::calcXcorrPrecursorCombinedCoelutionScore()
    {
      OPENSWATH_PRECONDITION(xcorr_precursor_combined_matrix_max_peak_.rows() > 1, "Expect cross-correlation matrix of at least 2x2");
      return meanStdevUpperTriangle(xcorr_precursor_combined_matrix_max_peak_);
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
      return meanUpperTriangle(xcorr_matrix_max_peak_sec_);
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
      OPENSWATH_PRECONDITION(xcorr_precursor_matrix_max_peak_sec_.rows() > 1, "Expect cross-correlation matrix of at least 2x2");
      return meanUpperTriangle(xcorr_precursor_matrix_max_peak_sec_);
    }

    double MRMScoring::calcXcorrPrecursorContrastSumFragShapeScore()
    {
      OPENSWATH_PRECONDITION(xcorr_precursor_contrast_matrix_max_peak_sec_.rows() > 0 && xcorr_precursor_contrast_matrix_max_peak_sec_.cols() > 0, "Expect cross-correlation matrix of at least 1x1");
      auto em = OpenMS::eigenView(xcorr_precursor_contrast_matrix_max_peak_sec_);
      return em.sum() / static_cast<double>(em.size());
    }

    double MRMScoring::calcXcorrPrecursorContrastShapeScore()
    {
      OPENSWATH_PRECONDITION(xcorr_precursor_contrast_matrix_max_peak_sec_.rows() > 0 && xcorr_precursor_contrast_matrix_max_peak_sec_.cols() > 1, "Expect cross-correlation matrix of at least 1x2");
      auto em = OpenMS::eigenView(xcorr_precursor_contrast_matrix_max_peak_sec_);
      return em.sum() / static_cast<double>(em.size());
    }

    double MRMScoring::calcXcorrPrecursorCombinedShapeScore()
    {
      OPENSWATH_PRECONDITION(xcorr_precursor_combined_matrix_max_peak_sec_.rows() > 1, "Expect cross-correlation matrix of at least 2x2");
      return meanUpperTriangle(xcorr_precursor_combined_matrix_max_peak_sec_);
    }

    void MRMScoring::calcLibraryScore(OpenSwath::IMRMFeature* mrmfeature, const std::vector<TransitionType>& transitions,
                                      double& correlation, double& norm_manhattan, double& manhattan, double& dotprod, double& spectral_angle, double& rmsd)
    {
      // reuse small 1-D temporaries from thread-local pool to avoid allocations
      mrm_poolA.ensure_1d_capacity(transitions.size());
      std::vector<double>& library_intensity = mrm_poolA.tmp_double_a;
      std::vector<double>& experimental_intensity = mrm_poolA.tmp_double_b;
      // ensure requested capacity
      library_intensity.clear();
      experimental_intensity.clear();
      library_intensity.reserve(transitions.size());
      experimental_intensity.reserve(transitions.size());
      const auto* openms_feature = dynamic_cast<const OpenMS::MRMFeatureOpenMS*>(mrmfeature);

      for (std::size_t k = 0; k < transitions.size(); k++)
      {
        double intensity = transitions[k].getLibraryIntensity();
        // the library intensity should never be below zero
        if (intensity < 0.0)
        {
          intensity = 0.0;
        }
        if (openms_feature != nullptr)
        {
          const std::string native_id = transitions[k].getNativeID();
          experimental_intensity.push_back(static_cast<double>(openms_feature->getFeatureIntensity(native_id, k)));
        }
        else
        {
          const std::string native_id = transitions[k].getNativeID();
          experimental_intensity.push_back(static_cast<double>(mrmfeature->getFeature(native_id)->getIntensity()));
        }
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

      const double feature_rt = mrmfeature->getRT();
      for (std::size_t k = 0; k < signal_noise_estimators.size(); k++)
      {
        sn_score += signal_noise_estimators[k]->getValueAtRT(feature_rt);
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
      sn_scores.reserve(signal_noise_estimators.size());

      const double feature_rt = mrmfeature->getRT();
      for (std::size_t k = 0; k < signal_noise_estimators.size(); k++)
      {
        const double sn_score = signal_noise_estimators[k]->getValueAtRT(feature_rt);
        if (sn_score < 1)
          // everything below S/N 1 can be set to zero (and the log safely applied)
        {
          sn_scores.push_back(0);
        }
        else
        {
          sn_scores.push_back(std::log(sn_score));
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
      mrm_poolA.ensure_size(native_ids.size());
      std::vector<std::vector<double>>& intensity = mrm_poolA.vv;
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
      mrm_poolA.ensure_size(native_ids_set1.size());
      mrm_poolB.ensure_size(native_ids_set2.size());
      std::vector<std::vector<double>>& intensityi = mrm_poolA.vv;
      std::vector<std::vector<double>>& intensityj = mrm_poolB.vv;
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
      mrm_poolA.ensure_size(precursor_ids.size());
      std::vector<std::vector<double>>& intensity = mrm_poolA.vv;
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
      mrm_poolA.ensure_size(precursor_ids.size());
      mrm_poolB.ensure_size(native_ids.size());
      std::vector<std::vector<double>>& intensityi = mrm_poolA.vv;
      std::vector<std::vector<double>>& intensityj = mrm_poolB.vv;
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
      // reuse pool sequentially for precursor then fragment
      mrm_poolA.ensure_size(precursor_ids.size());
      std::vector<std::vector<double>>& intensity = mrm_poolA.vv;
      fillIntensityFromPrecursorFeature(mrmfeature, precursor_ids, intensity);
      std::vector<unsigned int> max_rank_vec = Scoring::computeRankVector(intensity, rank_vec);
      // now reuse same buffer for native_ids
      intensity.clear();
      mrm_poolA.ensure_size(native_ids.size());
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
