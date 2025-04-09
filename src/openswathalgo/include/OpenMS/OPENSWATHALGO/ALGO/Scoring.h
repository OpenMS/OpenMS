// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest$
// $Authors: Hannes Roest$
// --------------------------------------------------------------------------

#pragma once

#include <numeric>
#include <map>
#include <vector>

#include <OpenMS/OPENSWATHALGO/OpenSwathAlgoConfig.h>

namespace OpenSwath
{

  /**
    @brief Scoring functions used by MRMScoring

    Many helper functions to calculate cross-correlations between data
  */
  namespace Scoring
  {
    /** @name Type defs and helper structures*/


    //@{
    typedef std::pair<unsigned int, unsigned int> pos2D;
    /// Simple hash function for Scoring::pos2D
    struct pair_hash 
    {
      template <class T1, class T2>
      std::size_t operator () (const std::pair<T1,T2> &p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ h2;  
      }
    };
    
    /// Cross Correlation array contains (lag,correlation) pairs
    typedef std::pair<int, double> XCorrEntry;
    struct XCorrArrayType 
    {
public:
    std::vector<XCorrEntry> data;

    // Access functions
    typedef std::vector<XCorrEntry>::iterator iterator;
    typedef std::vector<XCorrEntry>::const_iterator const_iterator;

    iterator begin() {return data.begin();}
    const_iterator begin() const {return data.begin();}
    iterator end() {return data.end();}
    const_iterator end() const {return data.end();}
    };
    //@}

    /** @name Helper functions */
    //@{
    /** @brief Calculate the normalized Manhattan distance between two arrays
     *
     * Equivalent to the function "delta_ratio_sum" from mQuest to calculate
     * similarity between library intensity and experimental ones.
     *
     * The delta_ratio_sum is calculated as follows:
     
       @f[
       d = \sqrt{\frac{1}{N}  \sum_{i=0}^N |\frac{x_i}{\mu_x} - \frac{y_i}{\mu_y}|) }
       @f]
    */
    OPENSWATHALGO_DLLAPI double NormalizedManhattanDist(double x[], double y[], int n);

    /** @brief Calculate the RMSD (root means square deviation)
     *
     * The RMSD is calculated as follows:
     
       @f[
       RMSD = \sqrt{\frac{1}{N}  \sum_{i=0}^N (x_i - y_i)^2 } 
       @f]
    */
    OPENSWATHALGO_DLLAPI double RootMeanSquareDeviation(double x[], double y[], int n);

    /** @brief Calculate the Spectral angle (acosine of the normalized dotproduct)
     *
     * The spectral angle is calculated as follows:
     
       @f[
       \theta = acos \left( \frac{\sum_{i=0}^N (x_i * y_i))}{\sqrt{\sum_{i=0}^N (x_i * x_i) \sum_{i=0}^N (y_i * y_i)} }  \right)
       @f]
    */
    OPENSWATHALGO_DLLAPI double SpectralAngle(double x[], double y[], int n);

    /// Calculate crosscorrelation on std::vector data - Deprecated!
    /// Legacy code, this is a 1:1 port of the function from mQuest
    OPENSWATHALGO_DLLAPI XCorrArrayType calcxcorr_legacy_mquest_(std::vector<double>& data1,
                                                  std::vector<double>& data2, bool normalize);

    /// Calculate crosscorrelation on std::vector data (which is first normalized)
    /// NOTE: this replaces calcxcorr 
    OPENSWATHALGO_DLLAPI XCorrArrayType normalizedCrossCorrelation(std::vector<double>& data1,
                                                                   std::vector<double>& data2, const int maxdelay, const int lag);

    /// Calculate crosscorrelation on std::vector data that is already normalized
    OPENSWATHALGO_DLLAPI XCorrArrayType normalizedCrossCorrelationPost(std::vector<double>& normalized_data1,
                                                                       std::vector<double>& normalized_data2, const int maxdelay, const int lag);                                                                   

    /// Calculate crosscorrelation on std::vector data without normalization
    OPENSWATHALGO_DLLAPI XCorrArrayType calculateCrossCorrelation(const std::vector<double>& data1,
                                                                  const std::vector<double>& data2, const int maxdelay, const int lag);

    /// Find best peak in an cross-correlation (highest apex)
    OPENSWATHALGO_DLLAPI XCorrArrayType::const_iterator xcorrArrayGetMaxPeak(const XCorrArrayType & array);

    /// Standardize a vector (subtract mean, divide by standard deviation)
    OPENSWATHALGO_DLLAPI void standardize_data(std::vector<double>& data);

    /// Divide each element of x by the sum of the vector
    OPENSWATHALGO_DLLAPI void normalize_sum(double x[], unsigned int n);

    // Compute rank of vector elements, append it to @p ranks and return the highest rank
    OPENSWATHALGO_DLLAPI unsigned int computeAndAppendRank(const std::vector<double>& v, std::vector<unsigned int>& ranks);

    // Compute rank of vector elements and its highest rank for each row in a 2D array
    OPENSWATHALGO_DLLAPI std::vector<unsigned int> computeRankVector(const std::vector<std::vector<double>>& intensity, std::vector<std::vector<unsigned int>>& ranks);

    // Estimate mutual information between two vectors of ranks
    OPENSWATHALGO_DLLAPI double rankedMutualInformation(std::vector<unsigned int>& ranked_data1, std::vector<unsigned int>& ranked_data2, const unsigned int max_rank1, const unsigned int max_rank2);

    //@}

  }
}

