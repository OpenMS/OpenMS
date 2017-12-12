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
// $Maintainer: Hannes Roest$
// $Authors: Hannes Roest$
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_OPENSWATH_OPENSWATHALGO_ALGO_SCORING_H
#define OPENMS_ANALYSIS_OPENSWATH_OPENSWATHALGO_ALGO_SCORING_H

#include <numeric>
#include <map>
#include <vector>

#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/OpenSwathAlgoConfig.h>

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
                                                                   std::vector<double>& data2, const int& maxdelay, const int& lag);

    /// Calculate crosscorrelation on std::vector data without normalization
    OPENSWATHALGO_DLLAPI XCorrArrayType calculateCrossCorrelation(std::vector<double>& data1,
                                                                  std::vector<double>& data2, const int& maxdelay, const int& lag);

    /// Find best peak in an cross-correlation (highest apex)
    OPENSWATHALGO_DLLAPI XCorrArrayType::const_iterator xcorrArrayGetMaxPeak(const XCorrArrayType & array);

    /// Standardize a vector (subtract mean, divide by standard deviation)
    OPENSWATHALGO_DLLAPI void standardize_data(std::vector<double>& data);

    /// divide each element of x by the sum of the vector
    OPENSWATHALGO_DLLAPI void normalize_sum(double x[], unsigned int n);
    //@}

  }
}

#endif // OPENMS_ANALYSIS_OPENSWATH_OPENSWATHALGO_ALGO_SCORING_H
