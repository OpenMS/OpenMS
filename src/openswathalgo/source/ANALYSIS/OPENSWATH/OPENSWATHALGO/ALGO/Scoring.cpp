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

#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/Scoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/Macros.h>
#include <cmath>

#include <boost/numeric/conversion/cast.hpp>

namespace OpenSwath
{
  namespace Scoring
  {

    void normalize_sum(double x[], unsigned int n)
    {
      double sumx = std::accumulate(&x[0], &x[0] + n, 0.0);
      if (sumx == 0.0)
      {
        return;
      } // do not divide by zero
      for (unsigned int i = 0; i < n; i++)
      {
        x[i] = x[i] / sumx;
      }
    }

    double NormalizedManhattanDist(double x[], double y[], int n)
    {
      OPENSWATH_PRECONDITION(n > 0, "Need at least one element");

      double delta_ratio_sum = 0;
      normalize_sum(x, n);
      normalize_sum(y, n);
      for (int i = 0; i < n; i++)
      {
        delta_ratio_sum += std::fabs(x[i] - y[i]);
      }
      return delta_ratio_sum / n;
    }

    double RootMeanSquareDeviation(double x[], double y[], int n)
    {
      OPENSWATH_PRECONDITION(n > 0, "Need at least one element");

      double result = 0;
      for (int i = 0; i < n; i++)
      {

        result += (x[i] - y[i]) * (x[i] - y[i]);
      }
      return std::sqrt(result / n);
    }

    double SpectralAngle(double x[], double y[], int n)
    {
      OPENSWATH_PRECONDITION(n > 0, "Need at least one element");

      double dotprod = 0;
      double x_len = 0;
      double y_len = 0;
      for (int i = 0; i < n; i++)
      {
        dotprod += x[i] * y[i];
        x_len += x[i] * x[i];
        y_len += y[i] * y[i];
      }
      x_len = std::sqrt(x_len);
      y_len = std::sqrt(y_len);

      return std::acos(dotprod / (x_len * y_len));
    }

    XCorrArrayType::const_iterator xcorrArrayGetMaxPeak(const XCorrArrayType& array)
    {
      OPENSWATH_PRECONDITION(array.data.size() > 0, "Cannot get highest apex from empty array.");

      XCorrArrayType::const_iterator max_it = array.begin();
      double max = array.begin()->second;
      for (XCorrArrayType::const_iterator it = array.begin(); it != array.end(); ++it)
      {
        if (it->second > max)
        {
          max = it->second;
          max_it = it;
        }
      }
      return max_it;
    }

    void standardize_data(std::vector<double>& data)
    {
      OPENSWATH_PRECONDITION(data.size() > 0, "Need non-empty array.");

      // subtract the mean and divide by the standard deviation
      double mean = std::accumulate(data.begin(), data.end(), 0.0) / (double) data.size();
      double sqsum = 0;
      for (std::vector<double>::iterator it = data.begin(); it != data.end(); ++it)
      {
        sqsum += (*it - mean) * (*it - mean);
      }
      double std = sqrt(sqsum / data.size()); // standard deviation

      for (std::size_t i = 0; i < data.size(); i++)
      {
        data[i] = (data[i] - mean) / std;
      }
    }

    XCorrArrayType normalizedCrossCorrelation(std::vector<double>& data1,
                                              std::vector<double>& data2, const int& maxdelay, const int& lag = 1)
    {
      OPENSWATH_PRECONDITION(data1.size() != 0 && data1.size() == data2.size(), "Both data vectors need to have the same length");

      // normalize the data
      standardize_data(data1);
      standardize_data(data2);
      XCorrArrayType result = calculateCrossCorrelation(data1, data2, maxdelay, lag);
      for (XCorrArrayType::iterator it = result.begin(); it != result.end(); ++it)
      {
        it->second = it->second / data1.size();
      }
      return result;
    }

    XCorrArrayType calculateCrossCorrelation(std::vector<double>& data1,
                                             std::vector<double>& data2, const int& maxdelay, const int& lag)
    {
      OPENSWATH_PRECONDITION(data1.size() != 0 && data1.size() == data2.size(), "Both data vectors need to have the same length");

      XCorrArrayType result;
      result.data.reserve( std::ceil((2*maxdelay + 1) / lag));
      int datasize = boost::numeric_cast<int>(data1.size());
      int i, j, delay;

      for (delay = -maxdelay; delay <= maxdelay; delay = delay + lag)
      {
        double sxy = 0;
        for (i = 0; i < datasize; ++i)
        {
          j = i + delay;
          if (j < 0 || j >= datasize)
          {
            continue;
          }
          sxy += (data1[i]) * (data2[j]);
        }
        result.data.push_back(std::make_pair(delay, sxy));
      }
      return result;
    }

    XCorrArrayType calcxcorr_legacy_mquest_(std::vector<double>& data1,
                                            std::vector<double>& data2, bool normalize)
    {
      OPENSWATH_PRECONDITION(!data1.empty() && data1.size() == data2.size(), "Both data vectors need to have the same length");
      int maxdelay = boost::numeric_cast<int>(data1.size());
      int lag = 1;

      double mean1 = std::accumulate(data1.begin(), data1.end(), 0.) / (double)data1.size();
      double mean2 = std::accumulate(data2.begin(), data2.end(), 0.) / (double)data2.size();
      double denominator = 1.0;
      int datasize = boost::numeric_cast<int>(data1.size());
      int i, j, delay;

      // Normalized cross-correlation = subtract the mean and divide by the standard deviation
      if (normalize)
      {
        double sqsum1 = 0;
        double sqsum2 = 0;
        for (std::vector<double>::iterator it = data1.begin(); it != data1.end(); ++it)
        {
          sqsum1 += (*it - mean1) * (*it - mean1);
        }

        for (std::vector<double>::iterator it = data2.begin(); it != data2.end(); ++it)
        {
          sqsum2 += (*it - mean2) * (*it - mean2);
        }
        // sigma_1 * sigma_2 * n
        denominator = sqrt(sqsum1 * sqsum2);
      }

      XCorrArrayType result;
      result.data.reserve( std::ceil((2*maxdelay + 1) / lag));
      int cnt = 0;
      for (delay = -maxdelay; delay <= maxdelay; delay = delay + lag, cnt++)
      {
        double sxy = 0;
        for (i = 0; i < datasize; i++)
        {
          j = i + delay;
          if (j < 0 || j >= datasize)
          {
            continue;
          }
          if (normalize)
          {
            sxy += (data1[i] - mean1) * (data2[j] - mean2);
          }
          else
          {
            sxy += (data1[i]) * (data2[j]);
          }
        }

        if (denominator > 0)
        {
          result.data.push_back(std::make_pair(delay, sxy/denominator));
        }
        else
        {
          // e.g. if all datapoints are zero
          result.data.push_back(std::make_pair(delay, 0));
        }
      }
      return result;
    }

  } //end namespace Scoring
}
