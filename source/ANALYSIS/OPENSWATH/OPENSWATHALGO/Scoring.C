// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/Scoring.h"
#include <cmath>

namespace OpenMS
{
  namespace Scoring
  {

    void normalize_sum(double x[], int n)
    {
      double sumx = std::accumulate(&x[0], &x[0] + n, 0.0);
      if (sumx == 0.0)
      {
        return;
      } // do not divide by zero
      for (int i = 0; i < n; i++)
      {
        x[i] = x[i] / sumx;
      }
    }

    double RMSD(double x[], double y[], int n)
    {
      //OPENMS_PRECONDITION(n > 0, "Need at least one element");

      double delta_ratio_sum = 0;
      normalize_sum(x, n);
      normalize_sum(y, n);
      for (int i = 0; i < n; i++)
      {
        delta_ratio_sum += std::fabs(x[i] - y[i]);
      }
      return delta_ratio_sum / n;
    }

    XCorrArrayType::iterator xcorrArrayGetMaxPeak(XCorrArrayType array)
    {
      //OPENMS_PRECONDITION(array.size() > 0, "Cannot get highest apex from empty array.");

      XCorrArrayType::iterator max_it = array.begin();
      double max = array.begin()->second;
      for (XCorrArrayType::iterator it = array.begin();
          it != array.end(); it++)
      {
        if (it->second > max)
        {
          max = it->second;
          max_it = it;
        }
      }
      return max_it;
    }

    void standardize_data(std::vector<double> & data)
    {
      // OPENMS_PRECONDITION(data.size() > 0, "Need non-empty array.");

      // subtract the mean and divide by the standard deviation
      double mean = std::accumulate(data.begin(), data.end(), 0.0)
          / (double) data.size();
      double sqsum = 0;
      for (std::vector<double>::iterator it = data.begin(); it != data.end();
          it++)
      {
        sqsum += (*it - mean) * (*it - mean);
      }
      double std = sqrt(sqsum / data.size()); // standard deviation

      for (std::size_t i = 0; i < data.size(); i++)
      {
        data[i] = (data[i] - mean) / std;
      }
    }

    XCorrArrayType normalizedCalcxcorr(std::vector<double> & data1,
        std::vector<double> & data2, int maxdelay, int lag = 1)
    {
      //OPENMS_PRECONDITION(data1.size() != 0 && data1.size() == data2.size(), "Both data vectors need to have the same length");

      // normalize the data
      standardize_data(data1);
      standardize_data(data2);
      std::map<int, double> result = calcxcorr_new(data1, data2, maxdelay, lag);
      for (std::map<int, double>::iterator it = result.begin();
          it != result.end(); it++)
      {
        it->second = it->second / data1.size();
      }
      return result;
    }

    XCorrArrayType calcxcorr_new(std::vector<double> & data1,
        std::vector<double> & data2, int maxdelay, int lag)
    {
      // OPENMS_PRECONDITION(data1.size() != 0 && data1.size() == data2.size(), "Both data vectors need to have the same length");

      XCorrArrayType result;
      int datasize = data1.size();
      int i, j, delay;
      double sxy;

      for (delay = -maxdelay; delay <= maxdelay; delay = delay + lag)
      {
        sxy = 0;
        for (i = 0; i < datasize; i++)
        {
          j = i + delay;
          if (j < 0 || j >= datasize)
          {
            continue;
          }
          sxy += (data1[i]) * (data2[j]);
        }
        result[delay] = sxy;
      }
      return result;
    }

    XCorrArrayType calcxcorr(std::vector<double> & data1,
        std::vector<double> & data2, bool normalize)
    {
      //OPENMS_PRECONDITION(data1.size() != 0 && data1.size() == data2.size(), "Both data vectors need to have the same length");
      int maxdelay = data1.size();
      int lag = 1;

      XCorrArrayType result;
      double mean1 = std::accumulate(data1.begin(), data1.end(), 0.)
          / (double)data1.size();
      double mean2 = std::accumulate(data2.begin(), data2.end(), 0.)
          / (double)data2.size();
      double denominator = 1;
      int datasize = data1.size();
      int i, j, delay;

      // Normalized cross-correlation = subtract the mean and divide by the standard deviation
      if (normalize)
      {
        double sqsum1 = 0;
        double sqsum2 = 0;
        for (std::vector<double>::iterator it = data1.begin();
            it != data1.end(); it++)
        {
          sqsum1 += (*it - mean1) * (*it - mean1);
        }

        for (std::vector<double>::iterator it = data2.begin();
            it != data2.end(); it++)
        {
          sqsum2 += (*it - mean2) * (*it - mean2);
        }
        // sigma_1 * sigma_2 * n
        denominator = sqrt(sqsum1 * sqsum2);
      }

      for (delay = -maxdelay; delay <= maxdelay; delay = delay + lag)
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
          result[delay] = sxy / denominator;
        }
        else
        {
          // e.g. if all datapoints are zero
          result[delay] = 0;
        }
      }
      return result;
    }

  } //end namespace Scoring
} //end namespace

