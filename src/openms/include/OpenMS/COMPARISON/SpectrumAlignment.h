// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/MatchedIterator.h>

#include <vector>
#include <map>
#include <utility>
#include <algorithm>

namespace OpenMS
{

  /**
      @brief Aligns the peaks of two sorted spectra
      Method 1: Using a banded (width via 'tolerance' parameter) alignment if absolute tolerances are given.
                Scoring function is the m/z distance between peaks. Intensity does not play a role!

      Method 2: If relative tolerance (ppm) is specified a simple matching of peaks is performed:
      Peaks from s1 (usually the theoretical spectrum) are assigned to the closest peak in s2 if it lies in the tolerance window
      @note: a peak in s2 can be matched to none, one or multiple peaks in s1. Peaks in s1 may be matched to none or one peak in s2.
      @note: intensity is ignored 
      TODO: improve time complexity, currently O(|s1|*log(|s2|))

      @htmlinclude OpenMS_SpectrumAlignment.parameters

      @ingroup SpectraComparison
  */

  class OPENMS_DLLAPI SpectrumAlignment :
    public DefaultParamHandler
  {
public:

    // @name Constructors and Destructors
    // @{
    /// default constructor
    SpectrumAlignment();

    /// copy constructor
    SpectrumAlignment(const SpectrumAlignment & source);

    /// destructor
    ~SpectrumAlignment() override;

    /// assignment operator
    SpectrumAlignment & operator=(const SpectrumAlignment & source);
    // @}

    template <typename SpectrumType1, typename SpectrumType2>
    void getSpectrumAlignment(std::vector<std::pair<Size, Size> >& alignment, const SpectrumType1& s1, const SpectrumType2& s2) const
    {
      if (!s1.isSorted() || !s2.isSorted())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Input to SpectrumAlignment is not sorted!");
      }

      // clear result
      alignment.clear();
      double tolerance = (double)param_.getValue("tolerance");

      if (!param_.getValue("is_relative_tolerance").toBool() )
      {
        std::map<Size, std::map<Size, std::pair<Size, Size> > > traceback;
        std::map<Size, std::map<Size, double> > matrix;

        // init the matrix with "gap costs" tolerance
        matrix[0][0] = 0;
        for (Size i = 1; i <= s1.size(); ++i)
        {
          matrix[i][0] = i * tolerance;
          traceback[i][0]  = std::make_pair(i - 1, 0);
        }
        for (Size j = 1; j <= s2.size(); ++j)
        {
          matrix[0][j] = j * tolerance;
          traceback[0][j] = std::make_pair(0, j - 1);
        }

        // fill in the matrix
        Size left_ptr(1);
        Size last_i(0), last_j(0);

        for (Size i = 1; i <= s1.size(); ++i)
        {
          double pos1(s1[i - 1].getMZ());

          for (Size j = left_ptr; j <= s2.size(); ++j)
          {
            bool off_band(false);
            // find min of the three possible directions
            double pos2(s2[j - 1].getMZ());
            double diff_align = fabs(pos1 - pos2);

            // running off the right border of the band?
            if (pos2 > pos1 && diff_align > tolerance)
            {
              if (i < s1.size() && j < s2.size() && s1[i].getMZ() < pos2)
              {
                off_band = true;
              }
            }

            // can we tighten the left border of the band?
            if (pos1 > pos2 && diff_align > tolerance && j > left_ptr + 1)
            {
              ++left_ptr;
            }

            double score_align = diff_align;

            if (matrix.contains(i - 1) && matrix[i - 1].contains(j - 1))
            {
              score_align += matrix[i - 1][j - 1];
            }
            else
            {
              score_align += (i - 1 + j - 1) * tolerance;
            }

            double score_up = tolerance;
            if (matrix.contains(i) && matrix[i].contains(j - 1))
            {
              score_up += matrix[i][j - 1];
            }
            else
            {
              score_up += (i + j - 1) * tolerance;
            }

            double score_left = tolerance;
            if (matrix.contains(i - 1) && matrix[i - 1].contains(j))
            {
              score_left += matrix[i - 1][j];
            }
            else
            {
              score_left += (i - 1 + j) * tolerance;
            }

          if (score_align <= score_up && score_align <= score_left && diff_align <= tolerance)
          {
             matrix[i][j] = score_align;
             traceback[i][j] = std::make_pair(i - 1, j - 1);
             last_i = i;
             last_j = j;
          }
          else
          {
            if (score_up <= score_left)
            {
              matrix[i][j] = score_up;
              traceback[i][j] = std::make_pair(i, j - 1);
            }
            else
            {
              matrix[i][j] = score_left;
              traceback[i][j] = std::make_pair(i - 1, j);
            }
          }

          if (off_band)
          {
            break;
          }
        }
      }

      // do traceback
      Size i = last_i;
      Size j = last_j;

      while (i >= 1 && j >= 1)
      {
        if (traceback[i][j].first == i - 1 && traceback[i][j].second == j - 1)
        {
          alignment.push_back(std::make_pair(i - 1, j - 1));
        }
        Size new_i = traceback[i][j].first;
        Size new_j = traceback[i][j].second;

        i = new_i;
        j = new_j;
      }

      std::reverse(alignment.begin(), alignment.end());

      }
      else  // relative alignment (ppm tolerance)
      {        
        // find  closest match of s1[i] in s2 for all i
        MatchedIterator<SpectrumType1, PpmTrait> it(s1, s2, tolerance);
        for (; it != it.end(); ++it) alignment.emplace_back(it.refIdx(), it.tgtIdx());
      }
    }
  };
}
