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

#define ALIGNMENT_DEBUG
#undef  ALIGNMENT_DEBUG

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
    /// traceback of a cell of the alignment matrix
    enum Move_ : unsigned char
    {
      MOVE_NONE_,  ///< not computed: leads to (0, 0)
      MOVE_ALIGN_, ///< from (i - 1, j - 1), aligning the two peaks
      MOVE_UP_,    ///< from (i, j - 1)
      MOVE_LEFT_   ///< from (i - 1, j)
    };

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
        // Banded alignment. Row i of the matrix is computed for the columns
        // [row_first[i], row_first[i] + row_length(i)), stored at row_offset[i] in the flat arrays.
        // Every other cell keeps its initial value: (i + j) * tolerance, which is also the gap cost
        // of the first row and column, and an empty traceback, which leads to (0, 0).
        const Size n1 = s1.size();
        std::vector<double> scores;
        std::vector<unsigned char> moves; // one of the Move_ values, for each computed cell
        std::vector<Size> row_first(n1 + 1, 0);
        std::vector<Size> row_offset(n1 + 2, 0);
        auto computed = [&](Size i, Size j) -> bool
        {
          return i >= 1 && i <= n1 && j >= row_first[i] && j - row_first[i] < row_offset[i + 1] - row_offset[i];
        };
        auto matrix = [&](Size i, Size j) -> double
        {
          return computed(i, j) ? scores[row_offset[i] + j - row_first[i]] : (i + j) * tolerance;
        };

        // fill in the matrix
        Size left_ptr(1);
        Size last_i(0), last_j(0);

        //Size off_band_counter(0);
        for (Size i = 1; i <= n1; ++i)
        {
          double pos1(s1[i - 1].getMZ());
          row_first[i] = left_ptr;
          row_offset[i + 1] = row_offset[i];

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

            double score_align = diff_align + matrix(i - 1, j - 1);
            double score_up = tolerance + matrix(i, j - 1);
            double score_left = tolerance + matrix(i - 1, j);

    #ifdef ALIGNMENT_DEBUG
          cerr << i << " " << j << " " << left_ptr << " " << pos1 << " " << pos2 << " " << score_align << " " << score_left << " " << score_up << endl;
    #endif

          // The cell (i, j) is appended to row i, which makes it visible to matrix() and computed().
          if (score_align <= score_up && score_align <= score_left && diff_align <= tolerance)
          {
             scores.push_back(score_align);
             moves.push_back(MOVE_ALIGN_);
             last_i = i;
             last_j = j;
          }
          else
          {
            if (score_up <= score_left)
            {
              scores.push_back(score_up);
              moves.push_back(MOVE_UP_);
            }
            else
            {
              scores.push_back(score_left);
              moves.push_back(MOVE_LEFT_);
            }
          }
          ++row_offset[i + 1];

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
        // a cell outside the band leads to (0, 0); seen from (1, 1) that is the aligning move
        const unsigned char move = computed(i, j) ? moves[row_offset[i] + j - row_first[i]]
                                                  : static_cast<unsigned char>(i == 1 && j == 1 ? MOVE_ALIGN_ : MOVE_NONE_);
        if (move == MOVE_ALIGN_)
        {
          alignment.push_back(std::make_pair(i - 1, j - 1));
        }
        switch (move)
        {
          case MOVE_ALIGN_: --i; --j; break;
          case MOVE_UP_: --j; break;
          case MOVE_LEFT_: --i; break;
          default: i = 0; j = 0; break;
        }
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
