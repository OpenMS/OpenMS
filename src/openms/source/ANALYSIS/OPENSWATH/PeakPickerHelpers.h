// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Types.h>

#include <cmath>
#include <vector>

namespace OpenMS
{
  namespace Internal
  {
    /**
      @brief Private helpers shared by PeakPickerChromatogram and PeakPickerMobilogram.

      Both pickers find the raw data point closest to a picked apex, resolve overlapping peak borders
      and integrate the peaks in the same way; only the container (MSChromatogram / Mobilogram), its
      position accessor (RT / ion mobility) and some logging details differ. The algorithms live here
      once; the classes' member functions forward to them. Not installed.

      All functions are templates so they inline into the callers.
    */
    namespace PeakPickerHelpers
    {
      /**
        @brief Finds the peak of @p container closest to @p target, starting the search at index @p current_peak

        The function assumes the closest peak is to the right of @p current_peak.

        @param[in] pos pos(peak): position (RT, ion mobility) of a peak
        @return The index of the closest peak (container.size() if @p target lies behind the last peak)
      */
      template <typename Container, typename PosOf>
      Size findClosestPeak(const Container& container, double target, Size current_peak, PosOf&& pos)
      {
        while (current_peak < container.size())
        {
          // check if we have walked past the position of the peak
          if (target < pos(container[current_peak]))
          {
            // see which one is closer, the current one or the one before
            if (current_peak > 0 &&
                std::fabs(target - pos(container[current_peak - 1])) <
                std::fabs(target - pos(container[current_peak])))
            {
              current_peak--;
            }

            return current_peak;
          }
          current_peak++;
        }
        return current_peak;
      }

      /**
        @brief Compute peak area (peak integration): sums the intensities of @p container within each [left, right] index pair
      */
      template <typename Container, typename WidthVector>
      void integratePeaks(const Container& container, const WidthVector& left_width, const WidthVector& right_width, std::vector<double>& integrated_intensities)
      {
        for (Size i = 0; i < left_width.size(); i++)
        {
          const int current_left_idx = left_width[i];
          const int current_right_idx = right_width[i];

          // Also integrate the intensities
          integrated_intensities[i] = 0;
          for (int k = current_left_idx; k <= current_right_idx; k++)
          {
            integrated_intensities[i] += container[k].getIntensity();
          }
        }
      }

      /**
        @brief Finds a new border between each pair of consecutive picked peaks whose borders overlap

        @param[in] container The (raw or smoothed) data the borders refer to
        @param[in] picked The picked peaks
        @param[in,out] left_width, right_width Left / right border (index into @p container) of each picked peak
        @param[in] integrated_intensities Intensity of each picked peak (only logged)
        @param[in] pos pos(peak): position (RT, ion mobility) of a peak
        @param[in] report_still_overlapping report_still_overlapping(new_left_border, new_right_border): called if the new borders still overlap
        @param[in] eol Line terminator streamed to the debug log (a character or a stream manipulator)
      */
      template <typename Container, typename WidthVector, typename PosOf, typename ReportStillOverlapping, typename Eol>
      void removeOverlappingPeaks(const Container& container, const Container& picked, WidthVector& left_width, WidthVector& right_width,
                                  const std::vector<double>& integrated_intensities, PosOf&& pos, ReportStillOverlapping&& report_still_overlapping, Eol eol)
      {
        if (picked.empty()) {return; }
        OPENMS_LOG_DEBUG << "Remove overlapping peaks now (size " << picked.size() << ")" << eol;
        Size current_peak = 0;
        // Find overlapping peaks
        for (Size i = 0; i < picked.size() - 1; i++)
        {
          // Check whether the current right overlaps with the next left
          // See whether we can correct this and find some border between the two
          // features ...
          if (right_width[i] > left_width[i + 1])
          {
            const int current_left_idx = left_width[i];
            const int current_right_idx = right_width[i];
            const int next_left_idx = left_width[i + 1];
            const int next_right_idx = right_width[i + 1];
            OPENMS_LOG_DEBUG << " Found overlapping " << i << " : " << current_left_idx << " " << current_right_idx << eol;
            OPENMS_LOG_DEBUG << "                   -- with  " << i + 1 << " : " << next_left_idx << " " << next_right_idx << eol;

            // Find the peak width and best position
            double central_peak_pos = pos(picked[i]);
            double next_peak_pos = pos(picked[i + 1]);
            current_peak = findClosestPeak(container, central_peak_pos, current_peak, pos);
            Size next_peak = findClosestPeak(container, next_peak_pos, current_peak, pos);

            // adjust the right border of the current and left border of next
            Size k = 1;
            while ((current_peak + k) < container.size()
                  && (container[current_peak + k].getIntensity() < container[current_peak + k - 1].getIntensity()))
            {
              ++k;
            }
            Size new_right_border = current_peak + k - 1;
            k = 1;
            while ((next_peak - k + 1) > 0
                  && (container[next_peak - k].getIntensity() < container[next_peak - k + 1].getIntensity()))
            {
              ++k;
            }
            Size new_left_border = next_peak - k + 1;

            // assert that the peaks are now not overlapping any more ...
            if (new_left_border < new_right_border)
            {
              report_still_overlapping(new_left_border, new_right_border);
              new_left_border = (new_left_border + new_right_border) / 2;
              new_right_border = (new_left_border + new_right_border) / 2;

            }

            OPENMS_LOG_DEBUG << "New peak l: " << pos(container[current_left_idx])
              << " " << pos(container[new_right_border])
              << " int " << integrated_intensities[i] << eol;
            OPENMS_LOG_DEBUG << "New peak r: " << pos(container[new_left_border])
              << " " << pos(container[next_right_idx])
              << " int " << integrated_intensities[i + 1] << eol;


            right_width[i] = new_right_border;
            left_width[i + 1] = new_left_border;

          }
        }
      }
    } // namespace PeakPickerHelpers
  } // namespace Internal
} // namespace OpenMS
