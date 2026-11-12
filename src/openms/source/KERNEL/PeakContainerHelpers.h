// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h> // for OPENMS_ASSERTIONS
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

namespace OpenMS
{
  namespace Internal
  {
    /**
      @brief Private helpers shared by the 1D peak containers MSSpectrum, MSChromatogram and Mobilogram.

      These containers hold a vector of peaks plus float/string/integer data arrays parallel to it, and
      implement the same index-based reordering and nearest-peak lookups. The algorithms live here once;
      the classes' public member functions forward to them, passing their peak storage, data arrays,
      position accessor (m/z, RT or mobility) and container name (for exception messages). Not installed.

      All functions are templates so they inline into the callers' hot paths.
    */
    namespace PeakContainerHelpers
    {
      /**
        @brief Throws if a non-empty data array's size differs from @p peaks.

        @param[in] peaks Number of peaks in the container
        @param[in] float_arrays, string_arrays, integer_arrays The container's data arrays
        @param[in] container_name Name used in the exception message ("spectrum", "chromatogram", ...)
        @param[in] function Name of the calling function, reported by the exception

        @exception Exception::Precondition if a non-empty data array's size differs from @p peaks
      */
      template<class FloatArrays, class StringArrays, class IntegerArrays>
      void checkDataArraySizes(Size peaks, const FloatArrays& float_arrays, const StringArrays& string_arrays,
                               const IntegerArrays& integer_arrays, const char* container_name, const char* function)
      {
        auto check = [peaks, container_name, function](const auto& arrays, const char* what)
        {
          for (Size i = 0; i < arrays.size(); ++i)
          {
            if (!arrays[i].empty() && arrays[i].size() != peaks)
            {
              throw Exception::Precondition(__FILE__, __LINE__, function,
                                            std::string(what) + "[" + StringUtils::toStr(i) + "] size (" + StringUtils::toStr(arrays[i].size()) +
                                              ") does not match " + container_name + " size (" + StringUtils::toStr(peaks) + ")");
            }
          }
        };
        check(float_arrays, "FloatDataArray");
        check(string_arrays, "StringDataArray");
        check(integer_arrays, "IntegerDataArray");
      }

      /**
        @brief Throws if any of @p indices is not smaller than @p peaks (the range check of select()).

        @exception Exception::Precondition if an index is out of range
      */
      inline void checkSelectIndices(const std::vector<Size>& indices, Size peaks, const char* container_name, const char* function)
      {
        for (Size i = 0; i < indices.size(); ++i)
        {
          if (indices[i] >= peaks)
          {
            throw Exception::Precondition(__FILE__, __LINE__, function,
                                          "Index " + StringUtils::toStr(indices[i]) + " is out of range for a " + container_name + " of size " +
                                            StringUtils::toStr(peaks));
          }
        }
      }

      /// Replaces each non-empty array in @p arrays by its elements at @p indices (in that order).
      template<class DataArrays>
      void reorderDataArrays(DataArrays& arrays, const std::vector<Size>& indices)
      {
        const Size snew = indices.size();
        std::vector<typename DataArrays::value_type::value_type> tmp;
        for (Size i = 0; i < arrays.size(); ++i)
        {
          if (arrays[i].empty())
          {
            continue;
          }
          tmp.clear();
          tmp.reserve(snew);
          for (Size j = 0; j < snew; ++j)
          {
            tmp.push_back(std::move(arrays[i][indices[j]]));
          }
          std::swap(arrays[i], tmp);
        }
      }

      /**
        @brief Keeps only the peaks (and data array entries) at @p indices, in that order (selectUnchecked()).

        Indices must be in range and duplicate-free (asserted in debug builds only); every element is
        moved, so a repeated index would read a moved-from value. The data array sizes must already
        have been checked against @p peaks.
      */
      template<class PeakVector, class FloatArrays, class StringArrays, class IntegerArrays>
      void selectByIndices(PeakVector& peaks, FloatArrays& float_arrays, StringArrays& string_arrays,
                           IntegerArrays& integer_arrays, const std::vector<Size>& indices)
      {
        const Size snew = indices.size();

#ifdef OPENMS_ASSERTIONS
        {
          std::vector<bool> seen(peaks.size(), false);
          for (Size i = 0; i < snew; ++i)
          {
            OPENMS_PRECONDITION(indices[i] < seen.size() && !seen[indices[i]],
                                "selectUnchecked(): indices must be in range and duplicate-free");
            seen[indices[i]] = true;
          }
        }
#endif

        PeakVector tmp;
        tmp.reserve(snew);
        for (Size i = 0; i < snew; ++i)
        {
          tmp.push_back(std::move(peaks[indices[i]]));
        }
        peaks.swap(tmp);

        reorderDataArrays(float_arrays, indices);
        reorderDataArrays(string_arrays, indices);
        reorderDataArrays(integer_arrays, indices);
      }

      /**
        @brief Index of the peak in [@p first, @p last) whose position is nearest to @p x (findNearest(x)).

        @p insert_pos is the lower bound of @p x in the (sorted) range; @p pos returns a peak's position.
        On a tie the left peak wins. The range must not be empty.
      */
      template<class ConstIterator, class CoordinateType, class PositionGetter>
      Size nearestIndex(ConstIterator first, ConstIterator last, ConstIterator insert_pos, CoordinateType x, PositionGetter pos)
      {
        // border cases
        if (insert_pos == first)
        {
          return 0;
        }
        if (insert_pos == last)
        {
          return Size(last - first) - 1;
        }
        // the peak before or the current peak are closest
        ConstIterator prev = insert_pos;
        --prev;
        if (std::fabs(pos(*insert_pos) - x) < std::fabs(pos(*prev) - x))
        {
          return Size(insert_pos - first);
        }
        else
        {
          return Size(prev - first);
        }
      }

      /// Index of the peak nearest to @p x if it lies within [x - tolerance, x + tolerance], -1 otherwise
      /// (findNearest(x, tolerance)). Uses @p container's findNearest(x).
      template<class Container, class CoordinateType, class PositionGetter>
      Int nearestIndexInTolerance(const Container& container, CoordinateType x, CoordinateType tolerance, PositionGetter pos)
      {
        if (container.empty())
        {
          return -1;
        }
        Size i = container.findNearest(x);
        const double found = pos(container[i]);
        if (found >= x - tolerance && found <= x + tolerance)
        {
          return static_cast<Int>(i);
        }
        else
        {
          return -1;
        }
      }

      /// Index of the peak nearest to @p x within the asymmetric window [x - tolerance_left, x + tolerance_right],
      /// -1 if there is none (findNearest(x, tolerance_left, tolerance_right)). Uses @p container's findNearest(x).
      template<class Container, class CoordinateType, class PositionGetter>
      Int nearestIndexInWindow(const Container& container, CoordinateType x, CoordinateType tolerance_left,
                               CoordinateType tolerance_right, PositionGetter pos)
      {
        if (container.empty())
        {
          return -1;
        }
        // do a binary search for nearest peak first
        Size i = container.findNearest(x);

        const double nearest = pos(container[i]);

        if (nearest < x)
        {
          if (nearest >= x - tolerance_left)
          {
            return i; // success: nearest peak is in left tolerance window
          }
          else
          {
            if (i == container.size() - 1)
            {
              return -1; // we are at the last peak which is too far left
            }
            // Nearest peak is too far left so there can't be a closer peak in the left window.
            // There still might be a peak to the right of x that falls in the right window
            ++i; // now we are at a peak exactly on or to the right of x
            const double next = pos(container[i]);
            if (next <= x + tolerance_right)
            {
              return i;
            }
          }
        }
        else
        {
          if (nearest <= x + tolerance_right)
          {
            return i; // success: nearest peak is in right tolerance window
          }
          else
          {
            if (i == 0)
            {
              return -1; // we are at the first peak which is too far right
            }
            --i; // now we are at a peak exactly on or to the left of x
            const double next = pos(container[i]);
            if (next >= x - tolerance_left)
            {
              return i;
            }
          }
        }

        // neither in the left nor the right tolerance window
        return -1;
      }

      /// Index (relative to @p first) of the most intense peak in [@p left, @p right), -1 if that range is
      /// empty (findHighestInWindow()). On a tie the leftmost peak wins.
      template<class ConstIterator>
      Int highestIndexInRange(ConstIterator first, ConstIterator left, ConstIterator right)
      {
        // no peak in the window
        if (left == right)
        {
          return -1;
        }
        using PeakType = typename std::iterator_traits<ConstIterator>::value_type;
        auto max_intensity_it = std::max_element(left, right, typename PeakType::IntensityLess());
        // find peak (index) with highest intensity to expected position
        return max_intensity_it - first;
      }

      /// The most intense peak in [@p first, @p last), or @p first if the range is empty (getBasePeak()).
      /// On a tie the leftmost peak wins.
      template<class ConstIterator>
      ConstIterator basePeak(ConstIterator first, ConstIterator last)
      {
        ConstIterator largest = first;
        if (first == last)
        {
          return largest;
        }
        ConstIterator current = first;
        ++current;
        for (; current != last; ++current)
        {
          if (largest->getIntensity() < current->getIntensity())
          {
            largest = current;
          }
        }
        return largest;
      }
    } // namespace PeakContainerHelpers
  } // namespace Internal
} // namespace OpenMS
