// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest, Luis Jacob Keller, Alen Saric$
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Macros.h>
#include <atomic>
#include <cmath>
#include <iterator>
#include <limits>
#include <vector>

namespace OpenMS
{
/// Shared process-wide flag controlling resampling-spacing warnings.
extern OPENMS_DLLAPI std::atomic<bool> suppress_resampling_spacing_warning;

namespace Internal
{
  /**
    @brief Shared linear resampling implementation without parameter or progress handling.

    Used by MSExperiment::calculateTIC() and LinearResamplerAlign. The latter
    retains the parameter interface and experiment-wide processing API.
  */
  class OPENMS_DLLAPI LinearResampling
  {
  public:
    /**
      @brief Construct a resampler with explicit spacing and units.
      @param[in] spacing Spacing of output points.
      @param[in] ppm Whether spacing is relative (ppm) instead of absolute.
    */
    explicit LinearResampling(double spacing, bool ppm = false): spacing_(spacing), ppm_(ppm)
    {
    }

    /**
      @brief Resample a peak container onto a grid spanning its first and last points.
      @param[in,out] container Container to resample, with the same requirements as LinearResamplerAlign::raster().
    */
    template<class PeakContainerT>
    void raster(PeakContainerT& container)
    {
      // return if nothing to do
      if (container.empty()) return;

      auto first = container.begin();
      auto last = container.end();

      double end_pos = (last - 1)->getPos();
      double start_pos = first->getPos();
      int number_resampled_points = (int)(ceil((end_pos - start_pos) / spacing_ + 1));

      std::vector<typename PeakContainerT::PeakType> resampled_peak_container;
      populateRaster(resampled_peak_container, start_pos, end_pos, number_resampled_points);

      raster(container.begin(), container.end(), resampled_peak_container.begin(), resampled_peak_container.end());

      container.swap(resampled_peak_container);
    }

    /**
      @brief Distribute input intensities onto an existing output grid.

      Intensities outside the grid accumulate at its endpoints. Existing
      output intensities are retained and incremented.
      @param[in] raw_it Start of the input range.
      @param[in] raw_end End of the input range.
      @param[in,out] resampled_begin Start of the nonempty output grid.
      @param[in,out] resampled_end End of the output grid.
    */
    template<typename PeakTypeIterator, typename ConstPeakTypeIterator>
    void raster(ConstPeakTypeIterator raw_it, ConstPeakTypeIterator raw_end, PeakTypeIterator resampled_begin, PeakTypeIterator resampled_end)
    {
      OPENMS_PRECONDITION(resampled_begin != resampled_end, "Output iterators cannot be identical") // as we use +1
      // OPENMS_PRECONDITION(raw_it != raw_end, "Input iterators cannot be identical")

      verifySpacing(raw_it, raw_end, [](auto x) { return x->getPos(); });

      PeakTypeIterator resample_start = resampled_begin;

      // need to get the raw iterator between two resampled iterators of the raw data
      while (raw_it != raw_end && raw_it->getPos() < resampled_begin->getPos())
      {
        resampled_begin->setIntensity(resampled_begin->getIntensity() + raw_it->getIntensity());
        raw_it++;
      }

      while (raw_it != raw_end)
      {
        // advance the resample iterator until our raw point is between two resampled iterators
        while (resampled_begin != resampled_end && resampled_begin->getPos() < raw_it->getPos())
        {
          resampled_begin++;
        }
        if (resampled_begin != resample_start) { resampled_begin--; }

        // if we have the last datapoint we break
        if ((resampled_begin + 1) == resampled_end) { break; }

        double dist_left = fabs(raw_it->getPos() - resampled_begin->getPos());
        double dist_right = fabs(raw_it->getPos() - (resampled_begin + 1)->getPos());

        // distribute the intensity of the raw point according to the distance to resample_it and resample_it+1
        resampled_begin->setIntensity(resampled_begin->getIntensity() + raw_it->getIntensity() * dist_right / (dist_left + dist_right));
        (resampled_begin + 1)->setIntensity((resampled_begin + 1)->getIntensity() + raw_it->getIntensity() * dist_left / (dist_left + dist_right));

        raw_it++;
      }

      // add the final intensity to the right
      while (raw_it != raw_end)
      {
        resampled_begin->setIntensity(resampled_begin->getIntensity() + raw_it->getIntensity());
        raw_it++;
      }
    }

    /**
      @brief Populate an absolute or ppm grid using the existing endpoint convention.
      @param[in,out] resampled_peak_container Output grid.
      @param[in] start_pos First grid position.
      @param[in] end_pos End position for ppm grids.
      @param[in] number_resampled_points Number of points for absolute grids.
    */
    template<typename PeakType>
    void populateRaster(std::vector<PeakType>& resampled_peak_container, double start_pos, double end_pos, int number_resampled_points)
    {
      if (! ppm_)
      {
        // generate the resampled peaks at positions origin+i*spacing_
        resampled_peak_container.resize(number_resampled_points);
        typename std::vector<PeakType>::iterator it = resampled_peak_container.begin();
        for (int i = 0; i < number_resampled_points; ++i)
        {
          it->setPos(start_pos + i * spacing_);
          ++it;
        }
      }
      else
      {
        // generate resampled peaks with ppm distance (not fixed)
        double current_mz = start_pos;
        while (current_mz < end_pos)
        {
          PeakType p;
          p.setIntensity(0);
          p.setPos(current_mz);
          resampled_peak_container.push_back(p);

          // increment current_mz
          current_mz += current_mz * (spacing_ / 1e6);
        }
      }
    }

    /**
      @brief Emit the existing spacing warning for absolute grids when enabled.
      @param[in] it Start of the input range.
      @param[in] end End of the input range.
      @param[in] access Function returning the position at an iterator.
    */
    template<typename PeakTypeIterator>
    void verifySpacing(PeakTypeIterator it, PeakTypeIterator end, auto access)
    {
      // ppm_ spacing is relative (parts-per-million) and cannot be compared
      // directly against the absolute neighbour distance computed below.
      if (ppm_) return;
      if (it == end || std::next(it) == end) return;
      double min_dist = std::numeric_limits<double>::infinity();
      double current_dist {};

      while (std::next(it) != end)
      {
        current_dist = (access(std::next(it)) - access(it));
        if (min_dist > current_dist) min_dist = current_dist;
        ++it;
      }

      if (spacing_ < min_dist && ! suppress_resampling_spacing_warning.load())
      {
        OPENMS_LOG_WARN << "Resampling spacing (" << spacing_ << ") is smaller than the smallest distance between data points (" << min_dist
                        << "). This approximates the detector dead time and may produce spurious peaks.\n";
      }
    }

  private:
    double spacing_;
    bool ppm_;
  };
} // namespace Internal
} // namespace OpenMS
