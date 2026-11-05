// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.h>

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMSInMemory.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <algorithm>
#include <cstddef>
#include <exception>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS
{

  namespace
  {
    /**
      @brief Advances @p mz_it to the first peak whose m/z is not below @p mz and returns the number of peaks skipped.

      The targets are visited in m/z order, so the peak is usually a few steps ahead: those are
      tried first. If the target is further away (e.g. the precursors of one isolation window in
      a full MS1 spectrum), a binary search over the m/z-sorted remainder replaces the walk.
    */
    inline std::ptrdiff_t advanceToMz_(std::vector<double>::const_iterator& mz_it,
                                       const std::vector<double>::const_iterator& mz_end,
                                       const double mz)
    {
      const std::vector<double>::const_iterator begin = mz_it;
      for (int step = 0; step < 8 && mz_it != mz_end && (*mz_it) < mz; ++step)
      {
        ++mz_it;
      }
      if (mz_it != mz_end && (*mz_it) < mz)
      {
        mz_it = std::lower_bound(mz_it, mz_end, mz);
      }
      return mz_it - begin;
    }
  }

  void ChromatogramExtractorAlgorithm::extract_value_tophat(
      const std::vector<double>::const_iterator& mz_start,
            std::vector<double>::const_iterator& mz_it,
      const std::vector<double>::const_iterator& mz_end,
            std::vector<double>::const_iterator& int_it,
      const double mz,
      double& integrated_intensity,
      const double mz_extraction_window,
      const bool ppm)
  {
    integrated_intensity = 0;
    if (mz_start == mz_end)
    {
      return;
    }

    // calculate extraction window
    double left, right;
    if (ppm)
    {
      left  = mz - mz * mz_extraction_window / 2.0 * 1.0e-6;
      right = mz + mz * mz_extraction_window / 2.0 * 1.0e-6;
    }
    else
    {
      left  = mz - mz_extraction_window / 2.0;
      right = mz + mz_extraction_window / 2.0;
    }

    std::vector<double>::const_iterator mz_walker;
    std::vector<double>::const_iterator int_walker;

    // advance the mz / int iterator until we hit the m/z value of the next transition
    const std::ptrdiff_t skipped = advanceToMz_(mz_it, mz_end, mz);
    int_it += skipped;

    // add the current peak if it is between right and left (past the end of
    // the spectrum there is none: the walk to the left starts with the last peak)
    if (mz_it != mz_end && (*mz_it) > left && (*mz_it) < right)
    {
      integrated_intensity += (*int_it);
    }

    // (i) Walk to the left until we go outside the window or have added the
    // first data point
    mz_walker  = mz_it;
    int_walker = int_it;
    while (mz_walker != mz_start)
    {
      --mz_walker;
      --int_walker;
      if (!((*mz_walker) > left && (*mz_walker) < right))
      {
        break;
      }
      integrated_intensity += (*int_walker);
    }

    // (ii) Walk to the right one step and then keep walking right until we are
    // outside the window
    mz_walker  = mz_it;
    int_walker = int_it;
    if (mz_it != mz_end)
    {
      ++mz_walker;
      ++int_walker;
    }
    while (mz_walker != mz_end && (*mz_walker) > left && (*mz_walker) < right)
    {
      integrated_intensity += (*int_walker);
      ++mz_walker;
      ++int_walker;
    }
  }

  void ChromatogramExtractorAlgorithm::extract_value_tophat(
      const std::vector<double>::const_iterator& mz_start,
            std::vector<double>::const_iterator& mz_it,
      const std::vector<double>::const_iterator& mz_end,
            std::vector<double>::const_iterator& int_it,
            std::vector<double>::const_iterator& im_it,
      const double mz,
      const double im,
      double& integrated_intensity,
      const double mz_extraction_window,
      const double im_extraction_window,
      const bool ppm)
  {
    // Note that we have a 3D spectrum with m/z, intensity and ion mobility.
    // The spectrum is sorted by m/z but we expect to have ion mobility
    // information for each m/z point as well. Right now we simply filter by
    // ion mobility and skip data that does not fall within the ion mobility
    // window.

    integrated_intensity = 0;
    if (mz_start == mz_end)
    {
      return;
    }

    // calculate extraction window
    double left, right;
    if (ppm)
    {
      left  = mz - mz * mz_extraction_window / 2.0 * 1.0e-6;
      right = mz + mz * mz_extraction_window / 2.0 * 1.0e-6;
    }
    else
    {
      left  = mz - mz_extraction_window / 2.0;
      right = mz + mz_extraction_window / 2.0;
    }
    double left_im  = im - im_extraction_window / 2.0;
    double right_im = im + im_extraction_window / 2.0;

    std::vector<double>::const_iterator mz_walker;
    std::vector<double>::const_iterator im_walker;
    std::vector<double>::const_iterator int_walker;

    // advance the mz / int iterator until we hit the m/z value of the next transition
    const std::ptrdiff_t skipped = advanceToMz_(mz_it, mz_end, mz);
    im_it += skipped;
    int_it += skipped;

    // add the current peak if it is between right and left (past the end of
    // the spectrum there is none: the walk to the left starts with the last peak)
    if (mz_it != mz_end && (*mz_it) > left && (*mz_it) < right && (*im_it) > left_im && (*im_it) < right_im)
    {
      integrated_intensity += (*int_it);
    }

    // (i) Walk to the left until we go outside the m/z window or have added
    // the first data point
    mz_walker  = mz_it;
    im_walker  = im_it;
    int_walker = int_it;
    while (mz_walker != mz_start)
    {
      --mz_walker;
      --im_walker;
      --int_walker;
      if (!((*mz_walker) > left && (*mz_walker) < right))
      {
        break;
      }
      if (*im_walker > left_im && *im_walker < right_im) integrated_intensity += (*int_walker);
    }

    // (ii) Walk to the right one step and then keep walking right until we are
    // outside the window
    mz_walker  = mz_it;
    im_walker  = im_it;
    int_walker = int_it;
    if (mz_it != mz_end)
    {
      ++im_walker;
      ++mz_walker;
      ++int_walker;
    }
    while (mz_walker != mz_end && (*mz_walker) > left && (*mz_walker) < right)
    {
      if (*im_walker > left_im && *im_walker < right_im) integrated_intensity += (*int_walker);
      ++mz_walker;
      ++im_walker;
      ++int_walker;
    }
  }

  void ChromatogramExtractorAlgorithm::extractChromatograms(const OpenSwath::SpectrumAccessPtr& input,
      std::vector< OpenSwath::ChromatogramPtr >& output,
      const std::vector<ExtractionCoordinates>& extraction_coordinates,
      double mz_extraction_window,
      bool ppm,
      double im_extraction_window,
      const std::string& filter)
  {
    Size input_size = input->getNrSpectra();
    if (input_size < 1)
    {
      return;
    }

    if (output.size() != extraction_coordinates.size())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Output and extraction coordinates need to have the same size: "+ StringUtils::toStr(output.size()) + " != " + StringUtils::toStr(extraction_coordinates.size()) );
    }

    int used_filter = getFilterNr_(filter);
    // assert that they are sorted!
    if (std::adjacent_find(extraction_coordinates.begin(), extraction_coordinates.end(),
          ExtractionCoordinates::SortExtractionCoordinatesReverseByMZ) != extraction_coordinates.end())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Input to extractChromatogram needs to be sorted by m/z");
    }

    // The RT of every spectrum, read once here and used by the extraction loop below (with an on-disk
    // access such as sqMass, every meta data read is a query).
    std::vector<double> spectrum_rts(input_size);
    for (Size scan_idx = 0; scan_idx < input_size; ++scan_idx)
    {
      spectrum_rts[scan_idx] = input->getSpectrumMetaById(scan_idx).RT;
    }

    // A chromatogram gets one point for every spectrum inside its RT window: reserve that many
    // up front instead of growing the arrays point by point (which leaves each of them with up
    // to twice the capacity it needs and fragments the heap with the outgrown buffers).
    {
      std::vector<double> sorted_rts;
      if (std::any_of(extraction_coordinates.begin(), extraction_coordinates.end(),
                      [](const ExtractionCoordinates& coord) { return coord.rt_end - coord.rt_start > 0; }))
      {
        sorted_rts = spectrum_rts;
        std::sort(sorted_rts.begin(), sorted_rts.end());
      }
      for (Size k = 0; k < extraction_coordinates.size(); ++k)
      {
        const ExtractionCoordinates& coord = extraction_coordinates[k];
        Size nr_points = input_size;
        if (coord.rt_end - coord.rt_start > 0)
        {
          nr_points = std::upper_bound(sorted_rts.begin(), sorted_rts.end(), coord.rt_end) -
                      std::lower_bound(sorted_rts.begin(), sorted_rts.end(), coord.rt_start);
        }
        std::vector<double>& times = output[k]->getTimeArray()->data;
        std::vector<double>& intensities = output[k]->getIntensityArray()->data;
        times.reserve(times.size() + nr_points);
        intensities.reserve(intensities.size() + nr_points);
      }
    }

    // Extracts the chromatograms [coord_begin, coord_end) from all spectra. Each chromatogram only
    // depends on its own coordinates, so disjoint ranges can be extracted concurrently.
    const bool has_im = (im_extraction_window > 0.0);
    auto extract_range = [&](const Size coord_begin, const Size coord_end, const bool report_progress)
    {
      for (Size scan_idx = 0; scan_idx < input_size; ++scan_idx)
      {
        if (report_progress)
        {
          setProgress(scan_idx);
        }

        OpenSwath::SpectrumPtr sptr = input->getSpectrumById(scan_idx);
        const double current_rt = spectrum_rts[scan_idx];

        OpenSwath::BinaryDataArrayPtr mz_arr = sptr->getMZArray();
        OpenSwath::BinaryDataArrayPtr int_arr = sptr->getIntensityArray();
        std::vector<double>::const_iterator mz_start = mz_arr->data.begin();
        std::vector<double>::const_iterator mz_end = mz_arr->data.end();
        std::vector<double>::const_iterator mz_it = mz_arr->data.begin();
        std::vector<double>::const_iterator int_it = int_arr->data.begin();
        std::vector<double>::const_iterator im_it;

        if (sptr->getMZArray()->data.empty())
        {
          continue;
        }

        // Look for ion mobility array
        if (has_im)
        {
          OpenSwath::BinaryDataArrayPtr im_arr = sptr->getDriftTimeArray();
          if (im_arr != nullptr)
          {
            im_it = im_arr->data.begin();
          }
          else
          {
            throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Requested ion mobility extraction but no ion mobility array found.");
          }
        }

        // go through all transitions / chromatograms which are sorted by
        // ProductMZ. We can use this to step through the spectrum and at the
        // same time step through the transitions. We increase the peak counter
        // until we hit the next transition and then extract the signal.
        for (Size k = coord_begin; k < coord_end; ++k)
        {
          double integrated_intensity = 0;
          if (extraction_coordinates[k].rt_end - extraction_coordinates[k].rt_start > 0 &&
               (current_rt < extraction_coordinates[k].rt_start ||
                current_rt > extraction_coordinates[k].rt_end) )
          {
            continue;
          }

          const bool use_im = (extraction_coordinates[k].ion_mobility >= 0.0 && has_im);
          if (!use_im && used_filter == 1)
          {
            const std::vector<double>::const_iterator mz_before = mz_it;
            extract_value_tophat(mz_start, mz_it, mz_end, int_it,
                                 extraction_coordinates[k].mz, integrated_intensity, mz_extraction_window, ppm);
            if (has_im)
            {
              // keep the ion mobility iterator at the same peak as the m/z iterator for the next
              // coordinate that uses ion mobility
              im_it += mz_it - mz_before;
            }
          }
          else if (use_im && used_filter == 1)
          {
            if (extraction_coordinates[k].ion_mobility < 0)
            {
              std::cerr << "WARNING : Drift time of ion is negative!" << std::endl;
            }
            extract_value_tophat(mz_start, mz_it, mz_end, int_it, im_it,
                                 extraction_coordinates[k].mz, extraction_coordinates[k].ion_mobility,
                                 integrated_intensity, mz_extraction_window, im_extraction_window, ppm);
          }
          else if (used_filter == 2)
          {
            throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
          }

          output[k]->getTimeArray()->data.push_back(current_rt);
          output[k]->getIntensityArray()->data.push_back(integrated_intensity);
        }
      }
    };

    // Called from a single thread (e.g. a SWATH workflow that extracts one SWATH at a time), the
    // otherwise idle threads extract contiguous m/z ranges of the coordinates. Every range visits
    // all spectra, so a range needs enough coordinates to pay for that, and the spectra need to be
    // in memory: other accesses read (and possibly transform) a spectrum on every call and are not
    // safe to share between threads.
    constexpr Size min_coordinates_per_range = 256;
    Size nr_ranges = 1;
#ifdef _OPENMP
    if (omp_in_parallel() == 0 && omp_get_max_threads() > 1 &&
        dynamic_cast<const SpectrumAccessOpenMSInMemory*>(input.get()) != nullptr)
    {
      // a few ranges per thread balance the uneven RT windows of the coordinates
      nr_ranges = std::min<Size>(static_cast<Size>(omp_get_max_threads()) * 4,
                                 extraction_coordinates.size() / min_coordinates_per_range);
    }
#endif

    //go through all spectra
    startProgress(0, input_size, "Extracting chromatograms");
    if (nr_ranges <= 1)
    {
      extract_range(0, extraction_coordinates.size(), true);
    }
    else
    {
      std::exception_ptr range_exception;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
      for (SignedSize range = 0; range < static_cast<SignedSize>(nr_ranges); ++range)
      {
        const Size coord_begin = extraction_coordinates.size() * range / nr_ranges;
        const Size coord_end = extraction_coordinates.size() * (range + 1) / nr_ranges;
        try
        {
          extract_range(coord_begin, coord_end, false);
        }
        catch (...)
        {
#ifdef _OPENMP
#pragma omp critical (chromatogram_extractor_exception)
#endif
          {
            if (range_exception == nullptr)
            {
              range_exception = std::current_exception();
            }
          }
        }
      }
      if (range_exception != nullptr)
      {
        std::rethrow_exception(range_exception);
      }
    }
    endProgress();
  }

  int ChromatogramExtractorAlgorithm::getFilterNr_(const std::string& filter)
  {
    if (filter == "tophat")
    {
      return 1;
    }
    else if (filter == "bartlett")
    {
      return 2;
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "Filter either needs to be tophat or bartlett");
    }
  }

}

