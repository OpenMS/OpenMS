// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.h>

#include <OpenMS/DATASTRUCTURES/String.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <iostream>

namespace OpenMS
{

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
    while (mz_it != mz_end && (*mz_it) < mz)
    {
      mz_it++;
      int_it++;
    }

    // walk right and left and add to our intensity
    mz_walker  = mz_it;
    int_walker = int_it;

    // if we moved past the end of the spectrum, we need to try the last peak
    // of the spectrum (it could still be within the window)
    if (mz_it == mz_end)
    {
      --mz_walker;
      --int_walker;
    }

    // add the current peak if it is between right and left
    if ((*mz_walker) > left && (*mz_walker) < right)
    {
      integrated_intensity += (*int_walker);
    }

    // (i) Walk to the left one step and then keep walking left until we go
    // outside the window. Note for the first step to the left we have to
    // check for the walker becoming equal to the first data point.
    mz_walker  = mz_it;
    int_walker = int_it;
    if (mz_it != mz_start)
    {
      --mz_walker;
      --int_walker;

      // Special case: target m/z is larger than first data point but the first
      // data point is inside the window.
      // Then, mz_it is the second data point, mz_walker now points to the very
      // first data point. If mz_it was the first data point, we already added
      // it above. We still need to add this point if it is inside the window
      // (while loop below will not catch it)
      if (mz_walker == mz_start && (*mz_walker) > left && (*mz_walker) < right)
      {
        integrated_intensity += (*int_walker);
      }
    }
    while (mz_walker != mz_start && (*mz_walker) > left && (*mz_walker) < right)
    {
      integrated_intensity += (*int_walker);
      --mz_walker;
      --int_walker;
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
    while (mz_it != mz_end && (*mz_it) < mz)
    {
      mz_it++;
      im_it++;
      int_it++;
    }

    // walk right and left and add to our intensity
    mz_walker  = mz_it;
    im_walker  = im_it;
    int_walker = int_it;

    // if we moved past the end of the spectrum, we need to try the last peak
    // of the spectrum (it could still be within the window)
    if (mz_it == mz_end)
    {
      --mz_walker;
      --im_walker;
      --int_walker;
    }

    // add the current peak if it is between right and left
    if ((*mz_walker) > left && (*mz_walker) < right && (*im_walker) > left_im && (*im_walker) < right_im)
    {
      integrated_intensity += (*int_walker);
    }

    // (i) Walk to the left one step and then keep walking left until we go
    // outside the window. Note for the first step to the left we have to
    // check for the walker becoming equal to the first data point.
    mz_walker  = mz_it;
    int_walker = int_it;
    im_walker = im_it;
    if (mz_it != mz_start)
    {
      --mz_walker;
      --im_walker;
      --int_walker;

      // Special case: target m/z is larger than first data point but the first
      // data point is inside the window.
      // Then, mz_it is the second data point, mz_walker now points to the very
      // first data point. If mz_it was the first data point, we already added
      // it above. We still need to add this point if it is inside the window
      // (while loop below will not catch it)
      if (mz_walker == mz_start && (*mz_walker) > left && (*mz_walker) < right && (*im_walker) > left_im && (*im_walker) < right_im)
      {
        integrated_intensity += (*int_walker);
      }
    }
    while (mz_walker != mz_start && (*mz_walker) > left && (*mz_walker) < right)
    {
      if (*im_walker > left_im && *im_walker < right_im) integrated_intensity += (*int_walker);
      --mz_walker;
      --im_walker;
      --int_walker;
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
      const String& filter)
  {
    Size input_size = input->getNrSpectra();
    if (input_size < 1)
    {
      return;
    }

    if (output.size() != extraction_coordinates.size())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Output and extraction coordinates need to have the same size: "+ String(output.size()) + " != " + String(extraction_coordinates.size()) );
    }

    int used_filter = getFilterNr_(filter);
    // assert that they are sorted!
    if (std::adjacent_find(extraction_coordinates.begin(), extraction_coordinates.end(),
          ExtractionCoordinates::SortExtractionCoordinatesReverseByMZ) != extraction_coordinates.end())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Input to extractChromatogram needs to be sorted by m/z");
    }

    //go through all spectra
    startProgress(0, input_size, "Extracting chromatograms");
    for (Size scan_idx = 0; scan_idx < input_size; ++scan_idx)
    {
      setProgress(scan_idx);

      OpenSwath::SpectrumPtr sptr = input->getSpectrumById(scan_idx);
      OpenSwath::SpectrumMeta s_meta = input->getSpectrumMetaById(scan_idx);

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
      bool has_im = (im_extraction_window > 0.0);
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
      for (Size k = 0; k < extraction_coordinates.size(); ++k)
      {
        double integrated_intensity = 0;
        double current_rt = s_meta.RT;
        if (extraction_coordinates[k].rt_end - extraction_coordinates[k].rt_start > 0 &&
             (current_rt < extraction_coordinates[k].rt_start ||
              current_rt > extraction_coordinates[k].rt_end) )
        {
          continue;
        }

        const bool use_im = (extraction_coordinates[k].ion_mobility >= 0.0 && has_im);
        if (!use_im && used_filter == 1)
        {
          extract_value_tophat(mz_start, mz_it, mz_end, int_it,
                               extraction_coordinates[k].mz, integrated_intensity, mz_extraction_window, ppm);
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
    endProgress();
  }

  int ChromatogramExtractorAlgorithm::getFilterNr_(const String& filter)
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

