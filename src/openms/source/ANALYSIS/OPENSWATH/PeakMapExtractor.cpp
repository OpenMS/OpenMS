// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/PeakMapExtractor.h>

#include <OpenMS/CONCEPT/Exception.h>

#include <algorithm>
#include <iterator>

namespace OpenMS
{
  void PeakMapExtractor::extractPeakMaps(const OpenSwath::SpectrumAccessPtr& input,
                                         std::vector<ExtractedPeakMap>& output,
                                         const std::vector<ExtractionCoordinates>& extraction_coordinates,
                                         double mz_extraction_window,
                                         bool ppm,
                                         double im_extraction_window,
                                         const String& filter)
  {
    const Size input_size = input->getNrSpectra();
    output.clear();
    output.resize(extraction_coordinates.size());

    for (Size k = 0; k < extraction_coordinates.size(); ++k)
    {
      output[k].native_id = extraction_coordinates[k].id;
      output[k].target_mz = extraction_coordinates[k].mz;
      output[k].target_ion_mobility = extraction_coordinates[k].ion_mobility;
      output[k].rt_start = extraction_coordinates[k].rt_start;
      output[k].rt_end = extraction_coordinates[k].rt_end;
      if (output[k].rt_end - output[k].rt_start > 0.0)
      {
        output[k].target_rt = (output[k].rt_start + output[k].rt_end) / 2.0;
      }
    }

    if (input_size < 1 || extraction_coordinates.empty())
    {
      return;
    }

    const int used_filter = getFilterNr_(filter);
    if (used_filter != 1)
    {
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    startProgress(0, input_size, "Extracting targeted peak maps");
    for (Size scan_idx = 0; scan_idx < input_size; ++scan_idx)
    {
      setProgress(scan_idx);

      OpenSwath::SpectrumPtr sptr = input->getSpectrumById(scan_idx);
      OpenSwath::SpectrumMeta s_meta = input->getSpectrumMetaById(scan_idx);

      OpenSwath::BinaryDataArrayPtr mz_arr = sptr->getMZArray();
      OpenSwath::BinaryDataArrayPtr int_arr = sptr->getIntensityArray();
      OpenSwath::BinaryDataArrayPtr im_arr = sptr->getDriftTimeArray();

      if (mz_arr->data.empty())
      {
        continue;
      }
      if (im_arr == nullptr)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "Peak map extraction requires ion mobility data on every input spectrum.");
      }

      const std::vector<double>& mz_data = mz_arr->data;
      const std::vector<double>& int_data = int_arr->data;
      const std::vector<double>& im_data = im_arr->data;

      if (mz_data.size() != int_data.size() || mz_data.size() != im_data.size())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "m/z, intensity, and ion mobility arrays need to have the same length.");
      }

      const double current_rt = s_meta.RT;
      for (Size k = 0; k < extraction_coordinates.size(); ++k)
      {
        const ExtractionCoordinates& coord = extraction_coordinates[k];
        if (coord.rt_end - coord.rt_start > 0 &&
            (current_rt < coord.rt_start || current_rt > coord.rt_end))
        {
          continue;
        }

        double left = 0.0;
        double right = 0.0;
        if (ppm)
        {
          left = coord.mz - coord.mz * mz_extraction_window / 2.0 * 1.0e-6;
          right = coord.mz + coord.mz * mz_extraction_window / 2.0 * 1.0e-6;
        }
        else
        {
          left = coord.mz - mz_extraction_window / 2.0;
          right = coord.mz + mz_extraction_window / 2.0;
        }

        const bool use_im_filter = (im_extraction_window > 0.0 && coord.ion_mobility >= 0.0);
        const double left_im = coord.ion_mobility - im_extraction_window / 2.0;
        const double right_im = coord.ion_mobility + im_extraction_window / 2.0;

        auto mz_it = std::lower_bound(mz_data.begin(), mz_data.end(), left);
        Size idx = static_cast<Size>(std::distance(mz_data.begin(), mz_it));
        while (mz_it != mz_data.end() && *mz_it < right)
        {
          const double point_mz = *mz_it;
          if (point_mz > left && point_mz < right)
          {
            const double point_im = im_data[idx];
            if (!use_im_filter || (point_im > left_im && point_im < right_im))
            {
              ExtractedPeakMap& peak_map = output[k];
              peak_map.mz.push_back(point_mz);
              peak_map.rt.push_back(current_rt);
              peak_map.ion_mobility.push_back(point_im);
              peak_map.intensity.push_back(int_data[idx]);
            }
          }
          ++mz_it;
          ++idx;
        }
      }
    }
    endProgress();
  }

  int PeakMapExtractor::getFilterNr_(const String& filter) const
  {
    if (filter == "tophat")
    {
      return 1;
    }
    if (filter == "bartlett")
    {
      return 2;
    }
    throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                     "Filter either needs to be tophat or bartlett");
  }
} // namespace OpenMS
