// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------
//

#pragma once

#include <OpenMS/DATASTRUCTURES/Param.h>

#include <cmath>

namespace OpenMS
{
  /**
    @brief The parameters of SignalToNoiseEstimatorMedian as plain values; the defaults are those of its parameters

    Kept apart from the estimator, so that classes that estimate many containers (e.g. PeakPickerHiRes) can hold the
    settings as a member and call SignalToNoiseEstimatorMedian::estimate() without building an estimator object.
    The type does not depend on the container type of the estimator.
  */
  struct SignalToNoiseEstimatorMedianSettings
  {
    double max_intensity = -1;
    double auto_max_stdev_factor = 3.0;
    double auto_max_percentile = 95;
    int auto_mode = 0; ///< -1: 'max_intensity', 0: mean + 'auto_max_stdev_factor' * stdev, 1: 'auto_max_percentile'
    double win_len = 200.0;
    int bin_count = 30;
    int min_required_elements = 10;
    double noise_for_empty_window = std::pow(10.0, 20);
    bool write_log_messages = true;

    /// Reads the settings from the parameters of a SignalToNoiseEstimatorMedian (e.g. a "SignalToNoise:" subsection)
    static SignalToNoiseEstimatorMedianSettings fromParam(const Param& p)
    {
      // the same conversions as SignalToNoiseEstimatorMedian::updateMembers_()
      SignalToNoiseEstimatorMedianSettings s;
      s.max_intensity = (double)p.getValue("max_intensity");
      s.auto_max_stdev_factor = (double)p.getValue("auto_max_stdev_factor");
      s.auto_max_percentile = p.getValue("auto_max_percentile");
      s.auto_mode = p.getValue("auto_mode");
      s.win_len = (double)p.getValue("win_len");
      s.bin_count = p.getValue("bin_count");
      s.min_required_elements = p.getValue("min_required_elements");
      s.noise_for_empty_window = (double)p.getValue("noise_for_empty_window");
      s.write_log_messages = (bool)p.getValue("write_log_messages").toBool();
      return s;
    }
  };
} // namespace OpenMS
