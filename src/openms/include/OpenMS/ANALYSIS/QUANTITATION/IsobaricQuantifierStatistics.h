// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>

#include <map>

namespace OpenMS
{
  class String;

  /**
    @brief Statistics for quantitation performance and comparison of NNLS vs. naive method (aka matrix inversion)
   */
  class OPENMS_DLLAPI IsobaricQuantifierStatistics
  {
  public:
    /**
     @brief Create stats object.
     */
    IsobaricQuantifierStatistics();

    /**
     @brief Reset statistics object.
     */
    void reset();

    Size channel_count; ///< 4plex, 6plex, or 8 plex?!
    Size iso_number_ms2_negative; ///< number of MS2 spectra where one or more channels had negative solution
    Size iso_number_reporter_negative; ///< number of channels where naive solution was negative
    Size iso_number_reporter_different; ///< number of channels >0 where naive solution was different; happens when naive solution is negative in other channels
    double iso_solution_different_intensity; ///< absolute intensity difference between both solutions (for channels > 0)
    double iso_total_intensity_negative; ///< only for spectra where naive solution is negative
    Size number_ms2_total; ///< total number of MS2 spectra
    Size number_ms2_empty; ///< number of empty MS2 (no reporters at all)
    std::map<String, Size> empty_channels; ///< Channel_ID -> Missing; indicating the number of empty channels from all MS2 scans, i.e., numbers are between number_ms2_empty and number_ms2_total

    /// Copy c'tor
    IsobaricQuantifierStatistics(const IsobaricQuantifierStatistics& other);

    /// Assignment operator
    IsobaricQuantifierStatistics& operator=(const IsobaricQuantifierStatistics& rhs);
  };
} // namespace
