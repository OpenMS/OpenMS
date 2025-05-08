// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantifierStatistics.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{

  IsobaricQuantifierStatistics::IsobaricQuantifierStatistics() :
    channel_count(0),
    iso_number_ms2_negative(0),
    iso_number_reporter_negative(0),
    iso_number_reporter_different(0),
    iso_solution_different_intensity(0),
    iso_total_intensity_negative(0),
    number_ms2_total(0),
    number_ms2_empty(0),
    empty_channels()
  {
  }

  void IsobaricQuantifierStatistics::reset()
  {
    channel_count = 0;
    iso_number_ms2_negative = 0;
    iso_number_reporter_negative = 0;
    iso_number_reporter_different = 0;
    iso_solution_different_intensity = 0;
    iso_total_intensity_negative = 0;
    number_ms2_total = 0;
    number_ms2_empty = 0;
    empty_channels.clear();
  }

  IsobaricQuantifierStatistics::IsobaricQuantifierStatistics(const IsobaricQuantifierStatistics& other)
  {
    channel_count = other.channel_count;
    iso_number_ms2_negative = other.iso_number_ms2_negative;
    iso_number_reporter_negative = other.iso_number_reporter_negative;
    iso_number_reporter_different = other.iso_number_reporter_different;
    iso_solution_different_intensity = other.iso_solution_different_intensity;
    iso_total_intensity_negative = other.iso_total_intensity_negative;
    number_ms2_total = other.number_ms2_total;
    number_ms2_empty = other.number_ms2_empty;
    empty_channels.clear();
    empty_channels.insert(other.empty_channels.begin(), other.empty_channels.end());
  }

  IsobaricQuantifierStatistics& IsobaricQuantifierStatistics::operator=(const IsobaricQuantifierStatistics& rhs)
  {
    if (this == &rhs) return *this;

    channel_count = rhs.channel_count;
    iso_number_ms2_negative = rhs.iso_number_ms2_negative;
    iso_number_reporter_negative = rhs.iso_number_reporter_negative;
    iso_number_reporter_different = rhs.iso_number_reporter_different;
    iso_solution_different_intensity = rhs.iso_solution_different_intensity;
    iso_total_intensity_negative = rhs.iso_total_intensity_negative;
    number_ms2_total = rhs.number_ms2_total;
    number_ms2_empty = rhs.number_ms2_empty;
    empty_channels.clear();
    empty_channels.insert(rhs.empty_channels.begin(), rhs.empty_channels.end());

    return *this;
  }

} // namespace
