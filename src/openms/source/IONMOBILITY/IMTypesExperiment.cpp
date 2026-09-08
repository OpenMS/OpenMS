// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <set>

namespace OpenMS
{
  IMFormat IMTypes::determineIMFormat(const MSExperiment& exp, int ms_level)
  {
    std::set<IMFormat> occs;
    for (const auto& spec : exp.getSpectra())
    {
      if (spec.getMSLevel() != ms_level) continue;
      occs.insert(determineIMFormat(spec));
    }
    occs.erase(IMFormat::NONE);

    if (occs.empty())
    {
      return IMFormat::NONE;
    }

    if (occs.size() == 1)
    {
      auto format = *occs.begin();
      if (format != IMFormat::IM_PEAK && format != IMFormat::IM_SPECTRUM)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "subfunction returned invalid value(s)", "Number of different values: " + StringUtils::toStr(occs.size()));
      }
      return format;
    }
    else
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "MSExperiment contains MS" + StringUtils::toStr(ms_level) + " spectra with different IM formats. "
        "Handle per-spectrum.", "Number of different formats: " + StringUtils::toStr(occs.size()));
    }
  }

  IMFormat IMTypes::determineIMFormat(const MSSpectrum& spec)
  {
    // First check if format is already set and not UNKNOWN
    IMFormat current_format = spec.getIMFormat();
    if (current_format != IMFormat::UNKNOWN)
    {
      return current_format;
    }
    
    // If format is UNKNOWN, determine it
    bool has_float_data = spec.containsIMData(); // cache value; query is 'expensive'
    bool has_drift_time = spec.getDriftTime() != DRIFTTIME_NOT_SET;

    if (has_float_data)
    {
      if (has_drift_time)
      {
        OPENMS_LOG_DEBUG << "both drift time and IM data array found in spectrum " << spec.getNativeID() << "\n. Support for both is experimental." << std::endl;
      }
      return IMFormat::IM_PEAK;
    }
    else if (has_drift_time)
    {
      if (spec.getDriftTimeUnit() == DriftTimeUnit::NONE)
      {
        OPENMS_LOG_WARN << "Warning: no drift time unit set for spectrum " << spec.getNativeID() << "\n";
      }
      return IMFormat::IM_SPECTRUM;
    }
    return IMFormat::NONE;
  }

} // namespace OpenMS
