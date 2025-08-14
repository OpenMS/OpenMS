// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/IONMOBILITY/IMTypes.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{

  const std::string NamesOfDriftTimeUnit[] = {"<NONE>", "ms", "1/K0", "FAIMS_CV"};
  const std::string NamesOfIMFormat[] = {"none", "concatenated", "multiple_spectra", "mixed", "centroided", "unknown"};


 DriftTimeUnit toDriftTimeUnit(const std::string& dtu_string)
  {
    auto first = &NamesOfDriftTimeUnit[0];
    auto last = &NamesOfDriftTimeUnit[(size_t) DriftTimeUnit::SIZE_OF_DRIFTTIMEUNIT];
    const auto it = std::find(first, last, dtu_string);
    if (it == last)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value unknown", dtu_string);
    }
    return DriftTimeUnit(it - first);
  }

  const std::string& toString(const DriftTimeUnit value)
  {
    if (value == DriftTimeUnit::SIZE_OF_DRIFTTIMEUNIT)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value not allowed", "SIZE_OF_DRIFTTIMEUNIT");
    }
    return NamesOfDriftTimeUnit[(size_t) value];
  }

  IMFormat toIMFormat(const std::string& IM_format)
  {
    auto first = &NamesOfIMFormat[0];
    auto last = &NamesOfIMFormat[(size_t) IMFormat::SIZE_OF_IMFORMAT];
    const auto it = std::find(first, last, IM_format);
    if (it == last)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value unknown", IM_format);
    }
    return IMFormat(it - first);
  }

  const std::string& toString(const IMFormat value)
  {
    if (value == IMFormat::SIZE_OF_IMFORMAT)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value not allowed", "SIZE_OF_IMFORMAT");
    }
    return NamesOfIMFormat[(size_t)value];
  }

  IMFormat IMTypes::determineIMFormat(const MSExperiment& exp)
  {
    std::set<IMFormat> occs;
    for (const auto& spec : exp.getSpectra())
    {
      occs.insert(determineIMFormat(spec));
    }
    occs.erase(IMFormat::NONE); // ignore NONE (i.e. normal spectra)

    if (occs.empty())
    {
      return IMFormat::NONE;
    }    

    if (occs.size() == 1) 
    {
      auto format = *occs.begin();
      if (format != IMFormat::CONCATENATED
          && format != IMFormat::MULTIPLE_SPECTRA
          && format != IMFormat::CENTROIDED)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "subfunction returned invalid value(s)", "Number of different values: " + String(occs.size()));
      }
      return format;
    }
    else    
    {
      return IMFormat::MIXED;
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
      return IMFormat::CONCATENATED; // TODO: or CENTROIDED. for now assume that no picking was done (otherwise we would have annotated it)
    }
    else if (has_drift_time)
    {
      if (spec.getDriftTimeUnit() == DriftTimeUnit::NONE)
      {
        OPENMS_LOG_WARN << "Warning: no drift time unit set for spectrum " << spec.getNativeID() << "\n";
      }
      return IMFormat::MULTIPLE_SPECTRA;
    }
    return IMFormat::NONE;
  }

  DIM_UNIT IMTypes::fromIMUnit(const DriftTimeUnit from)
  {
    switch (from)
    {
      case DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE:
        return DIM_UNIT::FAIMS_CV;
      case DriftTimeUnit::MILLISECOND:
        return DIM_UNIT::IM_MS;
      case DriftTimeUnit::VSSC:
        return DIM_UNIT::IM_VSSC;
      default:
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Cannot convert from " + toString(from) + " to a DIM_UNIT.");
    }
  }
}// namespace OpenMS
