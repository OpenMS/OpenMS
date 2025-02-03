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
  const std::string NamesOfIMFormat[] = {"none", "concatenated", "multiple_spectra", "mixed"};


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
    if (occs.size() == 1 && (occs.find(IMFormat::CONCATENATED) != occs.end() || occs.find(IMFormat::MULTIPLE_SPECTRA) != occs.end()))
    {
      return *occs.begin();
    }
    if (occs.size() == 2 && occs.find(IMFormat::CONCATENATED) != occs.end() && occs.find(IMFormat::MULTIPLE_SPECTRA) != occs.end())
    {
      return IMFormat::MIXED;
    }

    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "subfunction returned invalid value(s)", "Number of different values: " + String(occs.size()));
  }

  IMFormat IMTypes::determineIMFormat(const MSSpectrum& spec)
  {
    bool has_float_data = spec.containsIMData(); // cache value; query is 'expensive'
    bool has_drift_time = spec.getDriftTime() != DRIFTTIME_NOT_SET;
    if (has_float_data && has_drift_time)
    {
      const auto& fda = spec.getFloatDataArrays()[spec.getIMData().first];
      String array_val = fda.empty() ? "[empty]" : String(fda[0]);
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MSSpectrum contains both an float-data-array and a single drift time. At most one is allowed per spectrum!", String("Array: ") + array_val + ", ... <> Spec: " + spec.getDriftTime());
    }

    if (has_float_data)
    {
      return IMFormat::CONCATENATED;
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
