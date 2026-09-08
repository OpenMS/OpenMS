// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/IONMOBILITY/IMTypes.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <algorithm>
#include <iterator>

#include <cmath>
#include <cstdlib>

namespace OpenMS
{

  const std::string NamesOfDriftTimeUnit[] = {"<NONE>", "ms", "1/K0", "FAIMS_CV", "CCS"};
  const std::string NamesOfIMFormat[] = {"none", "im_peak", "im_spectrum", "unknown"};


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

  const std::string& driftTimeUnitToString(const DriftTimeUnit value)
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

  const std::string& imFormatToString(const IMFormat value)
  {
    if (value == IMFormat::SIZE_OF_IMFORMAT)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value not allowed", "SIZE_OF_IMFORMAT");
    }
    return NamesOfIMFormat[(size_t)value];
  }

  const std::string NamesOfIMPeakType[] = {"im_profile", "im_centroided", "unknown"};

  IMPeakType toIMPeakType(const std::string& im_peak_type)
  {
    auto idx = std::find(NamesOfIMPeakType, NamesOfIMPeakType + (int)IMPeakType::SIZE_OF_IMPEAKTYPE, im_peak_type);
    if (idx == NamesOfIMPeakType + (int)IMPeakType::SIZE_OF_IMPEAKTYPE)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid IMPeakType", im_peak_type);
    }
    return (IMPeakType)std::distance(NamesOfIMPeakType, idx);
  }

  const std::string& imPeakTypeToString(IMPeakType im_peak_type)
  {
    if ((size_t)im_peak_type >= (size_t)IMPeakType::SIZE_OF_IMPEAKTYPE)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid IMPeakType index", std::to_string((int)im_peak_type));
    }
    return NamesOfIMPeakType[(int)im_peak_type];
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
      case DriftTimeUnit::CCS:
        return DIM_UNIT::IM_CCS;
      default:
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Cannot convert from " + driftTimeUnitToString(from) + " to a DIM_UNIT.");
    }
  }

  namespace
  {
    /// Bruker Mason-Schamp calibration constant relating 1/K0 [V*s/cm^2] and CCS [Angstrom^2] for an N2
    /// drift gas; value confirmed against alphatims and MaxQuant CCS values (also used by OpenNuXL).
    /// CCS = (MASON_SCHAMP_CONSTANT * |z| / sqrt(reduced_mass)) * (1/K0).
    constexpr double MASON_SCHAMP_CONSTANT = 1059.62245;

    /// Ion-gas reduced mass [Da] with the ion mass approximated as mz * |charge|.
    double reducedMass_(double mz, int charge, double buffer_gas_mass)
    {
      const double ion_mass = mz * std::abs(charge);
      return (ion_mass * buffer_gas_mass) / (ion_mass + buffer_gas_mass);
    }
  }

  double IMTypes::oneOverK0ToCCS(double one_over_k0, double mz, int charge, double buffer_gas_mass)
  {
    if (one_over_k0 <= 0.0 || mz <= 0.0 || charge == 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "oneOverK0ToCCS requires one_over_k0 > 0, mz > 0 and charge != 0",
        "1/K0=" + StringUtils::toStr(one_over_k0) + ", mz=" + StringUtils::toStr(mz) + ", charge=" + StringUtils::toStr(charge));
    }
    const double mu = reducedMass_(mz, charge, buffer_gas_mass);
    return MASON_SCHAMP_CONSTANT * std::abs(charge) / std::sqrt(mu) * one_over_k0;
  }

  double IMTypes::ccsToOneOverK0(double ccs, double mz, int charge, double buffer_gas_mass)
  {
    if (ccs <= 0.0 || mz <= 0.0 || charge == 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "ccsToOneOverK0 requires ccs > 0, mz > 0 and charge != 0",
        "CCS=" + StringUtils::toStr(ccs) + ", mz=" + StringUtils::toStr(mz) + ", charge=" + StringUtils::toStr(charge));
    }
    const double mu = reducedMass_(mz, charge, buffer_gas_mass);
    return ccs * std::sqrt(mu) / (MASON_SCHAMP_CONSTANT * std::abs(charge));
  }
}// namespace OpenMS
