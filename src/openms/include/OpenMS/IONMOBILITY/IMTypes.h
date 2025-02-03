// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/CommonEnums.h>
#include <OpenMS/OpenMSConfig.h>

#include <string>

namespace OpenMS
{
  class MSExperiment;
  class MSSpectrum;

  /// Drift time unit for ion mobility
  enum class DriftTimeUnit
  {
    NONE,                      ///< No unit
    MILLISECOND,               ///< milliseconds
    VSSC,                      ///< volt-second per square centimeter (i.e. 1/K_0)
    FAIMS_COMPENSATION_VOLTAGE,///< compensation voltage
    SIZE_OF_DRIFTTIMEUNIT
  };

  /// Names of IM Units. Should be usable as axis annotation.
  OPENMS_DLLAPI extern const std::string NamesOfDriftTimeUnit[(size_t) DriftTimeUnit::SIZE_OF_DRIFTTIMEUNIT];

  /// convert an entry in NamesOfDriftTimeUnit[] to DriftTimeUnit enum
  /// @throws Exception::InvalidValue if @p dtu_string is not contained in NamesOfDriftTimeUnit[]
  OPENMS_DLLAPI DriftTimeUnit toDriftTimeUnit(const std::string& dtu_string);

  /// convert a DriftTimeUnit enum to String
  /// @throws Exception::InvalidValue if @p value is SIZE_OF_DRIFTTIMEUNIT
  OPENMS_DLLAPI const std::string& toString(const DriftTimeUnit value);

  /// Different ways to represent ion mobility data in a spectrum
  enum class IMFormat
  {
    NONE,            ///< no ion mobility
    CONCATENATED,    ///< ion mobility frame is stacked in a single spectrum (i.e. has an IM float data array)
    MULTIPLE_SPECTRA,///< ion mobility is recorded as multiple spectra per frame (i.e. has one IM annotation per spectrum)
    MIXED,           ///< an MSExperiment contains both CONCATENATED and MULTIPLE_SPECTRA
    SIZE_OF_IMFORMAT
  };
  /// Names of IMFormat
  OPENMS_DLLAPI extern const std::string NamesOfIMFormat[(size_t) IMFormat::SIZE_OF_IMFORMAT];
  
  /// convert an entry in NamesOfIMFormat[] to IMFormat enum
  /// @throws Exception::InvalidValue if @p IM_format is not contained in NamesOfIMFormat[]
  OPENMS_DLLAPI IMFormat toIMFormat(const std::string& IM_format);
  /// convert an IMFormat enum to String
  /// @throws Exception::InvalidValue if @p value is SIZE_OF_IMFORMAT
  OPENMS_DLLAPI const std::string& toString(const IMFormat value);

  class OPENMS_DLLAPI IMTypes
  {
  public:
    /// If drift time for a spectrum is unavailable (i.e. not an IM spectrum), it will have this value
    inline static constexpr double DRIFTTIME_NOT_SET = -1.0;

    /// Checks the all spectra for their type (see overload)
    /// and returns the common type (or IMFormat::MIXED if both CONCATENATED and MULTIPLE_SPECTRA are present)
    /// If @p exp is empty or contains no IM spectra at all, IMFormat::NONE is returned
    /// @throws Exception::InvalidValue if IM values are annotated as single drift time and float array for any single spectrum
    static IMFormat determineIMFormat(const MSExperiment& exp);

    /** 
        @brief Checks for existence of a single driftTime (using spec.getDriftTime()) or an ion-mobility float data array (using spec.hasIMData()) 
        
        If neither is found, IMFormat::NONE is returned.
        If a single drift time (== IMFormat::MULTIPLE_SPECTRA) is found, but no unit, a warning is issued.

        @throws Exception::InvalidValue if IM values are annotated as single drift time and float array in the given spectrum
    */
    static IMFormat determineIMFormat(const MSSpectrum& spec);

    /**
     * \brief 
     * \param from Drift unit to convert from
     * \return A more general DIM_UNIT (or exception)
     * \throws Exception::ConversionError if @p from has invalid value (e.g. 'NONE')
     */
    static DIM_UNIT fromIMUnit(const DriftTimeUnit from);
  };

};

