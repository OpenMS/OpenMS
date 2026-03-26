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
    CCS,                       ///< collisional cross section (square angstrom)
    SIZE_OF_DRIFTTIMEUNIT
  };

  /// Names of IM Units. Should be usable as axis annotation.
  OPENMS_DLLAPI extern const std::string NamesOfDriftTimeUnit[(size_t) DriftTimeUnit::SIZE_OF_DRIFTTIMEUNIT];

  /// convert an entry in NamesOfDriftTimeUnit[] to DriftTimeUnit enum
  /// @throws Exception::InvalidValue if @p dtu_string is not contained in NamesOfDriftTimeUnit[]
  OPENMS_DLLAPI DriftTimeUnit toDriftTimeUnit(const std::string& dtu_string);

  /// convert a DriftTimeUnit enum to String
  /// @throws Exception::InvalidValue if @p value is SIZE_OF_DRIFTTIMEUNIT
  OPENMS_DLLAPI const std::string& driftTimeUnitToString(const DriftTimeUnit value);

  /// Different ways to represent ion mobility data in a spectrum
  /// Note: UNKNOWN should be used if the format is not yet determined.
  /// FileHandler or e.g. IM peak picker should ideally set the format a known value.
  enum class IMFormat
  {
    NONE,            ///< no ion mobility
    IM_PEAK,         ///< full TIMS frame / per-scan IM-resolved data: ion mobility is annotated per peak in a float data array
    IM_SPECTRUM,     ///< conventional spectrum with one precursor IM value (i.e. has one IM annotation per spectrum via getDriftTime())
    UNKNOWN,         ///< ion mobility format not yet determined.
    SIZE_OF_IMFORMAT
  };
  /// Names of IMFormat
  OPENMS_DLLAPI extern const std::string NamesOfIMFormat[(size_t) IMFormat::SIZE_OF_IMFORMAT];
  
  /// convert an entry in NamesOfIMFormat[] to IMFormat enum
  /// @throws Exception::InvalidValue if @p IM_format is not contained in NamesOfIMFormat[]
  OPENMS_DLLAPI IMFormat toIMFormat(const std::string& IM_format);
  /// convert an IMFormat enum to String
  /// @throws Exception::InvalidValue if @p value is SIZE_OF_IMFORMAT
  OPENMS_DLLAPI const std::string& imFormatToString(const IMFormat value);

  /// Processing state of ion mobility data in the IM dimension.
  /// Analogous to SpectrumSettings::SpectrumType for the m/z dimension.
  enum class IMPeakType
  {
    IM_PROFILE,        ///< raw/unprocessed IM data (e.g. full TIMS frame before IM centroiding)
    IM_CENTROIDED,     ///< IM data after centroiding in the IM dimension
    UNKNOWN,           ///< IM peak type not yet determined
    SIZE_OF_IMPEAKTYPE
  };

  /// Names of IMPeakType entries
  OPENMS_DLLAPI extern const std::string NamesOfIMPeakType[(size_t) IMPeakType::SIZE_OF_IMPEAKTYPE];

  /// @brief converts a string to IMPeakType
  /// @param[in] im_peak_type string representation (e.g. "im_profile", "im_centroided", "unknown")
  /// @return the corresponding IMPeakType enum value
  /// @throws Exception::InvalidValue if the string does not match any known IMPeakType
  OPENMS_DLLAPI IMPeakType toIMPeakType(const std::string& im_peak_type);

  /// @brief converts an IMPeakType to a string
  /// @param[in] im_peak_type the enum value to convert
  /// @return reference to the string representation
  /// @throws Exception::InvalidValue if the enum value is out of range
  OPENMS_DLLAPI const std::string& imPeakTypeToString(IMPeakType im_peak_type);

  class OPENMS_DLLAPI IMTypes
  {
  public:
    /// If drift time for a spectrum is unavailable (i.e. not an IM spectrum), it will have this value
    inline static constexpr double DRIFTTIME_NOT_SET = -1.0;

    /// Checks only spectra of the given MS level for their IM format and returns the common type.
    /// If no spectra of @p ms_level exist or none have IM data, IMFormat::NONE is returned.
    /// @throws Exception::InvalidValue if spectra of the given MS level have different IM formats
    static IMFormat determineIMFormat(const MSExperiment& exp, int ms_level);

    /** 
        @brief Checks for existence of a single driftTime (using spec.getDriftTime()) or an ion-mobility float data array (using spec.hasIMData()) 
        
        If neither is found, IMFormat::NONE is returned.
        If a single drift time (== IMFormat::IM_SPECTRUM) is found, but no unit, a warning is issued.

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

