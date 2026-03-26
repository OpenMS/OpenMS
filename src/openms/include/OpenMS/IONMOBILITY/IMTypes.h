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
  /// Note: 
  /// 1. MIXED is only used for MSExperiment, not for MSSpectrum
  /// 2. UNKNOWN should be used if the format is not yet determined. 
  /// FileHandler or e.g. IM peak picker should ideally set the format a known value.
  enum class IMFormat
  {
    NONE,            ///< no ion mobility
    IM_PEAK,         ///< full TIMS frame / per-scan IM-resolved data: ion mobility is annotated per peak in a float data array
    IM_SPECTRUM,     ///< conventional spectrum with one precursor IM value (i.e. has one IM annotation per spectrum via getDriftTime())
    MIXED,           ///< an MSExperiment contains both IM_PEAK and IM_SPECTRUM
    CENTROIDED,      ///< @deprecated Use IMPeakType::IM_CENTROIDED instead. Will be removed.
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

  /// Convert string to IMPeakType
  OPENMS_DLLAPI IMPeakType toIMPeakType(const std::string& im_peak_type);

  /// Convert IMPeakType to string
  OPENMS_DLLAPI const std::string& imPeakTypeToString(IMPeakType im_peak_type);

  class OPENMS_DLLAPI IMTypes
  {
  public:
    /// If drift time for a spectrum is unavailable (i.e. not an IM spectrum), it will have this value
    inline static constexpr double DRIFTTIME_NOT_SET = -1.0;

    /// Checks the all spectra for their type (see overload)
    /// and returns the common type (or IMFormat::MIXED if both IM_PEAK and IM_SPECTRUM are present)
    /// If @p exp is empty or contains no IM spectra at all, IMFormat::NONE is returned
    /// @throws Exception::InvalidValue if IM values are annotated as single drift time and float array for any single spectrum
    static IMFormat determineIMFormat(const MSExperiment& exp);

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

