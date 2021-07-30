// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>

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
  };

};

