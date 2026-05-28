// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#ifdef WITH_THERMO_RAW

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

namespace OpenMS
{
  /**
    @brief Reader for Thermo RAW files via the openms-thermo-bridge.

    Uses the openms-thermo-bridge library to access Thermo Fisher
    RawFileReader through the .NET host runtime. Reads spectra (MS1 and MSn)
    including retention times, precursor information, and instrument metadata.

    Requires the openms-thermo-bridge managed runtime files
    (ThermoWrapperManaged.dll and its dependencies) to be installed alongside
    the OpenMS binaries.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI ThermoRawFile : public ProgressLogger
  {
  public:
    /**
      @brief Load a Thermo RAW file into an MSExperiment.

      Reads all scans from the RAW file and populates the experiment with
      spectra, retention times, MS levels, precursor information (for MSn),
      and source file metadata.

      @param[in] path Path to the .raw file
      @param[out] exp The experiment to populate

      @throws Exception::FileNotFound if the file does not exist
      @throws Exception::ParseError if the file cannot be read by the thermo bridge
    */
    void load(const String& path, MSExperiment& exp);
  };

} // namespace OpenMS

#endif // WITH_THERMO_RAW
