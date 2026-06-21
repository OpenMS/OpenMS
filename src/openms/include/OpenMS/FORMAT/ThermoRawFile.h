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
    RawFileReader through the .NET host runtime.

    Per-spectrum metadata populated (mapped to PSI-MS CV terms):
    - Retention time, MS level, precursor m/z and isolation window
    - Polarity (mapped from the Thermo PolarityType enum: Negative=0, Positive=1)
    - Total ion current, base-peak m/z and intensity (MS:1000285/MS:1000504/MS:1000505)
    - Lowest/highest observed m/z (MS:1000528/MS:1000527)
    - Scan filter string (stored as a UserParam)
    - Scan window (scanWindowLower/scanWindowUpper)
    - Ion injection time as an Acquisition element (MS:1000927)
    - FAIMS compensation voltage (driftTime + FAIMS_COMPENSATION_VOLTAGE user param)

    Experiment-level metadata:
    - Creation date (experiment.setDateTime())
    - Sample name

    Instrument block:
    - Ion source type (electrospray, nanoelectrospray, etc.)
    - Mass analyzer type(s) per scan (Orbitrap/FTMS disambiguated via model name, IT, Q, TOF, sector)
    - Ion detector type(s)
    - Software version

    Chromatograms:
    - A total ion current (TIC) chromatogram derived from the spectra

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

      Reads all scans and populates @p exp with the full set of spectra and
      metadata documented in the class description, including per-spectrum
      CV-term annotations (TIC, base peak, ion injection time, scan filter,
      FAIMS CV), a TIC chromatogram, and instrument/sample/software metadata.

      @param[in] path Path to the .raw file
      @param[out] exp The experiment to populate

      @throws Exception::FileNotFound if the file does not exist
      @throws Exception::ParseError if the file cannot be read by the thermo bridge
    */
    void load(const std::string& path, MSExperiment& exp);
  };

} // namespace OpenMS

#endif // WITH_THERMO_RAW
