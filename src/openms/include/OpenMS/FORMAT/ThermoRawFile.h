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
    RawFileReader through the .NET host runtime.  The output MSExperiment is
    populated with the same metadata that msconvert/ProteoWizard emits,
    mapped to PSI-MS CV terms:

    - **Per spectrum**: total ion current, base peak m/z and intensity,
      lowest/highest observed m/z, scan filter string, scan window (parsed
      from the filter), ion injection time (stored as MS:1000927 on an
      Acquisition so it serialises inside \c \<scan\>), and FAIMS
      compensation voltage (driftTime + FAIMS_COMPENSATION_VOLTAGE).
    - **Instrument**: ion source, mass analyser(s), ion detector(s), and
      software version.
    - **Experiment**: creation date and sample name.
    - **Chromatogram**: a total ion current (TIC) chromatogram.

    Spectrum polarity is mapped from the Thermo PolarityType enum
    (Negative=0, Positive=1, Any=2) to the corresponding OpenMS values.
    Activation method is mapped independently of collision energy, so
    zero-CE ETD scans retain their activation type; PQD is also handled.

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
      spectra, chromatograms, and the full set of mzML-equivalent metadata
      described in the class documentation above.

      @param[in] path Path to the .raw file
      @param[out] exp The experiment to populate

      @throws Exception::FileNotFound if the file does not exist
      @throws Exception::ParseError if the file cannot be read by the thermo bridge
    */
    void load(const std::string& path, MSExperiment& exp);
  };

} // namespace OpenMS

#endif // WITH_THERMO_RAW
