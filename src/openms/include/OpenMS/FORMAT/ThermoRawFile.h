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

    Requires a .NET 8 (or newer) runtime on the machine and the
    openms-thermo-bridge managed runtime files (ThermoWrapperManaged.dll, its
    runtimeconfig.json and the Thermo CommonCore assemblies). These are looked
    up in this order:
      1. the directory named by the OPENMS_THERMO_MANAGED_DIR environment variable
         (skipped with a warning if it does not hold the files),
      2. <OpenMS share dir>/openms_thermo_bridge/managed (installed by OpenMS,
         also shipped inside pyOpenMS wheels),
      3. the bridge's own default, i.e. a 'managed' directory next to the
         openms_thermo_bridge shared library.
    Set DOTNET_ROOT if the .NET runtime is installed in a non-standard location.

    See doc/thermo_raw_mzml.md for the metadata mapping to mzML and how to
    build against a local openms-thermo-bridge checkout.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI ThermoRawFile : public ProgressLogger
  {
  public:
    /// Optional data exported by ThermoRawFileParser, plus lossless vendor metadata.
    struct Options
    {
      bool centroid = false; ///< Preserve acquired representation by default; true matches TRFP peak picking.
      bool charge_data = false; ///< Instrument-assigned centroid peak charge array.
      bool noise_data = false; ///< Noise/baseline arrays and their independent mass coordinates.
      bool all_detectors = false; ///< UV/PDA/Analog/MSAnalog traces and PDA spectra.
      bool preserve_trailers = true; ///< Retain all original scan trailer label/value pairs.
      bool instrument_methods = true; ///< Retain embedded method text once per run.
      bool checksum = true; ///< Compute source SHA-1 for mzML provenance.
    };

    /// Returns current loading options.
    const Options& getOptions() const { return options_; }
    /// Set loading options.
    void setOptions(const Options& options) { options_ = options; }
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
    void load(const std::string& path, MSExperiment& exp);

  private:
    Options options_;
  };

} // namespace OpenMS

#endif // WITH_THERMO_RAW
