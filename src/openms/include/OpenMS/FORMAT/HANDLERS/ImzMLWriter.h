// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Aditya Sarna $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/HANDLERS/ImzMLHandlerHelper.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{
namespace Internal
{

  /**
    @brief Writer for imzML 1.1.0 files (.imzML + companion .ibd).

    @ingroup FileIO
  */
  class OPENMS_DLLAPI ImzMLWriter
  {
  public:
    /**
      @brief Store an MSExperiment as imzML + .ibd.

      Supports continuous (shared m/z array) and processed (per-spectrum m/z)
      modes. Pixel coordinates are read from @p imzml:x/y/z MetaValues on each
      spectrum (required for export); dataset metadata is taken from experiment
      MetaValues when present.

      Export writes external binary arrays with a 16-byte UUID header in the
      @c .ibd file. Binary precision follows @p PeakFileOptions (@p getMz32Bit,
      @p getIntensity32Bit). @p PeakFileOptions spectrum/peak filters (MS level,
      RT, precursor m/z, m/z and intensity ranges, metadata-only, sort-by-m/z)
      are applied to a temporary copy before export.

      @param[in] imzml_path Path to the output @c .imzML file.
      @param[in] exp        Experiment to store (must contain at least one spectrum).
      @param[in] options    Peak file options (filtering, sort, binary precision).
      @param[in] logger     Progress logger for status output.

      @throws Exception::MissingInformation if @p exp has no spectra or lacks imzml:x/y on any spectrum.
      @throws Exception::InvalidValue if pixel coordinates are invalid or duplicate.
      @throws Exception::InvalidParameter if continuous export is requested but spectra are incompatible.
      @throws Exception::UnableToCreateFile if output files cannot be written.
      @throws Exception::ParseError if binary array serialization fails.
    */
    static void store(const std::string& imzml_path,
                      const MSExperiment& exp,
                      const PeakFileOptions& options,
                      ProgressLogger& logger);
  };

} // namespace Internal
} // namespace OpenMS
