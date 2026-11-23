// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>

namespace OpenMS
{
/** @brief Load MS1 peaks and scan metadata for the MS1 labeling workflow. */
class MS1LabeledSpectra
{
public:
  /// Read mzML or Thermo RAW (requires WITH_THERMO_RAW), retaining all scan metadata but only MS1 peaks.
  static void load(const std::string& filename, MSExperiment& metadata, MSExperiment& ms1, ProgressLogger::LogType log_type);

  /// Repair absent/empty references from MS2 RTs. Call FAIMS CV annotation/repair first for FAIMS data.
  /// Existing references are preserved; missing matches throw MissingInformation.
  static void addMissingSpectrumReferences(const MSExperiment& metadata, PeptideIdentificationList& ids);
};
} // namespace OpenMS
