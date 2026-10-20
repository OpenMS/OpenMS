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
  /// Read mzML, Thermo RAW (requires WITH_THERMO_RAW) or Bruker timsTOF .d / .d.zip (requires WITH_OPENTIMS),
  /// retaining all scan metadata but only MS1 peaks. Bruker MS1 peaks are summed over ion mobility (see collapseIonMobility()).
  static void load(const std::string& filename, MSExperiment& metadata, MSExperiment& ms1, ProgressLogger::LogType log_type);

  /// Remove the ion mobility of the peaks of @p spectrum: peaks within @p mz_tolerance_ppm of the first peak of a run
  /// (in m/z order) are merged into one with the summed intensity and the intensity-weighted m/z. All data arrays are
  /// removed. A spectrum without ion mobility data is left unchanged.
  static void collapseIonMobility(MSSpectrum& spectrum, double mz_tolerance_ppm);

  /// Repair absent/empty references from MS2 RTs. Call FAIMS CV annotation/repair first for FAIMS data.
  /// Existing references are preserved; missing matches throw MissingInformation.
  static void addMissingSpectrumReferences(const MSExperiment& metadata, PeptideIdentificationList& ids);
};
} // namespace OpenMS
