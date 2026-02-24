// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/TargetedDataFileLoader.h>
#include <OpenMS/FORMAT/SwathFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MRMFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CONCEPT/LogStream.h>

using namespace std;

namespace OpenMS
{

std::vector<::OpenSwath::SwathMap> TargetedDataFileLoader::loadFile(const String& file,
                                                                  const String& tmp,
                                                                  std::shared_ptr<ExperimentalSettings>& exp_meta,
                                                                  const String& readoptions,
                                                                  Interfaces::IMSDataConsumer* plugin_consumer)
{
  // Quick probe: load metadata only (no data) to check for presence of spectra
  PeakMap probe;
  try
  {
    FileHandler fh;
    fh.getOptions().setFillData(false);
    fh.loadExperiment(file, probe, {FileTypes::MZML});
  }
  catch (...) {
    // If probing fails, fall back to SwathFile loader
    OPENMS_LOG_DEBUG << "TargetedDataFileLoader: probe failed for " << file << ", falling back to SwathFile::loadMzML" << std::endl;
    SwathFile sw;
    sw.setLogType(ProgressLogger::LogType::NONE);
    return sw.loadMzML(file, tmp, exp_meta, readoptions, plugin_consumer);
  }

  // If there are no spectra but chromatograms exist, treat as SRM/MRM/chrom-only
  if (probe.getSpectra().empty() && !probe.getChromatograms().empty())
  {
    OPENMS_LOG_DEBUG << "TargetedDataFileLoader: detected chromatogram-only mzML -> using SRM/MRM loader" << std::endl;
    return MRMFile::loadMzML(file, tmp, exp_meta);
  }

  // Otherwise use SwathFile loader
  SwathFile sw;
  sw.setLogType(ProgressLogger::LogType::NONE);
  return sw.loadMzML(file, tmp, exp_meta, readoptions, plugin_consumer);
}

} // namespace OpenMS
