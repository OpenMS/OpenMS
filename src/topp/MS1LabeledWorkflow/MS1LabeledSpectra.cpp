// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include "MS1LabeledSpectra.h"

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/METADATA/SpectrumLookup.h>
#include <algorithm>
#include <cmath>

namespace OpenMS
{
void MS1LabeledSpectra::load(const std::string& filename, MSExperiment& metadata, MSExperiment& ms1, ProgressLogger::LogType log_type)
{
  metadata.clear(true);
  ms1.clear(true);
  FileHandler reader;
  if (FileHandler::getType(filename) == FileTypes::RAW)
  {
    // The Thermo reader decodes all scans even with PeakFileOptions filters. Read once,
    // retain peak-free metadata for ID/FAIMS lookup, and discard MSn peaks before quantification.
    reader.loadExperiment(filename, ms1, {FileTypes::RAW}, log_type);
    static_cast<ExperimentalSettings&>(metadata) = ms1;
    metadata.reserve(ms1.size());
    for (const auto& spectrum : ms1)
    {
      MSSpectrum scan(spectrum);
      scan.clear(false);
      scan.shrink_to_fit();
      metadata.addSpectrum(std::move(scan));
    }
    std::erase_if(ms1.getSpectra(), [](const MSSpectrum& spectrum) { return spectrum.getMSLevel() != 1; });
    ms1.getChromatograms().clear();
    ms1.updateRanges();
  }
  else
  {
    // mzML can skip peak decoding and MSn scans while parsing.
    reader.getOptions().setFillData(false);
    reader.loadExperiment(filename, metadata, {FileTypes::MZML}, log_type);
    reader.getOptions().setFillData(true);
    reader.getOptions().setMSLevels({1});
    reader.loadExperiment(filename, ms1, {FileTypes::MZML}, log_type);
  }
}

void MS1LabeledSpectra::addMissingSpectrumReferences(const MSExperiment& metadata, PeptideIdentificationList& ids)
{
  if (std::none_of(ids.begin(), ids.end(), [](const auto& id) { return id.getSpectrumReference().empty(); })) { return; }
  std::vector<MSSpectrum> ms2;
  for (const auto& spectrum : metadata)
  {
    if (spectrum.getMSLevel() == 2) { ms2.push_back(spectrum); }
  }
  SpectrumLookup lookup;
  lookup.readSpectra(ms2);
  for (auto& id : ids)
  {
    if (! id.getSpectrumReference().empty()) { continue; }
    try
    {
      if (! std::isfinite(id.getRT()))
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Identification has no finite retention time.");
      }
      const auto& native_id = ms2.at(lookup.findByRT(id.getRT())).getNativeID();
      if (native_id.empty())
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Matching MS2 spectrum has no native ID.");
      }
      id.setSpectrumReference(native_id);
    }
    catch (const Exception::ElementNotFound&)
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Cannot repair spectrum reference: no MS2 spectrum matches identification RT " + std::to_string(id.getRT())
                                            + ".");
    }
  }
}
} // namespace OpenMS
