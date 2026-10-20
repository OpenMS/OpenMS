// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include "MS1LabeledSpectra.h"

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/METADATA/SpectrumLookup.h>
#ifdef WITH_OPENTIMS
#include <OpenMS/FORMAT/BrukerTimsFile.h>
#endif
#include <algorithm>
#include <cmath>

namespace OpenMS
{
void MS1LabeledSpectra::load(const std::string& filename, MSExperiment& metadata, MSExperiment& ms1, ProgressLogger::LogType log_type)
{
  metadata.clear(true);
  ms1.clear(true);
  FileHandler reader;
  const FileTypes::Type type = FileHandler::getType(filename);
  if (type == FileTypes::RAW || type == FileTypes::BRUKER_TDF)
  {
    // The Thermo and Bruker readers decode all scans even with PeakFileOptions filters. Read once,
    // retain peak-free metadata for ID/FAIMS lookup, and discard MSn peaks before quantification.
    if (type == FileTypes::RAW)
    {
      reader.loadExperiment(filename, ms1, {FileTypes::RAW}, log_type);
    }
    else
    {
#ifdef WITH_OPENTIMS
      // timsTOF MS1 frames hold one peak per (m/z, ion mobility scan). Centroid them in m/z and ion mobility
      // as ProteomicsLFQ does (Sage defaults), then sum each m/z over ion mobility: the multiplet detection
      // works on m/z and RT only.
      BrukerTimsFile tdf;
      tdf.setLogType(log_type);
      BrukerTimsFile::Config config;
      config.ms1_centroid_mz_ppm = 5.0f;
      config.ms1_centroid_im_pct = 3.0f;
      tdf.load(filename, ms1, config);
      for (MSSpectrum& spectrum : ms1)
      {
        if (spectrum.getMSLevel() == 1) { collapseIonMobility(spectrum, 5.0); }
      }
#else
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Bruker .d input requires an OpenMS build with OpenTIMS support (WITH_OPENTIMS).", filename);
#endif
    }
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

void MS1LabeledSpectra::collapseIonMobility(MSSpectrum& spectrum, double mz_tolerance_ppm)
{
  if (!spectrum.containsIMData()) { return; }
  spectrum.getFloatDataArrays().clear();
  spectrum.getStringDataArrays().clear();
  spectrum.getIntegerDataArrays().clear();
  spectrum.setIMFormat(IMFormat::NONE);
  spectrum.setIMPeakType(IMPeakType::UNKNOWN);
  spectrum.sortByPosition();
  // merge runs of peaks within the tolerance of the first peak of the run: summed intensity, intensity-weighted m/z
  Size merged = 0; // merged peaks are written in place: their number never exceeds the number of peaks read
  for (Size begin = 0; begin < spectrum.size();)
  {
    const double first_mz = spectrum[begin].getMZ();
    const double max_mz = first_mz + first_mz * mz_tolerance_ppm * 1e-6;
    double intensity = 0.0, weighted_mz = 0.0;
    Size end = begin;
    for (; end < spectrum.size() && spectrum[end].getMZ() <= max_mz; ++end)
    {
      intensity += spectrum[end].getIntensity();
      weighted_mz += spectrum[end].getMZ() * spectrum[end].getIntensity();
    }
    spectrum[merged].setMZ(intensity > 0.0 ? weighted_mz / intensity : first_mz);
    spectrum[merged].setIntensity(static_cast<float>(intensity));
    ++merged;
    begin = end;
  }
  spectrum.resize(merged);
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
