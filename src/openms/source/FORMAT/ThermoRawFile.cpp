// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU
// Berlin SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#ifdef WITH_THERMO_RAW

  #include <OpenMS/CONCEPT/Exception.h>
  #include <OpenMS/FORMAT/HANDLERS/ThermoRawScanReader.h>
  #include <OpenMS/FORMAT/ThermoRawFile.h>
  #include <OpenMS/KERNEL/MSChromatogram.h>
  #include <OpenMS/KERNEL/MSSpectrum.h>
  #include <algorithm>
  #include <cmath>
  #include <exception>
  #include <string>
  #include <vector>

namespace OpenMS
{

void ThermoRawFile::load(const std::string& path, MSExperiment& exp)
{
  exp = MSExperiment();
  const auto& levels = options_.ms_levels;
  auto level_selected = [&levels](int level) { return levels.empty() || std::find(levels.begin(), levels.end(), level) != levels.end(); };
  auto rt_selected = [this](double rt) { return rt >= options_.rt_min && rt < options_.rt_max; };
  try
  {
    Internal::ThermoRawScanReader reader(path, options_);
    reader.fillSettings(exp);
    if (reader.hasMSData())
    {
      const int first = reader.firstScan(), last = reader.lastScan();
      // Converted spectrum of each scan (for the precursor intensity estimate), -1 if not loaded.
      std::vector<SignedSize> index_by_scan(static_cast<Size>(std::max(0, last - first + 1)), -1);
      const Internal::ThermoRawScanReader::ParentLookup parent = [&](int scan) -> const MSSpectrum* {
        if (scan < first || scan > last) { return nullptr; }
        const SignedSize index = index_by_scan[scan - first];
        return index < 0 ? nullptr : &exp[index];
      };
      const bool filtered = ! levels.empty() || options_.first_scan > first || options_.precursor_mz_min > -std::numeric_limits<double>::infinity()
                            || options_.precursor_mz_max < std::numeric_limits<double>::infinity() || options_.rt_min > -std::numeric_limits<double>::infinity();
      if (! filtered) { exp.reserve(last - first + 1); }
      MSChromatogram tic;
      bool complete = true; // all headers read, so the TIC can be built from them
      startProgress(first, last, "Loading Thermo RAW file");
      for (int scan = first; scan <= last; ++scan)
      {
        // Scans of one controller are in acquisition order, so no later scan can be selected.
        if (options_.last_scan > 0 && scan > options_.last_scan)
        {
          complete = false;
          break;
        }
        setProgress(scan);
        const auto& info = reader.header(scan);
        if (! options_.skip_chromatograms) { tic.push_back(ChromatogramPeak(info.rt, info.tic)); }
        if (info.rt >= options_.rt_max)
        {
          complete = false;
          break;
        }
        if (scan < options_.first_scan || ! rt_selected(info.rt) || ! level_selected(info.ms_level)) { continue; }
        if (info.ms_level > 1 && ! std::isnan(info.precursor_mz)
            && (info.precursor_mz < options_.precursor_mz_min || info.precursor_mz >= options_.precursor_mz_max))
        {
          continue;
        }
        MSSpectrum spectrum;
        reader.convert(spectrum, exp, parent, ! options_.metadata_only);
        index_by_scan[scan - first] = exp.size();
        exp.addSpectrum(std::move(spectrum));
      }
      endProgress();
      reader.addInstrumentConfigurations(exp);
      if (! options_.skip_chromatograms)
      {
        if (! complete) { tic = reader.ticChromatogram(); }
        else
        {
          tic.setNativeID("TIC");
          tic.setChromatogramType(ChromatogramSettings::ChromatogramType::TOTAL_ION_CURRENT_CHROMATOGRAM);
          tic.getDataProcessing().push_back(reader.dataProcessing());
        }
        if (! tic.empty()) { exp.addChromatogram(std::move(tic)); }
      }
    }
    if (options_.all_detectors)
    {
      std::vector<MSChromatogram> chromatograms;
      std::vector<MSSpectrum> spectra;
      reader.readDetectors(chromatograms, spectra);
      if (! options_.skip_chromatograms)
      {
        for (auto& chromatogram : chromatograms)
        {
          exp.addChromatogram(std::move(chromatogram));
        }
      }
      if (level_selected(0))
      {
        for (auto& spectrum : spectra)
        {
          if (! rt_selected(spectrum.getRT())) { continue; }
          if (options_.metadata_only) { spectrum.clear(false); }
          exp.addSpectrum(std::move(spectrum));
        }
      }
    }
    // MS scans come in acquisition order; only interleaved PDA spectra need sorting by RT.
    const auto& spectra = exp.getSpectra();
    if (! std::is_sorted(spectra.begin(), spectra.end(), MSSpectrum::RTLess())) { exp.sortSpectra(false); }
    exp.updateRanges();
  }
  catch (const Exception::FileNotFound&)
  {
    throw;
  }
  catch (const std::exception& error)
  {
    throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path, std::string("Thermo RAW reader: ") + error.what());
  }
}

} // namespace OpenMS
#endif // WITH_THERMO_RAW
