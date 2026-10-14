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
  #include <OpenMS/INTERFACES/IMSDataConsumer.h>
  #include <algorithm>
  #include <cmath>
  #include <deque>
  #include <exception>
  #include <limits>
  #include <string>
  #include <vector>

namespace OpenMS
{
namespace
{
  using ScanInfo = Internal::ThermoRawScanReader::ScanInfo;

  /// Scan selection by ThermoRawFile::Options, shared by load() and transform().
  struct Selection
  {
    const ThermoRawFile::Options& options;

    bool level(int level) const
    {
      return options.ms_levels.empty() || std::find(options.ms_levels.begin(), options.ms_levels.end(), level) != options.ms_levels.end();
    }
    bool rt(double rt) const { return rt >= options.rt_min && rt < options.rt_max; }
    bool precursorFiltered() const
    {
      return options.precursor_mz_min > -std::numeric_limits<double>::infinity() || options.precursor_mz_max < std::numeric_limits<double>::infinity();
    }
    /// True if scans are selected by values of their header.
    bool headerFiltered() const
    {
      return ! options.ms_levels.empty() || precursorFiltered() || options.rt_min > -std::numeric_limits<double>::infinity()
             || options.rt_max < std::numeric_limits<double>::infinity();
    }
    /// Last scan to read: scans of one controller are in acquisition order, so no later scan can be selected.
    int lastScan(int last) const { return options.last_scan > 0 ? std::min(last, options.last_scan) : last; }
    /// True if no scan from @p info on can be selected (checked after reading its header).
    bool beyondRT(const ScanInfo& info) const { return info.rt >= options.rt_max; }
    bool selected(const ScanInfo& info) const
    {
      if (info.scan < options.first_scan || ! rt(info.rt) || ! level(info.ms_level)) { return false; }
      return info.ms_level <= 1 || std::isnan(info.precursor_mz)
             || (info.precursor_mz >= options.precursor_mz_min && info.precursor_mz < options.precursor_mz_max);
    }
  };

  /// Copies of the most recent converted spectra that later scans can reference as parent.
  class ParentCache
  {
  public:
    explicit ParentCache(Size capacity) : capacity_(capacity) {}

    const MSSpectrum* find(int scan) const
    {
      for (const auto& [number, spectrum] : entries_)
      {
        if (number == scan) { return &spectrum; }
      }
      return nullptr;
    }

    /// Keep the peaks of @p spectrum (the estimate needs nothing else); the oldest entry is evicted.
    void add(int scan, const MSSpectrum& spectrum)
    {
      if (capacity_ == 0) { return; }
      if (entries_.size() >= capacity_) { entries_.pop_front(); }
      entries_.emplace_back(scan, MSSpectrum());
      MSSpectrum& copy = entries_.back().second;
      copy.reserve(spectrum.size());
      copy.insert(copy.end(), spectrum.begin(), spectrum.end());
    }

  private:
    Size capacity_;
    std::deque<std::pair<int, MSSpectrum>> entries_;
  };
} // namespace

void ThermoRawFile::load(const std::string& path, MSExperiment& exp)
{
  exp = MSExperiment();
  const Selection select{options_};
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
      const bool header_filtered = select.headerFiltered();
      if (! header_filtered && options_.first_scan <= first) { exp.reserve(last - first + 1); }
      MSChromatogram tic;
      bool complete = true; // all headers read, so the TIC can be built from them
      const bool prefetch_peaks = ! header_filtered && ! options_.metadata_only;
      reader.readAhead(select.lastScan(last),
                       prefetch_peaks ? options_.first_scan : std::numeric_limits<int>::max());
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
        // Without header-based filters every scan from first_scan on is converted, so its peaks
        // are read together with the header.
        const bool prefetch = prefetch_peaks && scan >= options_.first_scan;
        const auto& info = reader.header(scan, prefetch);
        if (! options_.skip_chromatograms) { tic.push_back(ChromatogramPeak(info.rt, info.tic)); }
        if (select.beyondRT(info))
        {
          complete = false;
          break;
        }
        if (! select.selected(info)) { continue; }
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
      if (select.level(0))
      {
        for (auto& spectrum : spectra)
        {
          if (! select.rt(spectrum.getRT())) { continue; }
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

void ThermoRawFile::transform(const std::string& path, Interfaces::IMSDataConsumer* consumer)
{
  if (consumer == nullptr) { throw Exception::NullPointer(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); }
  const Selection select{options_};
  try
  {
    Internal::ThermoRawScanReader reader(path, options_);
    ExperimentalSettings settings;
    reader.fillSettings(settings);
    const bool has_ms = reader.hasMSData();
    const int first = has_ms ? reader.firstScan() : 0;
    const int last = has_ms ? select.lastScan(reader.lastScan()) : -1;

    // Header pass: count the selected scans and register the instrument configurations, both
    // of which the consumer needs before the first spectrum (e.g. for the mzML header).
    Size spectra_count = 0, chromatogram_count = 0;
    int max_level = 0;
    if (has_ms)
    {
      reader.scanHeaders(first, last, select.precursorFiltered(), [&](const ScanInfo& info) {
        max_level = std::max(max_level, info.ms_level);
        if (select.beyondRT(info)) { return false; }
        if (select.selected(info)) { ++spectra_count; }
        return true;
      });
      reader.addInstrumentConfigurations(settings);
      if (! options_.skip_chromatograms) { ++chromatogram_count; }
    }
    std::vector<MSChromatogram> detector_chromatograms;
    std::vector<MSSpectrum> pda_spectra;
    if (options_.all_detectors)
    {
      reader.readDetectors(detector_chromatograms, pda_spectra);
      if (options_.skip_chromatograms) { detector_chromatograms.clear(); }
      if (! select.level(0)) { pda_spectra.clear(); }
      pda_spectra.erase(std::remove_if(pda_spectra.begin(), pda_spectra.end(), [&](const MSSpectrum& spectrum) { return ! select.rt(spectrum.getRT()); }),
                        pda_spectra.end());
      std::stable_sort(pda_spectra.begin(), pda_spectra.end(), MSSpectrum::RTLess());
      if (options_.metadata_only)
      {
        for (auto& spectrum : pda_spectra) { spectrum.clear(false); }
      }
      spectra_count += pda_spectra.size();
      chromatogram_count += detector_chromatograms.size();
    }
    consumer->setExpectedSize(spectra_count, chromatogram_count);
    consumer->setExperimentalSettings(settings);

    if (has_ms)
    {
      // Recent parent candidates (all scans below the highest MS level); a parent that is not
      // kept is decoded again by the reader, so the estimate does not depend on the capacity.
      ParentCache parents(options_.metadata_only ? 0 : 32);
      const Internal::ThermoRawScanReader::ParentLookup parent = [&parents](int scan) { return parents.find(scan); };
      const bool header_filtered = select.headerFiltered();
      MSChromatogram tic;
      bool complete = last == reader.lastScan(); // all headers read, so the TIC can be built from them
      const bool prefetch_peaks = ! header_filtered && ! options_.metadata_only;
      reader.readAhead(last, prefetch_peaks ? options_.first_scan : std::numeric_limits<int>::max());
      startProgress(first, last, "Streaming Thermo RAW file");
      for (int scan = first; scan <= last; ++scan)
      {
        setProgress(scan);
        const bool prefetch = prefetch_peaks && scan >= options_.first_scan;
        const ScanInfo& info = reader.header(scan, prefetch);
        if (! options_.skip_chromatograms) { tic.push_back(ChromatogramPeak(info.rt, info.tic)); }
        if (select.beyondRT(info))
        {
          complete = false;
          break;
        }
        if (! select.selected(info)) { continue; }
        MSSpectrum spectrum;
        reader.convert(spectrum, settings, parent, ! options_.metadata_only);
        if (info.ms_level < max_level) { parents.add(scan, spectrum); }
        consumer->consumeSpectrum(spectrum);
      }
      endProgress();
      for (auto& spectrum : pda_spectra) { consumer->consumeSpectrum(spectrum); }
      pda_spectra.clear();
      if (! options_.skip_chromatograms)
      {
        if (! complete) { tic = reader.ticChromatogram(); }
        else
        {
          tic.setNativeID("TIC");
          tic.setChromatogramType(ChromatogramSettings::ChromatogramType::TOTAL_ION_CURRENT_CHROMATOGRAM);
          tic.getDataProcessing().push_back(reader.dataProcessing());
        }
        if (! tic.empty()) { consumer->consumeChromatogram(tic); }
      }
    }
    else
    {
      for (auto& spectrum : pda_spectra) { consumer->consumeSpectrum(spectrum); }
    }
    for (auto& chromatogram : detector_chromatograms) { consumer->consumeChromatogram(chromatogram); }
  }
  catch (const Exception::FileNotFound&)
  {
    throw;
  }
  catch (const Exception::BaseException&)
  {
    throw; // e.g. from the consumer
  }
  catch (const std::exception& error)
  {
    throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path, std::string("Thermo RAW reader: ") + error.what());
  }
}

} // namespace OpenMS
#endif // WITH_THERMO_RAW
