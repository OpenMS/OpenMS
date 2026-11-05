// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/OnDiscThermoRawExperiment.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/HANDLERS/ThermoRawScanReader.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>

#include <algorithm>
#include <exception>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <vector>

namespace OpenMS
{
  using ScanReader = Internal::ThermoRawScanReader;

  /// Immutable after openFile(); shared by all copies.
  struct OnDiscThermoRawExperiment::Index
  {
    struct Entry
    {
      double rt = 0.0;
      UInt ms_level = 0;
      bool pda = false;
      Size pos = 0; ///< into records or pda
    };

    std::string filename;
    Options options;
    std::shared_ptr<ExperimentalSettings> settings = std::make_shared<ExperimentalSettings>();
    std::vector<ScanReader::ScanRecord> records;
    std::vector<MSSpectrum> pda;
    std::vector<Entry> entries;
    std::vector<MSChromatogram> chromatograms;
    std::unordered_map<std::string, Size> by_native_id;
    int first_scan = 0;
    std::vector<SignedSize> by_scan; ///< scan - first_scan -> spectrum index (-1 if not indexed)
    bool sorted = true;
    /// Only used to clone accessors, which is thread-safe.
    std::unique_ptr<ScanReader::Accessor> prototype;
  };

  /// Per-instance file handle.
  struct OnDiscThermoRawExperiment::Handle
  {
    std::unique_ptr<ScanReader::Accessor> accessor;
  };

  namespace
  {
    /// Run @p f and report errors of the vendor library as Exception::ParseError, as ThermoRawFile does.
    template <typename F>
    auto guarded(const std::string& path, F&& f)
    {
      try
      {
        return f();
      }
      catch (const Exception::BaseException&)
      {
        throw;
      }
      catch (const std::exception& error)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path, std::string("Thermo RAW reader: ") + error.what());
      }
    }
  } // namespace

  OnDiscThermoRawExperiment::OnDiscThermoRawExperiment() = default;

  OnDiscThermoRawExperiment::OnDiscThermoRawExperiment(const OnDiscThermoRawExperiment& source) :
    IMSExperimentRandomAccess(source),
    index_(source.index_)
  {
    if (index_ && index_->prototype)
    {
      handle_ = std::make_unique<Handle>();
      handle_->accessor = guarded(index_->filename, [&] { return index_->prototype->clone(); });
    }
  }

  OnDiscThermoRawExperiment::OnDiscThermoRawExperiment(OnDiscThermoRawExperiment&& source) noexcept = default;

  OnDiscThermoRawExperiment& OnDiscThermoRawExperiment::operator=(const OnDiscThermoRawExperiment& source)
  {
    if (this != &source)
    {
      OnDiscThermoRawExperiment copy(source);
      *this = std::move(copy);
    }
    return *this;
  }

  OnDiscThermoRawExperiment& OnDiscThermoRawExperiment::operator=(OnDiscThermoRawExperiment&& source) noexcept = default;

  OnDiscThermoRawExperiment::~OnDiscThermoRawExperiment() = default;

  void OnDiscThermoRawExperiment::openFile(const std::string& path, const Options& options)
  {
    auto index = std::make_shared<Index>();
    index->filename = path;
    index->options = options;
    index->options.checksum = false;
    const Internal::ThermoRawScanSelection select{index->options};
    const Options& opt = index->options;
    guarded(path, [&] {
      ScanReader reader(path, opt);
      reader.fillSettings(*index->settings);
      if (reader.hasMSData())
      {
        const int first = reader.firstScan(), last = reader.lastScan();
        index->first_scan = first;
        index->by_scan.assign(static_cast<Size>(std::max(0, last - first + 1)), -1);
        MSChromatogram tic;
        bool complete = true; // all headers read, so the TIC can be built from them
        // Same scan loop as ThermoRawFile::load(), reading headers only.
        reader.readAhead(select.lastScan(last), std::numeric_limits<int>::max());
        for (int scan = first; scan <= last; ++scan)
        {
          if (opt.last_scan > 0 && scan > opt.last_scan)
          {
            complete = false;
            break;
          }
          const auto& info = reader.header(scan, false);
          if (! opt.skip_chromatograms) { tic.push_back(ChromatogramPeak(info.rt, info.tic)); }
          if (select.beyondRT(info))
          {
            complete = false;
            break;
          }
          if (! select.selected(info)) { continue; }
          index->entries.push_back(Index::Entry{info.rt, static_cast<UInt>(info.ms_level), false, index->records.size()});
          index->records.push_back(reader.record());
        }
        reader.addInstrumentConfigurations(*index->settings);
        if (! opt.skip_chromatograms)
        {
          if (! complete) { tic = reader.ticChromatogram(); }
          else
          {
            tic.setNativeID("TIC");
            tic.setChromatogramType(ChromatogramSettings::ChromatogramType::TOTAL_ION_CURRENT_CHROMATOGRAM);
            tic.getDataProcessing().push_back(reader.dataProcessing());
          }
          if (! tic.empty())
          {
            tic.updateRanges();
            index->chromatograms.push_back(std::move(tic));
          }
        }
        index->prototype = reader.accessor();
      }
      if (opt.all_detectors)
      {
        std::vector<MSChromatogram> chromatograms;
        std::vector<MSSpectrum> spectra;
        reader.readDetectors(chromatograms, spectra);
        if (! opt.skip_chromatograms)
        {
          for (auto& chromatogram : chromatograms)
          {
            chromatogram.updateRanges();
            index->chromatograms.push_back(std::move(chromatogram));
          }
        }
        if (select.level(0))
        {
          for (auto& spectrum : spectra)
          {
            if (! select.rt(spectrum.getRT())) { continue; }
            if (opt.metadata_only) { spectrum.clear(false); }
            spectrum.updateRanges();
            index->entries.push_back(Index::Entry{spectrum.getRT(), 0, true, index->pda.size()});
            index->pda.push_back(std::move(spectrum));
          }
        }
      }
      return 0;
    });
    // MS scans come in acquisition order; only interleaved PDA spectra need sorting (stable, as in load()).
    auto rt_less = [](const Index::Entry& a, const Index::Entry& b) { return a.rt < b.rt; };
    if (! std::is_sorted(index->entries.begin(), index->entries.end(), rt_less))
    {
      std::stable_sort(index->entries.begin(), index->entries.end(), rt_less);
    }
    index->sorted = true;
    index->by_native_id.reserve(index->entries.size());
    for (Size i = 0; i < index->entries.size(); ++i)
    {
      const Index::Entry& entry = index->entries[i];
      if (entry.pda)
      {
        index->by_native_id.emplace(index->pda[entry.pos].getNativeID(), i);
      }
      else
      {
        const auto& record = index->records[entry.pos];
        index->by_native_id.emplace(record.native_id, i);
        index->by_scan[record.info.scan - index->first_scan] = static_cast<SignedSize>(i);
      }
    }
    auto handle = std::make_unique<Handle>();
    if (index->prototype) { handle->accessor = guarded(path, [&] { return index->prototype->clone(); }); }
    index_ = std::move(index);
    handle_ = std::move(handle);
  }

  bool OnDiscThermoRawExperiment::isOpen() const
  {
    return index_ != nullptr;
  }

  const std::string& OnDiscThermoRawExperiment::getFilename() const
  {
    static const std::string empty;
    return index_ ? index_->filename : empty;
  }

  const OnDiscThermoRawExperiment::Options& OnDiscThermoRawExperiment::getOptions() const
  {
    static const Options defaults;
    return index_ ? index_->options : defaults;
  }

  Size OnDiscThermoRawExperiment::getNrSpectra() const
  {
    return index_ ? index_->entries.size() : 0;
  }

  Size OnDiscThermoRawExperiment::getNrChromatograms() const
  {
    return index_ ? index_->chromatograms.size() : 0;
  }

  void OnDiscThermoRawExperiment::checkIndex_(Size index) const
  {
    if (index >= getNrSpectra())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, static_cast<SignedSize>(index), getNrSpectra());
    }
  }

  MSSpectrum OnDiscThermoRawExperiment::getSpectrum(Size index)
  {
    checkIndex_(index);
    const Index::Entry& entry = index_->entries[index];
    if (entry.pda) { return index_->pda[entry.pos]; }
    MSSpectrum spectrum;
    guarded(index_->filename, [&] {
      handle_->accessor->convert(index_->records[entry.pos], spectrum, ! index_->options.metadata_only);
      return 0;
    });
    spectrum.updateRanges();
    return spectrum;
  }

  MSSpectrum OnDiscThermoRawExperiment::getSpectrumMetadata(Size index)
  {
    checkIndex_(index);
    const Index::Entry& entry = index_->entries[index];
    MSSpectrum spectrum;
    if (entry.pda)
    {
      spectrum = index_->pda[entry.pos];
      spectrum.clear(false);
    }
    else
    {
      guarded(index_->filename, [&] {
        handle_->accessor->convert(index_->records[entry.pos], spectrum, false);
        return 0;
      });
    }
    spectrum.updateRanges();
    return spectrum;
  }

  MSSpectrum OnDiscThermoRawExperiment::getSpectrumPeaks(Size index)
  {
    checkIndex_(index);
    const Index::Entry& entry = index_->entries[index];
    if (entry.pda) { return index_->pda[entry.pos]; }
    const auto& record = index_->records[entry.pos];
    MSSpectrum spectrum;
    if (! index_->options.metadata_only)
    {
      guarded(index_->filename, [&] {
        handle_->accessor->readPeaks(record.info.scan, spectrum);
        return 0;
      });
    }
    spectrum.setRT(record.info.rt);
    spectrum.setMSLevel(record.info.ms_level);
    spectrum.setNativeID(record.native_id);
    spectrum.updateRanges();
    return spectrum;
  }

  MSSpectrum OnDiscThermoRawExperiment::getSpectrumByNativeId(const std::string& native_id)
  {
    if (index_)
    {
      auto it = index_->by_native_id.find(native_id);
      if (it != index_->by_native_id.end()) { return getSpectrum(it->second); }
    }
    throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, native_id);
  }

  double OnDiscThermoRawExperiment::getSpectrumRT(Size index) const
  {
    checkIndex_(index);
    return index_->entries[index].rt;
  }

  UInt OnDiscThermoRawExperiment::getSpectrumMSLevel(Size index) const
  {
    checkIndex_(index);
    return index_->entries[index].ms_level;
  }

  const std::string& OnDiscThermoRawExperiment::getSpectrumNativeID(Size index) const
  {
    checkIndex_(index);
    const Index::Entry& entry = index_->entries[index];
    return entry.pda ? index_->pda[entry.pos].getNativeID() : index_->records[entry.pos].native_id;
  }

  int OnDiscThermoRawExperiment::getScanNumber(Size index) const
  {
    checkIndex_(index);
    const Index::Entry& entry = index_->entries[index];
    return entry.pda ? 0 : index_->records[entry.pos].info.scan;
  }

  Size OnDiscThermoRawExperiment::getSpectrumIndexByScanNumber(int scan) const
  {
    if (! index_ || scan < index_->first_scan || static_cast<Size>(scan - index_->first_scan) >= index_->by_scan.size()) { return getNrSpectra(); }
    const SignedSize i = index_->by_scan[scan - index_->first_scan];
    return i < 0 ? getNrSpectra() : static_cast<Size>(i);
  }

  bool OnDiscThermoRawExperiment::isSortedByRT() const
  {
    return true;
  }

  Size OnDiscThermoRawExperiment::findFirstRT(double rt) const
  {
    if (! index_) { return 0; }
    const auto& entries = index_->entries;
    auto it = std::lower_bound(entries.begin(), entries.end(), rt, [](const Index::Entry& entry, double value) { return entry.rt < value; });
    return static_cast<Size>(std::distance(entries.begin(), it));
  }

  MSChromatogram OnDiscThermoRawExperiment::getChromatogram(Size index)
  {
    if (index >= getNrChromatograms())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, static_cast<SignedSize>(index), getNrChromatograms());
    }
    return index_->chromatograms[index];
  }

  std::shared_ptr<const ExperimentalSettings> OnDiscThermoRawExperiment::getExperimentalSettings() const
  {
    if (! index_) { return std::make_shared<const ExperimentalSettings>(); }
    return index_->settings;
  }

  std::unique_ptr<IMSExperimentRandomAccess> OnDiscThermoRawExperiment::clone() const
  {
    return std::make_unique<OnDiscThermoRawExperiment>(*this);
  }
} // namespace OpenMS
