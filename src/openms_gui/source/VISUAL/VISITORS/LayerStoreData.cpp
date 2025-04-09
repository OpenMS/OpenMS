// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/VISITORS/LayerStoreData.h>

#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/MetaInfo.h>
#include <OpenMS/SYSTEM/File.h>

using namespace std;

namespace OpenMS
{
  FileTypes::Type LayerStoreData::getSupportedExtension_(const String& filename) const
  {
    auto type = FileHandler::getTypeByFileName(filename);
    if (type == FileTypes::UNKNOWN)
      return storage_formats_.getTypes().front();
    if (!storage_formats_.contains(type))
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "Format is not supported.");
    return type;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // helper for saving a peakmap to a file
  void savePeakMapToFile(const String& path, const PeakMap& pm, const ProgressLogger::LogType lt, const FileTypes::Type /*ext*/)
  {
    FileHandler().storeExperiment(path, pm, {}, lt);
  }

  void LayerStoreDataPeakMapVisible::saveToFile(const String& path, const ProgressLogger::LogType lt) const
  {
    savePeakMapToFile(path, pm_, lt, getSupportedExtension_(path));
  }

  // helper to filter a single MS1 spectrum
  // Returns true if filtered spectrum contains data
  bool filterSpectrum(const MSSpectrum& in, MSSpectrum& out, const RangeAllType& visible_range, const DataFilters& layer_filters)
  {
    out = in;
    out.clear(false); // keep metadata
    auto it_end = in.MZEnd(visible_range.getMaxMZ());
    for (auto it = in.MZBegin(visible_range.getMinMZ()); it != it_end; ++it)
    {
      if (layer_filters.passes(in, it - in.begin()))
      {
        out.push_back(*it);
      }
    }
    return !out.empty();
  }

  void LayerStoreDataPeakMapVisible::storeVisibleSpectrum(const MSSpectrum& spec, const RangeAllType& visible_range, const DataFilters& layer_filters)
  {
    pm_.clear(true);
    MSSpectrum filtered;
    if (filterSpectrum(spec, filtered, visible_range, layer_filters))
    {
      pm_.addSpectrum(filtered);
    }
  }

  // helper to filter a single MSChromatogram
  // Returns true if filtered chromatogram contains data
  bool filterChrom(const MSChromatogram& in, MSChromatogram& out, const RangeAllType& visible_range, const DataFilters& layer_filters)
  {
    out = in;
    out.clear(false); // keep metadata
    auto it_end = in.RTEnd(visible_range.getMaxRT());
    for (auto it = in.RTBegin(visible_range.getMinRT()); it != it_end; ++it)
    {
      if (layer_filters.passes(in, it - in.begin()))
      {
        out.push_back(*it);
      }
    }
    return !out.empty();
  }

  void LayerStoreDataPeakMapVisible::storeVisibleChromatogram(const MSChromatogram& chrom, const RangeAllType& visible_range, const DataFilters& layer_filters)
  {
    pm_.clear(true);
    MSChromatogram filtered;
    if (filterChrom(chrom, filtered, visible_range, layer_filters))
    {
      pm_.addChromatogram(filtered);
    }
  }

  void LayerStoreDataPeakMapVisible::storeVisibleExperiment(const PeakMap& exp, const RangeAllType& visible_range, const DataFilters& layer_filters)
  {
    pm_.clear(true);
    // copy experimental settings
    pm_.ExperimentalSettings::operator=(exp);
    // get begin / end of the range
    auto peak_start = exp.begin();
    auto begin = exp.RTBegin(visible_range.getMinRT());
    auto end = exp.RTEnd(visible_range.getMaxRT());
    Size begin_idx = std::distance(peak_start, begin);
    Size end_idx = std::distance(peak_start, end);

    // reserve space for the correct number of spectra in RT range
    pm_.reserve(end - begin);
    // copy spectra
    for (Size idx = begin_idx; idx < end_idx; ++idx)
    {
      const MSSpectrum& spectrum_ref = exp[idx];

      // MS^n (n>1) spectra are copied if their precursor is in the m/z range
      if (spectrum_ref.getMSLevel() > 1 && !spectrum_ref.getPrecursors().empty())
      {
        if (visible_range.containsMZ(spectrum_ref.getPrecursors()[0].getMZ()))
        {
          pm_.addSpectrum(spectrum_ref);
        }
      }
      else
      {
        // MS1 spectra are cropped to the m/z range
        MSSpectrum filtered;
        if (filterSpectrum(spectrum_ref, filtered, visible_range, layer_filters))
        {
          pm_.addSpectrum(filtered);
        }
      }
      // do not use map.addSpectrum() here, otherwise empty spectra which did not pass the filters above will be added
    }
  }

  void LayerStoreDataPeakMapAll::saveToFile(const String& path, const ProgressLogger::LogType lt) const
  {
    savePeakMapToFile(path, *full_exp_, lt, getSupportedExtension_(path));
  }

  void LayerStoreDataPeakMapAll::storeFullExperiment(const PeakMap& exp)
  {
    full_exp_ = &exp;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // helper for saving a FeatureMap to a file
  void saveFeatureMapToFile(const String& path, const FeatureMap& fm, const ProgressLogger::LogType lt, const FileTypes::Type /*ext*/)
  {
    FileHandler().storeFeatures(path, fm, {FileTypes::FEATUREXML}, lt);
  }

  void LayerStoreDataFeatureMapVisible::saveToFile(const String& path, const ProgressLogger::LogType lt) const
  {
    saveFeatureMapToFile(path, fm_, lt, this->getSupportedExtension_(path));
  }

  void LayerStoreDataFeatureMapVisible::storeVisibleFM(const FeatureMap& fm, const RangeAllType& visible_range, const DataFilters& layer_filters)
  {
    // clear output experiment
    fm_.clear(true);

    // copy meta data
    fm_.setIdentifier(fm.getIdentifier());
    fm_.setProteinIdentifications(fm.getProteinIdentifications());
    // copy features
    for (auto it = fm.begin(); it != fm.end(); ++it)
    {
      if (layer_filters.passes(*it) && visible_range.containsRT(it->getRT()) && visible_range.containsMZ(it->getMZ()))
      {
        fm_.push_back(*it);
      }
    }
  }

  void LayerStoreDataFeatureMapAll::saveToFile(const String& path, const ProgressLogger::LogType lt) const
  {
    saveFeatureMapToFile(path, *full_fm_, lt, this->getSupportedExtension_(path));
  }

  void LayerStoreDataFeatureMapAll::storeFullFM(const FeatureMap& fm)
  {
    full_fm_ = &fm;
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // helper for saving a ConsensusMap to a file
  void saveConsensusMapToFile(const String& path, const ConsensusMap& fm, const ProgressLogger::LogType lt, const FileTypes::Type /*ext*/)
  {
    FileHandler().storeConsensusFeatures(path, fm, {FileTypes::CONSENSUSXML}, lt);
  }

  void LayerStoreDataConsensusMapVisible::saveToFile(const String& path, const ProgressLogger::LogType lt) const
  {
    return saveConsensusMapToFile(path, cm_, lt, this->getSupportedExtension_(path));
  }

  void LayerStoreDataConsensusMapVisible::storeVisibleCM(const ConsensusMap& cm, const RangeAllType& visible_range, const DataFilters& layer_filters)
  {
    // clear output experiment
    cm_.clear(true);

    // copy file descriptions
    cm_.getColumnHeaders() = cm.getColumnHeaders();
    // copy features
    for (auto it = cm.begin(); it != cm.end(); ++it)
    {
      if (layer_filters.passes(*it) && visible_range.containsRT(it->getRT()) && visible_range.containsMZ(it->getMZ()))
      {
        cm_.push_back(*it);
      }
    }
  }

  void LayerStoreDataConsensusMapAll::saveToFile(const String& path, const ProgressLogger::LogType lt) const
  {
    return saveConsensusMapToFile(path, *full_cm_, lt, this->getSupportedExtension_(path));
  }

  void LayerStoreDataConsensusMapAll::storeFullCM(const ConsensusMap& cm)
  {
    full_cm_ = &cm;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // helper for saving a PepIDs to a file
  void savePepIdsToFile(const String& path, const IPeptideIds::PepIds& ids, const ProgressLogger::LogType lt, const FileTypes::Type /*ext*/)
  {
    FileHandler().storeIdentifications(path, {}, ids, {FileTypes::IDXML}, lt);
  }

  void LayerStoreDataIdentVisible::saveToFile(const String& path, const ProgressLogger::LogType lt) const
  {
    return savePepIdsToFile(path, ids_, lt, this->getSupportedExtension_(path));
  }

  void LayerStoreDataIdentVisible::storeVisibleIdent(const IPeptideIds::PepIds& ids, const RangeAllType& visible_range, const DataFilters& /*layer_filters*/)
  {
    ids_.clear();

    // copy peptides, if visible
    for (const auto& p : ids)
    {
      double rt = p.getRT();
      double mz = p.getMZ();
      // TODO: if (layer.filters.passes(*it) && ...)
      if (visible_range.containsRT(rt) && visible_range.containsMZ(mz))
      {
        ids_.push_back(p);
      }
    }
  }

  void LayerStoreDataIdentAll::saveToFile(const String& path, const ProgressLogger::LogType lt) const
  {
    return savePepIdsToFile(path, *full_ids_, lt, this->getSupportedExtension_(path));
  }

  void LayerStoreDataIdentAll::storeFullIdent(const IPeptideIds::PepIds& ids)
  {
    full_ids_ = &ids;
  }
} // namespace OpenMS