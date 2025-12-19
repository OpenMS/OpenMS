// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/OnDiscMSExperiment.h>

#include <OpenMS/FORMAT/FileHandler.h>

namespace OpenMS
{
  bool OnDiscMSExperiment::openFile(const String& filename, bool skipMetaData)
  {
    filename_ = filename;
    indexed_mzml_file_.openFile(filename);
    if (!filename.empty() && !skipMetaData)
    {
      loadMetaData_(filename);
    }
    return indexed_mzml_file_.getParsingSuccess();
  }

  void OnDiscMSExperiment::setSkipXMLChecks(bool skip)
  {
    indexed_mzml_file_.setSkipXMLChecks(skip);
  }

  PeakFileOptions& OnDiscMSExperiment::getOptions()
  {
    return options_;
  }

  const PeakFileOptions& OnDiscMSExperiment::getOptions() const
  {
    return options_;
  }

  void OnDiscMSExperiment::setOptions(const PeakFileOptions& options)
  {
    options_ = options;
  }

  OpenMS::Interfaces::ChromatogramPtr OnDiscMSExperiment::getChromatogramById(Size id)
  {
    return indexed_mzml_file_.getChromatogramById(id);
  }

  void OnDiscMSExperiment::loadMetaData_(const String& filename)
  {
    meta_ms_experiment_ = std::shared_ptr< PeakMap >(new PeakMap);

    FileHandler f;
    // Always load ALL metadata without RT/MS level filters to maintain proper index correspondence
    // with indexed_mzml_file_. Filtering is applied at retrieval time in getSpectrum()/getChromatogram().
    // We only set fillData=false to avoid loading peak data during metadata loading.
    PeakFileOptions options;
    options.setFillData(false);
    f.setOptions(options);
    f.loadExperiment(filename, *meta_ms_experiment_.get(), {FileTypes::MZML});
  }

  MSChromatogram OnDiscMSExperiment::getMetaChromatogramById_(const std::string& id)
  {
    if (chromatograms_native_ids_.empty())
    {
      for (Size k = 0; k < meta_ms_experiment_->getChromatograms().size(); k++)
      {
        chromatograms_native_ids_.emplace(meta_ms_experiment_->getChromatograms()[k].getNativeID(), k);
      }
    }

    if (chromatograms_native_ids_.find(id) == chromatograms_native_ids_.end())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          String("Could not find chromatogram with id '") + id + "'.");
    }
    return meta_ms_experiment_->getChromatogram(chromatograms_native_ids_[id]);
  }

  MSChromatogram OnDiscMSExperiment::getChromatogramByNativeId(const std::string& id)
  {
    if (!meta_ms_experiment_)
    {
      MSChromatogram chromatogram;
      indexed_mzml_file_.getMSChromatogramByNativeId(id, chromatogram);
      return chromatogram;
    }

    MSChromatogram chromatogram = getMetaChromatogramById_(id);
    indexed_mzml_file_.getMSChromatogramByNativeId(id, chromatogram);
    return chromatogram;
  }

  MSSpectrum OnDiscMSExperiment::getMetaSpectrumById_(const std::string& id)
  {
    if (spectra_native_ids_.empty())
    {
      for (Size k = 0; k < meta_ms_experiment_->getSpectra().size(); k++)
      {
        spectra_native_ids_.emplace(meta_ms_experiment_->getSpectra()[k].getNativeID(), k);
      }
    }

    if (spectra_native_ids_.find(id) == spectra_native_ids_.end())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          String("Could not find spectrum with id '") + id + "'.");
    }
    return meta_ms_experiment_->getSpectrum(spectra_native_ids_[id]);
  }

  MSSpectrum OnDiscMSExperiment::getSpectrumByNativeId(const std::string& id)
  {
    if (!meta_ms_experiment_)
    {
      MSSpectrum spec;
      indexed_mzml_file_.getMSSpectrumByNativeId(id, spec);
      return spec;
    }

    MSSpectrum spec = getMetaSpectrumById_(id);
    indexed_mzml_file_.getMSSpectrumByNativeId(id, spec);
    return spec;
  }

  MSSpectrum OnDiscMSExperiment::getSpectrum(Size id)
  {
    MSSpectrum spectrum;

    // If we have metadata, we can check RT range and MS level filters BEFORE loading peak data
    // This allows skipping expensive I/O for spectra that will be filtered out anyway
    if (meta_ms_experiment_)
    {
      spectrum = meta_ms_experiment_->operator[](id);

      // Check RT range filter - skip loading peaks if outside range
      if (options_.hasRTRange() && !options_.getRTRange().encloses(DPosition<1>(spectrum.getRT())))
      {
        // Return spectrum with metadata but no peaks (filtered out)
        return spectrum;
      }

      // Check MS level filter - skip loading peaks if MS level doesn't match
      if (options_.hasMSLevels() && !options_.containsMSLevel(spectrum.getMSLevel()))
      {
        // Return spectrum with metadata but no peaks (filtered out)
        return spectrum;
      }

      // Check precursor m/z range filter - skip loading peaks if precursor m/z is outside range
      // This only applies to MS2+ spectra that have precursor information
      if (options_.hasPrecursorMZRange() && spectrum.getMSLevel() >= 2 && !spectrum.getPrecursors().empty())
      {
        double precursor_mz = spectrum.getPrecursors()[0].getMZ();
        if (!options_.getPrecursorMZRange().encloses(DPosition<1>(precursor_mz)))
        {
          // Return spectrum with metadata but no peaks (filtered out)
          return spectrum;
        }
      }

      // Spectrum passes metadata filters, load peak data
      indexed_mzml_file_.getMSSpectrumById(int(id), spectrum);
    }
    else
    {
      // No metadata available, must load full spectrum
      spectrum = indexed_mzml_file_.getMSSpectrumById(int(id));
    }

    // Apply m/z and intensity range filters if set (requires peak data)
    if (options_.hasMZRange() || options_.hasIntensityRange())
    {
      MSSpectrum filtered;
      filtered.SpectrumSettings::operator=(spectrum);
      filtered.reserve(spectrum.size());

      for (const auto& peak : spectrum)
      {
        bool pass_mz = !options_.hasMZRange() || options_.getMZRange().encloses(DPosition<1>(peak.getMZ()));
        bool pass_int = !options_.hasIntensityRange() || options_.getIntensityRange().encloses(DPosition<1>(peak.getIntensity()));
        if (pass_mz && pass_int)
        {
          filtered.push_back(peak);
        }
      }
      return filtered;
    }

    return spectrum;
  }

  MSChromatogram OnDiscMSExperiment::getChromatogram(Size id)
  {
    MSChromatogram chromatogram;
    if (!meta_ms_experiment_)
    {
      chromatogram = indexed_mzml_file_.getMSChromatogramById(int(id));
    }
    else
    {
      chromatogram = meta_ms_experiment_->getChromatogram(id);
      indexed_mzml_file_.getMSChromatogramById(int(id), chromatogram);
    }

    // Apply RT and intensity range filters if set (RT range for chromatograms filters on RT dimension)
    if (options_.hasRTRange() || options_.hasIntensityRange())
    {
      MSChromatogram filtered;
      filtered.ChromatogramSettings::operator=(chromatogram);
      filtered.reserve(chromatogram.size());

      for (const auto& peak : chromatogram)
      {
        bool pass_rt = !options_.hasRTRange() || options_.getRTRange().encloses(DPosition<1>(peak.getRT()));
        bool pass_int = !options_.hasIntensityRange() || options_.getIntensityRange().encloses(DPosition<1>(peak.getIntensity()));
        if (pass_rt && pass_int)
        {
          filtered.push_back(peak);
        }
      }
      return filtered;
    }

    return chromatogram;
  }

} //namespace OpenMS

