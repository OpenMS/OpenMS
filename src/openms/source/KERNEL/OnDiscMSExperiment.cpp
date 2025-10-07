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

  OpenMS::Interfaces::ChromatogramPtr OnDiscMSExperiment::getChromatogramById(Size id)
  {
    return indexed_mzml_file_.getChromatogramById(id);
  }

  void OnDiscMSExperiment::loadMetaData_(const String& filename)
  {
    meta_ms_experiment_ = boost::shared_ptr< PeakMap >(new PeakMap);

    FileHandler f;
    PeakFileOptions options = f.getOptions();
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

} //namespace OpenMS

