// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/DATAACCESS/MSDataSqlConsumer.h>

#include <OpenMS/FORMAT/HANDLERS/MzMLSqliteHandler.h>

namespace OpenMS
{

  MSDataSqlConsumer::MSDataSqlConsumer(const String& filename, UInt64 run_id, int flush_after, bool full_meta, bool lossy_compression, double linear_mass_acc) :
        filename_(filename),
        handler_(new OpenMS::Internal::MzMLSqliteHandler(filename, run_id) ),
        flush_after_(flush_after),
        full_meta_(full_meta)
  {
    spectra_.reserve(flush_after_);
    chromatograms_.reserve(flush_after_);

    handler_->setConfig(full_meta, lossy_compression, linear_mass_acc, flush_after_);
    handler_->createTables();
  }

  MSDataSqlConsumer::~MSDataSqlConsumer()
  {
    flush();

    // Write run level information into the file (e.g. run id, run name and mzML structure)
    peak_meta_.setLoadedFilePath(filename_);
    handler_->writeRunLevelInformation(peak_meta_, full_meta_);

    delete handler_;
  }

  void MSDataSqlConsumer::flush()
  {
    if (!spectra_.empty() ) 
    {
      handler_->writeSpectra(spectra_);
      spectra_.clear();
      spectra_.reserve(flush_after_);
    }

    if (!chromatograms_.empty() ) 
    {
      handler_->writeChromatograms(chromatograms_);
      chromatograms_.clear();
      chromatograms_.reserve(flush_after_);
    }
  }

  void MSDataSqlConsumer::consumeSpectrum(SpectrumType & s)
  {
    spectra_.push_back(s);
    s.clear(false);
    if (full_meta_)
    {
      peak_meta_.addSpectrum(s);
    }
    if (spectra_.size() >= flush_after_)
    {
      flush();
    }
  }

  void MSDataSqlConsumer::consumeChromatogram(ChromatogramType & c)
  {
    chromatograms_.push_back(c);
    c.clear(false);
    if (full_meta_)
    {
      peak_meta_.addChromatogram(c);
    }
    if (chromatograms_.size() >= flush_after_)
    {
      flush();
    }
  }

  void MSDataSqlConsumer::setExpectedSize(Size /* expectedSpectra */, Size /* expectedChromatograms */) {;}

  void MSDataSqlConsumer::setExperimentalSettings(const ExperimentalSettings& /* exp */) {;}

} // namespace OpenMS

