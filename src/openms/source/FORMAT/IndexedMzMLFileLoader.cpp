// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IndexedMzMLFileLoader.h>

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>

namespace OpenMS
{

  IndexedMzMLFileLoader::IndexedMzMLFileLoader() = default;

  IndexedMzMLFileLoader::~IndexedMzMLFileLoader() = default;

  PeakFileOptions & IndexedMzMLFileLoader::getOptions()
  {
    return options_;
  }

  const PeakFileOptions & IndexedMzMLFileLoader::getOptions() const
  {
    return options_;
  }

  void IndexedMzMLFileLoader::setOptions(const PeakFileOptions & options)
  {
      options_ = options;
  }

  bool IndexedMzMLFileLoader::load(const String& filename, OnDiscPeakMap& exp)
  {
    return exp.openFile(filename);
  }

  void IndexedMzMLFileLoader::store(const String& filename, OnDiscPeakMap& exp)
  {
    // Create a writing data consumer which consumes the experiment (writes it to disk)
    PlainMSDataWritingConsumer consumer(filename);
    consumer.setExpectedSize(exp.getNrSpectra(), exp.getNrChromatograms());
    consumer.setExperimentalSettings(*exp.getExperimentalSettings().get());
    options_.setWriteIndex(true);  // ensure that we write the index
    consumer.setOptions(options_);
    for (Size i = 0; i < exp.getNrSpectra(); i++)
    {
      MSSpectrum s = exp.getSpectrum(i);
      consumer.consumeSpectrum(s);
    }
    for (Size i = 0; i < exp.getNrChromatograms(); i++)
    {
      MSChromatogram c = exp.getChromatogram(i);
      consumer.consumeChromatogram(c);
    }
  }

  void IndexedMzMLFileLoader::store(const String& filename, PeakMap& exp)
  {
    MzMLFile f;
    options_.setWriteIndex(true);  // ensure that we write the index
    f.setOptions(options_);
    f.store(filename, exp);
  }
}
