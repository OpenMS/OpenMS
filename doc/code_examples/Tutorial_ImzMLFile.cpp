// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// Demonstrates all imzML load modes in OpenMS.
//
// Usage:
//   Tutorial_ImzMLFile [file.imzML]
//
// When no path is given, uses the bundled example dataset
// doc/code_examples/data/Tutorial_ImzMLFile.imzML.

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/ImzMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/ImzMLHandlerHelper.h>
#include <OpenMS/IMAGING/MSImagingExperiment.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/OnDiscImzMLExperiment.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/openms_data_path.h> // OPENMS_DOC_PATH (location of bundled tutorial data)

#include <algorithm>
#include <iostream>

using namespace OpenMS;

namespace
{
  void printPeaks(const MSSpectrum& spec, Size max_peaks = 3)
  {
    const Size n = std::min(spec.size(), max_peaks);
    for (Size i = 0; i < n; ++i)
    {
      std::cout << "    m/z=" << spec[i].getMZ() << "  I=" << spec[i].getIntensity() << "\n";
    }
    if (spec.size() > n)
    {
      std::cout << "    ... (" << spec.size() << " peaks total)\n";
    }
  }

  void printMetaValueIfPresent(const MSExperiment& exp, const char* key)
  {
    if (exp.metaValueExists(key))
    {
      std::cout << "  " << key << ": " << exp.getMetaValue(key) << "\n";
    }
    else
    {
      std::cout << "  " << key << ": (not present)\n";
    }
  }

  void printDatasetMetadata(const MSExperiment& exp, const ImzMLMeta& meta, Size index_entries)
  {
    std::cout << "--- MSExperiment MetaValues ---\n";
    printMetaValueIfPresent(exp, "imzml:imaging_mode");
    printMetaValueIfPresent(exp, "imzml:max_count_x");
    printMetaValueIfPresent(exp, "imzml:max_count_y");
    printMetaValueIfPresent(exp, "imzml:max_count_z");
    printMetaValueIfPresent(exp, "imzml:pixel_size_x");
    printMetaValueIfPresent(exp, "imzml:pixel_size_y");
    printMetaValueIfPresent(exp, "imzml:max_dim_x");
    printMetaValueIfPresent(exp, "imzml:max_dim_y");
    printMetaValueIfPresent(exp, "imzml:uuid");
    printMetaValueIfPresent(exp, "imzml:ibd_md5");
    printMetaValueIfPresent(exp, "imzml:ibd_sha1");
    printMetaValueIfPresent(exp, "imzml:mz_data_type");
    printMetaValueIfPresent(exp, "imzml:int_data_type");
    printMetaValueIfPresent(exp, "imzml:scan_pattern");
    printMetaValueIfPresent(exp, "imzml:scan_direction");
    printMetaValueIfPresent(exp, "imzml:line_scan_direction");
    printMetaValueIfPresent(exp, "imzml:polarity");
    printMetaValueIfPresent(exp, "imzml:ibd_path");

    std::cout << "--- ImzMLMeta (loadSpectraIndex) ---\n";
    std::cout << "  imaging_mode: " << meta.imaging_mode << "\n";
    std::cout << "  grid: " << meta.max_count_x << " x " << meta.max_count_y
              << " x " << meta.max_count_z << "\n";
    std::cout << "  pixel_size (um): " << meta.pixel_size_x << " x " << meta.pixel_size_y << "\n";
    std::cout << "  index entries: " << index_entries << "\n";

    if (!exp.empty() && exp[0].metaValueExists("imzml:x"))
    {
      std::cout << "--- spectrum[0] pixel ---\n";
      std::cout << "  imzml:x/y/z: " << exp[0].getMetaValue("imzml:x") << ", "
                << exp[0].getMetaValue("imzml:y") << ", "
                << exp[0].getMetaValue("imzml:z") << "\n";
      if (!exp[0].empty())
      {
        std::cout << "  first peak m/z: " << exp[0][0].getMZ() << "\n";
      }
    }
  }
} // namespace

int main(int argc, const char** argv)
{
  // Default to the bundled example dataset; an explicit path may be passed as argv[1].
  const std::string path = (argc >= 2)
                             ? std::string(argv[1])
                             : OPENMS_DOC_PATH + std::string("/code_examples/data/Tutorial_ImzMLFile.imzML");
  ImzMLFile loader;
  loader.setLogType(ProgressLogger::CMD);

  // ------------------------------------------------------------------
  // Mode 1: full load into MSExperiment
  // ------------------------------------------------------------------
  std::cout << "=== Mode 1: ImzMLFile::load -> MSExperiment ===\n";
  MSExperiment exp;
  loader.load(path, exp);
  std::cout << "spectra: " << exp.getNrSpectra() << "\n";

  ImzMLMeta dataset_meta;
  std::vector<ImzMLSpectrumIndex> dataset_index;
  loader.loadSpectraIndex(path, dataset_meta, dataset_index);
  printDatasetMetadata(exp, dataset_meta, dataset_index.size());

  if (exp.getNrSpectra() > 0)
  {
    std::cout << "spectrum[0] peaks:\n";
    printPeaks(exp[0]);
  }

  // ------------------------------------------------------------------
  // Mode 2: FileHandler
  // ------------------------------------------------------------------
  std::cout << "\n=== Mode 2: FileHandler::loadExperiment ===\n";
  MSExperiment fh_exp;
  FileHandler().loadExperiment(path, fh_exp);
  std::cout << "spectra via FileHandler: " << fh_exp.getNrSpectra() << "\n";

  // ------------------------------------------------------------------
  // Mode 3: streaming consumer
  // ------------------------------------------------------------------
  std::cout << "\n=== Mode 3: ImzMLFile::load -> IMSDataConsumer ===\n";
  struct CountConsumer : public Interfaces::IMSDataConsumer
  {
    Size count {0};
    void setExpectedSize(Size, Size) override {}
    void setExperimentalSettings(const ExperimentalSettings&) override {}
    void consumeChromatogram(MSChromatogram&) override {}
    void consumeSpectrum(MSSpectrum&) override { ++count; }
  } consumer;
  loader.load(path, consumer);
  std::cout << "spectra streamed: " << consumer.count << "\n";

  // ------------------------------------------------------------------
  // Mode 4: index only (feeds OnDiscImzMLExperiment)
  // ------------------------------------------------------------------
  std::cout << "\n=== Mode 4: loadSpectraIndex ===\n";
  ImzMLMeta meta;
  std::vector<ImzMLSpectrumIndex> index;
  loader.loadSpectraIndex(path, meta, index);
  std::cout << "index entries: " << index.size() << ", mode: " << meta.imaging_mode << "\n";

  // ------------------------------------------------------------------
  // Mode 5: on-disc random access (1-based imzML coordinates)
  // ------------------------------------------------------------------
  std::cout << "\n=== Mode 5: OnDiscImzMLExperiment ===\n";
  OnDiscImzMLExperiment od;
  od.open(path);
  std::cout << "grid: " << od.gridWidth() << " x " << od.gridHeight()
            << ", spectra: " << od.getNrSpectra() << "\n";
  if (od.getNrSpectra() > 0)
  {
    const auto& e0 = od.getIndex(0);
    MSSpectrum disc = od.getSpectrumAtCoord(e0.x, e0.y, e0.z);
    std::cout << "pixel (" << e0.x << "," << e0.y << "," << e0.z << ") peaks:\n";
    printPeaks(disc);
  }

  // ------------------------------------------------------------------
  // Mode 6: in-memory pixel lookup (zero-based geometry coordinates)
  // ------------------------------------------------------------------
  std::cout << "\n=== Mode 6: MSImagingExperiment (RAM lookup) ===\n";
  MSImagingExperiment imaging;
  loader.load(path, imaging);
  std::cout << "geometry: " << imaging.getGeometry().getWidth()
            << " x " << imaging.getGeometry().getHeight()
            << ", pixels: " << imaging.getNumberOfPixels() << "\n";
  if (imaging.hasPixel(0, 0))
  {
    std::cout << "geometry pixel (0,0) peaks:\n";
    printPeaks(imaging.getSpectrum(0, 0));
  }

  // ------------------------------------------------------------------
  // Mode 7: export (store) and reload
  // ------------------------------------------------------------------
  std::cout << "\n=== Mode 7: ImzMLFile::store round-trip ===\n";
  std::string tmp_imzml = File::getTemporaryFile("Tutorial_ImzMLFile_export.imzML");
  loader.store(tmp_imzml, exp);
  MSExperiment reloaded;
  loader.load(tmp_imzml, reloaded);
  std::cout << "reloaded spectra: " << reloaded.getNrSpectra() << "\n";
  if (reloaded.metaValueExists("imzml:ibd_md5"))
  {
    std::cout << "ibd MD5: " << reloaded.getMetaValue("imzml:ibd_md5") << "\n";
  }

  return 0;
}
