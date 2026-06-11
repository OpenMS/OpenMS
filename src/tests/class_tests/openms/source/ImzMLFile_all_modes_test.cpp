// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Aditya Sarna $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/ImzMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/ImzMLHandlerHelper.h>
#include <OpenMS/IMAGING/MSImagingExperiment.h>
#include <OpenMS/IMAGING/MSImagingGeometry.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/OnDiscImzMLExperiment.h>

using namespace OpenMS;
using namespace std;

namespace
{
  const Size k_continuous_spectra = 9;
  const UInt k_grid = 3;

  Size findSpectrumIndexByImzMLCoord_(const MSExperiment& exp, Int x, Int y)
  {
    for (Size i = 0; i < exp.getNrSpectra(); ++i)
    {
      const MSSpectrum& s = exp[i];
      if (s.metaValueExists("imzml:x") && s.metaValueExists("imzml:y")
          && static_cast<Int>(s.getMetaValue("imzml:x")) == x
          && static_cast<Int>(s.getMetaValue("imzml:y")) == y)
      {
        return i;
      }
    }
    throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                     "(" + OpenMS::StringConversions::toString(x) + "," + OpenMS::StringConversions::toString(y) + ")");
  }
} // namespace

START_TEST(ImzMLFile_all_modes, "$Id$")

const std::string continuous_path = std::string(OPENMS_GET_TEST_DATA_PATH("ImzMLFile_1_Example_Continuous.imzML"));
const std::string processed_path = std::string(OPENMS_GET_TEST_DATA_PATH("ImzMLFile_2_Example_Processed.imzML"));

/////////////////////////////////////////////////////////////
START_SECTION(mode 1: full load into MSExperiment)
{
  MSExperiment exp;
  ImzMLFile().load(continuous_path, exp);

  TEST_EQUAL(exp.getNrSpectra(), k_continuous_spectra)
  TEST_EQUAL(exp.metaValueExists("imzml:imaging_mode"), true)
  TEST_EQUAL(OpenMS::StringConversions::toString(exp.getMetaValue("imzml:imaging_mode")), std::string("continuous"))
  TEST_EQUAL(static_cast<UInt>(exp.getMetaValue("imzml:max_count_x")), k_grid)
  TEST_EQUAL(static_cast<UInt>(exp.getMetaValue("imzml:max_count_y")), k_grid)
  if (!exp.empty())
  {
    TEST_EQUAL(exp[0].size() > 0, true)
    TEST_EQUAL(static_cast<Int>(exp[0].getMetaValue("imzml:x")), 1)
    TEST_EQUAL(static_cast<Int>(exp[0].getMetaValue("imzml:y")), 1)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
START_SECTION(mode 2: FileHandler and FileTypes)
{
  TEST_EQUAL(FileHandler::getTypeByFileName(continuous_path), FileTypes::IMZML)
  TEST_EQUAL(FileHandler::getTypeByContent(continuous_path), FileTypes::IMZML)

  MSExperiment exp;
  FileHandler().loadExperiment(continuous_path, exp);
  TEST_EQUAL(exp.getNrSpectra(), k_continuous_spectra)
  if (!exp.empty())
  {
    TEST_EQUAL(exp[0].size() > 0, true)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
START_SECTION(mode 3: streaming IMSDataConsumer)
{
  struct CollectConsumer : public Interfaces::IMSDataConsumer
  {
    Size count {0};
    Size first_size {0};
    void setExpectedSize(Size, Size) override {}
    void setExperimentalSettings(const ExperimentalSettings&) override {}
    void consumeChromatogram(MSChromatogram&) override {}
    void consumeSpectrum(MSSpectrum& s) override
    {
      if (count == 0)
      {
        first_size = s.size();
      }
      ++count;
    }
  } consumer;

  ImzMLFile().load(continuous_path, consumer);
  TEST_EQUAL(consumer.count, k_continuous_spectra)
  TEST_EQUAL(consumer.first_size > 0, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
START_SECTION(mode 4: loadSpectraIndex for on-disc access)
{
  ImzMLMeta meta;
  std::vector<ImzMLSpectrumIndex> index;
  ImzMLFile().loadSpectraIndex(continuous_path, meta, index);

  TEST_EQUAL(index.size(), k_continuous_spectra)
  TEST_EQUAL(meta.imaging_mode, "continuous")
  TEST_EQUAL(meta.max_count_x, k_grid)
  TEST_EQUAL(meta.max_count_y, k_grid)
  if (!index.empty())
  {
    TEST_EQUAL(index[0].x, 1)
    TEST_EQUAL(index[0].y, 1)
    TEST_EQUAL(index[0].mz_length > 0, true)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
START_SECTION(mode 5: OnDiscImzMLExperiment random access)
{
  OnDiscImzMLExperiment od;
  od.open(continuous_path);

  TEST_EQUAL(od.isOpen(), true)
  TEST_EQUAL(od.getNrSpectra(), k_continuous_spectra)
  TEST_EQUAL(od.gridWidth(), k_grid)
  TEST_EQUAL(od.gridHeight(), k_grid)
  TEST_EQUAL(od.getImzMLMeta().imaging_mode, "continuous")

  if (od.getNrSpectra() > 0)
  {
    MSSpectrum by_index = od.getSpectrum(0);
    TEST_EQUAL(by_index.size() > 0, true)

    const ImzMLSpectrumIndex& e0 = od.getIndex(0);
    MSSpectrum by_coord = od.getSpectrumAtCoord(e0.x, e0.y, e0.z);
    TEST_EQUAL(by_coord.size(), by_index.size())
    TEST_EQUAL(by_coord.size() > 0, true)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
START_SECTION(mode 6: MSImagingExperiment in-memory pixel lookup)
{
  MSImagingExperiment imaging;
  ImzMLFile().load(continuous_path, imaging);

  TEST_EQUAL(imaging.getNumberOfSpectra(), k_continuous_spectra)
  TEST_EQUAL(imaging.getNumberOfPixels(), k_continuous_spectra)
  TEST_EQUAL(imaging.getGeometry().getWidth(), k_grid)
  TEST_EQUAL(imaging.getGeometry().getHeight(), k_grid)
  TEST_EQUAL(imaging.hasPixel(0, 0), true)
  TEST_EQUAL(imaging.hasPixel(2, 2), true)
  if (imaging.hasPixel(0, 0))
  {
    TEST_EQUAL(imaging.getSpectrum(0, 0).size() > 0, true)
  }
  if (imaging.hasPixel(2, 2))
  {
    TEST_EQUAL(imaging.getSpectrum(2, 2).size() > 0, true)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
START_SECTION(mode 7: buildImagingGeometry from MSExperiment)
{
  MSExperiment exp;
  ImzMLFile().load(continuous_path, exp);

  MSImagingGeometry geom;
  ImzMLFile::buildImagingGeometry(exp, geom);
  TEST_EQUAL(geom.getNumberOfPixels(), k_continuous_spectra)
  if (geom.getNumberOfPixels() > 0)
  {
    TEST_EQUAL(geom.getSpectrumIndex(0, 0), 0)
    TEST_EQUAL(geom.getSpectrumIndex(2, 2), findSpectrumIndexByImzMLCoord_(exp, 3, 3))
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
START_SECTION(cross-mode consistency: full load vs on-disc vs RAM lookup)
{
  MSExperiment exp;
  ImzMLFile().load(continuous_path, exp);

  MSImagingExperiment imaging;
  ImzMLFile().load(continuous_path, imaging);

  OnDiscImzMLExperiment od;
  od.open(continuous_path);

  const Int px = 2;
  const Int py = 3;
  const Size idx = findSpectrumIndexByImzMLCoord_(exp, px, py);
  const MSSpectrum& ref = exp[idx];
  Int pz = 1;
  if (ref.metaValueExists("imzml:z"))
  {
    pz = static_cast<Int>(ref.getMetaValue("imzml:z"));
  }
  const Size full_size = ref.size();
  const Size ram_size = imaging.getSpectrum(static_cast<UInt>(px - 1), static_cast<UInt>(py - 1)).size();
  const Size disc_size = od.getSpectrumAtCoord(static_cast<uint32_t>(px), static_cast<uint32_t>(py),
                                                static_cast<uint32_t>(pz)).size();

  TEST_EQUAL(full_size > 0, true)
  TEST_EQUAL(ram_size, full_size)
  TEST_EQUAL(disc_size, full_size)
}
END_SECTION

/////////////////////////////////////////////////////////////
START_SECTION(processed imzML encoding)
{
  MSExperiment exp;
  ImzMLFile().load(processed_path, exp);

  TEST_EQUAL(exp.getNrSpectra() > 0, true)
  TEST_EQUAL(OpenMS::StringConversions::toString(exp.getMetaValue("imzml:imaging_mode")), std::string("processed"))
  if (!exp.empty())
  {
    TEST_EQUAL(exp[0].size() > 0, true)
  }

  OnDiscImzMLExperiment od;
  od.open(processed_path);
  TEST_EQUAL(od.getImzMLMeta().imaging_mode, "processed")
  if (od.getNrSpectra() > 0)
  {
    TEST_EQUAL(od.getSpectrum(0).size() > 0, true)
  }
}
END_SECTION

END_TEST
