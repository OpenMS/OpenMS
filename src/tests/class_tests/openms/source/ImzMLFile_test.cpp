// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Aditya Sarna, Patrick Boschmann $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/ImzMLFile.h>
#include <OpenMS/KERNEL/OnDiscImzMLExperiment.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/IMAGING/MSImagingExperiment.h>
#include <OpenMS/IMAGING/MSImagingGeometry.h>
#include <OpenMS/IMAGING/IonImage.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/METADATA/Instrument.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>

#include <fstream>
#include <cstdio>
#include <sstream>

using namespace OpenMS;
using namespace std;

namespace
{
  MSExperiment loadImzMLExperiment_(const std::string& path,
                                    const PeakFileOptions& opts = PeakFileOptions())
  {
    ImzMLFile f;
    f.setOptions(opts);
    MSImagingExperiment img;
    f.load(path, img);
    return img.getMSExperiment();
  }

  MSSpectrum makePixelSpectrum_(Int x, Int y, double mz, float intensity)
  {
    MSSpectrum spec;
    spec.push_back(Peak1D(mz, intensity));
    spec.setMetaValue("imzml:x", x);
    spec.setMetaValue("imzml:y", y);
    spec.setMetaValue("imzml:z", 1);
    return spec;
  }

  void writeMissingPixelCoordImzML_(const std::string& imzml_path)
  {
    std::string ibd_path = imzml_path;
    std::string lower = imzml_path;
    StringUtils::toLower(lower);
    if (StringUtils::hasSuffix(lower, ".imzml"))
    {
      ibd_path = imzml_path.substr(0, imzml_path.size() - 6) + ".ibd";
    }
    else
    {
      ibd_path = imzml_path + ".ibd";
    }

    std::ofstream ibd(ibd_path.c_str(), std::ios::binary);
    if (!ibd.good())
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, ibd_path);
    }

    std::ofstream xml(imzml_path.c_str());
    if (!xml.good())
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, imzml_path);
    }
    xml << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
        << "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" version=\"1.1.0\">\n"
        << "  <cvList count=\"2\">\n"
        << "    <cv id=\"MS\" fullName=\"MS\" version=\"1.0\" URI=\"\"/>\n"
        << "    <cv id=\"IMS\" fullName=\"Imaging MS Ontology\" version=\"1.1.0\" URI=\"\"/>\n"
        << "  </cvList>\n"
        << "  <fileDescription><fileContent>\n"
        << "    <cvParam accession=\"IMS:1000030\" name=\"continuous\" value=\"\"/>\n"
        << "  </fileContent></fileDescription>\n"
        << "  <referenceableParamGroupList count=\"1\">\n"
        << "    <referenceableParamGroup id=\"mzArray\">\n"
        << "      <cvParam accession=\"MS:1000514\" name=\"m/z array\" value=\"\"/>\n"
        << "    </referenceableParamGroup>\n"
        << "  </referenceableParamGroupList>\n"
        << "  <softwareList count=\"1\"><software id=\"sw1\" version=\"1\">\n"
        << "    <cvParam accession=\"MS:1000799\" name=\"custom unreleased software tool\" value=\"test\"/>\n"
        << "  </software></softwareList>\n"
        << "  <instrumentConfigurationList count=\"1\">\n"
        << "    <instrumentConfiguration id=\"IC1\">\n"
        << "      <cvParam accession=\"MS:1000031\" name=\"instrument model\" value=\"test\"/>\n"
        << "    </instrumentConfiguration>\n"
        << "  </instrumentConfigurationList>\n"
        << "  <dataProcessingList count=\"1\">\n"
        << "    <dataProcessing id=\"dp1\"><processingMethod order=\"1\" softwareRef=\"sw1\">\n"
        << "      <cvParam accession=\"MS:1000544\" name=\"Conversion to imzML\" value=\"\"/>\n"
        << "    </processingMethod></dataProcessing>\n"
        << "  </dataProcessingList>\n"
        << "  <run defaultInstrumentConfigurationRef=\"IC1\">\n"
        << "    <spectrumList count=\"1\" defaultDataProcessingRef=\"dp1\">\n"
        << "      <spectrum index=\"0\" id=\"spectrum=1\" defaultArrayLength=\"0\">\n"
        << "        <scanList count=\"1\">\n"
        << "          <scan instrumentConfigurationRef=\"IC1\"/>\n"
        << "        </scanList>\n"
        << "        <binaryDataArrayList count=\"0\"/>\n"
        << "      </spectrum>\n"
        << "    </spectrumList>\n"
        << "  </run>\n"
        << "</mzML>\n";
  }

  void writeDuplicatePixelImzML_(const std::string& imzml_path)
  {
    std::string ibd_path = imzml_path;
    std::string lower = imzml_path;
    StringUtils::toLower(lower);
    if (StringUtils::hasSuffix(lower, ".imzml"))
    {
      ibd_path = imzml_path.substr(0, imzml_path.size() - 6) + ".ibd";
    }
    else
    {
      ibd_path = imzml_path + ".ibd";
    }

    std::ofstream ibd(ibd_path.c_str(), std::ios::binary);
    if (!ibd.good())
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, ibd_path);
    }

    auto writeSpectrum = [](std::ofstream& xml, Size index, uint32_t x, uint32_t y)
    {
      xml << "      <spectrum index=\"" << index << "\" id=\"spectrum=" << (index + 1)
          << "\" defaultArrayLength=\"0\">\n"
          << "        <scanList count=\"1\">\n"
          << "          <scan instrumentConfigurationRef=\"IC1\">\n"
          << "            <cvParam accession=\"IMS:1000050\" name=\"position x\" value=\"" << x << "\"/>\n"
          << "            <cvParam accession=\"IMS:1000051\" name=\"position y\" value=\"" << y << "\"/>\n"
          << "            <cvParam accession=\"IMS:1000052\" name=\"position z\" value=\"1\"/>\n"
          << "          </scan>\n"
          << "        </scanList>\n"
          << "        <binaryDataArrayList count=\"0\"/>\n"
          << "      </spectrum>\n";
    };

    std::ofstream xml(imzml_path.c_str());
    if (!xml.good())
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, imzml_path);
    }
    xml << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
        << "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" version=\"1.1.0\">\n"
        << "  <cvList count=\"2\">\n"
        << "    <cv id=\"MS\" fullName=\"MS\" version=\"1.0\" URI=\"\"/>\n"
        << "    <cv id=\"IMS\" fullName=\"Imaging MS Ontology\" version=\"1.1.0\" URI=\"\"/>\n"
        << "  </cvList>\n"
        << "  <fileDescription><fileContent>\n"
        << "    <cvParam accession=\"IMS:1000030\" name=\"continuous\" value=\"\"/>\n"
        << "  </fileContent></fileDescription>\n"
        << "  <referenceableParamGroupList count=\"1\">\n"
        << "    <referenceableParamGroup id=\"mzArray\">\n"
        << "      <cvParam accession=\"MS:1000514\" name=\"m/z array\" value=\"\"/>\n"
        << "    </referenceableParamGroup>\n"
        << "  </referenceableParamGroupList>\n"
        << "  <softwareList count=\"1\"><software id=\"sw1\" version=\"1\">\n"
        << "    <cvParam accession=\"MS:1000799\" name=\"custom unreleased software tool\" value=\"test\"/>\n"
        << "  </software></softwareList>\n"
        << "  <instrumentConfigurationList count=\"1\">\n"
        << "    <instrumentConfiguration id=\"IC1\">\n"
        << "      <cvParam accession=\"MS:1000031\" name=\"instrument model\" value=\"test\"/>\n"
        << "    </instrumentConfiguration>\n"
        << "  </instrumentConfigurationList>\n"
        << "  <dataProcessingList count=\"1\">\n"
        << "    <dataProcessing id=\"dp1\"><processingMethod order=\"1\" softwareRef=\"sw1\">\n"
        << "      <cvParam accession=\"MS:1000544\" name=\"Conversion to imzML\" value=\"\"/>\n"
        << "    </processingMethod></dataProcessing>\n"
        << "  </dataProcessingList>\n"
        << "  <run defaultInstrumentConfigurationRef=\"IC1\">\n"
        << "    <spectrumList count=\"2\" defaultDataProcessingRef=\"dp1\">\n";
    writeSpectrum(xml, 0, 1, 1);
    writeSpectrum(xml, 1, 1, 1);
    xml << "    </spectrumList>\n"
        << "  </run>\n"
        << "</mzML>\n";
  }
} // namespace

START_TEST(ImzMLFile, "$Id$")

const std::string imzml_path = std::string(OPENMS_GET_TEST_DATA_PATH("ImzMLFile_1_Example_Continuous.imzML"));
const std::string imzml_processed_path = std::string(OPENMS_GET_TEST_DATA_PATH("ImzMLFile_2_Example_Processed.imzML"));


START_SECTION(void load(const std::string& filename, MSExperiment& exp))
{
  MSExperiment exp = loadImzMLExperiment_(imzml_path);

  TEST_EQUAL(exp.getNrSpectra() > 0, true)
  if (!exp.empty())
  {
    TEST_EQUAL(exp[0].size() > 0, true)

    if (exp[0].metaValueExists("imzml:x"))
    {
      TEST_EQUAL(exp[0].getMetaValue("imzml:x"), 1)
      TEST_EQUAL(exp[0].getMetaValue("imzml:y"), 1)
    }
  }
}
END_SECTION


START_SECTION(OnDiscImzMLExperiment random access)
{
  OnDiscImzMLExperiment od;
  od.open(imzml_path);

  TEST_EQUAL(od.getNrSpectra() > 0, true)

  if (od.getNrSpectra() > 0)
  {
    MSSpectrum s = od.getSpectrum(0);
    TEST_EQUAL(s.size() > 0, true)

    const auto& e0 = od.getIndex(0);
    MSSpectrum sc = od.getSpectrumAtCoord(e0.x, e0.y, e0.z);
    TEST_EQUAL(sc.size() > 0, true)
  }
}
END_SECTION


START_SECTION(const MSImagingGeometry& getGeometry() const)
{
  OnDiscImzMLExperiment od;
  od.open(imzml_path);

  const MSImagingGeometry& geom = od.getGeometry();
  // The on-disc geometry must match the dataset grid and (for an all-z==1 fixture)
  // map every spectrum to a pixel.
  TEST_EQUAL(geom.getWidth(), od.gridWidth())
  TEST_EQUAL(geom.getHeight(), od.gridHeight())
  TEST_EQUAL(geom.getNumberOfPixels(), od.getNrSpectra())

  // The 0-based geometry lookup must agree with the 1-based getSpectrumAtCoord().
  const auto& e0 = od.getIndex(0);
  TEST_EQUAL(geom.hasPixel(e0.x - 1, e0.y - 1), true)
  TEST_EQUAL(geom.getSpectrumIndex(e0.x - 1, e0.y - 1), 0)
}
END_SECTION


START_SECTION(IonImage extractIonImage(double mz, double tolerance_ppm) const)
{
  // The on-disc extraction must match the in-memory MSImagingExperiment path bit-for-bit:
  // both walk the same shared geometry and sum the same decoded peaks, the on-disc path just
  // reads each pixel's spectrum from the .ibd lazily instead of holding the dataset in memory.
  MSImagingExperiment img;
  ImzMLFile().load(imzml_path, img);

  OnDiscImzMLExperiment od;
  od.open(imzml_path);

  // Pick a real m/z present in the data (first peak of pixel 0) and a wide tolerance window.
  const MSSpectrum first = od.getSpectrum(0);
  TEST_EQUAL(first.empty(), false)
  const double mz = first[0].getMZ();
  const double tol_ppm = 1000.0;

  IonImage mem = img.extractIonImage(mz, tol_ppm);
  IonImage disc = od.extractIonImage(mz, tol_ppm);

  TEST_EQUAL(disc.getWidth(), mem.getWidth())
  TEST_EQUAL(disc.getHeight(), mem.getHeight())

  bool masks_match = true;
  for (UInt y = 0; y < mem.getHeight(); ++y)
  {
    for (UInt x = 0; x < mem.getWidth(); ++x)
    {
      if (mem.hasPixel(x, y) != disc.hasPixel(x, y))
      {
        masks_match = false;
      }
      if (mem.hasPixel(x, y) && disc.hasPixel(x, y))
      {
        TEST_REAL_SIMILAR(disc.getIntensity(x, y), mem.getIntensity(x, y))
      }
    }
  }
  TEST_EQUAL(masks_match, true)

  // Invalid arguments must throw.
  TEST_EXCEPTION(Exception::InvalidValue, od.extractIonImage(-1.0, tol_ppm))
  TEST_EXCEPTION(Exception::InvalidValue, od.extractIonImage(mz, -1.0))
}
END_SECTION


START_SECTION([EXTRA] OnDiscImzMLExperiment::open honours an explicit .ibd path override)
{
  // Copy the processed .imzML to a temp path whose *inferred* .ibd does NOT exist, then
  // open() pointing at the real .ibd via the override. This fails unless the override is
  // threaded through to the index load + UUID check (regression guard for the case where
  // the inferred sibling is missing/stale and the override is the file actually read).
  std::string tmp_imzml;
  NEW_TMP_FILE_EXT(tmp_imzml, ".imzML")
  {
    std::ifstream src(imzml_processed_path.c_str(), std::ios::binary);
    std::ofstream dst(tmp_imzml.c_str(), std::ios::binary);
    dst << src.rdbuf();
  }
  const std::string real_ibd = imzml_processed_path.substr(0, imzml_processed_path.size() - 6) + ".ibd";

  OnDiscImzMLExperiment od;
  od.open(tmp_imzml, real_ibd);
  TEST_EQUAL(od.getNrSpectra() > 0, true)
  if (od.getNrSpectra() > 0)
  {
    TEST_EQUAL(od.getSpectrum(0).size() > 0, true)
  }
}
END_SECTION


START_SECTION(void load(const std::string& filename, Interfaces::IMSDataConsumer& consumer))
{
  struct CountConsumer : public Interfaces::IMSDataConsumer
  {
    Size count {0};
    void setExpectedSize(Size, Size) override {}
    void setExperimentalSettings(const ExperimentalSettings&) override {}
    void consumeChromatogram(MSChromatogram&) override {}
    void consumeSpectrum(MSSpectrum&) override { ++count; }
  } consumer;

  ImzMLFile f;
  f.load(imzml_path, consumer);
  TEST_EQUAL(consumer.count > 0, true)
}
END_SECTION


START_SECTION(void load Example_Processed imzML)
{
  MSExperiment exp = loadImzMLExperiment_(imzml_processed_path);
  TEST_EQUAL(exp.getNrSpectra() > 0, true)
  if (!exp.empty())
  {
    TEST_EQUAL(exp[0].size() > 0, true)
  }
}
END_SECTION


START_SECTION(dataset metadata mirrored on MSExperiment after load)
{
  auto checkGridMeta_ = [](const MSExperiment& exp)
  {
    TEST_EQUAL(static_cast<UInt>(exp.getMetaValue("imzml:max_count_x")), 3)
    TEST_EQUAL(static_cast<UInt>(exp.getMetaValue("imzml:max_count_y")), 3)
    TEST_EQUAL(static_cast<UInt>(exp.getMetaValue("imzml:max_count_z")), 1)
    TEST_REAL_SIMILAR(static_cast<double>(exp.getMetaValue("imzml:pixel_size_x")), 100.0)
    TEST_REAL_SIMILAR(static_cast<double>(exp.getMetaValue("imzml:pixel_size_y")), 100.0)
    TEST_REAL_SIMILAR(static_cast<double>(exp.getMetaValue("imzml:max_dim_x")), 300.0)
    TEST_REAL_SIMILAR(static_cast<double>(exp.getMetaValue("imzml:max_dim_y")), 300.0)
    TEST_EQUAL(exp.metaValueExists("imzml:ibd_path"), true)
  };

  auto checkMatchesIndexMeta_ = [](const MSExperiment& exp, const ImzMLMeta& meta)
  {
    TEST_EQUAL(OpenMS::StringConversions::toString(exp.getMetaValue("imzml:imaging_mode")), meta.imaging_mode)
    TEST_EQUAL(static_cast<UInt>(exp.getMetaValue("imzml:max_count_x")), meta.max_count_x)
    TEST_EQUAL(static_cast<UInt>(exp.getMetaValue("imzml:max_count_y")), meta.max_count_y)
    TEST_EQUAL(static_cast<UInt>(exp.getMetaValue("imzml:max_count_z")), meta.max_count_z)
    TEST_EQUAL(OpenMS::StringConversions::toString(exp.getMetaValue("imzml:uuid")), meta.uuid)
    if (!meta.ibd_md5.empty())
    {
      TEST_EQUAL(OpenMS::StringConversions::toString(exp.getMetaValue("imzml:ibd_md5")), meta.ibd_md5)
    }
    if (!meta.ibd_sha1.empty())
    {
      TEST_EQUAL(OpenMS::StringConversions::toString(exp.getMetaValue("imzml:ibd_sha1")), meta.ibd_sha1)
    }
    if (!meta.scan_pattern.empty())
    {
      TEST_EQUAL(OpenMS::StringConversions::toString(exp.getMetaValue("imzml:scan_pattern")), meta.scan_pattern)
    }
    if (!meta.scan_direction.empty())
    {
      TEST_EQUAL(OpenMS::StringConversions::toString(exp.getMetaValue("imzml:scan_direction")), meta.scan_direction)
    }
    if (!meta.line_scan_direction.empty())
    {
      TEST_EQUAL(OpenMS::StringConversions::toString(exp.getMetaValue("imzml:line_scan_direction")), meta.line_scan_direction)
    }
    if (!meta.polarity.empty())
    {
      TEST_EQUAL(OpenMS::StringConversions::toString(exp.getMetaValue("imzml:polarity")), meta.polarity)
    }
  };

  // Continuous reference file: shared m/z axis, MD5 checksum, float32 arrays.
  {
    MSExperiment exp = loadImzMLExperiment_(imzml_path);

    TEST_EQUAL(OpenMS::StringConversions::toString(exp.getMetaValue("imzml:imaging_mode")), std::string("continuous"))
    checkGridMeta_(exp);
    TEST_EQUAL(OpenMS::StringConversions::toString(exp.getMetaValue("imzml:uuid")), std::string("12345678-1234-1234-1234-123456789012"))
    TEST_EQUAL(OpenMS::StringConversions::toString(exp.getMetaValue("imzml:ibd_md5")), std::string("4b5dd9fa84fafc955cfdd301f9ed55d7"))
    TEST_EQUAL(exp.metaValueExists("imzml:ibd_sha1"), false)
    TEST_EQUAL(OpenMS::StringConversions::toString(exp.getMetaValue("imzml:mz_data_type")), std::string("float32"))
    TEST_EQUAL(OpenMS::StringConversions::toString(exp.getMetaValue("imzml:int_data_type")), std::string("float32"))

    ImzMLMeta meta;
    std::vector<ImzMLSpectrumIndex> index;
    ImzMLFile().loadSpectraIndex(imzml_path, meta, index);
    checkMatchesIndexMeta_(exp, meta);
    TEST_EQUAL(index.size(), exp.getNrSpectra())

    TEST_EQUAL(exp[0].metaValueExists("imzml:x"), true)
    TEST_EQUAL(exp[0].getMetaValue("imzml:x"), 1)
    TEST_EQUAL(exp[0].getMetaValue("imzml:y"), 1)
    TEST_EQUAL(exp[0].getMetaValue("imzml:z"), 1)
    TEST_REAL_SIMILAR(exp[0][0].getMZ(), 100.0)
  }

  // Processed reference file: per-spectrum m/z, SHA-1, acquisition geometry terms.
  {
    MSExperiment exp = loadImzMLExperiment_(imzml_processed_path);

    TEST_EQUAL(OpenMS::StringConversions::toString(exp.getMetaValue("imzml:imaging_mode")), std::string("processed"))
    checkGridMeta_(exp);
    TEST_EQUAL(OpenMS::StringConversions::toString(exp.getMetaValue("imzml:uuid")), std::string("9d501bdc53444916b7e97e795b02c856"))
    TEST_EQUAL(OpenMS::StringConversions::toString(exp.getMetaValue("imzml:ibd_sha1")), std::string("7e8fdb93053915d3edb51b70aa0619ac209964df"))
    TEST_EQUAL(exp.metaValueExists("imzml:ibd_md5"), false)
    TEST_EQUAL(OpenMS::StringConversions::toString(exp.getMetaValue("imzml:scan_pattern")), std::string("top down"))
    TEST_EQUAL(OpenMS::StringConversions::toString(exp.getMetaValue("imzml:scan_direction")), std::string("horizontal"))
    TEST_EQUAL(OpenMS::StringConversions::toString(exp.getMetaValue("imzml:line_scan_direction")), std::string("left-right"))
    TEST_EQUAL(OpenMS::StringConversions::toString(exp.getMetaValue("imzml:polarity")), std::string("negative"))
    TEST_EQUAL(OpenMS::StringConversions::toString(exp.getMetaValue("imzml:mz_data_type")), std::string("float32"))
    TEST_EQUAL(OpenMS::StringConversions::toString(exp.getMetaValue("imzml:int_data_type")), std::string("float32"))

    ImzMLMeta meta;
    std::vector<ImzMLSpectrumIndex> index;
    ImzMLFile().loadSpectraIndex(imzml_processed_path, meta, index);
    checkMatchesIndexMeta_(exp, meta);
    TEST_EQUAL(index.size(), exp.getNrSpectra())

    TEST_EQUAL(exp[0].getMetaValue("imzml:x"), 1)
    TEST_EQUAL(exp[0].getMetaValue("imzml:y"), 1)
    TEST_REAL_SIMILAR(exp[0][0].getMZ(), 100.083336)
  }
}
END_SECTION


START_SECTION(void load(const std::string& filename, MSImagingExperiment& exp))
{
  MSImagingExperiment imaging;
  ImzMLFile f;
  f.load(imzml_path, imaging);

  TEST_EQUAL(imaging.getNumberOfSpectra() > 0, true)
  TEST_EQUAL(imaging.getGeometry().getWidth(), 3)
  TEST_EQUAL(imaging.getGeometry().getHeight(), 3)
  TEST_EQUAL(imaging.hasPixel(0, 0), true)
  if (imaging.hasPixel(0, 0))
  {
    TEST_EQUAL(imaging.getSpectrum(0, 0).size() > 0, true)
  }
  TEST_EQUAL(imaging.hasPixel(2, 2), true)
  if (imaging.hasPixel(2, 2))
  {
    TEST_EQUAL(imaging.getSpectrum(2, 2).size() > 0, true)
  }
}
END_SECTION


START_SECTION(void buildImagingGeometry(const MSExperiment& exp, MSImagingGeometry& geom))
{
  MSExperiment exp = loadImzMLExperiment_(imzml_path);

  MSImagingGeometry geom;
  ImzMLFile::buildImagingGeometry(exp, geom);
  TEST_EQUAL(geom.getNumberOfPixels(), exp.getNrSpectra())
  if (geom.getNumberOfPixels() > 0)
  {
    TEST_EQUAL(geom.getSpectrumIndex(0, 0), 0)
  }
}
END_SECTION


START_SECTION(static void buildImagingGeometry(const std::vector<ImzMLSpectrumIndex>& index, const ImzMLMeta& meta, MSImagingGeometry& geom))
{
  // Source-of-truth builder: geometry is derived straight from a parsed index with NO
  // MSExperiment and NO imzml:* MetaValues anywhere (as the loaders now use it).
  std::vector<ImzMLSpectrumIndex> index(4);
  index[0].x = 1; index[0].y = 1; // 1-based imzML coords
  index[1].x = 2; index[1].y = 1;
  index[2].x = 1; index[2].y = 2;
  index[3].x = 2; index[3].y = 2;

  ImzMLMeta meta;
  meta.max_count_x = 2;
  meta.max_count_y = 2;
  meta.pixel_size_x = 25.0;
  meta.pixel_size_y = 25.0;

  MSImagingGeometry geom;
  ImzMLFile::buildImagingGeometry(index, meta, geom);

  TEST_EQUAL(geom.getNumberOfPixels(), 4)
  TEST_EQUAL(geom.getWidth(), 2)
  TEST_EQUAL(geom.getHeight(), 2)
  TEST_REAL_SIMILAR(geom.getPixelSizeX(), 25.0)
  TEST_EQUAL(geom.hasPixel(0, 0), true) // (1,1) -> (0,0)
  TEST_EQUAL(geom.getSpectrumIndex(0, 0), 0)
  TEST_EQUAL(geom.hasPixel(1, 1), true) // (2,2) -> (1,1)
  TEST_EQUAL(geom.getSpectrumIndex(1, 1), 3)
}
END_SECTION


START_SECTION(void store round-trip continuous imzML)
{
  MSExperiment original = loadImzMLExperiment_(imzml_path);

  ImzMLFile f;
  std::string tmp_imzml;
  NEW_TMP_FILE_EXT(tmp_imzml, ".imzML");
  f.store(tmp_imzml, original);

  MSExperiment reloaded = loadImzMLExperiment_(tmp_imzml);

  TEST_EQUAL(reloaded.getNrSpectra(), original.getNrSpectra())
  if (!original.empty() && !reloaded.empty())
  {
    TEST_EQUAL(reloaded[0].size(), original[0].size())
    TEST_REAL_SIMILAR(reloaded[0][0].getMZ(), original[0][0].getMZ())
    TEST_REAL_SIMILAR(reloaded[0][0].getIntensity(), original[0][0].getIntensity())
  }
  TEST_EQUAL(reloaded.metaValueExists("imzml:imaging_mode"), true)
  TEST_EQUAL(reloaded.getMetaValue("imzml:imaging_mode"), "continuous")
  TEST_EQUAL(reloaded.metaValueExists("imzml:ibd_md5"), true)
  TEST_NOT_EQUAL(OpenMS::StringConversions::toString(reloaded.getMetaValue("imzml:ibd_md5")).empty(), true)
}
END_SECTION


START_SECTION(void store(const std::string& filename, const MSImagingExperiment& exp))
{
  // Geometry-driven store: build an MSImagingExperiment whose spectra carry NO imzml:x/y
  // MetaValues (as the BrukerTimsImagingFile path produces) — only the geometry holds the
  // coordinates — and verify store derives them from the geometry and they round-trip.
  MSExperiment ms;
  for (int i = 0; i < 4; ++i)
  {
    MSSpectrum s;
    Peak1D p;
    p.setMZ(100.0 + i);
    p.setIntensity(10.0 * (i + 1));
    s.push_back(p);
    ms.addSpectrum(s);
  }
  TEST_EQUAL(ms[0].metaValueExists("imzml:x"), false) // no per-spectrum coordinates anywhere

  MSImagingGeometry geom;
  geom.setDimensions(2, 2);
  geom.setPixelSize(25.0, 25.0, "micrometer");
  geom.addPixel(0, 0, 0);
  geom.addPixel(1, 0, 1);
  geom.addPixel(0, 1, 2);
  geom.addPixel(1, 1, 3);

  MSImagingExperiment img;
  img.setMSExperiment(ms);
  img.setGeometry(geom);

  std::string tmp_imzml;
  NEW_TMP_FILE_EXT(tmp_imzml, ".imzML");
  ImzMLFile().store(tmp_imzml, img); // must NOT throw despite missing imzml:x/y MetaValues

  MSImagingExperiment reloaded;
  ImzMLFile().load(tmp_imzml, reloaded);
  TEST_EQUAL(reloaded.getMSExperiment().getNrSpectra(), 4)
  TEST_EQUAL(reloaded.getGeometry().getNumberOfPixels(), 4)
  TEST_EQUAL(reloaded.getGeometry().getWidth(), 2)
  TEST_EQUAL(reloaded.getGeometry().getHeight(), 2)
  TEST_EQUAL(reloaded.getGeometry().hasPixel(0, 0), true)
  TEST_EQUAL(reloaded.getGeometry().hasPixel(1, 1), true)
}
END_SECTION


START_SECTION(void store round-trip processed imzML)
{
  MSExperiment original = loadImzMLExperiment_(imzml_processed_path);

  ImzMLFile f;
  std::string tmp_imzml;
  NEW_TMP_FILE_EXT(tmp_imzml, ".imzML");
  f.store(tmp_imzml, original);

  MSExperiment reloaded = loadImzMLExperiment_(tmp_imzml);

  TEST_EQUAL(reloaded.getNrSpectra(), original.getNrSpectra())
  if (!original.empty() && !reloaded.empty())
  {
    TEST_EQUAL(reloaded[0].size(), original[0].size())
    TEST_REAL_SIMILAR(reloaded[0][0].getMZ(), original[0][0].getMZ())
  }
  TEST_EQUAL(reloaded.getMetaValue("imzml:imaging_mode"), "processed")
}
END_SECTION


START_SECTION(bool isValid(const std::string& filename, std::ostream& os))
{
  MSExperiment exp = loadImzMLExperiment_(imzml_path);

  ImzMLFile f;
  std::string tmp_imzml;
  NEW_TMP_FILE_EXT(tmp_imzml, ".imzML");
  f.store(tmp_imzml, exp);

  std::ostringstream os;
  TEST_EQUAL(f.isValid(tmp_imzml, os), true)
}
END_SECTION


START_SECTION(void store rejects missing pixel coordinates)
{
  MSExperiment exp;
  MSSpectrum spec;
  spec.push_back(Peak1D(100.0, 1000.0));
  exp.addSpectrum(spec);

  ImzMLFile f;
  std::string tmp_imzml;
  NEW_TMP_FILE_EXT(tmp_imzml, ".imzML");
  TEST_EXCEPTION(Exception::MissingInformation, f.store(tmp_imzml, exp))
}
END_SECTION


START_SECTION(void store tolerates duplicate pixel coordinates by default)
{
  // The reader accepts duplicate coordinates by default, so the writer must too:
  // a dataset that loads has to be storable again without an unswitchable error.
  MSExperiment exp;
  exp.addSpectrum(makePixelSpectrum_(1, 1, 100.0, 1000.0));
  exp.addSpectrum(makePixelSpectrum_(1, 1, 101.0, 900.0));

  ImzMLFile f;
  std::string tmp_imzml;
  NEW_TMP_FILE_EXT(tmp_imzml, ".imzML");
  f.store(tmp_imzml, exp);

  // Both spectra are written out; only the first claims the shared pixel on re-read.
  MSImagingExperiment img;
  f.load(tmp_imzml, img);
  TEST_EQUAL(img.getMSExperiment().getNrSpectra(), 2u)
  TEST_EQUAL(img.getGeometry().getNumberOfPixels(), 1u)
  TEST_EQUAL(img.getGeometry().getSpectrumIndex(0, 0), 0u)

  std::string ibd_path = tmp_imzml;
  std::string lower = tmp_imzml;
  StringUtils::toLower(lower);
  if (StringUtils::hasSuffix(lower, ".imzml"))
  {
    ibd_path = tmp_imzml.substr(0, tmp_imzml.size() - 6) + ".ibd";
  }
  remove(ibd_path.c_str());
}
END_SECTION


START_SECTION(void store rejects incompatible continuous mode)
{
  MSExperiment exp;
  exp.setMetaValue("imzml:imaging_mode", "continuous");
  exp.addSpectrum(makePixelSpectrum_(1, 1, 100.0, 1000.0));
  exp.addSpectrum(makePixelSpectrum_(2, 1, 200.0, 800.0));

  ImzMLFile f;
  std::string tmp_imzml;
  NEW_TMP_FILE_EXT(tmp_imzml, ".imzML");
  TEST_EXCEPTION(Exception::InvalidParameter, f.store(tmp_imzml, exp))
}
END_SECTION


START_SECTION(void buildImagingGeometry tolerates duplicate pixels by default)
{
  // Must agree with the loaders: an experiment that loaded fine cannot fail here.
  MSExperiment exp;
  exp.addSpectrum(makePixelSpectrum_(1, 1, 100.0, 1000.0));
  exp.addSpectrum(makePixelSpectrum_(1, 1, 101.0, 900.0));
  exp.addSpectrum(makePixelSpectrum_(2, 1, 102.0, 800.0));

  MSImagingGeometry geom;
  ImzMLFile::buildImagingGeometry(exp, geom);
  TEST_EQUAL(geom.getNumberOfPixels(), 2u)
  TEST_EQUAL(geom.getSpectrumIndex(0, 0), 0u) // first spectrum keeps the shared pixel
  TEST_EQUAL(geom.getSpectrumIndex(1, 0), 2u)
}
END_SECTION


START_SECTION(void load rejects spectrum missing pixel coordinates)
{
  std::string tmp_imzml;
  NEW_TMP_FILE_EXT(tmp_imzml, ".imzML");
  writeMissingPixelCoordImzML_(tmp_imzml);

  ImzMLFile f;
  MSImagingExperiment img;
  TEST_EXCEPTION(Exception::ParseError, f.load(tmp_imzml, img))

  std::string ibd_path = tmp_imzml;
  std::string lower = tmp_imzml;
  StringUtils::toLower(lower);
  if (StringUtils::hasSuffix(lower, ".imzml"))
  {
    ibd_path = tmp_imzml.substr(0, tmp_imzml.size() - 6) + ".ibd";
  }
  remove(ibd_path.c_str());
}
END_SECTION


START_SECTION(void load tolerates duplicate pixel coordinates by default)
{
  std::string tmp_imzml;
  NEW_TMP_FILE_EXT(tmp_imzml, ".imzML");
  writeDuplicatePixelImzML_(tmp_imzml);

  ImzMLFile f;
  MSImagingExperiment img;
  f.load(tmp_imzml, img);

  // Both spectra are loaded and reachable by index, but the shared pixel (1,1) maps to
  // the first of them only.
  TEST_EQUAL(img.getMSExperiment().getNrSpectra(), 2u)
  TEST_EQUAL(img.getGeometry().getNumberOfPixels(), 1u)
  TEST_EQUAL(img.getGeometry().getSpectrumIndex(0, 0), 0u)

  std::string ibd_path = tmp_imzml;
  std::string lower = tmp_imzml;
  StringUtils::toLower(lower);
  if (StringUtils::hasSuffix(lower, ".imzml"))
  {
    ibd_path = tmp_imzml.substr(0, tmp_imzml.size() - 6) + ".ibd";
  }
  remove(ibd_path.c_str());
}
END_SECTION


START_SECTION(OnDiscImzMLExperiment tolerates duplicate pixel coordinates by default)
{
  std::string tmp_imzml;
  NEW_TMP_FILE_EXT(tmp_imzml, ".imzML");
  writeDuplicatePixelImzML_(tmp_imzml);

  OnDiscImzMLExperiment od;
  od.open(tmp_imzml);

  TEST_EQUAL(od.size(), 2u)
  TEST_EQUAL(od.getGeometry().getNumberOfPixels(), 1u)
  TEST_EQUAL(od.getGeometry().getSpectrumIndex(0, 0), 0u)

  std::string ibd_path = tmp_imzml;
  std::string lower = tmp_imzml;
  StringUtils::toLower(lower);
  if (StringUtils::hasSuffix(lower, ".imzml"))
  {
    ibd_path = tmp_imzml.substr(0, tmp_imzml.size() - 6) + ".ibd";
  }
  remove(ibd_path.c_str());
}
END_SECTION


START_SECTION(void store applies PeakFileOptions m/z range filter)
{
  MSExperiment original = loadImzMLExperiment_(imzml_path);
  ImzMLFile f;
  TEST_EQUAL(original.getNrSpectra() > 0, true)
  if (original.getNrSpectra() > 0 && !original[0].empty())
  {
    const Size full_peaks = original[0].size();
    const double lo_mz = original[0][0].getMZ();
    const double hi_mz = original[0][full_peaks - 1].getMZ();
    const double mid_mz = (lo_mz + hi_mz) / 2.0;

    PeakFileOptions opts;
    opts.setMZRange(DRange<1>(lo_mz, mid_mz));
    f.setOptions(opts);

    std::string tmp_imzml;
    NEW_TMP_FILE_EXT(tmp_imzml, ".imzML");
    f.store(tmp_imzml, original);

    MSExperiment filtered = loadImzMLExperiment_(tmp_imzml);
    TEST_EQUAL(filtered.getNrSpectra(), original.getNrSpectra())
    TEST_EQUAL(filtered[0].size() < full_peaks, true)
    TEST_EQUAL(filtered[0][0].getMZ() >= lo_mz, true)
    TEST_EQUAL(filtered[0].back().getMZ() <= mid_mz, true)
  }
}
END_SECTION


START_SECTION(void store honors PeakFileOptions binary precision)
{
  MSExperiment original = loadImzMLExperiment_(imzml_path);
  ImzMLFile f;

  PeakFileOptions opts;
  opts.setMz32Bit(false);
  opts.setIntensity32Bit(false);
  f.setOptions(opts);

  std::string tmp_imzml;
  NEW_TMP_FILE_EXT(tmp_imzml, ".imzML");
  f.store(tmp_imzml, original);

  ImzMLMeta meta;
  std::vector<ImzMLSpectrumIndex> index;
  ImzMLFile().loadSpectraIndex(tmp_imzml, meta, index);
  TEST_EQUAL(index.empty(), false)
  TEST_EQUAL(index[0].mz_type, ImzMLSpectrumIndex::DataType::FLOAT64)
  TEST_EQUAL(index[0].int_type, ImzMLSpectrumIndex::DataType::FLOAT64)
}
END_SECTION


START_SECTION(void store metadata round-trip)
{
  MSExperiment original = loadImzMLExperiment_(imzml_path);
  ImzMLFile f;

  original.setMetaValue("imzml:scan_pattern", "top down");
  original.setMetaValue("imzml:scan_direction", "flyback");
  original.setMetaValue("imzml:line_scan_direction", "left-right");
  original.setMetaValue("imzml:polarity", "positive");
  Instrument inst = original.getInstrument();
  inst.setModel("Test MSI Instrument");
  original.setInstrument(inst);
  if (!original.empty())
  {
    original[0].setRT(12.34);
    original[0].setMSLevel(1);
  }

  std::string tmp_imzml;
  NEW_TMP_FILE_EXT(tmp_imzml, ".imzML");
  f.store(tmp_imzml, original);

  MSExperiment reloaded = loadImzMLExperiment_(tmp_imzml);
  TEST_EQUAL(OpenMS::StringConversions::toString(reloaded.getMetaValue("imzml:scan_pattern")), std::string("top down"))
  TEST_EQUAL(OpenMS::StringConversions::toString(reloaded.getMetaValue("imzml:scan_direction")), std::string("flyback"))
  TEST_EQUAL(OpenMS::StringConversions::toString(reloaded.getMetaValue("imzml:line_scan_direction")), std::string("left-right"))
  TEST_EQUAL(OpenMS::StringConversions::toString(reloaded.getMetaValue("imzml:polarity")), std::string("positive"))
  TEST_EQUAL(reloaded.metaValueExists("imzml:ibd_sha1"), true)
  TEST_NOT_EQUAL(OpenMS::StringConversions::toString(reloaded.getMetaValue("imzml:ibd_sha1")).empty(), true)
  TEST_REAL_SIMILAR(static_cast<double>(reloaded.getMetaValue("imzml:max_dim_x")), 300.0)
  TEST_REAL_SIMILAR(static_cast<double>(reloaded.getMetaValue("imzml:max_dim_y")), 300.0)
  if (!reloaded.empty())
  {
    TEST_REAL_SIMILAR(reloaded[0].getRT(), 12.34)
    TEST_EQUAL(reloaded[0].getMSLevel(), 1)
  }
}
END_SECTION


START_SECTION(IonImage extractIonImage(double mz, double tolerance_ppm, Size region_id) const)
{
  // On-disc region extraction must produce the same result as in-memory for the same region.
  // This covers the 3-arg overload of OnDiscImzMLExperiment which was previously untested.
  MSImagingExperiment img;
  ImzMLFile().load(imzml_path, img);

  OnDiscImzMLExperiment od;
  od.open(imzml_path);

  const MSSpectrum first = od.getSpectrum(0);
  TEST_EQUAL(first.empty(), false)
  const double mz      = first[0].getMZ();
  const double tol_ppm = 1000.0;

  // Add the same region (single-pixel rectangle at origin) to both geometries,
  // then compare the extracted images pixel-for-pixel.
  const MSImagingRegion region = MSImagingRegion::rectangle(1, "roi", 0, 0, 0, 0);
  img.getGeometry().addRegion(region);
  od.getGeometry().addRegion(region);

  IonImage mem  = img.extractIonImage(mz, tol_ppm, 1);
  IonImage disc = od.extractIonImage(mz, tol_ppm, 1);

  TEST_EQUAL(disc.getWidth(),  mem.getWidth())
  TEST_EQUAL(disc.getHeight(), mem.getHeight())

  bool masks_match = true;
  bool compared_region_pixel = false;
  for (UInt y = 0; y < mem.getHeight(); ++y)
  {
    for (UInt x = 0; x < mem.getWidth(); ++x)
    {
      if (mem.hasPixel(x, y) != disc.hasPixel(x, y)) { masks_match = false; }
      if (mem.hasPixel(x, y) && disc.hasPixel(x, y))
      {
        compared_region_pixel = true;
        TEST_REAL_SIMILAR(disc.getIntensity(x, y), mem.getIntensity(x, y))
      }
    }
  }
  TEST_EQUAL(masks_match, true)
  TEST_TRUE(compared_region_pixel)

  // Unknown region id must throw.
  TEST_EXCEPTION(Exception::ElementNotFound, od.extractIonImage(mz, tol_ppm, 99))
}
END_SECTION

END_TEST
