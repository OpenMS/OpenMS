// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Aditya Sarna $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FORMAT/ImzMLFile.h>
#include <OpenMS/KERNEL/OnDiscImzMLExperiment.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/IMAGING/MSImagingExperiment.h>
#include <OpenMS/IMAGING/MSImagingGeometry.h>
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
  MSSpectrum makePixelSpectrum_(Int x, Int y, double mz, float intensity)
  {
    MSSpectrum spec;
    spec.push_back(Peak1D(mz, intensity));
    spec.setMetaValue("imzml:x", x);
    spec.setMetaValue("imzml:y", y);
    spec.setMetaValue("imzml:z", 1);
    return spec;
  }

  void writeMissingPixelCoordImzML_(const String& imzml_path)
  {
    String ibd_path = imzml_path;
    String lower = imzml_path;
    lower.toLower();
    if (lower.hasSuffix(".imzml"))
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

  void writeDuplicatePixelImzML_(const String& imzml_path)
  {
    String ibd_path = imzml_path;
    String lower = imzml_path;
    lower.toLower();
    if (lower.hasSuffix(".imzml"))
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

const String imzml_path = String(OPENMS_GET_TEST_DATA_PATH("ImzMLFile_1_Example_Continuous.imzML"));
const String imzml_processed_path = String(OPENMS_GET_TEST_DATA_PATH("ImzMLFile_2_Example_Processed.imzML"));


START_SECTION(void load(const String& filename, MSExperiment& exp))
{
  MSExperiment exp;
  ImzMLFile f;
  f.setLogType(ProgressLogger::CMD);
  f.load(imzml_path, exp);

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


START_SECTION(void load(const String& filename, Interfaces::IMSDataConsumer& consumer))
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
  MSExperiment exp;
  ImzMLFile f;
  f.load(imzml_processed_path, exp);
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
    TEST_EQUAL(String(exp.getMetaValue("imzml:imaging_mode")), meta.imaging_mode)
    TEST_EQUAL(static_cast<UInt>(exp.getMetaValue("imzml:max_count_x")), meta.max_count_x)
    TEST_EQUAL(static_cast<UInt>(exp.getMetaValue("imzml:max_count_y")), meta.max_count_y)
    TEST_EQUAL(static_cast<UInt>(exp.getMetaValue("imzml:max_count_z")), meta.max_count_z)
    TEST_EQUAL(String(exp.getMetaValue("imzml:uuid")), meta.uuid)
    if (!meta.ibd_md5.empty())
    {
      TEST_EQUAL(String(exp.getMetaValue("imzml:ibd_md5")), meta.ibd_md5)
    }
    if (!meta.ibd_sha1.empty())
    {
      TEST_EQUAL(String(exp.getMetaValue("imzml:ibd_sha1")), meta.ibd_sha1)
    }
    if (!meta.scan_pattern.empty())
    {
      TEST_EQUAL(String(exp.getMetaValue("imzml:scan_pattern")), meta.scan_pattern)
    }
    if (!meta.scan_direction.empty())
    {
      TEST_EQUAL(String(exp.getMetaValue("imzml:scan_direction")), meta.scan_direction)
    }
    if (!meta.line_scan_direction.empty())
    {
      TEST_EQUAL(String(exp.getMetaValue("imzml:line_scan_direction")), meta.line_scan_direction)
    }
    if (!meta.polarity.empty())
    {
      TEST_EQUAL(String(exp.getMetaValue("imzml:polarity")), meta.polarity)
    }
  };

  // Continuous reference file: shared m/z axis, MD5 checksum, float32 arrays.
  {
    MSExperiment exp;
    ImzMLFile().load(imzml_path, exp);

    TEST_EQUAL(String(exp.getMetaValue("imzml:imaging_mode")), String("continuous"))
    checkGridMeta_(exp);
    TEST_EQUAL(String(exp.getMetaValue("imzml:uuid")), String("12345678-1234-1234-1234-123456789012"))
    TEST_EQUAL(String(exp.getMetaValue("imzml:ibd_md5")), String("d8654952b0bb69e3f1e70c4a5fd3e06d"))
    TEST_EQUAL(exp.metaValueExists("imzml:ibd_sha1"), false)
    TEST_EQUAL(String(exp.getMetaValue("imzml:mz_data_type")), String("float32"))
    TEST_EQUAL(String(exp.getMetaValue("imzml:int_data_type")), String("float32"))

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
    MSExperiment exp;
    ImzMLFile().load(imzml_processed_path, exp);

    TEST_EQUAL(String(exp.getMetaValue("imzml:imaging_mode")), String("processed"))
    checkGridMeta_(exp);
    TEST_EQUAL(String(exp.getMetaValue("imzml:uuid")), String("9d501bdc53444916b7e97e795b02c856"))
    TEST_EQUAL(String(exp.getMetaValue("imzml:ibd_sha1")), String("7e8fdb93053915d3edb51b70aa0619ac209964df"))
    TEST_EQUAL(exp.metaValueExists("imzml:ibd_md5"), false)
    TEST_EQUAL(String(exp.getMetaValue("imzml:scan_pattern")), String("top down"))
    TEST_EQUAL(String(exp.getMetaValue("imzml:scan_direction")), String("horizontal"))
    TEST_EQUAL(String(exp.getMetaValue("imzml:line_scan_direction")), String("left-right"))
    TEST_EQUAL(String(exp.getMetaValue("imzml:polarity")), String("negative"))
    TEST_EQUAL(String(exp.getMetaValue("imzml:mz_data_type")), String("float32"))
    TEST_EQUAL(String(exp.getMetaValue("imzml:int_data_type")), String("float32"))

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


START_SECTION(void load(const String& filename, MSImagingExperiment& exp))
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
  MSExperiment exp;
  ImzMLFile f;
  f.load(imzml_path, exp);

  MSImagingGeometry geom;
  ImzMLFile::buildImagingGeometry(exp, geom);
  TEST_EQUAL(geom.getNumberOfPixels(), exp.getNrSpectra())
  if (geom.getNumberOfPixels() > 0)
  {
    TEST_EQUAL(geom.getSpectrumIndex(0, 0), 0)
  }
}
END_SECTION


START_SECTION(void store round-trip continuous imzML)
{
  MSExperiment original;
  ImzMLFile f;
  f.load(imzml_path, original);

  String tmp_imzml;
  NEW_TMP_FILE_EXT(tmp_imzml, ".imzML");
  f.store(tmp_imzml, original);

  MSExperiment reloaded;
  f.load(tmp_imzml, reloaded);

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
  TEST_NOT_EQUAL(String(reloaded.getMetaValue("imzml:ibd_md5")).empty(), true)
}
END_SECTION


START_SECTION(void store round-trip processed imzML)
{
  MSExperiment original;
  ImzMLFile f;
  f.load(imzml_processed_path, original);

  String tmp_imzml;
  NEW_TMP_FILE_EXT(tmp_imzml, ".imzML");
  f.store(tmp_imzml, original);

  MSExperiment reloaded;
  f.load(tmp_imzml, reloaded);

  TEST_EQUAL(reloaded.getNrSpectra(), original.getNrSpectra())
  if (!original.empty() && !reloaded.empty())
  {
    TEST_EQUAL(reloaded[0].size(), original[0].size())
    TEST_REAL_SIMILAR(reloaded[0][0].getMZ(), original[0][0].getMZ())
  }
  TEST_EQUAL(reloaded.getMetaValue("imzml:imaging_mode"), "processed")
}
END_SECTION


START_SECTION(bool isValid(const String& filename, std::ostream& os))
{
  ImzMLFile f;
  MSExperiment exp;
  f.load(imzml_path, exp);

  String tmp_imzml;
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
  String tmp_imzml;
  NEW_TMP_FILE_EXT(tmp_imzml, ".imzML");
  TEST_EXCEPTION(Exception::MissingInformation, f.store(tmp_imzml, exp))
}
END_SECTION


START_SECTION(void store rejects duplicate pixel coordinates)
{
  MSExperiment exp;
  exp.addSpectrum(makePixelSpectrum_(1, 1, 100.0, 1000.0));
  exp.addSpectrum(makePixelSpectrum_(1, 1, 101.0, 900.0));

  ImzMLFile f;
  String tmp_imzml;
  NEW_TMP_FILE_EXT(tmp_imzml, ".imzML");
  TEST_EXCEPTION(Exception::InvalidValue, f.store(tmp_imzml, exp))
}
END_SECTION


START_SECTION(void store rejects incompatible continuous mode)
{
  MSExperiment exp;
  exp.setMetaValue("imzml:imaging_mode", "continuous");
  exp.addSpectrum(makePixelSpectrum_(1, 1, 100.0, 1000.0));
  exp.addSpectrum(makePixelSpectrum_(2, 1, 200.0, 800.0));

  ImzMLFile f;
  String tmp_imzml;
  NEW_TMP_FILE_EXT(tmp_imzml, ".imzML");
  TEST_EXCEPTION(Exception::InvalidParameter, f.store(tmp_imzml, exp))
}
END_SECTION


START_SECTION(void buildImagingGeometry rejects duplicate pixels)
{
  MSExperiment exp;
  exp.addSpectrum(makePixelSpectrum_(1, 1, 100.0, 1000.0));
  exp.addSpectrum(makePixelSpectrum_(1, 1, 101.0, 900.0));

  MSImagingGeometry geom;
  TEST_EXCEPTION(Exception::InvalidValue, ImzMLFile::buildImagingGeometry(exp, geom))
}
END_SECTION


START_SECTION(void load rejects spectrum missing pixel coordinates)
{
  String tmp_imzml;
  NEW_TMP_FILE_EXT(tmp_imzml, ".imzML");
  writeMissingPixelCoordImzML_(tmp_imzml);

  ImzMLFile f;
  MSExperiment exp;
  TEST_EXCEPTION(Exception::ParseError, f.load(tmp_imzml, exp))

  String ibd_path = tmp_imzml;
  String lower = tmp_imzml;
  lower.toLower();
  if (lower.hasSuffix(".imzml"))
  {
    ibd_path = tmp_imzml.substr(0, tmp_imzml.size() - 6) + ".ibd";
  }
  remove(ibd_path.c_str());
}
END_SECTION


START_SECTION(OnDiscImzMLExperiment rejects duplicate pixel coordinates)
{
  String tmp_imzml;
  NEW_TMP_FILE_EXT(tmp_imzml, ".imzML");
  writeDuplicatePixelImzML_(tmp_imzml);

  OnDiscImzMLExperiment od;
  od.open(tmp_imzml);
  TEST_EXCEPTION(Exception::InvalidValue, od.getSpectrumAtCoord(1, 1, 1))

  String ibd_path = tmp_imzml;
  String lower = tmp_imzml;
  lower.toLower();
  if (lower.hasSuffix(".imzml"))
  {
    ibd_path = tmp_imzml.substr(0, tmp_imzml.size() - 6) + ".ibd";
  }
  remove(ibd_path.c_str());
}
END_SECTION


START_SECTION(void store applies PeakFileOptions m/z range filter)
{
  MSExperiment original;
  ImzMLFile f;
  f.load(imzml_path, original);
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

    String tmp_imzml;
    NEW_TMP_FILE_EXT(tmp_imzml, ".imzML");
    f.store(tmp_imzml, original);

    MSExperiment filtered;
    ImzMLFile f2;
    f2.load(tmp_imzml, filtered);
    TEST_EQUAL(filtered.getNrSpectra(), original.getNrSpectra())
    TEST_EQUAL(filtered[0].size() < full_peaks, true)
    TEST_EQUAL(filtered[0][0].getMZ() >= lo_mz, true)
    TEST_EQUAL(filtered[0].back().getMZ() <= mid_mz, true)
  }
}
END_SECTION


START_SECTION(void store honors PeakFileOptions binary precision)
{
  MSExperiment original;
  ImzMLFile f;
  f.load(imzml_path, original);

  PeakFileOptions opts;
  opts.setMz32Bit(false);
  opts.setIntensity32Bit(false);
  f.setOptions(opts);

  String tmp_imzml;
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
  MSExperiment original;
  ImzMLFile f;
  f.load(imzml_path, original);

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

  String tmp_imzml;
  NEW_TMP_FILE_EXT(tmp_imzml, ".imzML");
  f.store(tmp_imzml, original);

  MSExperiment reloaded;
  f.load(tmp_imzml, reloaded);
  TEST_EQUAL(String(reloaded.getMetaValue("imzml:scan_pattern")), String("top down"))
  TEST_EQUAL(String(reloaded.getMetaValue("imzml:scan_direction")), String("flyback"))
  TEST_EQUAL(String(reloaded.getMetaValue("imzml:line_scan_direction")), String("left-right"))
  TEST_EQUAL(String(reloaded.getMetaValue("imzml:polarity")), String("positive"))
  TEST_EQUAL(reloaded.metaValueExists("imzml:ibd_sha1"), true)
  TEST_NOT_EQUAL(String(reloaded.getMetaValue("imzml:ibd_sha1")).empty(), true)
  TEST_REAL_SIMILAR(static_cast<double>(reloaded.getMetaValue("imzml:max_dim_x")), 300.0)
  TEST_REAL_SIMILAR(static_cast<double>(reloaded.getMetaValue("imzml:max_dim_y")), 300.0)
  if (!reloaded.empty())
  {
    TEST_REAL_SIMILAR(reloaded[0].getRT(), 12.34)
    TEST_EQUAL(reloaded[0].getMSLevel(), 1)
  }
}
END_SECTION

END_TEST
