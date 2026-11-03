// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/XIPMParquetFile.h>

using namespace OpenMS;
using namespace std;

START_TEST(XIPMParquetFile, "$Id$")

START_SECTION(void load(std::vector<XIPMPeakMap>& output) const)
{
  const std::string file = OPENMS_GET_TEST_DATA_PATH("XIPMParquetFile_reader_input.xipm");
  XIPMParquetFile xipm(file);

  std::vector<XIPMParquetFile::XIPMPeakMap> peak_maps;
  xipm.load(peak_maps);
  TEST_EQUAL(peak_maps.size(), 2)
  TEST_EQUAL(peak_maps[0].mz.size(), peak_maps[0].rt.size())
  TEST_EQUAL(peak_maps[0].mz.size(), peak_maps[0].ion_mobility.size())
  TEST_EQUAL(peak_maps[0].mz.size(), peak_maps[0].intensity.size())
  TEST_EQUAL(peak_maps[0].has_target_rt, true)
  TEST_REAL_SIMILAR(peak_maps[0].target_rt, 100.0)
}
END_SECTION

START_SECTION(void getPeakMaps(...filters...) const)
{
  const std::string file = OPENMS_GET_TEST_DATA_PATH("XIPMParquetFile_reader_input.xipm");
  XIPMParquetFile xipm(file);

  std::vector<XIPMParquetFile::XIPMPeakMap> precursor_peak_maps;
  xipm.getPeakMaps(precursor_peak_maps, -1, -1, "", -1, -1, 1, 7, "precursor");
  TEST_EQUAL(precursor_peak_maps.size(), 1)
  TEST_STRING_EQUAL(precursor_peak_maps[0].peakmap_type, "precursor")

  std::vector<XIPMParquetFile::XIPMPeakMap> transition_peak_maps;
  xipm.getPeakMaps(transition_peak_maps, -1, 1, "", -1, -1, 2, 7, "transition");
  TEST_EQUAL(transition_peak_maps.size(), 1)
  TEST_STRING_EQUAL(transition_peak_maps[0].annotation, "y7^1")
  TEST_EQUAL(transition_peak_maps[0].has_target_rt, true)
  TEST_REAL_SIMILAR(transition_peak_maps[0].target_rt, 100.0)
}
END_SECTION

START_SECTION(void getPeakMaps_multi_file)
{
  const std::string file = OPENMS_GET_TEST_DATA_PATH("XIPMParquetFile_reader_input.xipm");
  std::vector<std::string> files = {file, file};
  XIPMParquetFile xipm(files);

  std::vector<XIPMParquetFile::XIPMPeakMap> peak_maps;
  xipm.getPeakMaps(peak_maps);
  TEST_EQUAL(peak_maps.size(), 4)
}
END_SECTION

START_SECTION(void getRuns(std::vector<XIPMRunInfo>& output) const)
{
  const std::string file = OPENMS_GET_TEST_DATA_PATH("XIPMParquetFile_reader_input.xipm");
  XIPMParquetFile xipm(file);

  std::vector<XIPMParquetFile::XIPMRunInfo> runs;
  xipm.getRuns(runs);
  TEST_EQUAL(runs.size(), 1)
  TEST_EQUAL(runs[0].run_id, 7)
}
END_SECTION

START_SECTION(void getPeakMaps_empty_file)
{
  const std::string file = OPENMS_GET_TEST_DATA_PATH("XIPMParquetFile_reader_empty.xipm");
  XIPMParquetFile xipm(file);

  std::vector<XIPMParquetFile::XIPMPeakMap> peak_maps;
  xipm.getPeakMaps(peak_maps);
  TEST_EQUAL(peak_maps.size(), 0)
}
END_SECTION

START_SECTION(void getRuns_empty_file)
{
  const std::string file = OPENMS_GET_TEST_DATA_PATH("XIPMParquetFile_reader_empty.xipm");
  XIPMParquetFile xipm(file);

  std::vector<XIPMParquetFile::XIPMRunInfo> runs;
  xipm.getRuns(runs);
  TEST_EQUAL(runs.size(), 0)
}
END_SECTION

START_SECTION(void getColumns(std::vector<std::string>& output) const)
{
  const std::string file = OPENMS_GET_TEST_DATA_PATH("XIPMParquetFile_reader_input.xipm");
  XIPMParquetFile xipm(file);

  std::vector<std::string> columns;
  xipm.getColumns(columns);
  TEST_EQUAL(columns.empty(), false)

  bool has_mz = false;
  bool has_rt = false;
  bool has_im = false;
  bool has_target_rt = false;
  for (const auto& col : columns)
  {
    if (col == "MZ_DATA") has_mz = true;
    if (col == "RT_DATA") has_rt = true;
    if (col == "MOBILITY_DATA") has_im = true;
    if (col == "TARGET_RT") has_target_rt = true;
  }
  TEST_EQUAL(has_mz, true)
  TEST_EQUAL(has_rt, true)
  TEST_EQUAL(has_im, true)
  TEST_EQUAL(has_target_rt, true)
}
END_SECTION

START_SECTION(void load_invalid_path)
{
  TEST_EXCEPTION(Exception::FileNotFound, XIPMParquetFile("no_such_file.xipm"))
}
END_SECTION

END_TEST
