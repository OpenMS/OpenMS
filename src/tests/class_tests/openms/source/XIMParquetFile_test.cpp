// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FORMAT/XIMParquetFile.h>
#include <OpenMS/FORMAT/ParquetFilter.h>

using namespace OpenMS;
using namespace std;

START_TEST(XIMParquetFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(void load(std::vector<XIMMobilogram>& output) const)
{
  XIMParquetFile xim(OPENMS_GET_TEST_DATA_PATH("XIMParquetFile_23_input.xim"));
  std::vector<XIMMobilogram> mobilograms;
  xim.load(mobilograms);

  TEST_EQUAL(mobilograms.empty(), false)
  TEST_EQUAL(mobilograms[0].mobility.empty(), false)
  TEST_EQUAL(mobilograms[0].intensity.empty(), false)

  for (const auto& m : mobilograms)
  {
    TEST_EQUAL(m.mobility.size(), m.intensity.size())
  }
}
END_SECTION

START_SECTION(void getMobilograms(std::vector<XIMMobilogram>&, Int64, Int64, const String&, Int64, Int64, Int64, Int64, const String&, Int64, double, const String&) const)
{
  XIMParquetFile xim(OPENMS_GET_TEST_DATA_PATH("XIMParquetFile_23_input.xim"));

  std::vector<XIMMobilogram> all;
  xim.getMobilograms(all);
  TEST_EQUAL(all.empty(), false)

  std::vector<XIMMobilogram> none;
  xim.getMobilograms(none, -1, -1, "", -1, -1, -1, -1, "", -1, -1.0, "feature_id=99999");
  TEST_EQUAL(none.size(), 0)

  std::vector<XIMMobilogram> filtered_string;
  xim.getMobilograms(filtered_string, -1, -1, "", -1, -1, -1, all[0].run_id, "", -1, -1.0, "");
  TEST_EQUAL(filtered_string.empty(), false)
  TEST_EQUAL(filtered_string.size() <= all.size(), true)

  std::vector<XIMMobilogram> filtered_args;
  xim.getMobilograms(filtered_args, -1, -1, "", -1, -1, -1, all[0].run_id, "", -1, -1.0, "");
  TEST_EQUAL(filtered_args.empty(), false)
  TEST_EQUAL(filtered_args.size() <= all.size(), true)

  ParquetFilter typed_filter;
  typed_filter.eq("RUN_ID", all[0].run_id);
  std::vector<XIMMobilogram> typed;
  xim.getMobilograms(typed, typed_filter);
  TEST_EQUAL(typed.empty(), false)
  TEST_EQUAL(typed.size() <= all.size(), true)

  std::vector<XIMMobilogram> invalid_filter;
  TEST_EXCEPTION(Exception::InvalidValue,
                 xim.getMobilograms(invalid_filter, -1, -1, "", -1, -1, -1, -1, "", -1, -1.0, "feature_id=="))
}
END_SECTION

START_SECTION(void getMobilograms_multi_file)
{
  const String file = OPENMS_GET_TEST_DATA_PATH("XIMParquetFile_23_input.xim");
  XIMParquetFile xim_single(file);
  std::vector<XIMMobilogram> mobilograms_single;
  xim_single.getMobilograms(mobilograms_single);

  std::vector<String> files;
  files.push_back(file);
  files.push_back(file);

  XIMParquetFile xim_multi(files);
  std::vector<XIMMobilogram> mobilograms_multi;
  xim_multi.getMobilograms(mobilograms_multi);
  TEST_EQUAL(mobilograms_multi.size(), 2 * mobilograms_single.size())
}
END_SECTION

START_SECTION(void getRuns(std::vector<XIMRunInfo>& output) const)
{
  XIMParquetFile xim(OPENMS_GET_TEST_DATA_PATH("XIMParquetFile_23_input.xim"));
  std::vector<XIMRunInfo> runs;
  xim.getRuns(runs);
  TEST_EQUAL(runs.empty(), false)
  TEST_NOT_EQUAL(runs[0].run_id, 0)
}
END_SECTION

START_SECTION(void load_invalid_path)
{
  TEST_EXCEPTION(Exception::FileNotFound, XIMParquetFile("no_such_file.xim"))
}
END_SECTION

START_SECTION(void getColumns(std::vector<String>&) const)
{
  XIMParquetFile xim(OPENMS_GET_TEST_DATA_PATH("XIMParquetFile_23_input.xim"));
  std::vector<String> columns;
  xim.getColumns(columns);

  TEST_EQUAL(columns.empty(), false)
  // Schema must contain at least the mandatory data columns
  bool has_run_id = false;
  bool has_mobility = false;
  bool has_intensity = false;
  for (const auto& col : columns)
  {
    if (col == "RUN_ID") has_run_id = true;
    if (col == "MOBILITY_DATA") has_mobility = true;
    if (col == "INTENSITY_DATA") has_intensity = true;
  }
  TEST_EQUAL(has_run_id, true)
  TEST_EQUAL(has_mobility, true)
  TEST_EQUAL(has_intensity, true)
}
END_SECTION

START_SECTION(void getAnalytes(std::vector<XIMAnalyte>&, const std::vector<String>&, bool) const)
{
  XIMParquetFile xim(OPENMS_GET_TEST_DATA_PATH("XIMParquetFile_23_input.xim"));

  // Default: all analytes with nested transitions
  std::vector<XIMAnalyte> analytes;
  xim.getAnalytes(analytes);
  TEST_EQUAL(analytes.empty(), false)

  // Non-nested: each row has scalar transition fields
  std::vector<XIMAnalyte> flat;
  xim.getAnalytes(flat, {}, false);
  TEST_EQUAL(flat.empty(), false)

  // With a valid analyte column filter
  std::vector<XIMAnalyte> filtered;
  xim.getAnalytes(filtered, {"PRECURSOR_ID", "MODIFIED_SEQUENCE", "PRECURSOR_CHARGE"}, false);
  TEST_EQUAL(filtered.empty(), false)

  // Passing a non-analyte column (e.g. a signal column) must throw
  std::vector<XIMAnalyte> bad;
  TEST_EXCEPTION(Exception::InvalidValue, xim.getAnalytes(bad, {"RUN_ID"}, false))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
