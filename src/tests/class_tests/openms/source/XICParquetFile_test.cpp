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

#include <OpenMS/FORMAT/XICParquetFile.h>

using namespace OpenMS;
using namespace std;

START_TEST(XICParquetFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(void load(std::vector<XICChromatogram>& output) const)
{
#ifdef WITH_PARQUET
  XICParquetFile xic(OPENMS_GET_TEST_DATA_PATH("XICParquetFile_1_input.xic"));
  std::vector<XICChromatogram> chroms;
  xic.load(chroms);
  TEST_EQUAL(chroms.size(), 18)
  TEST_EQUAL(chroms[0].rt.empty(), false)
  TEST_EQUAL(chroms[0].intensity.empty(), false)

  // RT/intensity lengths should match for all rows
  for (const auto& c : chroms)
  {
    TEST_EQUAL(c.rt.size(), c.intensity.size())
  }
#else
  XICParquetFile xic("dummy.xic");
  std::vector<XICChromatogram> chroms;
  TEST_EXCEPTION(Exception::MissingFeature, xic.load(chroms))
#endif
}
END_SECTION

START_SECTION(void getChromatograms(std::vector<XICChromatogram>&, Int64, Int64, const String&, Int64, Int64, Int64, Int64, const String&) const)
{
#ifdef WITH_PARQUET
  XICParquetFile xic(OPENMS_GET_TEST_DATA_PATH("XICParquetFile_1_input.xic"));
  std::vector<XICChromatogram> chroms_all;
  xic.getChromatograms(chroms_all);
  TEST_EQUAL(chroms_all.size(), 18)

  std::vector<XICChromatogram> chroms_none;
  xic.getChromatograms(chroms_none, -1, -1, "", -1, -1, -1, -1, "precursor_id=99999");
  TEST_EQUAL(chroms_none.size(), 0)

  std::vector<XICChromatogram> chroms_filtered;
  xic.getChromatograms(chroms_filtered, -1, -1, "", -1, -1, -1, -1, "precursor_id=2");
  TEST_EQUAL(chroms_filtered.size() > 0, true)
  TEST_EQUAL(chroms_filtered.size() <= chroms_all.size(), true)

  // Invalid filter syntax should throw a parse/bind error
  std::vector<XICChromatogram> chroms_invalid_filter;
  TEST_EXCEPTION(Exception::InvalidValue,
                 xic.getChromatograms(chroms_invalid_filter, -1, -1, "", -1, -1, -1, -1, "precursor_id=="))
#endif
}
END_SECTION

START_SECTION(void getChromatograms_multi_file)
{
#ifdef WITH_PARQUET
  std::vector<String> files;
  files.emplace_back(OPENMS_GET_TEST_DATA_PATH("XICParquetFile_1_input.xic"));
  files.emplace_back(OPENMS_GET_TEST_DATA_PATH("XICParquetFile_2_input.xic"));

  XICParquetFile xic_multi(files);
  std::vector<XICChromatogram> chroms_multi;
  xic_multi.getChromatograms(chroms_multi);
  TEST_EQUAL(chroms_multi.size(), 36)
#endif
}
END_SECTION

START_SECTION(void load_invalid_path)
{
#ifdef WITH_PARQUET
  TEST_EXCEPTION(Exception::FileNotFound, XICParquetFile("no_such_file.xic"))
#endif
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
