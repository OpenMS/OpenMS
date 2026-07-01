// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ANALYSIS/QUANTITATION/DDAWorkflowCommons.h>
#include <OpenMS/SYSTEM/File.h>

#include <map>
#include <string>

///////////////////////////

START_TEST(DDAWorkflowCommons, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

START_SECTION((static std::map<std::string, std::string> mapId2MzMLs(const std::map<std::string, std::string>& m2i)))
{
  // mapId2MzMLs simply inverts the mzML->id map into an id->mzML map
  std::map<std::string, std::string> m2i;
  m2i["/a/runA.mzML"] = "/x/runA.idXML";
  m2i["/b/runB.mzML"] = "/y/runB.idXML";

  std::map<std::string, std::string> i2m = DDAWorkflowCommons::mapId2MzMLs(m2i);

  TEST_EQUAL(i2m.size(), 2)
  TEST_STRING_EQUAL(i2m["/x/runA.idXML"], "/a/runA.mzML")
  TEST_STRING_EQUAL(i2m["/y/runB.idXML"], "/b/runB.mzML")

  // empty input -> empty output
  std::map<std::string, std::string> empty_in;
  TEST_EQUAL(DDAWorkflowCommons::mapId2MzMLs(empty_in).empty(), true)
}
END_SECTION

START_SECTION((static std::map<std::string, std::string> mapMzML2Ids(StringList& in, StringList& in_ids)))
{
  // Matching base names in the same order -> map of absolute mzML path -> absolute id path
  StringList in{"dirA/sample1.mzML", "dirB/sample2.mzML"};
  StringList ids{"idsA/sample1.idXML", "idsB/sample2.idXML"};

  std::map<std::string, std::string> mz2id = DDAWorkflowCommons::mapMzML2Ids(in, ids);

  TEST_EQUAL(mz2id.size(), 2)
  TEST_STRING_EQUAL(mz2id[File::absolutePath("dirA/sample1.mzML")], File::absolutePath("idsA/sample1.idXML"))
  TEST_STRING_EQUAL(mz2id[File::absolutePath("dirB/sample2.mzML")], File::absolutePath("idsB/sample2.idXML"))

  // Different number of files -> rejected
  StringList in_short{"a.mzML"};
  StringList ids_long{"a.idXML", "b.idXML"};
  TEST_EXCEPTION(Exception::IllegalArgument, DDAWorkflowCommons::mapMzML2Ids(in_short, ids_long))

  // Same base names but different order -> rejected
  StringList in_ord{"a.mzML", "b.mzML"};
  StringList ids_ord{"b.idXML", "a.idXML"};
  TEST_EXCEPTION(Exception::IllegalArgument, DDAWorkflowCommons::mapMzML2Ids(in_ord, ids_ord))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
