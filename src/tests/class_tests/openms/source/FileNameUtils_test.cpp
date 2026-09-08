// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileNameUtils.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/SYSTEM/PathUtils.h>
#include <OpenMS/CONCEPT/ClassTest.h>

using namespace OpenMS;

START_TEST(FileNameUtils, "$Id$")

START_SECTION((filename rules and legacy entry points))
  struct Case
  {
    const char* filename;
    FileTypes::Type type;
    const char* stem;
  };
  const Case cases[] = {
    {"", FileTypes::UNKNOWN, ""},
    {"fid", FileTypes::XMASS, "fid"},
    {"sample.mzML.gz", FileTypes::MZML, "sample"},
    {"sample.mzML.bz2", FileTypes::MZML, "sample"},
    {"sample.d.zip", FileTypes::BRUKER_TDF, "sample"},
    {"sample.pep.xml", FileTypes::PEPXML, "sample.pep"},
    {"sample.prot.xml", FileTypes::PROTXML, "sample.prot"},
    {"sample.pep.xml.gz", FileTypes::PEPXML, "sample.pep.xml"},
    {"sample.newEnding", FileTypes::UNKNOWN, "sample"},
    {"/dotted.directory/sample", FileTypes::UNKNOWN, "/dotted.directory/sample"},
    {R"(C:\dotted.directory\sample.mzML.gz)", FileTypes::MZML, R"(C:\dotted.directory\sample)"},
    {R"(C:\dotted.directory\sample)", FileTypes::UNKNOWN, R"(C:\dotted.directory\sample)"}
  };
  for (const auto& entry : cases)
  {
    TEST_EQUAL(FileNameUtils::getTypeByFileName(entry.filename), entry.type)
    TEST_EQUAL(FileHandler::getTypeByFileName(entry.filename), entry.type)
    TEST_EQUAL(FileNameUtils::stripExtension(entry.filename), entry.stem)
    TEST_EQUAL(FileHandler::stripExtension(entry.filename), entry.stem)
    const std::string replacement = std::string(entry.stem) + ".mzML";
    TEST_EQUAL(FileNameUtils::swapExtension(entry.filename, FileTypes::MZML), replacement)
    TEST_EQUAL(FileHandler::swapExtension(entry.filename, FileTypes::MZML), replacement)
  }
  TEST_TRUE(FileNameUtils::hasValidExtension("out.unknownSuffix", FileTypes::MZML))
  TEST_TRUE(FileHandler::hasValidExtension("out.unknownSuffix", FileTypes::MZML))
  TEST_FALSE(FileNameUtils::hasValidExtension("out.idXML", FileTypes::MZML))
  TEST_FALSE(FileHandler::hasValidExtension("out.idXML", FileTypes::MZML))
END_SECTION

START_SECTION((platform-independent basename))
  TEST_EQUAL(PathUtils::basename(""), "")
  TEST_EQUAL(PathUtils::basename("sample"), "sample")
  TEST_EQUAL(PathUtils::basename("/some/path/"), "")
  TEST_EQUAL(PathUtils::basename(R"(C:\some/path\sample)"), "sample")
END_SECTION

END_TEST
