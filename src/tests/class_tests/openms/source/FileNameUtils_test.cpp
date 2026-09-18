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
    {"/dotted.directory/fid", FileTypes::XMASS, "/dotted.directory/fid"},
    // 'fid' has no extension of its own, so a compression suffix is the whole extension
    {"fid.gz", FileTypes::XMASS, "fid"},
    {"fid.bz2", FileTypes::XMASS, "fid"},
    {"/dotted.directory/fid.zip", FileTypes::XMASS, "/dotted.directory/fid"},
    {"sample.mzML.gz", FileTypes::MZML, "sample"},
    {"sample.mzML.bz2", FileTypes::MZML, "sample"},
    {"sample.d.zip", FileTypes::BRUKER_TDF, "sample"},
    // compound extensions are stripped as a whole, and a compression suffix goes with them
    {"sample.pep.xml", FileTypes::PEPXML, "sample"},
    {"sample.prot.xml", FileTypes::PROTXML, "sample"},
    {"sample.pep.xml.gz", FileTypes::PEPXML, "sample"},
    {"sample.pepXML", FileTypes::PEPXML, "sample"},
    // the longest match wins, so these beat a plain '.xml'
    {"sample.xquest.xml", FileTypes::XQUESTXML, "sample"},
    {"sample.spec.xml", FileTypes::SPECXML, "sample"},
    {"sample.xml", FileTypes::XML, "sample"},
    // aliases are recognized alongside the preferred extension, case-insensitively
    {"database.fasta", FileTypes::FASTA, "database"},
    {"database.fa", FileTypes::FASTA, "database"},
    {"database.faa", FileTypes::FASTA, "database"},
    {"DATABASE.FA", FileTypes::FASTA, "DATABASE"},
    {"database.fa.gz", FileTypes::FASTA, "database"},
    {"table.pqt", FileTypes::PARQUET, "table"},
    {"table.parquet", FileTypes::PARQUET, "table"},
    // a dot inside the stem must not be mistaken for the extension
    {"HeLa_1.2ug.mzML", FileTypes::MZML, "HeLa_1.2ug"},
    {"sample.newEnding", FileTypes::UNKNOWN, "sample"},
    {"archive.gz", FileTypes::UNKNOWN, "archive"},
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
  // an alias is as valid as the preferred extension, and a genuine mismatch still fails
  TEST_TRUE(FileNameUtils::hasValidExtension("db.fa", FileTypes::FASTA))
  TEST_TRUE(FileNameUtils::hasValidExtension("db.fasta", FileTypes::FASTA))
  TEST_FALSE(FileNameUtils::hasValidExtension("db.fa", FileTypes::MZML))
  TEST_FALSE(FileNameUtils::hasValidExtension("out.csv", FileTypes::TSV))

  // swapExtension always writes the preferred extension, never an alias
  TEST_STRING_EQUAL(FileNameUtils::swapExtension("db.fa", FileTypes::FASTA), "db.fasta")
  TEST_STRING_EQUAL(FileNameUtils::swapExtension("sample.pep.xml", FileTypes::PEPXML), "sample.pepXML")
  TEST_STRING_EQUAL(FileNameUtils::swapExtension("table.pqt", FileTypes::PARQUET), "table.parquet")

  // a compression suffix is visible separately from the type it wraps
  TEST_TRUE(FileNameUtils::hasCompressionSuffix("sample.mzML.gz"))
  TEST_TRUE(FileNameUtils::hasCompressionSuffix("sample.mzML.bz2"))
  TEST_TRUE(FileNameUtils::hasCompressionSuffix("sample.d.zip"))
  TEST_TRUE(FileNameUtils::hasCompressionSuffix("sample.MZML.GZ"))
  TEST_FALSE(FileNameUtils::hasCompressionSuffix("sample.mzML"))
  TEST_FALSE(FileNameUtils::hasCompressionSuffix("fid"))
  TEST_FALSE(FileNameUtils::hasCompressionSuffix(""))
  // a dot in a directory name is not an extension, let alone a compression suffix
  TEST_FALSE(FileNameUtils::hasCompressionSuffix("/my.gz/sample"))

  // filename-only detection must work for output files that do not exist yet
  TEST_EQUAL(FileNameUtils::getTypeByFileName("does_not_exist_anywhere.fa"), FileTypes::FASTA)
END_SECTION

START_SECTION((platform-independent basename))
  TEST_EQUAL(PathUtils::basename(""), "")
  TEST_EQUAL(PathUtils::basename("sample"), "sample")
  TEST_EQUAL(PathUtils::basename("/some/path/"), "")
  TEST_EQUAL(PathUtils::basename(R"(C:\some/path\sample)"), "sample")
END_SECTION

END_TEST
