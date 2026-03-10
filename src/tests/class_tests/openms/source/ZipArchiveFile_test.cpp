// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/ZipArchiveFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <fstream>

using namespace OpenMS;

START_TEST(ZipArchiveFile, "$Id$")

START_SECTION(void addOrReplaceFromFile(const String&, const String&, const String&))
{
  // prepare temporary workspace
  File::TempDir tmp;
  const String base = tmp.getPath() + "/workspace";
  File::makeDir(base);
  File::makeDir(base + "/library");

  const String file1 = base + "/library/precursors.parquet";
  // write initial content
  {
    std::ofstream ofs(file1.c_str(), std::ios::binary);
  TEST_TRUE(ofs.is_open());
    ofs << "version1";
  }

  const String archive = tmp.getPath() + "/test.oswpq";

  // add file into archive
  ZipArchiveFile::addOrReplaceFromFile(archive, "library/precursors.parquet", file1);
  TEST_EQUAL(File::exists(archive), true)

  // verify listing contains the entry
  auto entries = ZipArchiveFile::listEntries(archive);
  bool found = false;
  for (const auto &e : entries) if (e == "library/precursors.parquet") found = true;
  TEST_EQUAL(found, true)

  // extract and verify content
  std::unique_ptr<File::TempDir> unpack_tmp;
  const String unpack_dir = ZipArchiveFile::unzipDirectory(archive, unpack_tmp);
  const String extracted = unpack_dir + "/library/precursors.parquet";
  TEST_EQUAL(File::exists(extracted), true)
  {
    std::ifstream ifs(extracted.c_str(), std::ios::binary);
  TEST_TRUE(ifs.is_open());
    std::string content;
    std::getline(ifs, content);
    TEST_EQUAL(content, "version1");
  }

  // replace source file content
  {
    std::ofstream ofs(file1.c_str(), std::ios::binary | std::ios::trunc);
  TEST_TRUE(ofs.is_open());
    ofs << "version2";
  }

  // perform replace in archive
  ZipArchiveFile::addOrReplaceFromFile(archive, "library/precursors.parquet", file1);

  // extract to new temp dir and verify replaced content
  std::unique_ptr<File::TempDir> unpack_tmp2;
  const String unpack_dir2 = ZipArchiveFile::unzipDirectory(archive, unpack_tmp2);
  const String extracted2 = unpack_dir2 + "/library/precursors.parquet";
  TEST_EQUAL(File::exists(extracted2), true)
  {
    std::ifstream ifs(extracted2.c_str(), std::ios::binary);
  TEST_TRUE(ifs.is_open());
    std::string content;
    std::getline(ifs, content);
    TEST_EQUAL(content, "version2");
  }
}
END_SECTION

END_TEST
