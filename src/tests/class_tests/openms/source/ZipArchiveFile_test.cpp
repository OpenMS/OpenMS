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
#include <OpenMS/SYSTEM/TempFiles.h>
#include <fstream>
#include <iterator>

using namespace OpenMS;

START_TEST(ZipArchiveFile, "$Id$")

START_SECTION(void addOrReplaceFromFile(const std::string&, const std::string&, const std::string&))
{
  // prepare temporary workspace
  TempDir tmp;
  const std::string base = tmp.getPath() + "/workspace";
  File::makeDir(base);
  File::makeDir(base + "/library");

  const std::string file1 = base + "/library/precursors.parquet";
  // write initial content
  {
    std::ofstream ofs(file1.c_str(), std::ios::binary);
    TEST_TRUE(ofs.is_open());
    ofs << "version1";
  }

  const std::string archive = tmp.getPath() + "/test.oswpq";

  // add file into archive
  ZipArchiveFile::addOrReplaceFromFile(archive, "library/precursors.parquet", file1);
  TEST_EQUAL(File::exists(archive), true)

  // verify listing contains the entry
  auto entries = ZipArchiveFile::listEntries(archive);
  bool found = false;
  for (const auto& e : entries)
    if (e == "library/precursors.parquet") found = true;
  TEST_EQUAL(found, true)

  // extract and verify content
  std::unique_ptr<TempDir> unpack_tmp;
  const std::string unpack_dir = ZipArchiveFile::unzipDirectory(archive, unpack_tmp);
  const std::string extracted = unpack_dir + "/library/precursors.parquet";
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
  std::unique_ptr<TempDir> unpack_tmp2;
  const std::string unpack_dir2 = ZipArchiveFile::unzipDirectory(archive, unpack_tmp2);
  const std::string extracted2 = unpack_dir2 + "/library/precursors.parquet";
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

START_SECTION([EXTRA] unzipDirectory error paths)
{
  TempDir tmp;

  // Corrupt archive: an existing, readable file that is not a valid ZIP (garbage
  // bytes, no central directory). libzip's zip_open() fails and unzipDirectory
  // must raise a clear exception rather than crash.
  std::string corrupt;
  NEW_TMP_FILE(corrupt);
  {
    std::ofstream ofs(corrupt.c_str(), std::ios::binary);
    TEST_TRUE(ofs.is_open())
    ofs << "this is not a valid ZIP archive -- just plain text bytes.";
  }
  TEST_EQUAL(File::exists(corrupt), true)
  {
    std::unique_ptr<TempDir> td;
    TEST_EXCEPTION(Exception::InvalidValue, ZipArchiveFile::unzipDirectory(corrupt, td))
  }

  // Missing / unreadable input: a clear FileNotFound.
  {
    std::unique_ptr<TempDir> td;
    TEST_EXCEPTION(Exception::FileNotFound, ZipArchiveFile::unzipDirectory("/nonexistent/path/to/archive.oswpq", td))
  }

  // Already-unpacked directory input: returned as-is (no extraction), supporting
  // the directory-or-zip bundle convention.
  {
    const std::string adir = tmp.getPath() + "/already_a_dir";
    File::makeDir(adir);
    std::unique_ptr<TempDir> td;
    TEST_STRING_EQUAL(ZipArchiveFile::unzipDirectory(adir, td), adir)
  }
}
END_SECTION

START_SECTION([EXTRA] unzipDirectory stops at an entry larger than it declares)
{
  // A crafted archive (zip bomb) can declare small entries that inflate to much more. libzip
  // does not stop an entry at its declared size, so unzipDirectory counts the bytes itself.
  // Here a 1 MiB entry claims 1000 bytes in its headers.
  TempDir tmp;
  const std::string src = tmp.getPath() + "/src";
  File::makeDir(src);
  {
    std::ofstream ofs((src + "/data.bin").c_str(), std::ios::binary);
    ofs << std::string(1 << 20, 'x');
  }
  const std::string archive = tmp.getPath() + "/lying.zip";
  ZipArchiveFile::zipDirectory(src, archive);

  std::string bytes;
  {
    std::ifstream ifs(archive.c_str(), std::ios::binary);
    bytes.assign(std::istreambuf_iterator<char>(ifs), std::istreambuf_iterator<char>());
  }
  const size_t local = bytes.find(std::string("PK\x03\x04", 4));
  const size_t central = bytes.find(std::string("PK\x01\x02", 4));
  TEST_NOT_EQUAL(local, std::string::npos)
  TEST_NOT_EQUAL(central, std::string::npos)
  const std::string declared("\xe8\x03\x00\x00", 4); // 1000, little endian
  bytes.replace(local + 22, 4, declared);   // uncompressed size in the local file header
  bytes.replace(central + 24, 4, declared); // and in the central directory header
  {
    std::ofstream ofs(archive.c_str(), std::ios::binary | std::ios::trunc);
    ofs << bytes;
  }

  std::unique_ptr<TempDir> td;
  std::string message;
  try
  {
    ZipArchiveFile::unzipDirectory(archive, td);
  }
  catch (const Exception::InvalidValue& e)
  {
    message = e.what();
  }
  TEST_EQUAL(message.find("larger than its declared size") != std::string::npos, true)
}
END_SECTION

START_SECTION([EXTRA] unzipDirectory rejects an entry that leaves the target directory)
{
  // The target directory is <TempDir>/parquet_unpacked. '../parquet_unpacked_x/f.txt' leads to
  // a sibling whose path starts with the target's path, so a string prefix check accepts it.
  TempDir tmp;
  const std::string file = tmp.getPath() + "/f.txt";
  {
    std::ofstream ofs(file.c_str());
    ofs << "payload";
  }
  const std::string archive = tmp.getPath() + "/escape.zip";
  ZipArchiveFile::addOrReplaceFromFile(archive, "../parquet_unpacked_x/f.txt", file);

  std::unique_ptr<TempDir> td;
  std::string message;
  try
  {
    ZipArchiveFile::unzipDirectory(archive, td);
  }
  catch (const Exception::InvalidValue& e)
  {
    message = e.what();
  }
  TEST_EQUAL(message.find("outside target directory") != std::string::npos, true)
}
END_SECTION

END_TEST
