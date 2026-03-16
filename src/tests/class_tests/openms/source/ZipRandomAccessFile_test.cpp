// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/ZipRandomAccessFile.h>
#include <OpenMS/FORMAT/ZipArchiveFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <fstream>
#include <vector>
#include <cstring>
#include <filesystem>

#if __has_include(<zip.h>)
#include <zip.h>
#endif

using namespace OpenMS;

START_TEST(ZipRandomAccessFile, "$Id$")

START_SECTION(static arrow::Result<std::shared_ptr<arrow::io::RandomAccessFile>> Open(const String&, const String&, std::unique_ptr<File::TempDir>&))
{
  // create temporary directory
  std::unique_ptr<File::TempDir> tmpdir(new File::TempDir(false));
  String tmp_path = tmpdir->getPath();

#if __has_include(<zip.h>)
  // create a zip archive with one entry
  std::string archive = (std::filesystem::path(std::string(tmp_path)) / "test.zip").string();
  int err = 0;
  zip_t* za = zip_open(archive.c_str(), ZIP_CREATE | ZIP_TRUNCATE, &err);
  TEST_NOT_EQUAL(za, nullptr)

  const char* entry_name = "hello.txt";
  const char* data = "Hello ZipRandomAccessFile";
  zip_source_t* src = zip_source_buffer(za, data, strlen(data), 0);
  TEST_NOT_EQUAL(src, nullptr)
  TEST_EQUAL(zip_file_add(za, entry_name, src, ZIP_FL_ENC_UTF_8) < 0, false)

  zip_close(za);

  // open via ZipRandomAccessFile
  auto raf_res = ZipRandomAccessFile::Open(archive, entry_name, tmpdir);
  TEST_TRUE(raf_res.ok())
  auto raf = raf_res.ValueOrDie();

  auto size_res = raf->GetSize();
  TEST_TRUE(size_res.ok())
  int64_t size = size_res.ValueOrDie();
  TEST_EQUAL(size, static_cast<int64_t>(strlen(data)));

  // read at offset 6 (should read "ZipRandomAccessFile" tail)
  std::vector<char> buf(10);
  auto read_res = raf->ReadAt(6, 10, buf.data());
  TEST_TRUE(read_res.ok())
  int64_t read_bytes = read_res.ValueOrDie();
  TEST_EQUAL(read_bytes, 10);

  // verify a few bytes
  std::string got(buf.begin(), buf.end());
  TEST_NOT_EQUAL(got.find("ZipRandom"), std::string::npos)

  // cleanup via Close
  TEST_TRUE(raf->Close().ok())

#else
  // If libzip not available, ensure the fallback path works (we'll create a directory layout)
  std::string dir = (std::filesystem::path(std::string(tmp_path)) / "unpacked").string();
  std::filesystem::create_directory(dir);
  std::string entry = (std::filesystem::path(dir) / "hello.txt").string();
  std::ofstream ofs(entry);
  ofs << "Hello ZipRandomAccessFile";
  ofs.close();

  auto raf_res = ZipRandomAccessFile::Open(dir, "hello.txt", tmpdir);
  TEST_TRUE(raf_res.ok())
  auto raf = raf_res.ValueOrDie();
  auto size_res2 = raf->GetSize();
  TEST_TRUE(size_res2.ok())
  int64_t size2 = size_res2.ValueOrDie();
  TEST_EQUAL(size2, static_cast<int64_t>(strlen("Hello ZipRandomAccessFile")));
  TEST_TRUE(raf->Close().ok());
#endif

}

END_SECTION

END_TEST
