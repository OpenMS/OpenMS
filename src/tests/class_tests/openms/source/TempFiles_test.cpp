// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/SYSTEM/TempFiles.h>
///////////////////////////

#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/SystemSettings.h>

#include <string>

using namespace OpenMS;
using namespace std;

START_TEST(TempFiles, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(TempDir::TempDir(bool keep_dir = false))
{
  TempDir* dir = new TempDir();
  TempDir* nullPointer = nullptr;
  TEST_NOT_EQUAL(dir, nullPointer)
  TEST_TRUE(File::isDirectory(dir->getPath()))
  // created below the configured temp directory
  TEST_EQUAL(File::absolutePath(dir->getPath()).find(File::absolutePath(SystemSettings::getTempDirectory())), 0)
  delete dir;
}
END_SECTION

START_SECTION(TempDir::TempDir(const std::string& base_dir, bool keep_dir = false))
{
  TempDir parent;
  std::string child_path;
  {
    TempDir child(parent.getPath(), false);
    child_path = child.getPath();
    TEST_TRUE(File::isDirectory(child_path))
    TEST_EQUAL(child_path.find(parent.getPath()), 0)
  }
  TEST_FALSE(File::exists(child_path))
  TEST_TRUE(File::exists(parent.getPath()))
}
END_SECTION

START_SECTION(TempDir::~TempDir())
{
  std::string path;
  {
    TempDir dir;
    path = dir.getPath();
    TEST_EQUAL(File::exists(path), 1)
  }
  TEST_EQUAL(File::exists(path), 0)
  if (File::exists(path)) File::removeDir(path);
  {
    TempDir dir2(true);
    path = dir2.getPath();
    TEST_EQUAL(File::exists(path), 1)
  }
  TEST_EQUAL(File::exists(path), 1)
  if (File::exists(path)) File::removeDir(path);
}
END_SECTION

START_SECTION(const std::string& TempDir::getPath() const)
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION(static std::string TempFiles::getTemporaryFile(const std::string& alternative_file = ""))
{
  const std::string first = TempFiles::getTemporaryFile();
  const std::string second = TempFiles::getTemporaryFile();
  TEST_FALSE(first.empty())
  TEST_NOT_EQUAL(first, second)
  // names are allocated below the configured temp directory; the file itself is not created
  TEST_EQUAL(File::absolutePath(first).find(File::absolutePath(SystemSettings::getTempDirectory())), 0)
  TEST_FALSE(File::exists(first))
  TEST_EQUAL(TempFiles::getTemporaryFile("retain-this-filename"), "retain-this-filename")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
