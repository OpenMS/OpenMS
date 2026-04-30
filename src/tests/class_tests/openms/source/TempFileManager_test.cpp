// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Tilman Aurich $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/TempFileManager.h>

#include <fstream>

///////////////////////////

using namespace OpenMS;

namespace
{
  // compute expected registry file path for a given registry id
  String computeRegistryPath(const String& registry_id)
  {
    String temp_dir = File::getTempDirectory();
    while (temp_dir.hasSuffix("/") || temp_dir.hasSuffix("\\"))
    {
      temp_dir = temp_dir.prefix(temp_dir.size() - 1);
    }
    return temp_dir + "/" + registry_id + ".list";
  }

  void writeTextFile(const String& file_path, const String& content)
  {
    std::ofstream os(file_path.c_str(), std::ios::out | std::ios::trunc);
    os << content;
  }

  Size countNonEmptyLines(const String& file_path)
  {
    std::ifstream is(file_path.c_str());
    Size count = 0;
    String line;
    while (std::getline(is, line))
    {
      if (!line.empty() && line.back() == '\r')
      {
        line.pop_back();
      }
      if (!line.empty())
      {
        ++count;
      }
    }
    return count;
  }
}

START_TEST(TempFileManager, "$Id$")

START_SECTION((TempFileManager(const String& registry_id)))
{
  // construction with unique registry id should succeed
  const String registry_id = "temp_file_manager_ctor_test" + File::getUniqueName(false);
  TempFileManager manager(registry_id);
  // pass marker
  NOT_TESTABLE
}
END_SECTION

// addFile lifecycle test
START_SECTION((void addFile(const String& file_path)))
{
  const String registry_id = "temp_file_manager_add_test_" + File::getUniqueName(false);
  const String registry_path = computeRegistryPath(registry_id);
  const String temp_file = File::getTempDirectory() + "/tfm_add_" + File::getUniqueName(false) + ".tmp";

  writeTextFile(temp_file, "tmp");
  TEST_TRUE(File::exists(temp_file))

  {
    TempFileManager manager(registry_id);
    manager.addFile(temp_file);
    TEST_TRUE(File::exists(registry_path))
  }

  // After manager destruction, the temp file and registry should be cleaned up
  TEST_FALSE(File::exists(temp_file))
  TEST_FALSE(File::exists(registry_path))
}
END_SECTION

// duplicate add should not create duplicate registry entries
START_SECTION((void addFile(const String& file_path)))
{
  const String registry_id = "temp_file_manager_add_idempotence_test_" + File::getUniqueName(false);
  const String registry_path = computeRegistryPath(registry_id);
  const String temp_file = File::getTempDirectory() + "/tfm_add_idempotent_" + File::getUniqueName(false) + ".tmp";

  writeTextFile(temp_file, "tmp");

  {
    TempFileManager manager(registry_id);
    manager.addFile(temp_file);
    manager.addFile(temp_file);

    TEST_TRUE(File::exists(registry_path))
    TEST_EQUAL(countNonEmptyLines(registry_path), 1)
  }

  TEST_FALSE(File::exists(temp_file))
  TEST_FALSE(File::exists(registry_path))
}
END_SECTION

// releaseFile lifecycle test
START_SECTION((void releaseFile(const String& file_path)))
{
  const String registry_id = "temp_file_manager_release_test_" + File::getUniqueName(false);
  const String registry_path = computeRegistryPath(registry_id);
  const String temp_file = File::getTempDirectory() + "/tfm_release_" + File::getUniqueName(false) + ".tmp";

  writeTextFile(temp_file, "tmp");
  TEST_TRUE(File::exists(temp_file))

  {
    TempFileManager manager(registry_id);
    manager.addFile(temp_file);
    manager.releaseFile(temp_file);
  }

  // temp_file should still exist, but registry should be cleaned up
  TEST_TRUE(File::exists(temp_file))
  TEST_TRUE(File::remove(temp_file))
  TEST_FALSE(File::exists(registry_path))
}
END_SECTION

// releasing a non-registered path should not affect registered files
START_SECTION((void releaseFile(const String& file_path)))
{
  const String registry_id = "temp_file_manager_release_non_registered_test_" + File::getUniqueName(false);
  const String registry_path = computeRegistryPath(registry_id);
  const String managed_file = File::getTempDirectory() + "/tfm_managed_" + File::getUniqueName(false) + ".tmp";
  const String unmanaged_file = File::getTempDirectory() + "/tfm_unmanaged_" + File::getUniqueName(false) + ".tmp";

  writeTextFile(managed_file, "managed");
  writeTextFile(unmanaged_file, "unmanaged");

  {
    TempFileManager manager(registry_id);
    manager.addFile(managed_file);
    manager.releaseFile(unmanaged_file);
  }

  TEST_FALSE(File::exists(managed_file))
  TEST_TRUE(File::exists(unmanaged_file))
  TEST_TRUE(File::remove(unmanaged_file))
  TEST_FALSE(File::exists(registry_path))
}
END_SECTION

// test that releaseFile throws on empty input
START_SECTION((void releaseFile(const String& file_path)))
{
  const String registry_id = "temp_file_manager_release_empty_input_test_" + File::getUniqueName(false);
  TempFileManager manager(registry_id);
  TEST_EXCEPTION(Exception::InvalidParameter, manager.releaseFile(""))
}
END_SECTION

START_SECTION((bool removeFileNow(const String& file_path) noexcept))
{
  const String registry_id = "temp_file_manager_remove_now_test_" + File::getUniqueName(false);
  const String temp_file = File::getTempDirectory() + "/tfm_remove_now_" + File::getUniqueName(false) + ".tmp";

  writeTextFile(temp_file, "tmp");
  TEST_EQUAL(File::exists(temp_file), true)

  TempFileManager manager(registry_id);
  manager.addFile(temp_file);

  // check if removal works and returns true
  TEST_TRUE(manager.removeFileNow(temp_file))
  TEST_FALSE(File::exists(temp_file))
}
END_SECTION

// test if cleanupNow removes all files and registry correctly
START_SECTION((void cleanupNow() noexcept))
{
  const String registry_id = "temp_file_manager_cleanup_now_test_" + File::getUniqueName(false);
  const String registry_path = computeRegistryPath(registry_id);
  const String temp_file = File::getTempDirectory() + "/tfm_cleanup_now_" + File::getUniqueName(false) + ".tmp";

  writeTextFile(temp_file, "tmp");
  writeTextFile(registry_path, temp_file + "\n");

  TEST_TRUE(File::exists(temp_file))
  TEST_TRUE(File::exists(registry_path))

  {
    TempFileManager manager(registry_id);
    manager.cleanupNow();
  }

  TEST_FALSE(File::exists(temp_file))
  TEST_FALSE(File::exists(registry_path))
}
END_SECTION

// cleanupNow should merge persisted and in-memory entries and remove both
START_SECTION((void cleanupNow() noexcept))
{
  const String registry_id = "temp_file_manager_cleanup_merge_test_" + File::getUniqueName(false);
  const String registry_path = computeRegistryPath(registry_id);
  const String persisted_file = File::getTempDirectory() + "/tfm_persisted_" + File::getUniqueName(false) + ".tmp";
  const String in_memory_file = File::getTempDirectory() + "/tfm_in_memory_" + File::getUniqueName(false) + ".tmp";

  writeTextFile(persisted_file, "persisted");
  writeTextFile(in_memory_file, "in-memory");
  writeTextFile(registry_path, persisted_file + "\n");

  {
    TempFileManager manager(registry_id);
    manager.addFile(in_memory_file);
    manager.cleanupNow();
  }

  TEST_FALSE(File::exists(persisted_file))
  TEST_FALSE(File::exists(in_memory_file))
  TEST_FALSE(File::exists(registry_path))
}
END_SECTION

// test that addFile throws on empty input
START_SECTION((void addFile(const String& file_path)))
{
  const String registry_id = "temp_file_manager_empty_input_test_" + File::getUniqueName(false);
  TempFileManager manager(registry_id);
  TEST_EXCEPTION(Exception::InvalidParameter, manager.addFile(""))
}
END_SECTION

END_TEST
