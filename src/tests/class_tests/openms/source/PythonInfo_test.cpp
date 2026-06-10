// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

/////////////////////////////////////////////////////////////

#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/SYSTEM/PythonInfo.h>
#include <OpenMS/SYSTEM/File.h>
            
#include <fstream>
#include <filesystem>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(TextFile, "$Id$")

/////////////////////////////////////////////////////////////

START_SECTION((static bool canRun(std::string& python_executable, std::string& error_msg)))
  // test for missing python executable
  std::string py = "does_not_exist_@@";
  std::string error_msg;
  TEST_EQUAL(PythonInfo::canRun(py, error_msg), false)
  TEST_EQUAL(StringUtils::hasSubstring(error_msg, "Python not found at"), true)

  auto tmp_file = File::getTemporaryFile();
  ofstream f(tmp_file); // create the file
  f.close(); 
  TEST_EQUAL(PythonInfo::canRun(tmp_file, error_msg), false)
  TEST_EQUAL(StringUtils::hasSubstring(error_msg, "failed to run"), true)  

  py = "python";
  if (PythonInfo::canRun(py, error_msg))
  { 
    TEST_EQUAL(File::exists(py), true)
    TEST_EQUAL(std::filesystem::path(std::string(py)).is_relative(), false)
  }

END_SECTION

START_SECTION(bool PythonInfo::isPackageInstalled(const std::string& python_executable, const std::string& package_name))
  std::string error_msg;
  std::string py = "python";
  if (PythonInfo::canRun(py, error_msg))
  {
    TEST_EQUAL(PythonInfo::isPackageInstalled(py, "veryWeirdPackage___@@__@"), false)
    TEST_EQUAL(PythonInfo::isPackageInstalled(py, "math"), true)
  }
END_SECTION

START_SECTION(static std::string getVersion(const std::string& python_executable))
  
  std::string py = "python";
  std::string error_msg;
  if (PythonInfo::canRun(py, error_msg))
  {
    std::string version = PythonInfo::getVersion(py);
    TEST_EQUAL(version.empty(), false)
  }
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
