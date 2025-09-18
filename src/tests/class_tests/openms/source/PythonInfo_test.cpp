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

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/SYSTEM/PythonInfo.h>
#include <OpenMS/SYSTEM/File.h>
            
#include <fstream>

#include <QDir>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(TextFile, "$Id$")

/////////////////////////////////////////////////////////////

START_SECTION((static bool canRun(String& python_executable, String& error_msg)))
  // test for missing python executable
  String py = "does_not_exist_@@";
  String error_msg;
  TEST_EQUAL(PythonInfo::canRun(py, error_msg), false)
  TEST_EQUAL(error_msg.hasSubstring("Python not found at"), true)

  auto tmp_file = File::getTemporaryFile();
  ofstream f(tmp_file); // create the file
  f.close(); 
  TEST_EQUAL(PythonInfo::canRun(tmp_file, error_msg), false)
  TEST_EQUAL(error_msg.hasSubstring("failed to run"), true)  

  py = "python";
  if (PythonInfo::canRun(py, error_msg))
  { 
    TEST_EQUAL(File::exists(py), true)
    TEST_EQUAL(QDir::isRelativePath(py.toQString()), false)
  }

END_SECTION

START_SECTION(bool PythonInfo::isPackageInstalled(const String& python_executable, const String& package_name))
  String error_msg;
  String py = "python";
  if (PythonInfo::canRun(py, error_msg))
  {
    TEST_EQUAL(PythonInfo::isPackageInstalled(py, "veryWeirdPackage___@@__@"), false)
    TEST_EQUAL(PythonInfo::isPackageInstalled(py, "math"), true)
  }
END_SECTION

START_SECTION(static String getVersion(const String& python_executable))
  
  String py = "python";
  String error_msg;
  if (PythonInfo::canRun(py, error_msg))
  {
    String version = PythonInfo::getVersion(py);
    TEST_EQUAL(version.empty(), false)
  }
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
