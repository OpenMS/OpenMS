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

#include <OpenMS/SYSTEM/Network.h>
#include <OpenMS/SYSTEM/File.h>

#include <filesystem>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(Network, "$Id$")

/////////////////////////////////////////////////////////////

START_SECTION(static void downloadFile(const std::string& url, const std::string& download_folder))
{
  // Use a local file:// URL to avoid external network dependency
  std::string fixture_path = OPENMS_GET_TEST_DATA_PATH("Network_test_fixture.txt");
  std::string url = "file://" + fixture_path;

  std::string folder = File::getTempDirectory();

  Network::downloadFile(url, folder);
  std::string output_file_path = folder + "/Network_test_fixture.txt";

  TEST_EQUAL(File::exists(output_file_path), true);
  if (File::exists(output_file_path))
  {
    std::filesystem::remove(output_file_path);
  }
}
END_SECTION

/////////////////////////////////////////////////////////////

END_TEST
