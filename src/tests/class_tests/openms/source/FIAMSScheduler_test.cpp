// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Svetlana Kutuzova, Douglas McCloskey $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/FIAMSScheduler.h>

#include <OpenMS/SYSTEM/File.h>

#include <filesystem>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FIAMSScheduler, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FIAMSScheduler* ptr_1 = nullptr;
FIAMSScheduler* null_ptr_2 = nullptr;
START_SECTION(FIAMSScheduler())
{
    ptr_1 = new FIAMSScheduler(
        StringUtils::toStr(OPENMS_GET_TEST_DATA_PATH("FIAMS_input/params_test.csv"))
    );
    TEST_NOT_EQUAL(ptr_1, null_ptr_2);
    TEST_EQUAL(ptr_1->getBaseDir(), "/");
}
END_SECTION

START_SECTION(virtual ~FIAMSScheduler())
{
    delete ptr_1;
}
END_SECTION

START_SECTION(FIAMSScheduler)
{
  std::string tmp_dir = std::filesystem::current_path().generic_string()  + "/"; // write output to current directory
  FIAMSScheduler fia_scheduler(
      StringUtils::toStr(OPENMS_GET_TEST_DATA_PATH("FIAMS_input/params_test.csv")),
      StringUtils::toStr(OPENMS_GET_TEST_DATA_PATH("")),
      tmp_dir
  );
  const vector<map<std::string, std::string>> samples = fia_scheduler.getSamples();
  TEST_EQUAL(samples[0].at("time"), "10");
  fia_scheduler.run();
  std::string outfile =StringUtils::toStr(OPENMS_GET_TEST_DATA_PATH("FIAMS_output/SerumTest_10.mzTab"));
  std::string outfile2 = tmp_dir + "FIAMS_output/SerumTest_10.mzTab";
  TEST_FILE_EQUAL(outfile2.c_str(), outfile.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST