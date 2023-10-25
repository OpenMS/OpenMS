// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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
        String(OPENMS_GET_TEST_DATA_PATH("FIAMS_input/params_test.csv"))
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
    FIAMSScheduler fia_scheduler(
        String(OPENMS_GET_TEST_DATA_PATH("FIAMS_input/params_test.csv")),
        String(OPENMS_GET_TEST_DATA_PATH(""))
    );
    const vector<map<String, String>> samples = fia_scheduler.getSamples();
    TEST_EQUAL(samples[0].at("time"), "10");
    fia_scheduler.run();
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST