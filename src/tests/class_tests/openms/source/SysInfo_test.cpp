// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/SYSTEM/SysInfo.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <iostream>

///////////////////////////

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>


using namespace OpenMS;


START_TEST(SysInfo, "$Id$")


START_SECTION(std::string bytesToHumanReadable(UInt64 bytes))
{
  TEST_EQUAL(bytesToHumanReadable(   2ull << 00), "2 byte")
  TEST_EQUAL(bytesToHumanReadable(2048ull << 00), "2 KiB")
  TEST_EQUAL(bytesToHumanReadable(2048ull << 10), "2 MiB")
  TEST_EQUAL(bytesToHumanReadable(2048ull << 20), "2 GiB")
  TEST_EQUAL(bytesToHumanReadable(2048ull << 30), "2 TiB")
  TEST_EQUAL(bytesToHumanReadable(2048ull << 40), "2 PiB")
}
END_SECTION

START_SECTION(static bool getProcessMemoryConsumption(size_t& mem_virtual))
{
  size_t first, after, final;
  TEST_EQUAL(SysInfo::getProcessMemoryConsumption(first), true);
  std::cout << "Memory consumed initally: " << first << " KB" << std::endl;

  {
    PeakMap exp;
    MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_5_long.mzML"), exp);

    TEST_EQUAL(SysInfo::getProcessMemoryConsumption(after), true);
    std::cout << "Memory consumed after reading 20 MB mzML : " << after << " KB" << std::endl;

    TEST_EQUAL(after - first > 10000, true)
  }

  TEST_EQUAL(SysInfo::getProcessMemoryConsumption(final), true);
  std::cout << "Memory consumed after release of MSExperiment: " << final << " KB" << std::endl;
  // just for fun. There is probably no guarantee that we get the whole mem back by the memory manager
  // (and indeed, it does not work on all OS's; e.g. on Linux, the page tables will remain in RAM, unless mem pressure is high)
  //TEST_EQUAL(after > final, true)

}
END_SECTION

END_TEST
