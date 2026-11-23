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
#include <thread>
#ifdef _OPENMP
  #include <omp.h>
#endif

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

START_SECTION(static UInt getAvailablePhysicalCores())
{
  const UInt cores = SysInfo::getAvailablePhysicalCores();
  std::cout << "Physical cores available: " << cores << " (hardware threads: " << std::thread::hardware_concurrency() << ")" << std::endl;
  TEST_EQUAL(cores >= 1, true)
  // no more cores than hardware threads, where these are known
  TEST_EQUAL(std::thread::hardware_concurrency() == 0 || cores <= std::thread::hardware_concurrency(), true)
  TEST_EQUAL(SysInfo::getAvailablePhysicalCores(), cores) // the same again
}
END_SECTION

START_SECTION(static int ompThreadsOfThisThread())
{
  const int threads = SysInfo::ompThreadsOfThisThread();
  TEST_EQUAL(threads >= 1, true)
  // a number set on another thread stays there
  int other = 0;
  std::thread([&other]() { SysInfo::setOmpThreadsOfThisThread(3); other = SysInfo::ompThreadsOfThisThread(); }).join();
  TEST_EQUAL(other, 3)
  TEST_EQUAL(SysInfo::ompThreadsOfThisThread(), threads)
}
END_SECTION

START_SECTION(static void setOmpThreadsOfThisThread(int threads))
{
  const int threads = SysInfo::ompThreadsOfThisThread();
  SysInfo::setOmpThreadsOfThisThread(2);
  TEST_EQUAL(SysInfo::ompThreadsOfThisThread(), 2)
  SysInfo::setOmpThreadsOfThisThread(0); // back to the runtime's number
  TEST_EQUAL(SysInfo::ompThreadsOfThisThread(), threads)
  SysInfo::setOmpThreadsOfThisThread(-5); // like 0
  TEST_EQUAL(SysInfo::ompThreadsOfThisThread(), threads)
}
END_SECTION

START_SECTION(static bool ompThreadCountIsGlobal())
{
#ifdef _OPENMP
  const int before = omp_get_max_threads();
#endif
  const bool global = SysInfo::ompThreadCountIsGlobal();
  TEST_EQUAL(SysInfo::ompThreadCountIsGlobal(), global) // probed once
#ifdef _OPENMP
  TEST_EQUAL(omp_get_max_threads(), before) // the probe leaves the caller's number as it was
  // what it says is what the runtime does: another thread's number reaches this thread only where the number is global
  std::thread([before]() { omp_set_num_threads(before + 2); }).join();
  TEST_EQUAL(omp_get_max_threads() != before, global)
  omp_set_num_threads(before);
  #if _OPENMP >= 200805
  TEST_EQUAL(global, false) // OpenMP 3.0 runtimes keep the number per thread
  #endif
#else
  TEST_EQUAL(global, false)
#endif
}
END_SECTION

END_TEST
