// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#ifdef _OPENMP
#include <omp.h>
#endif

#include <OpenMS/METADATA/MetaInfoRegistry.h>

///////////////////////////

START_TEST(MetaInfoRegistry, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

MetaInfoRegistry* test = nullptr;
MetaInfoRegistry* nullPointer = nullptr;
START_SECTION((MetaInfoRegistry()))
	test = new MetaInfoRegistry;
  TEST_NOT_EQUAL(test, nullPointer)
END_SECTION


START_SECTION((~MetaInfoRegistry()))
	delete test;
END_SECTION

MetaInfoRegistry mir;

START_SECTION((UInt registerName(const String& name, const String& description = "", const String& unit = "")))
{
	UInt testname = mir.registerName("testname", "this is just a test");
  TEST_EQUAL(testname, 1024);
	UInt retention_time = mir.registerName("retention time", "this is just another test", "sec");
  TEST_EQUAL(retention_time, 1025);
	UInt another_testname = mir.registerName("another testname", "i will be set later", "me too");
  TEST_EQUAL(another_testname, 1026);
}
END_SECTION

START_SECTION((void setDescription(UInt index, const String& description)))
	mir.setDescription(1026, "foo");
	TEST_STRING_EQUAL(mir.getDescription(1026), "foo")
END_SECTION

START_SECTION((void setDescription(const String& name, const String& description)))
	mir.setDescription("another testname", "bar");
	TEST_STRING_EQUAL(mir.getDescription(1026), "bar")
END_SECTION

START_SECTION((void setUnit(UInt index, const String& unit)))
	mir.setUnit(1026, "foo");
	TEST_STRING_EQUAL(mir.getUnit(1026), "foo")
END_SECTION

START_SECTION((void setUnit(const String& name, const String& unit)))
	mir.setUnit("another testname", "bar");
	TEST_STRING_EQUAL(mir.getUnit(1026), "bar")
END_SECTION

START_SECTION((UInt getIndex(const String& name) const))
	TEST_EQUAL(mir.getIndex("testname"), 1024)
	TEST_EQUAL(mir.getIndex("retention time"), 1025)
	TEST_EQUAL(mir.getIndex("isotopic_range"), 1)
	TEST_EQUAL(mir.getIndex("cluster_id"), 2)
  TEST_EQUAL(mir.getIndex("unregistered name"), UInt(-1))
END_SECTION

START_SECTION((String getName(UInt index) const))
	TEST_STRING_EQUAL(mir.getName(1), "isotopic_range")
	TEST_STRING_EQUAL(mir.getName(2), "cluster_id")
	TEST_STRING_EQUAL(mir.getName(3), "label")
	TEST_STRING_EQUAL(mir.getName(4), "icon")
	TEST_STRING_EQUAL(mir.getName(1024), "testname")
	TEST_STRING_EQUAL(mir.getName(1025), "retention time")
END_SECTION

START_SECTION((String getDescription(UInt index) const))
	TEST_STRING_EQUAL(mir.getDescription(1024), "this is just a test")
	TEST_STRING_EQUAL(mir.getDescription(1025), "this is just another test")
	TEST_STRING_EQUAL(mir.getDescription(1), "consecutive numbering of the peaks in an isotope pattern. 0 is the monoisotopic peak")
	TEST_STRING_EQUAL(mir.getDescription(2), "consecutive numbering of isotope clusters in a spectrum")
END_SECTION

START_SECTION((String getDescription(const String& name) const))
	TEST_STRING_EQUAL(mir.getDescription("testname"), "this is just a test")
	TEST_STRING_EQUAL(mir.getDescription("retention time"), "this is just another test")
	TEST_STRING_EQUAL(mir.getDescription("isotopic_range"), "consecutive numbering of the peaks in an isotope pattern. 0 is the monoisotopic peak")
	TEST_STRING_EQUAL(mir.getDescription("cluster_id"), "consecutive numbering of isotope clusters in a spectrum")
END_SECTION

START_SECTION((String getUnit(UInt index) const))
	TEST_STRING_EQUAL(mir.getUnit(1024), "")
	TEST_STRING_EQUAL(mir.getUnit(1025), "sec")
	TEST_STRING_EQUAL(mir.getUnit(1), "")
	TEST_STRING_EQUAL(mir.getUnit(2), "")
END_SECTION

START_SECTION((String getUnit(const String& name) const))
	TEST_STRING_EQUAL(mir.getUnit("testname"), "")
  TEST_STRING_EQUAL(mir.getUnit("retention time"), "sec")
	TEST_STRING_EQUAL(mir.getUnit("isotopic_range"), "")
	TEST_STRING_EQUAL(mir.getUnit("cluster_id"), "")
END_SECTION

START_SECTION((MetaInfoRegistry(const MetaInfoRegistry& rhs)))
	MetaInfoRegistry mir2(mir);
  TEST_EQUAL(mir2.getIndex("testname"), 1024)
  TEST_EQUAL(mir2.getIndex("retention time"), 1025)
	TEST_STRING_EQUAL(mir2.getName(1), "isotopic_range")
	TEST_STRING_EQUAL(mir2.getName(1024), "testname")
	TEST_STRING_EQUAL(mir2.getName(1025), "retention time")
	TEST_STRING_EQUAL(mir2.getDescription(1024), "this is just a test")
	TEST_STRING_EQUAL(mir2.getDescription(1025), "this is just another test")
	TEST_STRING_EQUAL(mir2.getDescription("testname"), "this is just a test")
	TEST_STRING_EQUAL(mir2.getDescription("retention time"), "this is just another test")
	TEST_STRING_EQUAL(mir2.getUnit(1024), "")
	TEST_STRING_EQUAL(mir2.getUnit(1025), "sec")
	TEST_STRING_EQUAL(mir2.getUnit("testname"), "")
	TEST_STRING_EQUAL(mir2.getUnit("retention time"), "sec")
END_SECTION

START_SECTION((MetaInfoRegistry& operator=(const MetaInfoRegistry& rhs)))
	MetaInfoRegistry mir2;
	mir2 = mir;
  TEST_EQUAL(mir2.getIndex("testname"), 1024)
  TEST_EQUAL(mir2.getIndex("retention time"), 1025)
	TEST_STRING_EQUAL(mir2.getName(1), "isotopic_range")
	TEST_STRING_EQUAL(mir2.getName(1024), "testname")
	TEST_STRING_EQUAL(mir2.getName(1025), "retention time")
	TEST_STRING_EQUAL(mir2.getDescription(1024), "this is just a test")
	TEST_STRING_EQUAL(mir2.getDescription(1025), "this is just another test")
	TEST_STRING_EQUAL(mir2.getDescription("testname"), "this is just a test")
	TEST_STRING_EQUAL(mir2.getDescription("retention time"), "this is just another test")
	TEST_STRING_EQUAL(mir2.getUnit(1024), "")
	TEST_STRING_EQUAL(mir2.getUnit(1025), "sec")
	TEST_STRING_EQUAL(mir2.getUnit("testname"), "")
	TEST_STRING_EQUAL(mir2.getUnit("retention time"), "sec")
END_SECTION

START_SECTION([EXTRA] multithreaded example)
{
  // All measurements are best of three (wall time, Linux, 8 threads)
  //
  // Serial execution of code:
  // 1e7 iterations -> 1.10 seconds
  // 1e6 iterations -> 0.31 seconds
  //
  // Parallel execution of code:
  // 1e7 iterations -> 3.41 seconds with omp critical
  // 1e7 iterations -> 4.30 seconds with std::mutex
  // 1e7 iterations -> ???  seconds with boost::shared_mutex
  //
  // 1e6 iterations -> 0.36 seconds with omp critical
  // 1e6 iterations -> 0.43 seconds with std::mutex
  // 1e6 iterations -> 5.96 seconds with boost::shared_mutex


   // 1e7 iterations with 1000 different values
  // Serial execution of code:
  // 1e7 iterations -> 6.84 seconds
  // 1e6 iterations -> 0.91 seconds
  //
  // Parallel execution of code:
  // 1e7 iterations -> 5.20 seconds with omp critical
  // 1e7 iterations -> 7.44 seconds with std::mutex
  // 1e7 iterations -> ???  seconds with boost::shared_mutex
  //
  // 1e6 iterations -> 0.55 seconds with omp critical
  // 1e6 iterations -> 0.78 seconds with std::mutex
  // 1e6 iterations -> 7.63 seconds with boost::shared_mutex

   // Note how the omp critical is 3x slower for the "same value" code than
  // serial code whereas it is slightly faster for the "1000 different values"
  // code.
  // The shared_mutex code is an order of magnitude slower (10-20x) for both
  // cases. It seems for this case, most case is spent in locking / unlocking
  // the mutex.
  int nr_iterations (1e6);
  int test = 0;
#pragma omp parallel for reduction(+: test)
  for (int k = 1; k < nr_iterations + 1; k++)
  {
#if 0
    std::string val = "newValue" + String(k%1000);
#else
    std::string val = "newValue";
#endif
    UInt index = mir.registerName(val);
    mir.setDescription(index, val);
    test += index;
  }
  TEST_EQUAL(test, 1027 * nr_iterations)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
