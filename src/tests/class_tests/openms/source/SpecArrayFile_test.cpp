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
#include <OpenMS/FORMAT/SpecArrayFile.h>
///////////////////////////

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

using namespace OpenMS;
using namespace std;

START_TEST(SpecArrayFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SpecArrayFile* ptr = nullptr;
SpecArrayFile* null_ptr = nullptr;
START_SECTION(SpecArrayFile())
{
	ptr = new SpecArrayFile();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(virtual ~SpecArrayFile())
{
	delete ptr;
}
END_SECTION


START_SECTION((template < typename FeatureMapType > void load(const String &filename, FeatureMapType &feature_map)))
{
  SpecArrayFile f;
  FeatureMap fm;
  f.load(OPENMS_GET_TEST_DATA_PATH("SpecArrayFile_test_1.peplist"), fm);
  TEST_EQUAL(fm.size(),2)
  ABORT_IF(fm.size()!=2)
  TEST_EQUAL(fm[0].getRT(), 60.1*60)
  TEST_REAL_SIMILAR(fm[0].getMZ(), 500.1)
  TEST_EQUAL(fm[0].getIntensity(), 4343534)
  TEST_EQUAL(fm[0].getCharge(), 5)
  TEST_EQUAL(double(fm[0].getMetaValue("s/n")), 3.2)
  TEST_EQUAL(fm[1].getRT(),  40.1*60)
  TEST_REAL_SIMILAR(fm[1].getMZ(), 700.1	)
  TEST_EQUAL(fm[1].getIntensity(), 222432)
  TEST_EQUAL(fm[1].getCharge(), 3)
  TEST_EQUAL(double(fm[1].getMetaValue("s/n")), 2.2)

  

  TEST_EXCEPTION(Exception::ParseError, f.load(OPENMS_GET_TEST_DATA_PATH("SpecArrayFile_test_2.peplist"), fm));
  
  TEST_EXCEPTION(Exception::FileNotFound, f.load(OPENMS_GET_TEST_DATA_PATH("SpecArrayFile_test_2_doesnotexist.peplist"), fm));
}
END_SECTION

START_SECTION((template < typename SpectrumType > void store(const String &filename, const SpectrumType &spectrum) const ))
{
  SpecArrayFile f;
  MSSpectrum spec;
  TEST_EXCEPTION(Exception::NotImplemented, f.store("bla", spec))
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



