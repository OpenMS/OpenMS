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
#include <OpenMS/FORMAT/MsInspectFile.h>
///////////////////////////

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

using namespace OpenMS;
using namespace std;

START_TEST(MsInspectFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MsInspectFile* ptr = nullptr;
MsInspectFile* null_ptr = nullptr;
START_SECTION(MsInspectFile())
{
	ptr = new MsInspectFile();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(virtual ~MsInspectFile())
{
	delete ptr;
}
END_SECTION

START_SECTION((template < typename FeatureMapType > void load(const String &filename, FeatureMapType &feature_map)))
{
  MsInspectFile f;
  FeatureMap fm;
  f.load(OPENMS_GET_TEST_DATA_PATH("MSInspectFile_test_1.msi"), fm);
  TEST_EQUAL(fm.size(), 2)
  ABORT_IF(fm.size()!=2)
  
  TEST_REAL_SIMILAR(fm[0].getRT(), 12.92)
  TEST_REAL_SIMILAR(fm[0].getMZ(), 501.51)
  TEST_REAL_SIMILAR(fm[0].getIntensity(), 45677)
  TEST_REAL_SIMILAR(fm[0].getOverallQuality(), 0.98)
  TEST_EQUAL(double(fm[0].getMetaValue("background")), 0.11)
  
  TEST_REAL_SIMILAR(fm[1].getRT(), 22.92)
  TEST_REAL_SIMILAR(fm[1].getMZ(), 601.51)
  TEST_REAL_SIMILAR(fm[1].getIntensity(), 245677)
  TEST_REAL_SIMILAR(fm[1].getOverallQuality(), 0.99)
  TEST_EQUAL(double(fm[1].getMetaValue("background")), 0.22)
  
}
END_SECTION

START_SECTION((template < typename SpectrumType > void store(const String &filename, const SpectrumType &spectrum) const ))
{
  MsInspectFile f;
  MSSpectrum spec;
  TEST_EXCEPTION(Exception::NotImplemented, f.store("bla", spec))
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



