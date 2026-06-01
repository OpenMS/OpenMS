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
#include <OpenMS/FORMAT/KroenikFile.h>
///////////////////////////

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

using namespace OpenMS;
using namespace std;

START_TEST(KroenikFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

KroenikFile* ptr = nullptr;
KroenikFile* null_ptr = nullptr;
START_SECTION(KroenikFile())
{
	ptr = new KroenikFile();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(virtual ~KroenikFile())
{
	delete ptr;
}
END_SECTION


START_SECTION((template < typename FeatureMapType > void load(const String &filename, FeatureMapType &feature_map)))
{
  KroenikFile f;
  FeatureMap fm;
  f.load(OPENMS_GET_TEST_DATA_PATH("KroenikFile_test_1.krf"), fm);
  TEST_EQUAL(fm.size(),3)
  ABORT_IF(fm.size()!=3)
  TEST_EQUAL(fm[0].getRT(), 63.2)
  TEST_REAL_SIMILAR(fm[0].getMZ(), 1002.11)
  TEST_EQUAL(fm[0].getIntensity(), 999999)
  TEST_EQUAL(fm[0].getCharge(), 1)
  TEST_EQUAL(String(fm[0].getMetaValue("AveragineModifications")), String("Carbamido"))
  TEST_EQUAL(fm[1].getRT(),  62.2)
  TEST_REAL_SIMILAR(fm[1].getMZ(), 252.057	)
  TEST_EQUAL(fm[1].getIntensity(), 9999)
  TEST_EQUAL(fm[1].getCharge(), 2)
  TEST_EQUAL(String(fm[1].getMetaValue("AveragineModifications")), String("Carbamido2"))

  

  TEST_EXCEPTION(Exception::ParseError, f.load(OPENMS_GET_TEST_DATA_PATH("KroenikFile_test_2.krf"), fm));
  
  TEST_EXCEPTION(Exception::FileNotFound, f.load(OPENMS_GET_TEST_DATA_PATH("KroenikFile_test_2_doesnotexist.edta"), fm));
}
END_SECTION

START_SECTION((template < typename SpectrumType > void store(const String &filename, const SpectrumType &spectrum) const ))
{
  KroenikFile f;
  MSSpectrum spec;
  TEST_EXCEPTION(Exception::NotImplemented, f.store("bla", spec))
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



