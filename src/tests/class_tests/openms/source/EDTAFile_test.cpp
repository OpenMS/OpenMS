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
#include <OpenMS/FORMAT/EDTAFile.h>
///////////////////////////

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>

using namespace OpenMS;
using namespace std;

START_TEST(EDTAFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

EDTAFile* ptr = nullptr;
EDTAFile* null_ptr = nullptr;
START_SECTION(EDTAFile())
{
	ptr = new EDTAFile();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(virtual ~EDTAFile())
{
	delete ptr;
}
END_SECTION


START_SECTION(void load(const String &filename, ConsensusMap &consensus_map))
{
  EDTAFile f;
  ConsensusMap fm;
  f.load(OPENMS_GET_TEST_DATA_PATH("EDTAFile_test_1.edta"), fm);
  TEST_EQUAL(fm.size(),2)
  ABORT_IF(fm.size()!=2)
  TEST_EQUAL(fm[0].getRT(), 321)
  TEST_EQUAL(fm[0].getMZ(), 405.233)
  TEST_EQUAL(fm[0].getIntensity(), 24543534)
  TEST_EQUAL(fm[0].getCharge(), 2)
  TEST_EQUAL(String(fm[0].getMetaValue("mymeta")), String("lala"))
  TEST_EQUAL(fm[1].getRT(), 322)
  TEST_EQUAL(fm[1].getMZ(), 406.207)
  TEST_EQUAL(fm[1].getIntensity(), 4343344)
  TEST_EQUAL(fm[1].getCharge(), 3)
  TEST_EQUAL(String(fm[1].getMetaValue("mymeta")), String("blubb"))

  
  f.load(OPENMS_GET_TEST_DATA_PATH("EDTAFile_test_3.edta"), fm);
  TEST_EQUAL(fm.size(),3)

  TEST_EXCEPTION(Exception::ParseError, f.load(OPENMS_GET_TEST_DATA_PATH("EDTAFile_test_2.edta"), fm));
  
  TEST_EXCEPTION(Exception::FileNotFound, f.load(OPENMS_GET_TEST_DATA_PATH("EDTAFile_test_3_doesnotexist.edta"), fm));
      
}
END_SECTION

START_SECTION((void store(const String& filename, const ConsensusMap& map) const))
{
  EDTAFile f;
  ConsensusMap cm;
  f.load(OPENMS_GET_TEST_DATA_PATH("EDTAFile_test_4.edta"), cm);
  
  String outfile;
  NEW_TMP_FILE(outfile)
  f.store(outfile, cm);

  ConsensusMap cm2;
  f.load(outfile, cm2);

  TEST_EQUAL(cm.size(),cm2.size())
  ABORT_IF(cm.size()!=cm2.size())
  for (Size i=0; i< cm.size(); ++i)
  {
    TEST_REAL_SIMILAR(cm[i].getRT(), cm2[i].getRT())
    TEST_REAL_SIMILAR(cm[i].getMZ(), cm2[i].getMZ())
    TEST_REAL_SIMILAR(cm[i].getIntensity(), cm2[i].getIntensity())
    TEST_EQUAL(cm[i].getCharge(), cm2[i].getCharge())
    TEST_EQUAL(cm[i].getFeatures().size(), cm2[i].getFeatures().size())
    // cannot test for metavalues, since they are not written to EDTA (yet)
  }

}
END_SECTION

 START_SECTION((void store(const String& filename, const FeatureMap& map) const))
{
  FeatureMap fm;
  FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("EDTAFile_test_out_1.featureXML"), fm);

  EDTAFile f;

  String outfile;
  NEW_TMP_FILE(outfile)

  f.store(outfile, fm);

  ConsensusMap cm;
  f.load(outfile, cm);


  TEST_EQUAL(fm.size(), cm.size());
  ABORT_IF(fm.size() != cm.size());
  for (Size i=0; i< fm.size(); ++i)
  {
    TEST_REAL_SIMILAR(fm[i].getRT(), cm[i].getRT())
    TEST_REAL_SIMILAR(fm[i].getMZ(), cm[i].getMZ())
    TEST_REAL_SIMILAR(fm[i].getIntensity(), cm[i].getIntensity())
  }

}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



