// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Markus Apel, Nora Heese $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/HydrophobicityProfile.h>


///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(HydrophobicityProfile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

HydrophobicityProfile* ptr = nullptr;
HydrophobicityProfile* null_ptr = nullptr;
START_SECTION(HydrophobicityProfile())
{
  ptr = new HydrophobicityProfile();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~HydrophobicityProfile())
{
  delete ptr;
}
END_SECTION

START_SECTION(int testfunktion())
{
    TEST_EQUAL(ptr->testfunktion(),1);
}
END_SECTION

START_SECTION(double computeGRAVY(const AASequence& seq))
{
  AASequence seq("ACDE");
  HydrophobicityProfile profile;
  double gravy = profile.computeGRAVY(seq);
  TEST_REAL_SIMILAR(gravy, -0.675); 
  AASequence seq_2;
  double gravy_2 = profile.computeGRAVY(seq_2);
  TEST_REAL_SIMILAR(gravy_2, 0); 
  AASequence seq_3("XXX");
  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidValue,profile.computeGRAVY(seq_3),
  "the value 'X' was used but is not valid; No hydrophobicity value known for this residue");
}
END_SECTION

START_SECTION(std::vector<double> computeProfile(const AASequence& seq, const HydrophobicityScaleMethod scale))
{
  HydrophobicityProfile profile;
  AASequence seq("ACDE");
  AASequence seq_2;
  AASequence seq_3("XXX");
  TEST_REAL_SIMILAR(profile.computeProfile(seq,HydrophobicityScaleMethod::KYTE_DOOLITTLE)[0],1.8);
  TEST_REAL_SIMILAR(profile.computeProfile(seq,HydrophobicityScaleMethod::KYTE_DOOLITTLE)[1],2.5);
  TEST_REAL_SIMILAR(profile.computeProfile(seq,HydrophobicityScaleMethod::KYTE_DOOLITTLE)[2],-3.5);
  TEST_REAL_SIMILAR(profile.computeProfile(seq,HydrophobicityScaleMethod::KYTE_DOOLITTLE)[3],-3.5); 
  TEST_REAL_SIMILAR(profile.computeProfile(seq,HydrophobicityScaleMethod::EISENBERG)[0],0.62);
  TEST_EXCEPTION(Exception::InvalidValue,profile.computeProfile(seq_3,HydrophobicityScaleMethod::EISENBERG))
}
END_SECTION

START_SECTION(computeWindowedProfile)
{
  HydrophobicityProfile profile;
  AASequence seq("ACDEF");
  TEST_REAL_SIMILAR(profile.computeWindowedProfile(seq,3,HydrophobicityScaleMethod::KYTE_DOOLITTLE)[0],0.266666666666667);
  TEST_REAL_SIMILAR(profile.computeWindowedProfile(seq,3,HydrophobicityScaleMethod::KYTE_DOOLITTLE)[1],-1.5);
  TEST_REAL_SIMILAR(profile.computeWindowedProfile(seq,3,HydrophobicityScaleMethod::KYTE_DOOLITTLE)[2],-1.4);
  TEST_EXCEPTION(Exception::InvalidSize,profile.computeWindowedProfile(seq,6,HydrophobicityScaleMethod::KYTE_DOOLITTLE))
}
END_SECTION

START_SECTION(computeHydrophobicMoment)
{
  HydrophobicityProfile profile;
  AASequence seq("A");
  TEST_REAL_SIMILAR(profile.computeHydrophobicMoment(seq,1,100)[0],0.62);
  AASequence seq_2("ACDEF");
  TEST_REAL_SIMILAR(profile.computeHydrophobicMoment(seq_2,3,100)[0],0.512);
  TEST_REAL_SIMILAR(profile.computeHydrophobicMoment(seq_2,3,100)[1],0.435);
  TEST_REAL_SIMILAR(profile.computeHydrophobicMoment(seq_2,3,100)[2],0.735);
}
END_SECTION

// Add more test sections for each method...

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST