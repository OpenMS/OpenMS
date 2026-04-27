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

START_SECTION(double computeGRAVY(const AASequence& seq))
{
  AASequence seq("ACDE");
  HydrophobicityProfile profile;
  double gravy = profile.computeGRAVY(seq);
  TEST_REAL_SIMILAR(gravy, -0.675); 
  AASequence seq_2;
  TEST_EXCEPTION(Exception::InvalidValue,profile.computeGRAVY(seq_2)); 
  AASequence seq_3("XXX");
  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidValue,profile.computeGRAVY(seq_3),
  "the value 'X' was used but is not valid; No hydrophobicity value known for this residue");
}
END_SECTION

START_SECTION(std::vector<double> computeProfile(const AASequence& seq, const HydrophobicityScaleMethod scale) const)
{
  HydrophobicityProfile profile;
  AASequence seq("ACDE");
  AASequence seq_2;
  AASequence seq_3("XXX");
  std::vector<double> vec = profile.computeProfile(seq,HydrophobicityScaleMethod::KYTE_DOOLITTLE);
  TEST_REAL_SIMILAR(vec[0],1.8);
  TEST_REAL_SIMILAR(vec[1],2.5);
  TEST_REAL_SIMILAR(vec[2],-3.5);
  TEST_REAL_SIMILAR(vec[3],-3.5); 
  TEST_REAL_SIMILAR(profile.computeProfile(seq,HydrophobicityScaleMethod::EISENBERG)[0],0.62);
  TEST_EXCEPTION(Exception::InvalidValue,profile.computeProfile(seq_3,HydrophobicityScaleMethod::EISENBERG))
}
END_SECTION

START_SECTION(std::vector<double> computeWindowedProfile(const AASequence& seq, Size window_size, const HydrophobicityScaleMethod scale) const)
{
  HydrophobicityProfile profile;
  AASequence seq("ACDEF");
  AASequence seq_2;
  std::vector<double> vec = profile.computeWindowedProfile(seq,3);
  TEST_REAL_SIMILAR(vec[0],0.266666666666667);
  TEST_REAL_SIMILAR(vec[1],-1.5);
  TEST_REAL_SIMILAR(vec[2],-1.4);
  TEST_REAL_SIMILAR(profile.computeWindowedProfile(seq,6)[0],0.02);
  TEST_EXCEPTION(Exception::InvalidSize,profile.computeWindowedProfile(seq,0));
  TEST_EXCEPTION(Exception::InvalidValue,profile.computeWindowedProfile(seq_2,3));
}
END_SECTION

START_SECTION(std::vector<double> computeHydrophobicMoment(const AASequence& seq, Size window_size, double angle) const)
{
  HydrophobicityProfile profile;
  AASequence seq_1("ACDEF");
  AASequence seq_2;
  std::vector<double> vec = profile.computeHydrophobicMoment(seq_1,3,100);
  TEST_REAL_SIMILAR(vec[0],0.511576803);
  TEST_REAL_SIMILAR(vec[1],0.435170599);
  TEST_REAL_SIMILAR(vec[2],0.734926405);
  TEST_EXCEPTION(Exception::InvalidSize,profile.computeHydrophobicMoment(seq_1,0));
  TEST_EXCEPTION(Exception::InvalidValue,profile.computeHydrophobicMoment(seq_2,3));
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST