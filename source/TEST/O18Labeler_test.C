// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/SIMULATION/LABELING/O18Labeler.h>
#include <OpenMS/SIMULATION/DigestSimulation.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(O18Labeler, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

O18Labeler* ptr = 0;
START_SECTION(O18Labeler())
{
	ptr = new O18Labeler();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~O18Labeler())
{
	delete ptr;
}
END_SECTION

START_SECTION((virtual ~O18Labeler()))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void preCheck(Param &param) const ))
{
  Param p;
  p.setValue("Digestion:enzyme","Trypsin","Test Param", StringList::create(""));

  O18Labeler labeler;
  labeler.preCheck(p);

  Param p_Exception;
  p_Exception.setValue("Digestion:enzyme","not-Trypsin","Test Param", StringList::create(""));
  TEST_EXCEPTION(Exception::InvalidParameter, labeler.preCheck(p_Exception))
}
END_SECTION

START_SECTION((void setUpHook(FeatureMapSimVector &)))
{
  O18Labeler labeler;

  FeatureMapSim fm1,fm2,fm3;
  FeatureMapSimVector fm_vec;

  fm_vec.push_back(fm1);
  TEST_EXCEPTION(Exception::IllegalArgument, labeler.setUpHook(fm_vec))
  fm_vec.push_back(fm2);
  labeler.setUpHook(fm_vec);
  fm_vec.push_back(fm3);
  TEST_EXCEPTION(Exception::IllegalArgument, labeler.setUpHook(fm_vec))
}
END_SECTION

START_SECTION((void postDigestHook(FeatureMapSimVector &)))
{
  FeatureMapSimVector feature_maps;

  // first feature map TVQMENQFVAFVDK,ACHKKKKHHACAC,AAAAHTKLRTTIPPEFG,RYCNHKTUIKL
  FeatureMapSim fm1,fm2;
  ProteinHit prothit1,prothit2,prothit3,prothit4,prothit5;

  // create first map
  prothit1.setSequence("AAAAAAAKHHHHHHHHHHH");
  prothit1.setMetaValue("description", "test sequence 1");
  prothit1.setAccession("ACC1");
  prothit1.setMetaValue("intensity", 200.0);

  prothit2.setSequence("CNHAAAAAAAAA");
  prothit2.setMetaValue("description", "test sequence 2");
  prothit2.setAccession("ACC2");
  prothit2.setMetaValue("intensity", 100.0);

  ProteinIdentification protIdent1;
  protIdent1.insertHit(prothit1);
  protIdent1.insertHit(prothit2);
  vector<ProteinIdentification> protIdents_vec1;
  protIdents_vec1.push_back(protIdent1);
  fm1.setProteinIdentifications(protIdents_vec1);

  // create second map
  prothit3.setSequence("AAAAAAAKHHHHHHHHHHH"); // same as protein 1 from first map
  prothit3.setMetaValue("description", "test sequence 3");
  prothit3.setAccession("ACC3");
  prothit3.setMetaValue("intensity", 50.0);

  prothit5.setSequence("CNHAAAAAAAAA");
  prothit5.setMetaValue("description", "test sequence 5");
  prothit5.setAccession("ACC5");
  prothit5.setMetaValue("intensity", 100.0);

  ProteinIdentification protIdent2;
  protIdent2.insertHit(prothit3);
  protIdent2.insertHit(prothit5);
  vector<ProteinIdentification> protIdents_vec2;
  protIdents_vec2.push_back(protIdent2);
  fm2.setProteinIdentifications(protIdents_vec2);

  feature_maps.push_back(fm1);
  feature_maps.push_back(fm2);

  // digest here
  DigestSimulation digest_sim;
  Param p;
  p.setValue("model", "naive");
  p.setValue("model_naive|missed_cleavages", 0);
  digest_sim.setParameters(p);
  std::cout << digest_sim.getParameters() << std::endl;
  for(FeatureMapSimVector::iterator iter = feature_maps.begin() ; iter != feature_maps.end() ; ++iter)
  {
    digest_sim.digest((*iter));
  }

  // maps are digested by now
  O18Labeler labeler;
  labeler.postDigestHook(feature_maps);

  TEST_EQUAL(feature_maps.size(), 1)
  ABORT_IF(feature_maps.size() != 1)

  TEST_EQUAL(feature_maps[0].size(), 5)
  ABORT_IF(feature_maps[0].size() != 5)
  TEST_EQUAL(feature_maps[0][0].getIntensity(), 38)
  TEST_EQUAL(feature_maps[0][0].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "AAAAAAAK(Label:18O(2))")
  TEST_EQUAL(feature_maps[0][1].getIntensity(), 150)
  TEST_EQUAL(feature_maps[0][1].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "AAAAAAAK")
  TEST_EQUAL(feature_maps[0][2].getIntensity(), 188)
  TEST_EQUAL(feature_maps[0][2].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "AAAAAAAKHHHHHHHHHHH")
  TEST_EQUAL(feature_maps[0][3].getIntensity(), 200)
  TEST_EQUAL(feature_maps[0][3].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "CNHAAAAAAAAA")
  TEST_EQUAL(feature_maps[0][4].getIntensity(), 188)
  TEST_EQUAL(feature_maps[0][4].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "HHHHHHHHHHH")

}
END_SECTION

START_SECTION((void postRTHook(FeatureMapSimVector &)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void postDetectabilityHook(FeatureMapSimVector &)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void postIonizationHook(FeatureMapSimVector &)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void postRawMSHook(FeatureMapSimVector &)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void postRawTandemMSHook(FeatureMapSimVector &, MSSimExperiment &)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((static BaseLabeler* create()))
{
  BaseLabeler* labeler = O18Labeler::create();
  TEST_NOT_EQUAL(labeler, 0)
  delete labeler;
}
END_SECTION

START_SECTION((static const String getProductName()))
{
  TEST_EQUAL(O18Labeler::getProductName(), "o18")
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



