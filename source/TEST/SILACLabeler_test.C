// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Stephan Aiche, Fabian Kriegel, Frederic Lehnert $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/SIMULATION/DigestSimulation.h>
///////////////////////////
#include <OpenMS/SIMULATION/LABELING/SILACLabeler.h>
///////////////////////////

using namespace OpenMS;
using namespace std;




void createTestFeatureMapSimVector_(FeatureMapSimVector& feature_maps, bool add3rd)
{
  feature_maps.clear();
  
  FeatureMapSim fm1,fm2,fm3;
  ProteinHit prothit1,prothit2,prothit3,prothit4,prothit5,prothit6,prothit7,prothit8,prothit9,prothit10, prothit11, prothit12;

  // create first map
  prothit1.setSequence("AAAAAAAKAAAAA"); // 2 Fragmente AAAAAAAK und AAAAA und kommt in allen Channels vor
  prothit1.setMetaValue("description", "test sequence 1");
  prothit1.setAccession("ACC1");
  prothit1.setMetaValue("intensity", 200.0);

  prothit2.setSequence("CNARCNCNCN"); // 2 Fragmente CNAR und CNCNCN und kommt in allen Channels vor
  prothit2.setMetaValue("description", "test sequence 2");
  prothit2.setAccession("ACC2");
  prothit2.setMetaValue("intensity", 80.0);

  prothit3.setSequence("CNHAADDAAAAA"); // ungelabelt, einzelnes Fragment
  prothit3.setMetaValue("description", "test sequence 3");
  prothit3.setAccession("ACC3");
  prothit3.setMetaValue("intensity", 100.0);

  prothit12.setSequence("VNAAAAAARVNCNCNAAAA"); // Ergebniss: CNAAAAAAR(Label Medium_R) , CNCNCNAAAA (einmal kommt in allen Channels vor)
  prothit12.setMetaValue("description", "test sequence 12");
  prothit12.setAccession("ACC5");
  prothit12.setMetaValue("intensity", 115.0);

  ProteinIdentification protIdent1;
  protIdent1.insertHit(prothit1);
  protIdent1.insertHit(prothit2);
  protIdent1.insertHit(prothit3);
  protIdent1.insertHit(prothit12);
  vector<ProteinIdentification> protIdents_vec1;
  protIdents_vec1.push_back(protIdent1);
  fm1.setProteinIdentifications(protIdents_vec1);

  // create labeled map
  prothit4.setSequence("AAAAAAAKAAAAA"); // Ergbeniss: AAAAAAAK(Label Medium_K) , AAAAA ( einmal kommt in allen Channels vor)
  prothit4.setMetaValue("description", "test sequence 4");
  prothit4.setAccession("ACC4");
  prothit4.setMetaValue("intensity", 50.0);

  prothit5.setSequence("CNARCNCNCN"); // Ergebniss: CNAR(Label Medium_R) , CNCNCN (einmal kommt in allen Channels vor)
  prothit5.setMetaValue("description", "test sequence 5");
  prothit5.setAccession("ACC5");
  prothit5.setMetaValue("intensity", 100.0);

  prothit6.setSequence("LDRCEL"); // Ergbeniss : LDR(label Medium_R) , CEL (einmal kommt in channel 2 und 3 vor)
  prothit6.setMetaValue("description", "test sequence 6");
  prothit6.setAccession("ACC6");
  prothit6.setMetaValue("intensity", 120.0);

  prothit11.setSequence("VNAAAAAARVNCNCNAAAA"); // Ergebniss: CNAAAAAAR(Label Medium_R) , CNCNCNAAAA (einmal kommt in allen Channels vor)
  prothit11.setMetaValue("description", "test sequence 11");
  prothit11.setAccession("ACC5");
  prothit11.setMetaValue("intensity", 110.0);


  ProteinIdentification protIdent2;
  protIdent2.insertHit(prothit4);
  protIdent2.insertHit(prothit5);
  protIdent2.insertHit(prothit6);
  protIdent2.insertHit(prothit11);
  vector<ProteinIdentification> protIdents_vec2;
  protIdents_vec2.push_back(protIdent2);
  fm2.setProteinIdentifications(protIdents_vec2);


  feature_maps.push_back(fm1);
  feature_maps.push_back(fm2);

  if (add3rd)
  {
	  prothit7.setSequence("AAAAAAAKAAAAA"); // Ergebniss : AAAAAAAK(Label Heavy_K) , AAAAA ( einmal kommt in allen Channels vor )
    prothit7.setMetaValue("description", "test sequence 7");
	  prothit7.setAccession("ACC7");
    prothit7.setMetaValue("intensity", 30.0);

	  prothit8.setSequence("CNARCNCNCN"); // Ergebniss: CNAR(Label Heavy_R) , CNCNCN (einmal kommt in allen Channels vor)
	  prothit8.setMetaValue("description", "test sequence 8");
	  prothit8.setAccession("ACC8");
    prothit8.setMetaValue("intensity", 130.0);

	  prothit9.setSequence("LDRCEL"); //Ergebniss: LDR(label Heavy_R) , CEL (einmal kommt in channel 2 und 3 vor)
	  prothit9.setMetaValue("description", "test sequence 9");
	  prothit9.setAccession("ACC9");
    prothit9.setMetaValue("intensity", 70.0);

	  prothit10.setSequence("YCYCY"); //Ergebniss: YCYCY kommt nur in diesem Channel vor
	  prothit10.setMetaValue("description", "test sequence 10");
	  prothit10.setAccession("ACC10");
	  prothit10.setMetaValue("intensity", 80.0);

	  ProteinIdentification protIdent3;
	  protIdent3.insertHit(prothit7);
	  protIdent3.insertHit(prothit8);
	  protIdent3.insertHit(prothit9);
	  protIdent3.insertHit(prothit10);
	  vector<ProteinIdentification> protIdents_vec3;
	  protIdents_vec3.push_back(protIdent3);
    fm3.setProteinIdentifications(protIdents_vec3);
    feature_maps.push_back(fm3);
  }
  
}


void digestFeaturesMapSimVector_(FeatureMapSimVector& feature_maps)
{
  // digest here
  DigestSimulation digest_sim;
  Param p;
  p.setValue("model", "naive");
  p.setValue("model_naive:missed_cleavages", 0);
  digest_sim.setParameters(p);
  std::cout << digest_sim.getParameters() << std::endl;
  for(FeatureMapSimVector::iterator iter = feature_maps.begin() ; iter != feature_maps.end() ; ++iter)
  {
    digest_sim.digest((*iter));
  }
}



START_TEST(SILACLabeler, "$Id$")
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SILACLabeler* ptr = 0;
START_SECTION(SILACLabeler())
{
	ptr = new SILACLabeler();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~SILACLabeler())
{
	delete ptr;
}
END_SECTION

START_SECTION((void preCheck(Param &param) const ))
{
/*
  SILACLabeler labeler;
  Param p;
  p.setValue("medium_channel:modification_lysine","UniMod:695");
	  TEST_EXCEPTION_WITH_MESSAGE(Exception::ElementNotFound, labeler.preCheck(p),"The label"<< std::endl << "UniMod:695" << std::endl << "has not the modification for "<< std::endl << "K"<< std::endl << ". P please check your label. The current version of the labeler needs R or K as modification sites");
  */
 SILACLabeler labeler;
 Param p;
 labeler.preCheck(p); // Aufruf mit den Defaults
}
END_SECTION

START_SECTION((void setUpHook(FeatureMapSimVector & )))
{
  SILACLabeler labeler;
  
  FeatureMapSimVector feature_maps;
  FeatureMapSim fm1,fm2,fm3,fm4;
  
  feature_maps.push_back(fm1);
  TEST_EXCEPTION_WITH_MESSAGE(Exception::IllegalArgument, labeler.setUpHook(feature_maps),"1 channel(s) given. We currently support only 2-channel SILAC. Please provide two FASTA files!")
  feature_maps.push_back(fm2);
  labeler.setUpHook(feature_maps);
  feature_maps.push_back(fm3);
  labeler.setUpHook(feature_maps);
  feature_maps.push_back(fm4);
  TEST_EXCEPTION_WITH_MESSAGE(Exception::IllegalArgument, labeler.setUpHook(feature_maps),"4 channel(s) given. We currently support only 2-channel SILAC. Please provide two FASTA files!")
}
END_SECTION

START_SECTION((void postDigestHook(FeatureMapSimVector & )))
{
  
  FeatureMapSimVector feature_maps;
  createTestFeatureMapSimVector_(feature_maps, false);

  SILACLabeler labeler;
  labeler.setUpHook(feature_maps);
  digestFeaturesMapSimVector_(feature_maps);

  // maps are digested by now
  labeler.postDigestHook(feature_maps);
  
  TEST_EQUAL(feature_maps.size(), 1)
  ABORT_IF(feature_maps.size() != 1)

  TEST_EQUAL(feature_maps[0].size(), 12)
  ABORT_IF(feature_maps[0].size() != 12)
  
  
  TEST_EQUAL(feature_maps[0][0].getIntensity(), 250)
  TEST_EQUAL(feature_maps[0][0].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "AAAAA")

  TEST_EQUAL(feature_maps[0][1].getIntensity(), 50)
  TEST_EQUAL(feature_maps[0][1].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "AAAAAAAK(Label:2H(4))")

  TEST_EQUAL(feature_maps[0][2].getIntensity(), 200)
  TEST_EQUAL(feature_maps[0][2].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "AAAAAAAK")
  
  TEST_EQUAL(feature_maps[0][3].getIntensity(), 120)
  TEST_EQUAL(feature_maps[0][3].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "CEL")

  TEST_EQUAL(feature_maps[0][4].getIntensity(), 100)
  TEST_EQUAL(feature_maps[0][4].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "CNAR(Label:13C(6))")

  TEST_EQUAL(feature_maps[0][5].getIntensity(), 80)
  TEST_EQUAL(feature_maps[0][5].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "CNAR")

  TEST_EQUAL(feature_maps[0][6].getIntensity(), 180)
  TEST_EQUAL(feature_maps[0][6].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "CNCNCN")

  TEST_EQUAL(feature_maps[0][7].getIntensity(), 120)
  TEST_EQUAL(feature_maps[0][7].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "LDR(Label:13C(6))")

  TEST_EQUAL(feature_maps[0][8].getIntensity(), 110)
  TEST_EQUAL(feature_maps[0][8].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "VNAAAAAAR(Label:13C(6))")

  TEST_EQUAL(feature_maps[0][9].getIntensity(), 115)
  TEST_EQUAL(feature_maps[0][9].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "VNAAAAAAR")

  TEST_EQUAL(feature_maps[0][10].getIntensity(), 225)
  TEST_EQUAL(feature_maps[0][10].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "VNCNCNAAAA")

  TEST_EQUAL(feature_maps[0][11].getIntensity(), 100)
  TEST_EQUAL(feature_maps[0][11].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "CNHAADDAAAAA")

 
  createTestFeatureMapSimVector_(feature_maps, true);

  SILACLabeler three_channel_labeler;
  three_channel_labeler.setUpHook(feature_maps);
  digestFeaturesMapSimVector_(feature_maps);

  // maps are digested by now
  three_channel_labeler.postDigestHook(feature_maps);

  TEST_EQUAL(feature_maps.size(), 1)
  ABORT_IF(feature_maps.size() != 1)

  TEST_EQUAL(feature_maps[0].size(), 16)
  ABORT_IF(feature_maps[0].size() != 16)

  TEST_EQUAL(feature_maps[0][0].getIntensity(), 280)
  TEST_EQUAL(feature_maps[0][0].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "AAAAA")

  TEST_EQUAL(feature_maps[0][1].getIntensity(), 30)
  TEST_EQUAL(feature_maps[0][1].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "AAAAAAAK(Label:13C(6)15N(2))")

  TEST_EQUAL(feature_maps[0][2].getIntensity(), 50)
  TEST_EQUAL(feature_maps[0][2].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "AAAAAAAK(Label:2H(4))")

  TEST_EQUAL(feature_maps[0][3].getIntensity(), 200)
  TEST_EQUAL(feature_maps[0][3].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "AAAAAAAK")

  TEST_EQUAL(feature_maps[0][4].getIntensity(), 190)
  TEST_EQUAL(feature_maps[0][4].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "CEL")

  TEST_EQUAL(feature_maps[0][5].getIntensity(), 130)
  TEST_EQUAL(feature_maps[0][5].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "CNAR(Label:13C(6)15N(4))")

  TEST_EQUAL(feature_maps[0][6].getIntensity(), 100)
  TEST_EQUAL(feature_maps[0][6].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "CNAR(Label:13C(6))")

  TEST_EQUAL(feature_maps[0][7].getIntensity(), 80)
  TEST_EQUAL(feature_maps[0][7].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "CNAR")

  TEST_EQUAL(feature_maps[0][8].getIntensity(), 310)
  TEST_EQUAL(feature_maps[0][8].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "CNCNCN")

  TEST_EQUAL(feature_maps[0][9].getIntensity(), 70)
  TEST_EQUAL(feature_maps[0][9].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "LDR(Label:13C(6)15N(4))")

  TEST_EQUAL(feature_maps[0][10].getIntensity(), 120)
  TEST_EQUAL(feature_maps[0][10].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "LDR(Label:13C(6))")

  TEST_EQUAL(feature_maps[0][11].getIntensity(), 80)
  TEST_EQUAL(feature_maps[0][11].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "YCYCY")

  TEST_EQUAL(feature_maps[0][12].getIntensity(), 110)
  TEST_EQUAL(feature_maps[0][12].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "VNAAAAAAR(Label:13C(6))")

  TEST_EQUAL(feature_maps[0][13].getIntensity(), 115)
  TEST_EQUAL(feature_maps[0][13].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "VNAAAAAAR")

  TEST_EQUAL(feature_maps[0][14].getIntensity(), 225)
  TEST_EQUAL(feature_maps[0][14].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "VNCNCNAAAA")

  TEST_EQUAL(feature_maps[0][15].getIntensity(), 100)
  TEST_EQUAL(feature_maps[0][15].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "CNHAADDAAAAA")

}
END_SECTION

START_SECTION((void postRawMSHook(FeatureMapSimVector &)))
{
  // TODO
}
END_SECTION

// just to call the methods once
SILACLabeler dummyLabeler;
FeatureMapSimVector empty;

START_SECTION((void preCheck(Param &param) const ))
{
  //Param p;
  //dummyLabeler.preCheck(p);

  // preCheck has no content
  NOT_TESTABLE
}
END_SECTION


START_SECTION((void postRTHook(FeatureMapSimVector &)))
{
  // we do not modify the map in this step
  //dummyLabeler.postRTHook(empty);
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void postDetectabilityHook(FeatureMapSimVector &)))
{
  // we do not modify the map in this step
  //dummyLabeler.postDetectabilityHook(empty);
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void postIonizationHook(FeatureMapSimVector &)))
{
  // we do not modify the map in this step
  //dummyLabeler.postIonizationHook(empty);
  NOT_TESTABLE
}
END_SECTION

MSSimExperiment exp;
START_SECTION((void postRawTandemMSHook(FeatureMapSimVector &, MSSimExperiment &)))
{
  // we do not modify the map in this step
  //dummyLabeler.postRawTandemMSHook(empty,exp);
  NOT_TESTABLE
}
END_SECTION

START_SECTION((static BaseLabeler* create()))
{
  BaseLabeler* labeler = SILACLabeler::create();
  TEST_NOT_EQUAL(labeler, 0)
  delete labeler;
}
END_SECTION

START_SECTION((static const String getProductName()))
{
  TEST_EQUAL(SILACLabeler::getProductName(), "SILAC")
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
