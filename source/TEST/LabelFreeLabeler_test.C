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
#include <OpenMS/SIMULATION/LABELING/LabelFreeLabeler.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(LabelFreeLabeler, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

LabelFreeLabeler* ptr = 0;
START_SECTION(LabelFreeLabeler())
{
	ptr = new LabelFreeLabeler();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~LabelFreeLabeler())
{
	delete ptr;
}
END_SECTION


START_SECTION((void preCheck(Param &param) const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void setUpHook(FeatureMapSimVector &)))
{
  FeatureMapSimVector feature_maps;

  // first feature map TVQMENQFVAFVDK,ACHKKKKHHACAC,AAAAHTKLRTTIPPEFG,RYCNHKTUIKL
  FeatureMapSim fm1,fm2;
  ProteinHit prothit1,prothit2,prothit3,prothit4,prothit5;

  // create first map
  prothit1.setSequence("TVQMENQFVAFVDK");
  prothit1.setMetaValue("description", "test sequence 1");
  prothit1.setAccession("ACC1");
  prothit1.setMetaValue("intensity", 100.0);

  prothit2.setSequence("ACHKKKKHHACAC");
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
  prothit3.setSequence("TVQMENQFVAFVDK"); // same as protein 1 from first map
  prothit3.setMetaValue("description", "test sequence 3");
  prothit3.setAccession("ACC3");
  prothit3.setMetaValue("intensity", 10.0);

  prothit4.setSequence("AAAAHTKLRTTIPPEFG");
  prothit4.setMetaValue("description", "test sequence 4");
  prothit4.setAccession("ACC4");
  prothit4.setMetaValue("intensity", 100.0);

  prothit5.setSequence("RYCNHKTUIKL");
  prothit5.setMetaValue("description", "test sequence 5");
  prothit5.setAccession("ACC5");
  prothit5.setMetaValue("intensity", 100.0);

  ProteinIdentification protIdent2;
  protIdent2.insertHit(prothit3);
  protIdent2.insertHit(prothit4);
  protIdent2.insertHit(prothit5);
  vector<ProteinIdentification> protIdents_vec2;
  protIdents_vec2.push_back(protIdent2);
  fm2.setProteinIdentifications(protIdents_vec2);

  feature_maps.push_back(fm1);
  feature_maps.push_back(fm2);

  LabelFreeLabeler labeler;
  labeler.setUpHook(feature_maps);

  TEST_EQUAL(feature_maps.size(), 1)
  ABORT_IF(feature_maps.size() != 1)

  TEST_EQUAL(feature_maps[0].getProteinIdentifications().size(), 1)
  TEST_EQUAL(feature_maps[0].getProteinIdentifications()[0].getHits().size(), 4)
  ABORT_IF(feature_maps[0].getProteinIdentifications()[0].getHits().size() != 4)


  TEST_EQUAL( feature_maps[0].getProteinIdentifications()[0].getHits()[0].getSequence(), "AAAAHTKLRTTIPPEFG")
  TEST_REAL_SIMILAR( feature_maps[0].getProteinIdentifications()[0].getHits()[0].getMetaValue("intensity"), 100.0)
  TEST_EQUAL( feature_maps[0].getProteinIdentifications()[0].getHits()[1].getSequence(), "ACHKKKKHHACAC")
  TEST_REAL_SIMILAR( feature_maps[0].getProteinIdentifications()[0].getHits()[1].getMetaValue("intensity"), 100.0)
  TEST_EQUAL( feature_maps[0].getProteinIdentifications()[0].getHits()[2].getSequence(), "RYCNHKTUIKL")
  TEST_REAL_SIMILAR( feature_maps[0].getProteinIdentifications()[0].getHits()[2].getMetaValue("intensity"), 100.0)
  TEST_EQUAL( feature_maps[0].getProteinIdentifications()[0].getHits()[3].getSequence(), "TVQMENQFVAFVDK")
  TEST_REAL_SIMILAR( feature_maps[0].getProteinIdentifications()[0].getHits()[3].getMetaValue("intensity"), 110.0) // merge happened
  TEST_EQUAL( feature_maps[0].getProteinIdentifications()[0].getHits()[3].getAccession(), "ACC1")
}
END_SECTION

START_SECTION((void postDigestHook(FeatureMapSimVector &)))
{
  NOT_TESTABLE
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
  BaseLabeler* labeler = LabelFreeLabeler::create();
  TEST_NOT_EQUAL(labeler, 0)
  delete labeler;
}
END_SECTION

START_SECTION((static const String getProductName()))
{
  TEST_EQUAL(LabelFreeLabeler::getProductName(), "labelfree")
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
