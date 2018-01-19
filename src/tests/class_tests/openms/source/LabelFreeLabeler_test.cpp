// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/SIMULATION/LABELING/LabelFreeLabeler.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(LabelFreeLabeler, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

LabelFreeLabeler* ptr = nullptr;
LabelFreeLabeler* nullPointer = nullptr;
BaseLabeler*      base_nullPointer = nullptr;
START_SECTION(LabelFreeLabeler())
{
	ptr = new LabelFreeLabeler();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~LabelFreeLabeler())
{
	delete ptr;
}
END_SECTION

START_SECTION((void setUpHook(SimTypes::FeatureMapSimVector &)))
{
  SimTypes::FeatureMapSimVector feature_maps;

  // first feature map TVQMENQFVAFVDK,ACHKKKKHHACAC,AAAAHTKLRTTIPPEFG,RYCNHKTUIKL
  SimTypes::FeatureMapSim fm1,fm2;
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

// just to call the methods once
LabelFreeLabeler dummyLabeler;
SimTypes::FeatureMapSimVector empty;

START_SECTION((void preCheck(Param &param) const ))
{
  Param p;
  dummyLabeler.preCheck(p);

  // preCheck has no content
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void postDigestHook(SimTypes::FeatureMapSimVector &)))
{
  // we do not modify the map in this step
  dummyLabeler.postDigestHook(empty);
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void postRTHook(SimTypes::FeatureMapSimVector &)))
{
  // we do not modify the map in this step
  dummyLabeler.postRTHook(empty);
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void postDetectabilityHook(SimTypes::FeatureMapSimVector &)))
{
  // we do not modify the map in this step
  dummyLabeler.postDetectabilityHook(empty);
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void postIonizationHook(SimTypes::FeatureMapSimVector &)))
{
  // we do not modify the map in this step
  dummyLabeler.postIonizationHook(empty);
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void postRawMSHook(SimTypes::FeatureMapSimVector &)))
{
  // we do not modify the map in this step
  dummyLabeler.postRawMSHook(empty);
  NOT_TESTABLE
}
END_SECTION

SimTypes::MSSimExperiment exp;
START_SECTION((void postRawTandemMSHook(SimTypes::FeatureMapSimVector &, SimTypes::MSSimExperiment &)))
{
  // we do not modify the map in this step
  dummyLabeler.postRawTandemMSHook(empty,exp);
  NOT_TESTABLE
}
END_SECTION

START_SECTION((static BaseLabeler* create()))
{
  BaseLabeler* labeler = LabelFreeLabeler::create();
  TEST_NOT_EQUAL(labeler, base_nullPointer)
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
