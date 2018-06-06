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
#include <OpenMS/SIMULATION/LABELING/O18Labeler.h>
#include <OpenMS/SIMULATION/DigestSimulation.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

void createTestFeatureMapSimVector_(SimTypes::FeatureMapSimVector& feature_maps)
{
  feature_maps.clear();

  // first feature map TVQMENQFVAFVDK,ACHKKKKHHACAC,AAAAHTKLRTTIPPEFG,RYCNHKTUIKL
  SimTypes::FeatureMapSim fm1,fm2;
  ProteinHit prothit1,prothit2,prothit3,prothit4,prothit5,prothit6;

  // create first map
  prothit1.setSequence("AAAAAAAKHHHHHHHHHHH");
  prothit1.setMetaValue("description", "test sequence 1");
  prothit1.setAccession("ACC1");
  prothit1.setMetaValue("intensity", 200.0);

  prothit2.setSequence("CNHAAAAAAAAA");
  prothit2.setMetaValue("description", "test sequence 2");
  prothit2.setAccession("ACC2");
  prothit2.setMetaValue("intensity", 100.0);

  prothit3.setSequence("LDCELR");
  prothit3.setMetaValue("description", "test sequence 3");
  prothit3.setAccession("ACC3");
  prothit3.setMetaValue("intensity", 100.0);

  ProteinIdentification protIdent1;
  protIdent1.insertHit(prothit1);
  protIdent1.insertHit(prothit2);
  protIdent1.insertHit(prothit3);
  vector<ProteinIdentification> protIdents_vec1;
  protIdents_vec1.push_back(protIdent1);
  fm1.setProteinIdentifications(protIdents_vec1);

  // create labeled map
  prothit4.setSequence("AAAAAAAKHHHHHHHHHHH"); // same as protein 1 from first map
  prothit4.setMetaValue("description", "test sequence 4");
  prothit4.setAccession("ACC4");
  prothit4.setMetaValue("intensity", 50.0);

  prothit5.setSequence("CNHAAAAAAAAA");
  prothit5.setMetaValue("description", "test sequence 5");
  prothit5.setAccession("ACC5");
  prothit5.setMetaValue("intensity", 100.0);

  prothit6.setSequence("CNHAADDAAAAA");
  prothit6.setMetaValue("description", "test sequence 6");
  prothit6.setAccession("ACC6");
  prothit6.setMetaValue("intensity", 120.0);

  ProteinIdentification protIdent2;
  protIdent2.insertHit(prothit4);
  protIdent2.insertHit(prothit5);
  protIdent2.insertHit(prothit6);
  vector<ProteinIdentification> protIdents_vec2;
  protIdents_vec2.push_back(protIdent2);
  fm2.setProteinIdentifications(protIdents_vec2);

  feature_maps.push_back(fm1);
  feature_maps.push_back(fm2);
}

void digestFeaturesMapSimVector_(SimTypes::FeatureMapSimVector& feature_maps)
{
  // digest here
  DigestSimulation digest_sim;
  Param p;
  p.setValue("model", "naive");
  p.setValue("model_naive:missed_cleavages", 0);
  digest_sim.setParameters(p);
  std::cout << digest_sim.getParameters() << std::endl;
  for(SimTypes::FeatureMapSimVector::iterator iter = feature_maps.begin() ; iter != feature_maps.end() ; ++iter)
  {
    digest_sim.digest((*iter));
  }
}


START_TEST(O18Labeler, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

O18Labeler*  ptr = nullptr;
O18Labeler*  nullPointer = nullptr;
BaseLabeler* base_nullPointer = nullptr;

START_SECTION(O18Labeler())
{
  ptr = new O18Labeler();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~O18Labeler())
{
  delete ptr;
}
END_SECTION

START_SECTION((void preCheck(Param &param) const ))
{
  Param p;
  p.setValue("Digestion:enzyme","Trypsin","Test Param", ListUtils::create<String>(""));

  O18Labeler labeler;
  labeler.preCheck(p);

  Param p_Exception;
  p_Exception.setValue("Digestion:enzyme","not-Trypsin","Test Param", ListUtils::create<String>(""));
  TEST_EXCEPTION(Exception::InvalidParameter, labeler.preCheck(p_Exception))
}
END_SECTION

START_SECTION((void setUpHook(SimTypes::FeatureMapSimVector &)))
{
  O18Labeler labeler;

  SimTypes::FeatureMapSim fm1,fm2,fm3;
  SimTypes::FeatureMapSimVector fm_vec;

  fm_vec.push_back(fm1);
  TEST_EXCEPTION_WITH_MESSAGE(Exception::IllegalArgument, labeler.setUpHook(fm_vec), "1 channel(s) given. 18O Labeling only works with 2 channels. Please provide two FASTA files!")
  fm_vec.push_back(fm2);
  labeler.setUpHook(fm_vec);
  fm_vec.push_back(fm3);
  TEST_EXCEPTION_WITH_MESSAGE(Exception::IllegalArgument, labeler.setUpHook(fm_vec), "3 channel(s) given. 18O Labeling only works with 2 channels. Please provide two FASTA files!")
}
END_SECTION

START_SECTION((void postDigestHook(SimTypes::FeatureMapSimVector &)))
{
  SimTypes::FeatureMapSimVector feature_maps;

  createTestFeatureMapSimVector_(feature_maps);
  digestFeaturesMapSimVector_(feature_maps);

  // maps are digested by now
  O18Labeler labeler;
  labeler.postDigestHook(feature_maps);

  TEST_EQUAL(feature_maps.size(), 1)
  ABORT_IF(feature_maps.size() != 1)

  TEST_EQUAL(feature_maps[0].size(), 6)
  ABORT_IF(feature_maps[0].size() != 6)
  TEST_EQUAL(feature_maps[0][0].getIntensity(), 50)
  TEST_EQUAL(feature_maps[0][0].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "AAAAAAAK.(Label:18O(2))")

  TEST_EQUAL(feature_maps[0][1].getIntensity(), 200)
  TEST_EQUAL(feature_maps[0][1].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "AAAAAAAK")

  TEST_EQUAL(feature_maps[0][2].getIntensity(), 250)
  TEST_EQUAL(feature_maps[0][2].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "HHHHHHHHHHH")

  TEST_EQUAL(feature_maps[0][3].getIntensity(), 200)
  TEST_EQUAL(feature_maps[0][3].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "CNHAAAAAAAAA")

  TEST_EQUAL(feature_maps[0][4].getIntensity(), 120)
  TEST_EQUAL(feature_maps[0][4].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "CNHAADDAAAAA")

  TEST_EQUAL(feature_maps[0][5].getIntensity(), 100)
  TEST_EQUAL(feature_maps[0][5].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "LDCELR")

  // Test ConsensusMap association
  ConsensusMap cm = labeler.getConsensus();
  TEST_EQUAL(cm.size(), 1)
  ABORT_IF(cm.size() != 1)
  TEST_EQUAL(cm[0].getFeatures().size(),2)

  ConsensusFeature::HandleSetType::const_iterator fhIt = cm[0].getFeatures().begin();
  TEST_EQUAL(feature_maps[0][1].getUniqueId(), fhIt->getUniqueId())
  ++fhIt;
  TEST_EQUAL(feature_maps[0][0].getUniqueId(), fhIt->getUniqueId())

  // now test the incomplete variant
  createTestFeatureMapSimVector_(feature_maps);
  digestFeaturesMapSimVector_(feature_maps);

  O18Labeler incomplete_labeler;
  Param p;
  p.setValue("labeling_efficiency", 0.7);
  incomplete_labeler.setParameters(p);

  incomplete_labeler.postDigestHook(feature_maps);

  TEST_EQUAL(feature_maps.size(), 1)
  ABORT_IF(feature_maps.size() != 1)

  TEST_EQUAL(feature_maps[0].size(), 7)
  ABORT_IF(feature_maps[0].size() != 7)

  TEST_EQUAL(feature_maps[0][0].getIntensity(), 24.5)
  TEST_EQUAL(feature_maps[0][0].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "AAAAAAAK.(Label:18O(2))")

  TEST_EQUAL(feature_maps[0][1].getIntensity(), 21)
  TEST_EQUAL(feature_maps[0][1].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "AAAAAAAK.(Label:18O(1))")

  TEST_EQUAL(feature_maps[0][2].getIntensity(), 204.5)
  TEST_EQUAL(feature_maps[0][2].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "AAAAAAAK")

  TEST_EQUAL(feature_maps[0][3].getIntensity(), 250)
  TEST_EQUAL(feature_maps[0][3].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "HHHHHHHHHHH")

  TEST_EQUAL(feature_maps[0][4].getIntensity(), 200)
  TEST_EQUAL(feature_maps[0][4].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "CNHAAAAAAAAA")

  TEST_EQUAL(feature_maps[0][5].getIntensity(), 120)
  TEST_EQUAL(feature_maps[0][5].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "CNHAADDAAAAA")

  TEST_EQUAL(feature_maps[0][6].getIntensity(), 100)
  TEST_EQUAL(feature_maps[0][6].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "LDCELR")

  // Test ConsensusMap association
  ConsensusMap incomplete_cm = incomplete_labeler.getConsensus();
  TEST_EQUAL(incomplete_cm.size(), 1)
  ABORT_IF(incomplete_cm.size() != 1)
  TEST_EQUAL(incomplete_cm[0].getFeatures().size(),3)

  ConsensusFeature::HandleSetType::const_iterator incomplete_fhIt = incomplete_cm[0].getFeatures().begin();
  TEST_EQUAL(feature_maps[0][2].getUniqueId(), incomplete_fhIt->getUniqueId())
  ++incomplete_fhIt;
  TEST_EQUAL(feature_maps[0][1].getUniqueId(), incomplete_fhIt->getUniqueId())
  ++incomplete_fhIt;
  TEST_EQUAL(feature_maps[0][0].getUniqueId(), incomplete_fhIt->getUniqueId())
}
END_SECTION

START_SECTION((void postRTHook(SimTypes::FeatureMapSimVector &)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void postDetectabilityHook(SimTypes::FeatureMapSimVector &)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void postIonizationHook(SimTypes::FeatureMapSimVector &)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void postRawMSHook(SimTypes::FeatureMapSimVector &)))
{
  NOT_TESTABLE
  // TODO Test ConsensusMap association
}
END_SECTION

START_SECTION((void postRawTandemMSHook(SimTypes::FeatureMapSimVector &, SimTypes::MSSimExperiment &)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((static BaseLabeler* create()))
{
  BaseLabeler* labeler = O18Labeler::create();
  TEST_NOT_EQUAL(labeler, base_nullPointer)
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



