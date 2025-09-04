// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Immanuel Luhn$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/FORMAT/IdXMLFile.h>

#include <vector>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(IDRipper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


///load input data
std::vector< ProteinIdentification > protein_identifications;
PeptideIdentificationList identifications;
String document_id;
IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IDScoreSwitcherAlgorithm_test_input.idXML"), protein_identifications, identifications, document_id);
PeptideIdentification identification = identifications[0];
ProteinIdentification protein_identification = protein_identifications[0];

IDScoreSwitcherAlgorithm* ptr = nullptr;
IDScoreSwitcherAlgorithm* null_ptr = nullptr;
START_SECTION(IDScoreSwitcherAlgorithm())
{
  ptr = new IDScoreSwitcherAlgorithm();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~IDScoreSwitcherAlgorithm())
{
  delete ptr;
}
END_SECTION

START_SECTION(switchToGeneralScoreType)
{
  IDScoreSwitcherAlgorithm switcher{};
  Size c(0);
  switcher.switchToGeneralScoreType(identifications, IDScoreSwitcherAlgorithm::ScoreType::PEP, c);
  TEST_EQUAL(identifications[0].getScoreType(), "Posterior Error Probability");
}
END_SECTION

START_SECTION(findScoreType)
{
  IDScoreSwitcherAlgorithm switcher{};
  
  // Test case 1: When main score is already a PEP score
  PeptideIdentification pep_id_with_pep_main;
  pep_id_with_pep_main.setScoreType("Posterior Error Probability");
  pep_id_with_pep_main.setHigherScoreBetter(false);
  
  PeptideHit hit1;
  hit1.setScore(0.05);
  pep_id_with_pep_main.insertHit(hit1);
  
  IDScoreSwitcherAlgorithm::ScoreSearchResult result1 = switcher.findScoreType(pep_id_with_pep_main, IDScoreSwitcherAlgorithm::ScoreType::PEP);
  TEST_EQUAL(result1.is_main_score_type, true);
  TEST_EQUAL(result1.score_name, "Posterior Error Probability");
  
  // Test case 2: When main score is not PEP but PEP is available in meta values
  PeptideIdentification pep_id_with_pep_meta;
  pep_id_with_pep_meta.setScoreType("XTandem");
  pep_id_with_pep_meta.setHigherScoreBetter(true);
  
  PeptideHit hit2;
  hit2.setScore(100.0);
  hit2.setMetaValue("pep", 0.01);
  hit2.setMetaValue("other_score", 0.5);
  pep_id_with_pep_meta.insertHit(hit2);
  
  auto result2 = switcher.findScoreType(pep_id_with_pep_meta, IDScoreSwitcherAlgorithm::ScoreType::PEP);
  TEST_EQUAL(result2.is_main_score_type, false);
  TEST_EQUAL(result2.score_name, "pep");
  
  // Test case 3: When main score is not PEP and no PEP available in meta values
  PeptideIdentification pep_id_no_pep;
  pep_id_no_pep.setScoreType("Mascot");
  pep_id_no_pep.setHigherScoreBetter(true);
  
  PeptideHit hit3;
  hit3.setScore(50.0);
  hit3.setMetaValue("e_value", 0.001);
  pep_id_no_pep.insertHit(hit3);
  
  auto result3 = switcher.findScoreType(pep_id_no_pep, IDScoreSwitcherAlgorithm::ScoreType::PEP);
  TEST_EQUAL(result3.is_main_score_type, false);
  TEST_EQUAL(result3.score_name.empty(), true);
  
  // Test case 4: Check various PEP score name variants from the enum collection
  PeptideIdentification pep_id_uppercase;
  pep_id_uppercase.setScoreType("Mascot");
  
  PeptideHit hit4;
  hit4.setScore(75.0);
  hit4.setMetaValue("PEP", 0.02);  // Uppercase variant
  pep_id_uppercase.insertHit(hit4);
  
  auto result4 = switcher.findScoreType(pep_id_uppercase, IDScoreSwitcherAlgorithm::ScoreType::PEP);
  TEST_EQUAL(result4.is_main_score_type, false);
  TEST_EQUAL(result4.score_name, "PEP");
  
  // Test case 5: Check _score suffix variant
  PeptideIdentification pep_id_suffix;
  pep_id_suffix.setScoreType("SEQUEST:xcorr");
  
  PeptideHit hit5;
  hit5.setScore(2.5);
  hit5.setMetaValue("pep_score", 0.03);  // With _score suffix
  pep_id_suffix.insertHit(hit5);
  
  auto result5 = switcher.findScoreType(pep_id_suffix, IDScoreSwitcherAlgorithm::ScoreType::PEP);
  TEST_EQUAL(result5.is_main_score_type, false);
  TEST_EQUAL(result5.score_name, "pep_score");
  
  // Test case 6: Test with Q-value score type
  PeptideIdentification qval_id_main;
  qval_id_main.setScoreType("q-value");
  qval_id_main.setHigherScoreBetter(false);
  
  PeptideHit hit6;
  hit6.setScore(0.02);
  qval_id_main.insertHit(hit6);
  
  auto result6 = switcher.findScoreType(qval_id_main, IDScoreSwitcherAlgorithm::ScoreType::QVAL);
  TEST_EQUAL(result6.is_main_score_type, true);
  TEST_EQUAL(result6.score_name, "q-value");
  
  // Test case 7: Test with Q-value in meta values
  PeptideIdentification qval_id_meta;
  qval_id_meta.setScoreType("Mascot");
  qval_id_meta.setHigherScoreBetter(true);
  
  PeptideHit hit7;
  hit7.setScore(60.0);
  hit7.setMetaValue("qvalue", 0.05);
  qval_id_meta.insertHit(hit7);
  
  auto result7 = switcher.findScoreType(qval_id_meta, IDScoreSwitcherAlgorithm::ScoreType::QVAL);
  TEST_EQUAL(result7.is_main_score_type, false);
  TEST_EQUAL(result7.score_name, "qvalue");
  
  // Test case 8: Test with FDR score type
  PeptideIdentification fdr_id_meta;
  fdr_id_meta.setScoreType("XTandem");
  fdr_id_meta.setHigherScoreBetter(true);
  
  PeptideHit hit8;
  hit8.setScore(120.0);
  hit8.setMetaValue("FDR", 0.01);
  fdr_id_meta.insertHit(hit8);
  
  auto result8 = switcher.findScoreType(fdr_id_meta, IDScoreSwitcherAlgorithm::ScoreType::FDR);
  TEST_EQUAL(result8.is_main_score_type, false);
  TEST_EQUAL(result8.score_name, "FDR");
}
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
