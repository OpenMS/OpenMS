// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/Scores.h>
///////////////////////////

#include <set>

using namespace OpenMS;
using namespace std;

using IDType = Scores::IDType;

START_TEST(Scores, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static Scores::IDType parseIDType(const std::string& score_type)))
{
  // type-name strings -> enum (case-insensitive; '-', '_', ' ' ignored)
  TEST_EQUAL(Scores::parseIDType("raw") == IDType::RAW, true)
  TEST_EQUAL(Scores::parseIDType("RAW") == IDType::RAW, true)
  TEST_EQUAL(Scores::parseIDType("raw evalue") == IDType::RAW_EVAL, true)
  TEST_EQUAL(Scores::parseIDType("q-value") == IDType::QVAL, true)
  TEST_EQUAL(Scores::parseIDType("qvalue") == IDType::QVAL, true)
  TEST_EQUAL(Scores::parseIDType("FDR") == IDType::FDR, true)
  TEST_EQUAL(Scores::parseIDType("false discovery rate") == IDType::FDR, true)
  TEST_EQUAL(Scores::parseIDType("Posterior Error Probability") == IDType::PEP, true)
  TEST_EQUAL(Scores::parseIDType("pep") == IDType::PEP, true)
  TEST_EQUAL(Scores::parseIDType("Posterior Probability") == IDType::PP, true)
  TEST_EQUAL(Scores::parseIDType("pp") == IDType::PP, true)

  // a trailing "_score" suffix is stripped before matching
  TEST_EQUAL(Scores::parseIDType("q-value_score") == IDType::QVAL, true)

  // an unknown type string throws
  TEST_EXCEPTION(Exception::MissingInformation, Scores::parseIDType("not-a-score-type"))
}
END_SECTION

START_SECTION((static bool isScoreType(const std::string& score_name, IDType type)))
{
  // exact (case-sensitive) membership in the type's name set, after the
  // "_score" suffix is stripped
  TEST_EQUAL(Scores::isScoreType("hyperscore", IDType::RAW), true)
  TEST_EQUAL(Scores::isScoreType("ln(hyperscore)", IDType::RAW), true)
  TEST_EQUAL(Scores::isScoreType("q-value", IDType::QVAL), true)
  TEST_EQUAL(Scores::isScoreType("q-value_score", IDType::QVAL), true) // suffix stripped
  TEST_EQUAL(Scores::isScoreType("expect", IDType::RAW_EVAL), true)
  TEST_EQUAL(Scores::isScoreType("FDR", IDType::FDR), true)

  // name not in the requested type's set
  TEST_EQUAL(Scores::isScoreType("hyperscore", IDType::QVAL), false)
  TEST_EQUAL(Scores::isScoreType("unknown_name", IDType::RAW), false)
}
END_SECTION

START_SECTION((static bool isHigherBetter(IDType type)))
{
  TEST_EQUAL(Scores::isHigherBetter(IDType::RAW), true)
  TEST_EQUAL(Scores::isHigherBetter(IDType::PP), true)
  TEST_EQUAL(Scores::isHigherBetter(IDType::RAW_EVAL), false)
  TEST_EQUAL(Scores::isHigherBetter(IDType::PEP), false)
  TEST_EQUAL(Scores::isHigherBetter(IDType::FDR), false)
  TEST_EQUAL(Scores::isHigherBetter(IDType::QVAL), false)
}
END_SECTION

START_SECTION((static std::vector<std::string> getAllIDScoreNames()))
{
  std::vector<std::string> names = Scores::getAllIDScoreNames();
  // pins the current registry: 9 RAW + 11 RAW_EVAL + 1 PP + 6 PEP + 4 FDR + 7 QVAL
  TEST_EQUAL(names.size(), 38)
  std::set<std::string> s(names.begin(), names.end());
  TEST_EQUAL(s.count("q-value"), 1)
  TEST_EQUAL(s.count("expect"), 1)
  TEST_EQUAL(s.count("hyperscore"), 1)
  TEST_EQUAL(s.count("FDR"), 1)
  TEST_EQUAL(s.count("Posterior Probability"), 1)
}
END_SECTION

START_SECTION((static const std::set<std::string>& getIDNamesForType(IDType type)))
{
  TEST_EQUAL(Scores::getIDNamesForType(IDType::RAW).size(), 9)
  TEST_EQUAL(Scores::getIDNamesForType(IDType::RAW).count("hyperscore"), 1)
  TEST_EQUAL(Scores::getIDNamesForType(IDType::RAW_EVAL).size(), 11)
  TEST_EQUAL(Scores::getIDNamesForType(IDType::RAW_EVAL).count("expect"), 1)
  TEST_EQUAL(Scores::getIDNamesForType(IDType::PP).size(), 1)
  TEST_EQUAL(Scores::getIDNamesForType(IDType::PEP).size(), 6)
  TEST_EQUAL(Scores::getIDNamesForType(IDType::FDR).size(), 4)
  TEST_EQUAL(Scores::getIDNamesForType(IDType::FDR).count("fdr"), 1)
  TEST_EQUAL(Scores::getIDNamesForType(IDType::QVAL).size(), 7)
  TEST_EQUAL(Scores::getIDNamesForType(IDType::QVAL).count("q-value"), 1)
}
END_SECTION

START_SECTION((static bool findIDTypeByName(const std::string& name, IDType& type)))
{
  IDType t = IDType::QVAL;
  TEST_EQUAL(Scores::findIDTypeByName("hyperscore", t), true)
  TEST_EQUAL(t == IDType::RAW, true)
  TEST_EQUAL(Scores::findIDTypeByName("expect", t), true)
  TEST_EQUAL(t == IDType::RAW_EVAL, true)
  TEST_EQUAL(Scores::findIDTypeByName("q-value", t), true)
  TEST_EQUAL(t == IDType::QVAL, true)
  TEST_EQUAL(Scores::findIDTypeByName("Posterior Error Probability", t), true)
  TEST_EQUAL(t == IDType::PEP, true)
  // search engine scores by PSI-MS accession, all lower-is-better
  const std::vector<std::pair<std::string, IDType>> accessions{
    {"MS:1001330", IDType::RAW_EVAL}, // X!Tandem:expect
    {"MS:1001328", IDType::RAW_EVAL}, // OMSSA:evalue
    {"MS:1001329", IDType::RAW_EVAL}, // OMSSA:pvalue
    {"MS:1001172", IDType::RAW_EVAL}, // Mascot:expectation value
    {"MS:1002931", IDType::RAW_EVAL}, // TopPIC:spectral p-value
    {"MS:1002056", IDType::PEP},      // MS-GF:PEP
    {"MS:1002929", IDType::FDR},      // TopPIC:spectral FDR
    {"MS:1002054", IDType::QVAL},     // MS-GF:QValue
    {"MS:1002055", IDType::QVAL}      // MS-GF:PepQValue
  };
  for (const auto& [accession, type] : accessions)
  {
    t = IDType::RAW;
    TEST_EQUAL(Scores::findIDTypeByName(accession, t), true)
    TEST_EQUAL(t == type, true)
    TEST_EQUAL(Scores::isHigherBetter(t), false)
  }

  // unknown name -> false
  TEST_EQUAL(Scores::findIDTypeByName("definitely_not_a_score", t), false)
}
END_SECTION

START_SECTION((static std::string normalizeScoreName(const std::string& score_name)))
{
  TEST_EQUAL(Scores::normalizeScoreName("q-value_score"), "q-value")
  TEST_EQUAL(Scores::normalizeScoreName("q-value"), "q-value")
  // "hyperscore" ends in "score" but not in the "_score" suffix -> unchanged
  TEST_EQUAL(Scores::normalizeScoreName("hyperscore"), "hyperscore")
}
END_SECTION

START_SECTION((static bool isKnownScoreType(const std::string& score_name)))
{
  TEST_EQUAL(Scores::isKnownScoreType("q-value"), true)
  TEST_EQUAL(Scores::isKnownScoreType("q-value_score"), true)  // normalized before lookup
  TEST_EQUAL(Scores::isKnownScoreType("hyperscore"), true)
  TEST_EQUAL(Scores::isKnownScoreType("expect"), true)
  TEST_EQUAL(Scores::isKnownScoreType("definitely_not_a_score"), false)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
