// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Axel Walter $
// $Authors: Axel Walter $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ANALYSIS/ID/IonIdentityMolecularNetworking.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <fstream>
#include <vector>
#include <string>

///////////////////////////

START_TEST(IonIdentityMolecularNetworking, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

// Helper: build a ConsensusFeature with a list of IIMN linked (adduct) groups.
auto mkCF = [](double mz, const StringList& groups, const std::string& best_ion = "") -> ConsensusFeature
{
  ConsensusFeature cf;
  cf.setMZ(mz);
  cf.setRT(100.0);
  if (!groups.empty())
  {
    cf.setMetaValue(Constants::UserParam::IIMN_LINKED_GROUPS, DataValue(groups));
  }
  if (!best_ion.empty())
  {
    cf.setMetaValue(Constants::UserParam::IIMN_BEST_ION, best_ion);
  }
  return cf;
};

START_SECTION((static void annotateConsensusMap(ConsensusMap& consensus_map)))
{
  // --- Scenario from the class documentation: two adduct partners + one separate network ---
  // CF0 groups {1,2}, CF1 groups {2,3} -> share group 2 -> network 1
  // CF2 groups {4,5} -> isolated -> network 2
  {
    ConsensusMap cm;
    cm.push_back(mkCF(100.0, StringList{"1", "2"}));
    cm.push_back(mkCF(122.0, StringList{"2", "3"}));
    cm.push_back(mkCF(200.0, StringList{"4", "5"}));

    IonIdentityMolecularNetworking::annotateConsensusMap(cm);

    // 1-based row IDs assigned to every feature
    TEST_EQUAL((int)cm[0].getMetaValue(Constants::UserParam::IIMN_ROW_ID), 1)
    TEST_EQUAL((int)cm[1].getMetaValue(Constants::UserParam::IIMN_ROW_ID), 2)
    TEST_EQUAL((int)cm[2].getMetaValue(Constants::UserParam::IIMN_ROW_ID), 3)

    // network (connected-component) numbers
    TEST_EQUAL((int)cm[0].getMetaValue(Constants::UserParam::IIMN_ANNOTATION_NETWORK_NUMBER), 1)
    TEST_EQUAL((int)cm[1].getMetaValue(Constants::UserParam::IIMN_ANNOTATION_NETWORK_NUMBER), 1)
    TEST_EQUAL((int)cm[2].getMetaValue(Constants::UserParam::IIMN_ANNOTATION_NETWORK_NUMBER), 2)

    // adduct partners (semicolon-separated partner row IDs)
    TEST_STRING_EQUAL(cm[0].getMetaValue(Constants::UserParam::IIMN_ADDUCT_PARTNERS).toString(), "2")
    TEST_STRING_EQUAL(cm[1].getMetaValue(Constants::UserParam::IIMN_ADDUCT_PARTNERS).toString(), "1")
    // CF2 has groups but no partner -> no partner annotation
    TEST_EQUAL(cm[2].metaValueExists(Constants::UserParam::IIMN_ADDUCT_PARTNERS), false)

    // linked-groups meta value is consumed/removed from every feature
    TEST_EQUAL(cm[0].metaValueExists(Constants::UserParam::IIMN_LINKED_GROUPS), false)
    TEST_EQUAL(cm[1].metaValueExists(Constants::UserParam::IIMN_LINKED_GROUPS), false)
    TEST_EQUAL(cm[2].metaValueExists(Constants::UserParam::IIMN_LINKED_GROUPS), false)
  }

  // --- Chain: CF0-CF1-CF2 transitively connected -> single network, multi-partner join ---
  {
    ConsensusMap cm;
    cm.push_back(mkCF(100.0, StringList{"1", "2"}));
    cm.push_back(mkCF(120.0, StringList{"2", "3"}));
    cm.push_back(mkCF(140.0, StringList{"3", "4"}));

    IonIdentityMolecularNetworking::annotateConsensusMap(cm);

    TEST_EQUAL((int)cm[0].getMetaValue(Constants::UserParam::IIMN_ANNOTATION_NETWORK_NUMBER), 1)
    TEST_EQUAL((int)cm[1].getMetaValue(Constants::UserParam::IIMN_ANNOTATION_NETWORK_NUMBER), 1)
    TEST_EQUAL((int)cm[2].getMetaValue(Constants::UserParam::IIMN_ANNOTATION_NETWORK_NUMBER), 1)

    TEST_STRING_EQUAL(cm[0].getMetaValue(Constants::UserParam::IIMN_ADDUCT_PARTNERS).toString(), "2")
    TEST_STRING_EQUAL(cm[1].getMetaValue(Constants::UserParam::IIMN_ADDUCT_PARTNERS).toString(), "1;3")
    TEST_STRING_EQUAL(cm[2].getMetaValue(Constants::UserParam::IIMN_ADDUCT_PARTNERS).toString(), "2")
  }

  // --- A feature without linked groups receives only a row ID ---
  {
    ConsensusMap cm;
    cm.push_back(mkCF(100.0, StringList{"1", "2"}));
    cm.push_back(mkCF(122.0, StringList{"2", "3"}));
    cm.push_back(mkCF(300.0, StringList{})); // no IIMN_LINKED_GROUPS

    IonIdentityMolecularNetworking::annotateConsensusMap(cm);

    TEST_EQUAL((int)cm[2].getMetaValue(Constants::UserParam::IIMN_ROW_ID), 3)
    TEST_EQUAL(cm[2].metaValueExists(Constants::UserParam::IIMN_ANNOTATION_NETWORK_NUMBER), false)
    TEST_EQUAL(cm[2].metaValueExists(Constants::UserParam::IIMN_ADDUCT_PARTNERS), false)
  }
}
END_SECTION

START_SECTION((static void writeSupplementaryPairTable(const ConsensusMap& consensus_map, const std::string& output_file)))
{
  ConsensusMap cm;
  cm.push_back(mkCF(100.0, StringList{"1", "2"}, "[M+H]+"));
  cm.push_back(mkCF(122.0, StringList{"2", "3"}, "[M+Na]+"));
  cm.push_back(mkCF(200.0, StringList{"4", "5"}, "[M+H]+"));

  IonIdentityMolecularNetworking::annotateConsensusMap(cm);

  std::string out_file;
  NEW_TMP_FILE(out_file)
  IonIdentityMolecularNetworking::writeSupplementaryPairTable(cm, out_file);

  std::vector<std::string> lines;
  std::ifstream is(out_file.c_str());
  std::string line;
  while (std::getline(is, line))
  {
    if (!line.empty()) lines.push_back(line);
  }
  is.close();

  // comma-separated table: header + exactly one edge (CF0 <-> CF1)
  TEST_EQUAL(lines.size(), 2)
  TEST_STRING_EQUAL(lines[0], "ID1,ID2,EdgeType,Score,Annotation")
  // Score = (partners of CF0) + (partners of CF1) - 2 = 1 + 1 - 2 = 0
  // Annotation = "<best ion 1> <best ion 2> dm/z=<|mz1 - mz2|>"
  TEST_STRING_EQUAL(lines[1], "1,2,MS1 annotation,0,[M+H]+ [M+Na]+ dm/z=22.0")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
