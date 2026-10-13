// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Swenja Wagner, Patricia Scheil $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/QC/MissedCleavages.h>
#include <OpenMS/QC/QCBase.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
#include <string>
#include <vector>

//////////////////////////

using namespace OpenMS;

START_TEST(MissedCleavages, "$Id$")
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


// construct PeptideHits
PeptideHit pep_hit_no_cut;
pep_hit_no_cut.setSequence(AASequence::fromString("AAAAAAAAAAAAAK"));
PeptideHit pep_hit_one_cut;
pep_hit_one_cut.setSequence(AASequence::fromString("AAAAKAAAAAR"));
PeptideHit pep_hit_three_cuts;
pep_hit_three_cuts.setSequence(AASequence::fromString("AAAKAAARAAAKAAAR"));

// construct vectors of PeptideHits
std::vector<PeptideHit> pep_hits_0 = {pep_hit_no_cut};
std::vector<PeptideHit> pep_hits_1 = {pep_hit_one_cut};
std::vector<PeptideHit> pep_hits_3 = {pep_hit_three_cuts};
std::vector<PeptideHit> pep_hits_empty = {};

// construct PeptideIdentification with PeptideHits
PeptideIdentification pep_id_0;
pep_id_0.setHits(pep_hits_0);
PeptideIdentification pep_id_1;
pep_id_1.setHits(pep_hits_1);
PeptideIdentification pep_id_empty;
pep_id_empty.setHits(pep_hits_empty);
PeptideIdentification pep_id_3;
pep_id_3.setHits(pep_hits_3);
// identification data needs a score type
for (auto* id : {&pep_id_0, &pep_id_1, &pep_id_empty, &pep_id_3}) id->setScoreType("score");

// construct vectors of PeptideIdentifications
PeptideIdentificationList pep_ids = {pep_id_0, pep_id_1, pep_id_empty};
PeptideIdentificationList pep_ids_1 = {pep_id_1, pep_id_1};
PeptideIdentificationList pep_ids_empty {};
PeptideIdentificationList pep_ids_3 = {pep_id_3};

// construct features (their peptide identifications are given with the identifications of the maps)
Feature feat_empty_pi;
feat_empty_pi.setMetaValue("FWHM", 32.21);

Feature feat;
feat.setMetaValue("FWHM", 32.21);

Feature feat_empty;
feat_empty.setMetaValue("FWHM", 32.21);

Feature feat_3;
feat_3.setMetaValue("FWHM", 32.21);

// FeatureMap
FeatureMap feature_map;
FeatureMap feature_map_3;
FeatureMap feature_map_empty;
FeatureMap feature_map_no_protein;
FeatureMap feature_map_no_enzyme;

// stores data in the FeatureMaps
feature_map.push_back(feat_empty_pi);
feature_map.push_back(feat);
feature_map.push_back(feat_empty);
MapPeptideIdentifications ids;
ids.features[{1}] = pep_ids;
ids.unassigned = pep_ids_1;
ids.proteins.resize(1);
ids.proteins[0].getSearchParameters().digestion_enzyme = *ProteaseDB::getInstance()->getEnzyme("trypsin");
ids.proteins[0].getSearchParameters().missed_cleavages = 2;
IdentificationDataConverter::importMapIdentifications(feature_map, ids);

// FeatureMap with more missed cleavages than allowed
feature_map_3.push_back(feat_3);
MapPeptideIdentifications ids_3;
ids_3.features[{0}] = pep_ids_3;
ids_3.proteins.resize(1);
ids_3.proteins[0].getSearchParameters().digestion_enzyme.setName("trypsin");
ids_3.proteins[0].getSearchParameters().missed_cleavages = 2;
IdentificationDataConverter::importMapIdentifications(feature_map_3, ids_3);

// FeatureMap without ProteinIdentifications
feature_map_no_protein.push_back(feat);
MapPeptideIdentifications ids_no_protein;
ids_no_protein.features[{0}] = pep_ids;

// FeatureMap without given enzyme
feature_map_no_enzyme.push_back(feat);
MapPeptideIdentifications ids_no_enzyme;
ids_no_enzyme.features[{0}] = pep_ids;
ids_no_enzyme.proteins.resize(1);
IdentificationDataConverter::importMapIdentifications(feature_map_no_enzyme, ids_no_enzyme);

//////////////////////////////////////////////////////////////////
// start Section
/////////////////////////////////////////////////////////////////

MissedCleavages* ptr = nullptr;
MissedCleavages* nulpt = nullptr;
START_SECTION(MissedCleavages())
{
  ptr = new MissedCleavages();
  TEST_NOT_EQUAL(ptr, nulpt)
}
END_SECTION


START_SECTION(~MissedCleavages())
{
  delete ptr;
}
END_SECTION


// tests compute function
START_SECTION(void compute(FeatureMap& fmap))
{
  // test with valid input
  MissedCleavages mc;
  mc.compute(feature_map);
  auto result = mc.getResults();

  TEST_EQUAL(result[0].size(), 2)
  TEST_EQUAL(result[0][0], 1)
  TEST_EQUAL(result[0][1], 3)

  std::vector<UInt32> frequ;

  // test if result is stored as MetaInformation in PeptideHits in FeatureMap
  auto lam = [&frequ](PeptideIdentification& pep_id) {
    if (pep_id.getHits().empty())
    {
      return;
    }
    if (pep_id.getHits()[0].metaValueExists("missed_cleavages"))
    {
      frequ.push_back(pep_id.getHits()[0].getMetaValue("missed_cleavages"));
    }
  };

  const auto exported = IdentificationDataConverter::exportMapIdentifications(feature_map);
  for (const auto* pep_id : exported.peptides())
  {
    PeptideIdentification copy = *pep_id;
    lam(copy);
  }
  TEST_EQUAL(frequ.size(), 4)
  TEST_EQUAL(frequ[0], 0)
  TEST_EQUAL(frequ[1], 1)
  TEST_EQUAL(frequ[2], 1)
  TEST_EQUAL(frequ[3], 1)

  // empty feature map
  MissedCleavages mc_empty;
  mc_empty.compute(feature_map_empty);
  auto result_empty = mc_empty.getResults();

  TEST_EQUAL(result_empty[0].empty(), true)

  // Missing informations in ProteinIdentifications
  // fmap.getProteinIdentifications().empty()
  // peptide identifications without their search run (ProteinIdentification) cannot be the identifications of a map
  TEST_EXCEPTION(Exception::InvalidParameter, IdentificationDataConverter::importMapIdentifications(feature_map_no_protein, ids_no_protein))

  // no given enzyme
  // enzyme == "unknown_enzyme"
  MissedCleavages mc_no_enzyme;
  TEST_EXCEPTION_WITH_MESSAGE(Exception::MissingInformation, mc_no_enzyme.compute(feature_map_no_enzyme), "No digestion enzyme in FeatureMap detected. No computation possible.")

  // Number of missed cleavages is greater than the allowed maximum number of missed cleavages.
  MissedCleavages mc_3;
  mc_3.compute(feature_map_3);
  auto result_3 = mc_3.getResults();

  TEST_EQUAL(result_3[0].size(), 1)
  TEST_EQUAL(result_3[0][3], 1)
}
END_SECTION


MissedCleavages mc;

START_SECTION(const std::string& getName() const override) {TEST_EQUAL(mc.getName(), "MissedCleavages")} END_SECTION


  START_SECTION(QCBase::Status requirements() const override) {TEST_EQUAL(QCBase::Status(QCBase::Requires::POSTFDRFEAT) == mc.requirements(), true)} END_SECTION

  /////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////
  END_TEST
