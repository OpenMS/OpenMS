// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/FuzzyStringComparator.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/OMSFile.h>
#include <OpenMS/SYSTEM/File.h>

///////////////////////////

using namespace OpenMS;
using namespace std;


START_TEST(OMSFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

String oms_tmp;
IdentificationData ids;

START_SECTION(void store(const String& filename, const IdentificationData& id_data))
{
  vector<ProteinIdentification> proteins_in;
  vector<PeptideIdentification> peptides_in;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IdXMLFile_whole.idXML"), proteins_in, peptides_in);
  // IdentificationData doesn't allow score types with the same name, but different orientations:
  peptides_in[0].setHigherScoreBetter(true);

  IdentificationDataConverter::importIDs(ids, proteins_in, peptides_in);
  // add an adduct (not supported by idXML):
  AdductInfo adduct("Cl-", EmpiricalFormula("Cl"), -1);
  auto adduct_ref = ids.registerAdduct(adduct);
  IdentificationData::ObservationMatch match = *ids.getObservationMatches().begin();
  match.adduct_opt = adduct_ref;
  ids.registerObservationMatch(match);

  NEW_TMP_FILE(oms_tmp);
  OMSFile().store(oms_tmp, ids);
  TEST_EQUAL(File::empty(oms_tmp), false);
}
END_SECTION

START_SECTION(void load(const String& filename, IdentificationData& id_data))
{
  IdentificationData out;
  OMSFile().load(oms_tmp, out);

  TEST_EQUAL(ids.getInputFiles().size(), out.getInputFiles().size());
  TEST_EQUAL(ids.getScoreTypes().size(), out.getScoreTypes().size());
  TEST_EQUAL(ids.getProcessingSoftwares().size(),
             out.getProcessingSoftwares().size());
  TEST_EQUAL(ids.getDBSearchParams().size(), out.getDBSearchParams().size());
  TEST_EQUAL(ids.getProcessingSteps().size(),
             out.getProcessingSteps().size());
  TEST_EQUAL(ids.getObservations().size(), out.getObservations().size());
  TEST_EQUAL(ids.getParentSequences().size(),
             out.getParentSequences().size());
  TEST_EQUAL(ids.getParentGroupSets().size(),
             out.getParentGroupSets().size());
  TEST_EQUAL(ids.getIdentifiedPeptides().size(),
             out.getIdentifiedPeptides().size());
  TEST_EQUAL(ids.getIdentifiedOligos().size(),
             out.getIdentifiedOligos().size());
  TEST_EQUAL(ids.getIdentifiedCompounds().size(),
             out.getIdentifiedCompounds().size());
  TEST_EQUAL(ids.getAdducts().size(), out.getAdducts().size());
  TEST_EQUAL(ids.getObservationMatches().size(),
             out.getObservationMatches().size());
  auto it1 = ids.getObservationMatches().begin();
  auto it2 = out.getObservationMatches().begin();
  auto adduct_it = out.getObservationMatches().end();
  for (; (it1 != ids.getObservationMatches().end()) &&
         (it2 != out.getObservationMatches().end()); ++it1, ++it2)
  {
    TEST_EQUAL(it1->steps_and_scores.size(),
               it2->steps_and_scores.size());
    if (it2->adduct_opt) adduct_it = it2; // found PSM with adduct
  }
  // check PSM with adduct:
  TEST_EQUAL(adduct_it != out.getObservationMatches().end(), true);
  ABORT_IF(adduct_it == out.getObservationMatches().end());
  TEST_EQUAL(adduct_it->observation_ref->data_id,
             ids.getObservationMatches().begin()->observation_ref->data_id);
  TEST_EQUAL(adduct_it->identified_molecule_var.toString(),
             ids.getObservationMatches().begin()->identified_molecule_var.toString());
  TEST_EQUAL((*adduct_it->adduct_opt)->getName(), "Cl-");
}
END_SECTION

START_SECTION(void store(const String& filename, const FeatureMap& features))
{
  FeatureMap features;
  FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFileOMStest_1.featureXML"), features);
  // protein and peptide IDs use same score type (name) with different orientations;
  // IdentificationData doesn't allow this, so change it here:
  for (auto& run : features.getProteinIdentifications())
  {
    run.setScoreType(run.getScoreType() + "_protein");
  }
  IdentificationDataConverter::importFeatureIDs(features);

  NEW_TMP_FILE(oms_tmp);
  OMSFile().store(oms_tmp, features);
  TEST_EQUAL(File::empty(oms_tmp), false);
}
END_SECTION

START_SECTION(void load(const String& filename, FeatureMap& features))
{
  FeatureMap features;
  OMSFile().load(oms_tmp, features);

  TEST_EQUAL(features.size(), 2);
  TEST_EQUAL(features.at(0).getSubordinates().size(), 2);

  IdentificationDataConverter::exportFeatureIDs(features);
  // sort for reproducibility
  auto& proteins = features.getProteinIdentifications();
  for (auto& protein : proteins)
  {
    protein.sort();
  }
  auto& un_peptides = features.getUnassignedPeptideIdentifications();
  for (auto& un_pep : un_peptides)
  {
    un_pep.sort();
  }
  //features.setProteinIdentifications(proteins);
  //features.setUnassignedPeptideIdentifications(un_peptides);
  features.sortByPosition();

  String fxml_tmp;
  NEW_TMP_FILE(fxml_tmp);
  FeatureXMLFile().store(fxml_tmp, features);

  FuzzyStringComparator fsc;
  fsc.setAcceptableRelative(1.001);
  fsc.setAcceptableAbsolute(1);
  StringList sl;
  sl.push_back("xml-stylesheet");
  sl.push_back("UnassignedPeptideIdentification");
  fsc.setWhitelist(sl);

  TEST_EQUAL(fsc.compareFiles(fxml_tmp, OPENMS_GET_TEST_DATA_PATH("OMSFile_test_2.featureXML")), true);
}
END_SECTION

START_SECTION(void store(const String& filename, const ConsensusMap& consensus))
{
  ConsensusMap consensus;
  ConsensusXMLFile().load(OPENMS_GET_TEST_DATA_PATH("ConsensusXMLFile_1.consensusXML"), consensus);
  // protein and peptide IDs use same score type (name) with different orientations;
  // IdentificationData doesn't allow this, so change it here:
  for (auto& run : consensus.getProteinIdentifications())
  {
    run.setScoreType(run.getScoreType() + "_protein");
  }
  IdentificationDataConverter::importConsensusIDs(consensus);

  NEW_TMP_FILE(oms_tmp);
  OMSFile().store(oms_tmp, consensus);
  TEST_EQUAL(File::empty(oms_tmp), false);
}
END_SECTION

START_SECTION(void load(const String& filename, ConsensusMap& consensus))
{
  ConsensusMap consensus;
  OMSFile().load(oms_tmp, consensus);

  TEST_EQUAL(consensus.size(), 6);
  TEST_EQUAL(consensus.at(0).getFeatures().size(), 1);
  TEST_EQUAL(consensus.at(1).getFeatures().size(), 2);

  IdentificationDataConverter::exportConsensusIDs(consensus);
  // sort for reproducibility
  auto& proteins = consensus.getProteinIdentifications();
  for (auto& protein : proteins)
  {
    protein.sort();
  }
  auto& un_peptides = consensus.getUnassignedPeptideIdentifications();
  for (auto& un_pep : un_peptides)
  {
    un_pep.sort();
  }
  consensus.sortByPosition();

  String cxml_tmp;
  NEW_TMP_FILE(cxml_tmp);
  ConsensusXMLFile().store(cxml_tmp, consensus);
  TEST_EQUAL(File::empty(cxml_tmp), false);

  /*
  FuzzyStringComparator fsc;
  fsc.setAcceptableRelative(1.001);
  fsc.setAcceptableAbsolute(1);
  StringList sl;
  sl.push_back("xml-stylesheet");
  sl.push_back("UnassignedPeptideIdentification");
  fsc.setWhitelist(sl);

  TEST_EQUAL(fsc.compareFiles(cxml_tmp, OPENMS_GET_TEST_DATA_PATH("OMSFile_test_2.consensusXML")), true);
  */
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
