// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
using namespace OpenMS;
using ID = IdentificationData;
namespace
{
ID fixture(ID::MoleculeKind kind = ID::MoleculeKind::PEPTIDE)
{
  ID data;
  auto& run = data.addRun("search", kind);
  ID::ScoreDefinition score;
  score.name = "score";
  score.software = "engine";
  run.setPrimaryScore(run.addScore(score));
  ID::RunSettings settings;
  settings.software = "engine";
  settings.search.db = "database";
  run.setSettings(settings);
  ID::SourceFile source;
  source.path = "input.mzML";
  ID::Observation observation;
  observation.data_id = "scan=42";
  observation.rt = 12;
  observation.mz = 500;
  auto query = run.addIdentification(run.addSource(source), observation);
  ID::MatchData match;
  match.representation = kind == ID::MoleculeKind::OLIGONUCLEOTIDE ? "ACUGp" : "PEPM(Oxidation)K";
  match.encoding = kind == ID::MoleculeKind::OLIGONUCLEOTIDE ? ID::Encoding::NA_SEQUENCE : ID::Encoding::AA_SEQUENCE;
  match.charge = 2;
  match.target_decoy = ID::TargetDecoy::TARGET;
  match.parent_evidence.push_back({{"database", "parent"}, 3, 8, "K", "R"});
  match.setMetaValue("numbers", IntList {1, 2, 3});
  run.addMatch(query, match, {99.0});
  ID::ParentRecord parent;
  parent.identity = {"database", "parent"};
  parent.sequence = "KKKPEPMK";
  parent.target_decoy = ID::TargetDecoy::TARGET;
  run.setParents(std::vector {parent});
  return data;
}
} // namespace
START_TEST(IdentificationDataConverter, "$Id$")
START_SECTION((peptide conversion preserves modifications, evidence, scores and typed metadata))
{
  auto data = fixture();
  std::vector<ProteinIdentification> proteins;
  PeptideIdentificationList peptides;
  IdentificationDataConverter::exportIDs(data, proteins, peptides);
  TEST_EQUAL(proteins.size(), 1);
  TEST_EQUAL(peptides.size(), 1);
  const auto& hit = peptides[0].getHits()[0];
  TEST_EQUAL(hit.getSequence().toString(), "PEPM(Oxidation)K");
  TEST_REAL_SIMILAR(hit.getScore(), 99);
  TEST_EQUAL(hit.getPeptideEvidences()[0].getStart(), 3);
  TEST_EQUAL(hit.getMetaValue("numbers"), DataValue(IntList {1, 2, 3}));
  ID imported;
  IdentificationDataConverter::importIDs(imported, proteins, peptides);
  TEST_EQUAL(imported.getRuns().front().getNumberOfMatches(), 1);
  TEST_EQUAL(imported.getRuns().front().getSources()[0].identifications[0].getMatches()[0].representation, "PEPM(Oxidation)K");
}
END_SECTION
START_SECTION((RNA idXML convention and mzTab export retain molecular identity))
{
  auto data = fixture(ID::MoleculeKind::OLIGONUCLEOTIDE);
  std::vector<ProteinIdentification> proteins;
  PeptideIdentificationList peptides;
  IdentificationDataConverter::exportIDs(data, proteins, peptides);
  TEST_EQUAL(peptides[0].getHits()[0].getSequence().empty(), true);
  TEST_EQUAL(peptides[0].getHits()[0].getMetaValue("label"), "ACUGp");
  TEST_EQUAL(peptides[0].getHits()[0].getMetaValue("molecule_type"), "RNA");
  ID restored;
  IdentificationDataConverter::importIDs(restored, proteins, peptides);
  TEST_EQUAL(restored.getRuns().front().getMoleculeKind() == ID::MoleculeKind::OLIGONUCLEOTIDE, true);
  TEST_EQUAL(restored.getRuns().front().getSources()[0].identifications[0].getMatches()[0].representation, "ACUGp");
  const auto tab = IdentificationDataConverter::exportMzTab(data);
  TEST_EQUAL(tab.getOSMSectionRows().size(), 1);
  TEST_EQUAL(tab.getOSMSectionRows()[0].sequence.get(), "ACUGp");
}
END_SECTION
START_SECTION((feature ownership survives copies and reduced - data export including empty queries))
{
  FeatureMap map;
  map.resize(1);
  map[0].getSubordinates().resize(1);
  map.getIdentificationData() = fixture();
  auto& run = map.getIdentificationData().getRun("search");
  const auto& query = run.getSources()[0].identifications[0];
  ID::MatchReference link {run.getUuid(), query.getMatches()[0].getId()};
  map[0].getSubordinates()[0].addIDMatch(link);
  ID::Observation empty;
  empty.data_id = "scan=43";
  auto empty_id = run.addIdentification(run.getSourceId(0), empty);
  map[0].addIDQuery({run.getUuid(), empty_id});
  auto copy = map;
  TEST_EQUAL(copy == map, true);
  TEST_EQUAL(copy.getUnassignedIDMatches().empty(), true);
  TEST_EQUAL(copy[0].getSubordinates()[0].getAnnotationState(copy.getIdentificationData()) == BaseFeature::AnnotationState::FEATURE_ID_SINGLE, true);
  IdentificationDataConverter::exportFeatureIDs(copy, false);
  TEST_EQUAL(copy[0].getPeptideIdentifications().size(), 1);
  TEST_EQUAL(copy[0].getPeptideIdentifications()[0].getHits().empty(), true);
  TEST_EQUAL(copy[0].getSubordinates()[0].getPeptideIdentifications()[0].getHits().size(), 1);
  IdentificationDataConverter::exportFeatureIDs(copy, true);
  TEST_EQUAL(copy.getIdentificationData().empty(), true);
  TEST_EQUAL(copy[0].getIDQueries().empty(), true);
  TEST_EQUAL(copy[0].getSubordinates()[0].getIDMatches().empty(), true);
  IdentificationDataConverter::importFeatureIDs(copy, true);
  TEST_EQUAL(copy[0].getPeptideIdentifications().empty(), true);
  TEST_EQUAL(copy.getIdentificationData().getRuns().front().getNumberOfIdentifications(), 2);
}
END_SECTION
START_SECTION((RNA feature conversion retains subordinate and empty - query associations))
{
  FeatureMap map;
  map.resize(1);
  map[0].getSubordinates().resize(1);
  map.getIdentificationData() = fixture(ID::MoleculeKind::OLIGONUCLEOTIDE);
  auto& run = map.getIdentificationData().getRun("search");
  const auto& query = run.getSources()[0].identifications[0];
  map[0].getSubordinates()[0].addIDMatch({run.getUuid(), query.getMatches()[0].getId()});
  ID::Observation empty;
  empty.data_id = "scan=empty";
  map[0].addIDQuery({run.getUuid(), run.addIdentification(run.getSourceId(0), empty)});
  IdentificationDataConverter::exportFeatureIDs(map, true);
  TEST_EQUAL(map[0].getPeptideIdentifications().size(), 1);
  TEST_EQUAL(map[0].getPeptideIdentifications()[0].getHits().empty(), true);
  TEST_EQUAL(map[0].getSubordinates()[0].getPeptideIdentifications()[0].getHits()[0].getMetaValue("label"), "ACUGp");
  IdentificationDataConverter::importFeatureIDs(map, true);
  TEST_EQUAL(map.getIdentificationData().getRuns().front().getMoleculeKind() == ID::MoleculeKind::OLIGONUCLEOTIDE, true);
  TEST_EQUAL(map[0].getIDQueries().size(), 1);
  TEST_EQUAL(map[0].getSubordinates()[0].getIDMatches().size(), 1);
  TEST_EQUAL(map.getUnassignedIDMatches().empty(), true);
}
END_SECTION
START_SECTION((merging independent identically named runs preserves stable associations))
{
  auto left = fixture();
  auto right = fixture();
  const auto uuid = right.getRuns().front().getUuid();
  left.merge(right);
  TEST_EQUAL(left.getRuns().size(), 2);
  TEST_EQUAL(left.findRunByUuid(uuid)->getIdentifier(), "search#2");
  TEST_EQUAL(left.findRunByUuid(uuid)->getNumberOfMatches(), 1);
  auto copy = left;
  left.merge(copy);
  TEST_EQUAL(left == copy, true);
}
END_SECTION
END_TEST
