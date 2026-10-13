// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
#include <OpenMS/CHEMISTRY/AdductInfo.h>
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/test_config.h>
#include <limits>
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
  run.setSettings(settings);
  ID::Database database;
  database.path = "database";
  const auto database_id = run.addDatabase(database);
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
  match.sequence_evidence.push_back({database_id, "parent", 3, 8, 'K', 'R'});
  match.setMetaValue("numbers", IntList {1, 2, 3});
  run.addMatch(query, match, {99.0});
  ID::DatabaseSequence parent;
  parent.database = database_id;
  parent.accession = "parent";
  parent.sequence = "KKKPEPMK";
  parent.target_decoy = ID::TargetDecoy::TARGET;
  run.setDatabaseSequences(std::vector {parent});
  return data;
}
} // namespace
namespace
{
  // the sequences of hits (or matches) of identifications, as one string per identification
  std::vector<std::string> sequences(const PeptideIdentificationList& ids)
  {
    std::vector<std::string> result;
    for (const auto& id : ids)
    {
      std::string sequences = id.getIdentifier() + ":" + std::to_string(id.getRT()) + ":";
      for (const auto& hit : id.getHits())
        sequences += hit.getSequence().toString() + "/" + std::to_string(hit.getScore()) + ",";
      result.push_back(sequences);
    }
    return result;
  }
  std::vector<std::string> sequences(const std::vector<ID::QueryMatches>& entries)
  {
    std::vector<std::string> result;
    for (const auto& entry : entries)
    {
      std::string sequences = entry.run->getIdentifier() + ":" + std::to_string(entry.query->rt.value_or(std::numeric_limits<double>::quiet_NaN())) + ":";
      for (const auto* match : entry.matches)
        sequences += match->representation + "/" + std::to_string(*entry.run->getScore(match->getId(), *entry.run->getPrimaryScore())) + ",";
      result.push_back(sequences);
    }
    return result;
  }
  /// The links of the features of a map against the identifications of the map as peptide identifications
  template<class Map>
  void checkLinks(const Map& map)
  {
    const auto ids = IdentificationDataConverter::exportMapIdentifications(map);
    const auto& data = map.getIdentificationData();
    Size best_matches = 0;
    for (Size i = 0; i < map.size(); ++i)
    {
      const auto& peptides = ids.at({i});
      TEST_TRUE(sequences(peptides) == sequences(map[i].getLinkedIdentifications(data)))
      // the best hit of the peptide identifications (the first of equally good ones)
      const PeptideHit* best_hit = nullptr;
      for (const auto& peptide : peptides)
        for (const auto& hit : peptide.getHits())
          if (! best_hit || (peptide.isHigherScoreBetter() ? hit.getScore() > best_hit->getScore() : hit.getScore() < best_hit->getScore()))
            best_hit = &hit;
      const auto best = map[i].getBestLinkedMatch(data);
      TEST_EQUAL(best.has_value(), best_hit != nullptr)
      if (! best || ! best_hit) continue;
      ++best_matches;
      TEST_EQUAL(best->matches.size(), 1)
      TEST_REAL_SIMILAR(*best->run->getScore(best->matches[0]->getId(), *best->run->getPrimaryScore()), best_hit->getScore())
    }
    TEST_TRUE(best_matches > 0)
    TEST_TRUE(! ids.unassigned.empty())
    TEST_TRUE(sequences(ids.unassigned) == sequences(map.getUnassignedIdentifications()))
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
START_SECTION(([EXTRA] RNA mzTab calculated m/z includes the assigned adduct))
{
  auto data = fixture(ID::MoleculeKind::OLIGONUCLEOTIDE);
  auto& run = data.getRun("search");
  const auto& match = run.getSources()[0].identifications[0].getMatches()[0];
  ID::MatchData sodium = match.getData();
  sodium.charge = -2;
  // the adduct of the ion, as NucleicAcidSearchEngine stores "Na:+" for charge -2: [M+Na-3H]2-
  sodium.details.emplace().adduct = AdductInfo("Na:+", EmpiricalFormula("Na") + EmpiricalFormula("H") * -3, -2);
  run.replaceMatch(match.getId(), sodium, run.getScores(match));
  const auto rows = IdentificationDataConverter::exportMzTab(data).getOSMSectionRows();
  TEST_EQUAL(rows.size(), 1);
  TEST_REAL_SIMILAR(rows[0].calc_mass_to_charge.get(),
                    (NASequence::fromString("ACUGp").getMonoWeight(NASequence::Full, -2) + sodium.details->adduct->getMassShift()) / 2);
  // sodium in place of a proton
  TEST_REAL_SIMILAR(rows[0].calc_mass_to_charge.get(), (NASequence::fromString("ACUGp").getMonoWeight(NASequence::Full, -2)
                                                        + EmpiricalFormula("Na").getMonoWeight() - EmpiricalFormula("H").getMonoWeight()) / 2);
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
  const auto ids = IdentificationDataConverter::exportMapIdentifications(copy);
  TEST_EQUAL(ids.at({0}).size(), 1);
  TEST_EQUAL(ids.at({0})[0].getHits().empty(), true);
  TEST_EQUAL(ids.at({0, 0})[0].getHits().size(), 1);
  copy.clearIdentifications();
  TEST_EQUAL(copy.getIdentificationData().empty(), true);
  TEST_EQUAL(copy[0].getIDQueries().empty(), true);
  TEST_EQUAL(copy[0].getSubordinates()[0].getIDMatches().empty(), true);
  IdentificationDataConverter::importMapIdentifications(copy, ids);
  TEST_EQUAL(copy.getIdentificationData().getRuns().front().getNumberOfIdentifications(), 2);
  TEST_EQUAL(copy[0].getIDQueries().size(), 1);
  TEST_EQUAL(copy[0].getSubordinates()[0].getIDMatches().size(), 1);
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
  const auto ids = IdentificationDataConverter::exportMapIdentifications(map);
  TEST_EQUAL(ids.at({0}).size(), 1);
  TEST_EQUAL(ids.at({0})[0].getHits().empty(), true);
  TEST_EQUAL(ids.at({0, 0})[0].getHits()[0].getMetaValue("label"), "ACUGp");
  map.clearIdentifications();
  IdentificationDataConverter::importMapIdentifications(map, ids);
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
START_SECTION(([EXTRA] maps: linked and unassigned identifications follow their peptide identifications))
{
  FeatureMap features;
  FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MQEvidence_3.featureXML"), features);
  checkLinks(features);
  ConsensusMap consensus;
  ConsensusXMLFile().load(OPENMS_GET_TEST_DATA_PATH("ExperimentalDesign_ProteomicsLFQ_1_subset_out.consensusXML"), consensus);
  checkLinks(consensus);
}
END_SECTION

START_SECTION((static MapPeptideIdentifications exportMapIdentifications(const FeatureMap& map) and importMapIdentifications()))
{
  // The view of a map is its identifications as peptide identifications; imported into the map without
  // identifications, it gives identification data and links whose export gives the view again.
  FeatureMap features;
  FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MQEvidence_3.featureXML"), features);
  const auto view = IdentificationDataConverter::exportMapIdentifications(features);
  Size assigned = 0;
  for (Size i = 0; i < features.size(); ++i)
  {
    TEST_EQUAL(view.at({i}).size(), features[i].getLinkedIdentifications(features.getIdentificationData()).size())
    assigned += view.at({i}).size();
  }
  TEST_NOT_EQUAL(assigned, 0)
  TEST_EQUAL(view.peptides().size(), assigned + view.unassigned.size())
  TEST_EQUAL(view.unassigned.size(), features.getUnassignedIdentifications().size())

  FeatureMap native = features;
  native.clearIdentifications();
  IdentificationDataConverter::importMapIdentifications(native, view);
  TEST_EQUAL(native.getIdentificationData().getRuns().size(), features.getIdentificationData().getRuns().size())
  for (Size i = 0; i < native.size(); ++i)
  {
    TEST_EQUAL(native[i].getLinkedIdentifications(native.getIdentificationData()).size(),
               features[i].getLinkedIdentifications(features.getIdentificationData()).size())
  }
  const auto exported = IdentificationDataConverter::exportMapIdentifications(native);
  TEST_EQUAL(exported.proteins == view.proteins, true)
  TEST_EQUAL(exported.features == view.features, true)
  TEST_EQUAL(exported.unassigned == view.unassigned, true)

  // a view refers to features by path; one that names no feature is rejected
  MapPeptideIdentifications dangling;
  dangling.proteins = view.proteins;
  dangling.features[{features.size()}] = view.at({0});
  TEST_EXCEPTION(Exception::InvalidParameter, IdentificationDataConverter::importMapIdentifications(native, dangling))

  // consensus maps
  ConsensusMap consensus;
  ConsensusXMLFile().load(OPENMS_GET_TEST_DATA_PATH("ExperimentalDesign_ProteomicsLFQ_1_subset_out.consensusXML"), consensus);
  const auto consensus_view = IdentificationDataConverter::exportMapIdentifications(consensus);
  ConsensusMap consensus_native = consensus;
  consensus_native.clearIdentifications();
  IdentificationDataConverter::importMapIdentifications(consensus_native, consensus_view);
  const auto consensus_exported = IdentificationDataConverter::exportMapIdentifications(consensus_native);
  TEST_EQUAL(consensus_exported.proteins == consensus_view.proteins, true)
  TEST_EQUAL(consensus_exported.features == consensus_view.features, true)
  TEST_EQUAL(consensus_exported.unassigned == consensus_view.unassigned, true)
}
END_SECTION
START_SECTION((MapPeptideIdentifications: at(), peptides() and eraseFeature()))
{
  MapPeptideIdentifications ids;
  PeptideIdentification a, b, c, d;
  a.setIdentifier("a");
  b.setIdentifier("b");
  c.setIdentifier("c");
  d.setIdentifier("d");
  ids.features[{0}] = {a};
  ids.features[{0, 1}] = {b};
  ids.features[{1}] = {c};
  ids.unassigned = {d};
  TEST_EQUAL(ids.at({0, 1})[0].getIdentifier(), "b")
  TEST_EQUAL(ids.at({2}).empty(), true)
  const auto all = ids.peptides();
  TEST_EQUAL(all.size(), 4)
  TEST_EQUAL(all[0]->getIdentifier() + all[1]->getIdentifier() + all[2]->getIdentifier() + all[3]->getIdentifier(), "abcd")
  ids.eraseFeature({0}); // with its subordinates
  TEST_EQUAL(ids.features.size(), 1)
  TEST_EQUAL(ids.at({1})[0].getIdentifier(), "c")
  TEST_EQUAL(ids.empty(), false)
}
END_SECTION
START_SECTION((static bool makeRunsDistinct(FeatureMap& map, std::set<std::string>& taken)))
{
  FeatureMap map;
  map.getIdentificationData() = fixture();
  auto& run = map.getIdentificationData().getRun("search");
  const auto old_uuid = run.getUuid();
  const auto query = run.getSources()[0].identifications[0].getId();
  const auto match = run.getSources()[0].identifications[0].getMatches()[0].getId();
  Feature feature, subordinate;
  feature.addIDMatch({old_uuid, match});
  subordinate.addIDQuery({old_uuid, query});
  feature.getSubordinates().push_back(subordinate);
  map.push_back(feature);
  ID::InferenceResult inference;
  inference.identifier = "inference";
  inference.inputs.push_back({"search", old_uuid, std::nullopt, ""});
  map.getIdentificationData().addInferenceResult(inference);
  const auto original = map.getIdentificationData();

  // runs that are not taken keep their UUIDs (and are taken then):
  std::set<std::string> taken;
  TEST_EQUAL(IdentificationDataConverter::makeRunsDistinct(map, taken), false)
  TEST_EQUAL(taken.size(), 1)
  TEST_EQUAL(taken.contains(old_uuid), true)
  TEST_EQUAL(map.getIdentificationData() == original, true)

  // a taken run gets a new UUID; links and inference inputs follow:
  TEST_EQUAL(IdentificationDataConverter::makeRunsDistinct(map, taken), true)
  const auto& renewed = map.getIdentificationData().getRun("search");
  TEST_NOT_EQUAL(renewed.getUuid(), old_uuid)
  TEST_EQUAL(taken.size(), 2)
  TEST_EQUAL(taken.contains(renewed.getUuid()), true)
  TEST_EQUAL(renewed.getIdentification(query).getMatches().size(), 1)
  TEST_EQUAL(*map[0].getIDMatches().begin() == (ID::MatchReference {renewed.getUuid(), match}), true)
  TEST_EQUAL(*map[0].getSubordinates()[0].getIDQueries().begin() == (ID::QueryReference {renewed.getUuid(), query}), true)
  TEST_EQUAL(map.getIdentificationData().getInferenceResults()[0].inputs[0].run_uuid, renewed.getUuid())
  TEST_EQUAL(map[0].getLinkedIdentifications(map.getIdentificationData()).size(), 1)
}
END_SECTION

START_SECTION((static bool makeRunsDistinct(ConsensusMap& map, std::set<std::string>& taken)))
{
  ConsensusMap map;
  map.getIdentificationData() = fixture();
  const auto& run = map.getIdentificationData().getRun("search");
  const auto old_uuid = run.getUuid();
  const auto query = run.getSources()[0].identifications[0].getId();
  ConsensusFeature feature;
  feature.addIDQuery({old_uuid, query});
  map.push_back(feature);
  std::set<std::string> taken {old_uuid};
  TEST_EQUAL(IdentificationDataConverter::makeRunsDistinct(map, taken), true)
  const auto& uuid = map.getIdentificationData().getRun("search").getUuid();
  TEST_NOT_EQUAL(uuid, old_uuid)
  TEST_EQUAL(*map[0].getIDQueries().begin() == (ID::QueryReference {uuid, query}), true)
}
END_SECTION

START_SECTION((static std::vector<ProteinIdentification> proteinIdentifications(const ConsensusMap& map)))
{
  ConsensusMap map;
  TEST_EQUAL(IdentificationDataConverter::proteinIdentifications(map).empty(), true)
  ProteinIdentification run;
  run.setIdentifier("search");
  run.setSearchEngine("engine");
  map.resize(1);
  PeptideIdentification peptide;
  peptide.setIdentifier("search");
  peptide.setScoreType("score");
  peptide.insertHit(PeptideHit(1.0, 1, 2, AASequence::fromString("PEPTIDE")));
  MapPeptideIdentifications ids;
  ids.proteins = {run};
  ids.features[{0}] = {peptide};
  IdentificationDataConverter::importMapIdentifications(map, ids);

  // the protein runs that export writes
  const auto proteins = IdentificationDataConverter::proteinIdentifications(map);
  ABORT_IF(proteins.size() != 1)
  TEST_EQUAL(proteins[0].getIdentifier(), "search")
  TEST_EQUAL(proteins[0].getSearchEngine(), "engine")
  TEST_EQUAL(proteins == IdentificationDataConverter::exportMapIdentifications(map).proteins, true)
}
END_SECTION

START_SECTION((static const std::vector<FeatureMap>& withDistinctRuns(const std::vector<FeatureMap>& maps, std::vector<FeatureMap>& converted)))
{
  // maps with distinct runs are used as they are:
  std::vector<FeatureMap> maps(2);
  maps[0].getIdentificationData() = fixture();
  maps[1].getIdentificationData() = fixture();
  std::vector<FeatureMap> converted;
  TEST_EQUAL(&IdentificationDataConverter::withDistinctRuns(maps, converted) == &maps, true)
  TEST_EQUAL(converted.empty(), true)
  // a run that is also in an earlier map gets a new UUID in a copy:
  maps[1] = maps[0];
  const auto& distinct = IdentificationDataConverter::withDistinctRuns(maps, converted);
  TEST_EQUAL(&distinct == &converted, true)
  ABORT_IF(converted.size() != 2)
  TEST_EQUAL(converted[0] == maps[0], true)
  TEST_NOT_EQUAL(converted[1].getIdentificationData().getRuns()[0].getUuid(), maps[0].getIdentificationData().getRuns()[0].getUuid())
  TEST_EQUAL(maps[1].getIdentificationData() == maps[0].getIdentificationData(), true)
}
END_SECTION

START_SECTION((static const std::vector<ConsensusMap>& withDistinctRuns(const std::vector<ConsensusMap>& maps, std::vector<ConsensusMap>& converted)))
{
  // the same map twice: the copy of the second gets new UUIDs for its runs
  std::vector<ConsensusMap> maps(2);
  ConsensusXMLFile().load(OPENMS_GET_TEST_DATA_PATH("ExperimentalDesign_ProteomicsLFQ_1_subset_out.consensusXML"), maps[0]);
  maps[1] = maps[0];
  std::vector<ConsensusMap> converted;
  const auto& prepared = IdentificationDataConverter::withDistinctRuns(maps, converted);
  TEST_EQUAL(&prepared == &converted, true)
  ABORT_IF(converted.size() != 2)
  TEST_EQUAL(converted[0] == maps[0], true)
  TEST_NOT_EQUAL(converted[1].getIdentificationData().getRuns()[0].getUuid(), maps[0].getIdentificationData().getRuns()[0].getUuid())
  // with the same identifications
  TEST_EQUAL(IdentificationDataConverter::exportMapIdentifications(converted[1]).features == IdentificationDataConverter::exportMapIdentifications(maps[1]).features, true)
  // maps with distinct runs are used as they are
  std::vector<ConsensusMap> again;
  TEST_EQUAL(&IdentificationDataConverter::withDistinctRuns(converted, again) == &converted, true)
  TEST_EQUAL(again.empty(), true)
}
END_SECTION

END_TEST
