// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
  template<class Map>
  void checkLinks(const Map& legacy)
  {
    Map native = legacy;
    TEST_TRUE(IdentificationDataConverter::moveToIdentificationData(native))
    TEST_TRUE(! IdentificationDataConverter::hasPeptideIdentifications(native))
    TEST_TRUE(! IdentificationDataConverter::moveToIdentificationData(native))
    const auto& data = native.getIdentificationData();
    Size best_matches = 0;
    for (Size i = 0; i < legacy.size(); ++i)
    {
      TEST_TRUE(sequences(legacy[i].getPeptideIdentifications()) == sequences(native[i].getLinkedIdentifications(data)))
      auto sorted = legacy[i];
      sorted.sortPeptideIdentifications();
      const auto best = native[i].getBestLinkedMatch(data);
      const bool has_hit = ! sorted.getPeptideIdentifications().empty() && ! sorted.getPeptideIdentifications()[0].getHits().empty();
      TEST_EQUAL(best.has_value(), has_hit)
      if (! best || ! has_hit) continue;
      ++best_matches;
      TEST_EQUAL(best->matches.size(), 1)
      TEST_EQUAL(best->matches[0]->representation, sorted.getPeptideIdentifications()[0].getHits()[0].getSequence().toString())
      TEST_REAL_SIMILAR(*best->run->getScore(best->matches[0]->getId(), *best->run->getPrimaryScore()), sorted.getPeptideIdentifications()[0].getHits()[0].getScore())
    }
    TEST_TRUE(best_matches > 0)
    TEST_TRUE(! legacy.getUnassignedPeptideIdentifications().empty())
    TEST_TRUE(sequences(legacy.getUnassignedPeptideIdentifications()) == sequences(native.getUnassignedIdentifications()))
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
START_SECTION(([EXTRA] maps: linked and unassigned identifications follow their peptide identifications))
{
  FeatureMap features;
  FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MQEvidence_3.featureXML"), features);
  checkLinks(features);
  ConsensusMap consensus;
  ConsensusXMLFile().load(OPENMS_GET_TEST_DATA_PATH("ExperimentalDesign_ProteomicsLFQ_1_subset_out.consensusXML"), consensus);
  checkLinks(consensus);

  // read access converts a copy, maps without peptide identifications are used as they are
  std::optional<ConsensusMap> converted;
  const auto& view = IdentificationDataConverter::withIdentificationData(consensus, converted);
  TEST_TRUE(converted.has_value() && &view == &*converted)
  TEST_TRUE(IdentificationDataConverter::hasPeptideIdentifications(consensus))
  TEST_TRUE(! view.getIdentificationData().empty())
  std::optional<ConsensusMap> again;
  TEST_TRUE(&IdentificationDataConverter::withIdentificationData(view, again) == &view && ! again)
  // both models at once are rejected
  ConsensusMap both = *converted;
  both.getUnassignedPeptideIdentifications() = consensus.getUnassignedPeptideIdentifications();
  TEST_EXCEPTION(Exception::InvalidParameter, IdentificationDataConverter::withIdentificationData(both, again))
  TEST_EXCEPTION(Exception::InvalidParameter, IdentificationDataConverter::moveToIdentificationData(both))
}
END_SECTION

START_SECTION((static void editAsIdentificationData(ConsensusMap& map, const std::function<void(ConsensusMap&)>& edit)))
{
  ConsensusMap consensus;
  ConsensusXMLFile().load(OPENMS_GET_TEST_DATA_PATH("ExperimentalDesign_ProteomicsLFQ_1_subset_out.consensusXML"), consensus);
  const auto original = consensus;
  // Callers keep references to protein runs (e.g. that of an inference result) across edits.
  const auto* run = &consensus.getProteinIdentifications()[0];
  bool native = false;
  IdentificationDataConverter::editAsIdentificationData(consensus, [&](ConsensusMap& map) {
    native = map.getUnassignedPeptideIdentifications().empty() && ! map.getIdentificationData().empty();
  });
  TEST_TRUE(native)
  TEST_EQUAL(&consensus.getProteinIdentifications()[0] == run, true)
  TEST_TRUE(consensus == original)
  // An edit that throws leaves the map with peptide identifications.
  TEST_EXCEPTION(Exception::InvalidValue, IdentificationDataConverter::editAsIdentificationData(consensus, [](ConsensusMap&) {
                   throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "failed", "");
                 }))
  TEST_EQUAL(&consensus.getProteinIdentifications()[0] == run, true)
  TEST_TRUE(consensus == original)
  // A map with identification data is edited as it is.
  auto native_map = original;
  IdentificationDataConverter::importConsensusIDs(native_map);
  const auto* data = &native_map.getIdentificationData();
  IdentificationDataConverter::editAsIdentificationData(native_map, [&](ConsensusMap& map) { TEST_EQUAL(&map.getIdentificationData() == data, true) });
  TEST_EQUAL(native_map.getUnassignedPeptideIdentifications().empty(), true)
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

START_SECTION((static const ConsensusMap& withPeptideIdentifications(const ConsensusMap& map, std::optional<ConsensusMap>& exported) and proteinIdentifications()))
{
  ConsensusMap legacy;
  ProteinIdentification run;
  run.setIdentifier("search");
  run.setSearchEngine("engine");
  legacy.setProteinIdentifications({run});
  legacy.resize(1);
  PeptideIdentification peptide;
  peptide.setIdentifier("search");
  peptide.setScoreType("score");
  peptide.insertHit(PeptideHit(1.0, 1, 2, AASequence::fromString("PEPTIDE")));
  legacy[0].getPeptideIdentifications().push_back(peptide);

  // a map with peptide identifications is its own view
  std::optional<ConsensusMap> exported;
  TEST_EQUAL(&IdentificationDataConverter::withPeptideIdentifications(legacy, exported), &legacy)
  TEST_EQUAL(exported.has_value(), false)
  TEST_EQUAL(IdentificationDataConverter::proteinIdentifications(legacy).size(), 1)

  // a map with identification data has a view with its identifications exported, and the protein runs export writes
  ConsensusMap native = legacy;
  IdentificationDataConverter::importConsensusIDs(native);
  const ConsensusMap& view = IdentificationDataConverter::withPeptideIdentifications(native, exported);
  TEST_EQUAL(exported.has_value(), true)
  TEST_EQUAL(view[0].getPeptideIdentifications().size(), 1)
  TEST_EQUAL(view.getIdentificationData().empty(), true)
  TEST_EQUAL(native[0].getPeptideIdentifications().empty(), true)
  const auto proteins = IdentificationDataConverter::proteinIdentifications(native);
  ABORT_IF(proteins.size() != 1)
  TEST_EQUAL(proteins[0].getIdentifier(), "search")
  TEST_EQUAL(proteins[0].getSearchEngine(), "engine")
}
END_SECTION

START_SECTION((static bool moveToPeptideIdentifications(FeatureMap& map) and (ConsensusMap& map)))
{
  FeatureMap legacy;
  ProteinIdentification run;
  run.setIdentifier("search");
  legacy.setProteinIdentifications({run});
  legacy.resize(1);
  PeptideIdentification peptide;
  peptide.setIdentifier("search");
  peptide.setScoreType("score");
  peptide.insertHit(PeptideHit(1.0, 1, 2, AASequence::fromString("PEPTIDE")));
  legacy[0].getPeptideIdentifications().push_back(peptide);
  PeptideIdentification unassigned = peptide;
  unassigned.setRT(10.0);
  legacy.getUnassignedPeptideIdentifications().push_back(unassigned);

  // a map with peptide identifications (or none at all) stays as it is
  FeatureMap unchanged = legacy;
  TEST_EQUAL(IdentificationDataConverter::moveToPeptideIdentifications(unchanged), false)
  TEST_EQUAL(unchanged == legacy, true)
  FeatureMap empty;
  TEST_EQUAL(IdentificationDataConverter::moveToPeptideIdentifications(empty), false)

  // a map with identification data gets them as peptide identifications again
  FeatureMap native = legacy;
  TEST_EQUAL(IdentificationDataConverter::moveToIdentificationData(native), true)
  TEST_EQUAL(IdentificationDataConverter::moveToPeptideIdentifications(native), true)
  TEST_EQUAL(native.getIdentificationData().empty(), true)
  TEST_EQUAL(native[0].getPeptideIdentifications().size(), 1)
  TEST_EQUAL(native.getUnassignedPeptideIdentifications().size(), 1)
  TEST_EQUAL(native.getProteinIdentifications().size(), 1)

  ConsensusMap consensus;
  consensus.setProteinIdentifications({run});
  consensus.resize(1);
  consensus[0].getPeptideIdentifications().push_back(peptide);
  TEST_EQUAL(IdentificationDataConverter::moveToPeptideIdentifications(consensus), false)
  TEST_EQUAL(IdentificationDataConverter::moveToIdentificationData(consensus), true)
  TEST_EQUAL(IdentificationDataConverter::moveToPeptideIdentifications(consensus), true)
  TEST_EQUAL(consensus[0].getPeptideIdentifications().size(), 1)
  TEST_EQUAL(consensus.getIdentificationData().empty(), true)
}
END_SECTION

START_SECTION((static const std::vector<FeatureMap>& withIdentificationData(const std::vector<FeatureMap>& maps, std::vector<FeatureMap>& converted)))
{
  // maps with distinct runs are used as they are:
  std::vector<FeatureMap> maps(2);
  maps[0].getIdentificationData() = fixture();
  maps[1].getIdentificationData() = fixture();
  std::vector<FeatureMap> converted;
  TEST_EQUAL(&IdentificationDataConverter::withIdentificationData(maps, converted) == &maps, true)
  TEST_EQUAL(converted.empty(), true)
  // a run that is also in an earlier map gets a new UUID in a copy:
  maps[1] = maps[0];
  const auto& distinct = IdentificationDataConverter::withIdentificationData(maps, converted);
  TEST_EQUAL(&distinct == &converted, true)
  ABORT_IF(converted.size() != 2)
  TEST_EQUAL(converted[0] == maps[0], true)
  TEST_NOT_EQUAL(converted[1].getIdentificationData().getRuns()[0].getUuid(), maps[0].getIdentificationData().getRuns()[0].getUuid())
  TEST_EQUAL(maps[1].getIdentificationData() == maps[0].getIdentificationData(), true)
}
END_SECTION

START_SECTION((static const std::vector<ConsensusMap>& withIdentificationData(const std::vector<ConsensusMap>& maps, std::vector<ConsensusMap>& converted)))
{
  // peptide identifications are converted:
  std::vector<ConsensusMap> maps(2);
  ConsensusXMLFile().load(OPENMS_GET_TEST_DATA_PATH("ExperimentalDesign_ProteomicsLFQ_1_subset_out.consensusXML"), maps[0]);
  std::vector<ConsensusMap> converted;
  const auto& prepared = IdentificationDataConverter::withIdentificationData(maps, converted);
  TEST_EQUAL(&prepared == &converted, true)
  ABORT_IF(converted.size() != 2)
  TEST_EQUAL(IdentificationDataConverter::hasPeptideIdentifications(converted[0]), false)
  TEST_EQUAL(converted[0].getIdentificationData().empty(), false)
  TEST_EQUAL(converted[1] == maps[1], true)
  TEST_EQUAL(IdentificationDataConverter::hasPeptideIdentifications(maps[0]), true)
}
END_SECTION

END_TEST
