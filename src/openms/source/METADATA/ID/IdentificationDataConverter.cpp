// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
#include <algorithm>
#include <set>
#include <cmath>
#include <limits>
#include <tuple>
#include <type_traits>
namespace OpenMS
{
namespace
{
  using ID = IdentificationData;
  using Adapter = IdentificationDataAdapter;
  [[noreturn]] void invalid(const std::string& message)
  { throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message); }
  Adapter::ExportOptions legacyOptions()
  {
    Adapter::ExportOptions options;
    options.loss_policy = Adapter::LossPolicy::ALLOW;
    return options;
  }
  void reportLegacyLosses(const std::vector<std::string>& losses)
  {
    for (const auto& message : losses)
      OPENMS_LOG_WARN << "IdentificationDataConverter: " << message << "; use native persistence to retain this information." << std::endl;
  }
  /// The database of a legacy run (which has at most one); an unknown database if its search names none.
  ID::DatabaseId legacyDatabase(ID::Run& run)
  { return run.getDatabases().empty() ? run.addDatabase(ID::Database {}) : run.getDatabaseId(0); }
  /// Remove a legacy meta value that a field holds and export writes back as @p restored (as IdentificationDataAdapter
  /// does); a value that export would write differently (other spelling or value type) is kept.
  void dropRestoredMetaValue(MetaInfoInterface& metadata, const std::string& name, const std::string& restored)
  {
    if (restored.empty() || ! metadata.metaValueExists(name)) return;
    const auto& value = metadata.getMetaValue(name);
    if (value.valueType() == DataValue::STRING_VALUE && value.toString() == restored) metadata.removeMetaValue(name);
  }
  std::string targetDecoyText(ID::TargetDecoy value)
  {
    return value == ID::TargetDecoy::TARGET ? "target" : value == ID::TargetDecoy::DECOY ? "decoy" : value == ID::TargetDecoy::BOTH ? "target+decoy" : "";
  }

  Adapter::ImportResult importGeneric(const std::vector<ProteinIdentification>& proteins, const PeptideIdentificationList& peptides)
  {
    ID imported;
    std::vector<ID::QueryReference> imported_queries;
    std::map<std::string, const ProteinIdentification*> processing;
    for (const auto& protein : proteins)
      if (! processing.emplace(protein.getIdentifier(), &protein).second) invalid("Duplicate input run identifier");
    std::map<std::pair<std::string, ID::MoleculeKind>, ID::Run*> runs;
    const auto hit_kind = [](const PeptideHit& hit) {
      if (! hit.metaValueExists("molecule_type")) return ID::MoleculeKind::PEPTIDE;
      const auto value = hit.getMetaValue("molecule_type").toString();
      if (value == "RNA") return ID::MoleculeKind::OLIGONUCLEOTIDE;
      if (value == "compound") return ID::MoleculeKind::COMPOUND;
      invalid("Unsupported idXML molecule_type: " + value);
    };
    for (const auto& item : peptides)
    {
      const auto original = processing.find(item.getIdentifier());
      if (original == processing.end()) invalid("Missing protein run for identification");
      auto kind = item.getHits().empty() ? ID::MoleculeKind::PEPTIDE : hit_kind(item.getHits().front());
      if (item.getHits().empty())
        for (const auto& candidate : peptides)
          if (candidate.getIdentifier() == item.getIdentifier() && ! candidate.getHits().empty())
          {
            kind = hit_kind(candidate.getHits().front());
            break;
          }
      if (std::any_of(item.getHits().begin(), item.getHits().end(), [&](const auto& hit) { return hit_kind(hit) != kind; }))
        invalid("Mixed molecular kinds in one idXML identification");
      const auto key = std::make_pair(item.getIdentifier(), kind);
      auto found = runs.find(key);
      if (found == runs.end())
      {
        auto name = item.getIdentifier();
        if (std::any_of(imported.getRuns().begin(), imported.getRuns().end(), [&](const auto& run) { return run.getIdentifier() == name; }))
          name += ":kind_" + std::to_string(static_cast<int>(kind));
        auto& run = imported.addRun(name, kind);
        auto metadata = *original->second;
        metadata.setIdentifier(name);
        run.setSettings(Adapter::settingsFromLegacy(metadata));
        // The files of the legacy run become the sources of the run.
        StringList files;
        metadata.getPrimaryMSRunPath(files);
        Adapter::addLegacySources(run, files);
        // The database of the legacy run's search becomes the database of the run.
        const auto database = Adapter::databaseFromLegacy(metadata.getSearchParameters());
        if (database != ID::Database {}) run.addDatabase(database);
        std::vector<ID::DatabaseSequence> sequences;
        // Database sequences need distinct accessions. A protein list with empty or repeated accessions (which legacy
        // runs allow) has none, as in IdentificationDataAdapter; its inference result keeps the protein hits.
        std::set<std::string> accessions;
        const bool catalogue = std::all_of(metadata.getHits().begin(), metadata.getHits().end(), [&](const ProteinHit& hit) {
          return ! hit.getAccession().empty() && accessions.insert(hit.getAccession()).second;
        });
        for (const auto& hit : catalogue ? metadata.getHits() : std::vector<ProteinHit> {})
        {
          ID::DatabaseSequence sequence;
          static_cast<MetaInfoInterface&>(sequence) = hit;
          sequence.database = legacyDatabase(run);
          sequence.accession = hit.getAccession();
          sequence.sequence = hit.getSequence();
          sequence.description = hit.getDescription();
          dropRestoredMetaValue(sequence, "Description", sequence.description);
          if (hit.getCoverage() >= 0) sequence.setMetaValue("coverage", hit.getCoverage() / 100.0);
          sequence.target_decoy = hit.getTargetDecoyType() == ProteinHit::TargetDecoyType::DECOY    ? ID::TargetDecoy::DECOY
                                  : hit.getTargetDecoyType() == ProteinHit::TargetDecoyType::TARGET ? ID::TargetDecoy::TARGET
                                                                                                    : ID::TargetDecoy::UNKNOWN;
          dropRestoredMetaValue(sequence, "target_decoy", targetDecoyText(sequence.target_decoy));
          sequences.push_back(std::move(sequence));
        }
        if (! sequences.empty()) run.setDatabaseSequences(std::move(sequences));
        if (! item.getScoreType().empty())
        {
          ID::ScoreDefinition definition;
          definition.name = item.getScoreType();
          definition.higher_better = item.isHigherScoreBetter();
          std::tie(definition.software, definition.software_version) = metadata.getScoreSoftware(definition.name);
          run.setPrimaryScore(run.addScore(definition));
        }
        found = runs.emplace(key, &run).first;
      }
      auto& run = *found->second;
      if (run.getPrimaryScore()
          && (run.getScoreDefinitions()[0].name != item.getScoreType() || run.getScoreDefinitions()[0].higher_better != item.isHigherScoreBetter()))
        invalid("Inconsistent idXML PSM score contract");
      ID::Observation observation;
      static_cast<MetaInfoInterface&>(observation) = item;
      // The source is the file, so the index into the legacy file list is not kept.
      observation.removeMetaValue("id_merge_index");
      observation.data_id = item.getSpectrumReference();
      dropRestoredMetaValue(observation, Constants::UserParam::SPECTRUM_REFERENCE, observation.data_id);
      if (item.hasRT()) observation.rt = item.getRT();
      if (item.hasMZ()) observation.mz = item.getMZ();
      const auto source = Adapter::legacySource(run, original->second->nrPrimaryMSRunPaths(), item);
      const auto query = run.addIdentification(source, observation);
      imported_queries.push_back({run.getUuid(), query});
      for (const auto& hit : item.getHits())
      {
        ID::MatchData match;
        static_cast<MetaInfoInterface&>(match) = hit;
        match.charge = hit.getCharge();
        match.peak_annotations = hit.getPeakAnnotations();
        if (kind == ID::MoleculeKind::PEPTIDE) match.representation = hit.getSequence().toString();
        else
        {
          if (! hit.metaValueExists("label")) invalid("RNA/compound idXML identification has no label");
          match.representation = hit.getMetaValue("label").toString();
          match.encoding = kind == ID::MoleculeKind::OLIGONUCLEOTIDE ? ID::Encoding::NA_SEQUENCE : ID::Encoding::DATABASE_ID;
          if (hit.metaValueExists("identification:encoding"))
            match.encoding = static_cast<ID::Encoding>(static_cast<int>(hit.getMetaValue("identification:encoding")));
        }
        if (hit.metaValueExists("identification:identifier_databases"))
        {
          const auto databases = hit.getMetaValue("identification:identifier_databases").toStringList();
          const auto accessions = hit.getMetaValue("identification:identifier_accessions").toStringList();
          if (databases.size() != accessions.size()) invalid("Mismatched compound identifier columns");
          for (Size i = 0; i < databases.size(); ++i)
            match.details.emplace().identifiers.push_back({databases[i], accessions[i]});
        }
        if (hit.metaValueExists("identification:formula")) match.details.emplace().formula = hit.getMetaValue("identification:formula").toString();
        if (hit.metaValueExists("identification:name")) match.details.emplace().name = hit.getMetaValue("identification:name").toString();
        if (hit.metaValueExists("identification:calculated_mz"))
          match.calculated_mz = static_cast<double>(hit.getMetaValue("identification:calculated_mz"));
        if (hit.metaValueExists("identification:adduct_formula"))
          match.details.emplace().adduct
            = AdductInfo(hit.getMetaValue("adduct").toString(), EmpiricalFormula(hit.getMetaValue("identification:adduct_formula").toString()),
                         match.charge, static_cast<int>(hit.getMetaValue("identification:adduct_multiplier")));
        if (hit.getTargetDecoyType() == PeptideHit::TargetDecoyType::TARGET) match.target_decoy = ID::TargetDecoy::TARGET;
        else if (hit.getTargetDecoyType() == PeptideHit::TargetDecoyType::DECOY)
          match.target_decoy = ID::TargetDecoy::DECOY;
        else if (hit.getTargetDecoyType() == PeptideHit::TargetDecoyType::TARGET_DECOY)
          match.target_decoy = ID::TargetDecoy::BOTH;
        dropRestoredMetaValue(match, "target_decoy", targetDecoyText(match.target_decoy));
        for (const auto& old : hit.getPeptideEvidences())
        {
          ID::SequenceEvidence evidence;
          evidence.database = legacyDatabase(run);
          evidence.accession = old.getProteinAccession();
          if (old.getStart() >= 0) evidence.start = old.getStart();
          if (old.getEnd() >= 0) evidence.end = old.getEnd();
          evidence.before = old.getAABefore();
          evidence.after = old.getAAAfter();
          match.sequence_evidence.push_back(std::move(evidence));
        }
        run.addMatch(query, match,
                     run.getPrimaryScore() ? std::vector<std::optional<double>> {hit.getScore()} : std::vector<std::optional<double>> {});
      }
    }
    // Reuse the checked protein adapter for inference and protein-only input runs.
    auto protein_data = Adapter::fromLegacy(proteins, {});
    for (const auto& protein : proteins)
      if (std::none_of(runs.begin(), runs.end(), [&](const auto& entry) { return entry.first.first == protein.getIdentifier(); }))
        imported.merge(Adapter::fromLegacy({protein}, {}));
    for (auto inference : protein_data.getInferenceResults())
    {
      // These runs carry no database sequences: keep the complete protein hits in the inference result.
      inference.proteins.setHits(Adapter::proteinHits(protein_data, inference));
      std::vector<ID::InferenceInput> inputs;
      for (const auto& input : inference.inputs)
      {
        bool found = false;
        for (const auto& [key, run] : runs)
          if (key.first == input.run_identifier)
          {
            auto replacement = input;
            replacement.run_identifier = run->getIdentifier();
            replacement.run_uuid = run->getUuid();
            inputs.push_back(std::move(replacement));
            found = true;
          }
        if (! found) continue; // Already included by the protein-only merge above.
      }
      if (! inputs.empty())
      {
        inference.inputs = std::move(inputs);
        imported.addInferenceResult(std::move(inference));
      }
    }
    // Queries and candidates were appended one by one.
    for (const auto& [key, run] : runs)
      run->shrinkToFit();
    imported.validate();
    return {std::move(imported), std::move(imported_queries)};
  }
  template<class Map>
  void importMap(Map& map, bool clear_original)
  {
    auto converted = [&]() {
      if (map.getIdentificationData().empty())
      {
        PeptideIdentificationList queries;
        std::vector<Adapter::FeatureAssociation> associations;
        const auto collect = [&](const auto& self, const auto& feature, std::vector<Size> path) -> void {
          for (const auto& query : feature.getPeptideIdentifications())
          {
            queries.push_back(query);
            Adapter::FeatureAssociation association;
            association.feature_path = path;
            associations.push_back(std::move(association));
          }
          if constexpr (std::is_same_v<std::remove_cvref_t<decltype(feature)>, Feature>)
            for (Size i = 0; i < feature.getSubordinates().size(); ++i)
            {
              auto child = path;
              child.push_back(i);
              self(self, feature.getSubordinates()[i], std::move(child));
            }
        };
        for (Size i = 0; i < map.size(); ++i)
          collect(collect, map[i], {i});
        for (const auto& query : map.getUnassignedPeptideIdentifications())
        {
          queries.push_back(query);
          Adapter::FeatureAssociation association;
          association.unassigned = true;
          associations.push_back(std::move(association));
        }
        const bool generic = std::any_of(queries.begin(), queries.end(), [](const auto& query) {
          return std::any_of(query.getHits().begin(), query.getHits().end(), [](const auto& hit) { return hit.metaValueExists("molecule_type"); });
        });
        if (generic)
        {
          auto imported = importGeneric(map.getProteinIdentifications(), queries);
          for (Size i = 0; i < associations.size(); ++i)
          {
            auto& association = associations[i];
            association.query = imported.queries.at(i);
            const auto& run = *imported.data.findRunByUuid(association.query.run_uuid);
            for (const auto& match : run.getIdentification(association.query.query).getMatches())
              association.matches.push_back(match.getId());
          }
          return Adapter::FeatureImportResult {std::move(imported.data), std::move(associations)};
        }
      }
      if constexpr (std::is_same_v<Map, FeatureMap>) return Adapter::fromFeatureMap(map);
      else
        return Adapter::fromConsensusMap(map);
    }();
    auto updated = map;
    updated.getIdentificationData() = std::move(converted.data);
    const auto clear = [&](const auto& self, auto& feature) -> void {
      feature.getIDMatches().clear();
      feature.getIDQueries().clear();
      if (clear_original) feature.getPeptideIdentifications().clear();
      if constexpr (std::is_same_v<std::remove_cvref_t<decltype(feature)>, Feature>)
        for (auto& child : feature.getSubordinates())
          self(self, child);
    };
    for (auto& feature : updated)
      clear(clear, feature);
    for (const auto& association : converted.associations)
    {
      if (association.unassigned) continue;
      if (association.feature_path.empty() || association.feature_path[0] >= updated.size()) invalid("Invalid feature association path");
      auto* feature = &updated[association.feature_path[0]];
      if constexpr (std::is_same_v<Map, FeatureMap>)
        for (Size i = 1; i < association.feature_path.size(); ++i)
          feature = &feature->getSubordinates().at(association.feature_path[i]);
      else if (association.feature_path.size() != 1)
        invalid("Invalid consensus association path");
      feature->addIDQuery(association.query);
      for (auto match : association.matches)
        feature->addIDMatch({association.query.run_uuid, match});
    }
    if (clear_original)
    {
      updated.getProteinIdentifications().clear();
      updated.getUnassignedPeptideIdentifications().clear();
    }
    map = std::move(updated);
  }
  template<class Map>
  void exportMap(Map& map, bool clear_original)
  {
    auto converted = [&]() {
      if constexpr (std::is_same_v<Map, FeatureMap>) return Adapter::fromFeatureMap(map);
      else
        return Adapter::fromConsensusMap(map);
    }();
    auto updated = map;
    const bool generic = std::any_of(converted.data.getRuns().begin(), converted.data.getRuns().end(),
                                     [](const auto& run) { return run.getMoleculeKind() != ID::MoleculeKind::PEPTIDE; });
    if (generic)
    {
      std::vector<ProteinIdentification> proteins;
      PeptideIdentificationList queries;
      IdentificationDataConverter::exportIDs(converted.data, proteins, queries, true);
      std::map<ID::QueryReference, Size> indices;
      Size index = 0;
      for (const auto& run : converted.data.getRuns())
        for (const auto& source : run.getSources())
          for (const auto& query : source.identifications)
            indices[{run.getUuid(), query.getId()}] = index++;
      const auto clear = [&](const auto& self, auto& feature) -> void {
        feature.getPeptideIdentifications().clear();
        if constexpr (std::is_same_v<std::remove_cvref_t<decltype(feature)>, Feature>)
          for (auto& child : feature.getSubordinates())
            self(self, child);
      };
      for (auto& feature : updated)
        clear(clear, feature);
      updated.getProteinIdentifications() = std::move(proteins);
      updated.getUnassignedPeptideIdentifications().clear();
      for (const auto& association : converted.associations)
      {
        auto query = queries.at(indices.at(association.query));
        const auto& run = *converted.data.findRunByUuid(association.query.run_uuid);
        const auto& matches = run.getIdentification(association.query.query).getMatches();
        std::vector<PeptideHit> hits;
        for (Size i = 0; i < matches.size(); ++i)
          if (std::find(association.matches.begin(), association.matches.end(), matches[i].getId()) != association.matches.end())
            hits.push_back(query.getHits().at(i));
        query.setHits(hits);
        if (association.unassigned) updated.getUnassignedPeptideIdentifications().push_back(std::move(query));
        else
        {
          auto* feature = &updated.at(association.feature_path.at(0));
          if constexpr (std::is_same_v<Map, FeatureMap>)
            for (Size i = 1; i < association.feature_path.size(); ++i)
              feature = &feature->getSubordinates().at(association.feature_path[i]);
          feature->getPeptideIdentifications().push_back(std::move(query));
        }
      }
    }
    else if constexpr (std::is_same_v<Map, FeatureMap>)
      reportLegacyLosses(
        Adapter::applyToFeatureMap(converted.data, converted.associations, updated, legacyOptions(), Adapter::MissingLinkPolicy::REJECT));
    else
      reportLegacyLosses(
        Adapter::applyToConsensusMap(converted.data, converted.associations, updated, legacyOptions(), Adapter::MissingLinkPolicy::REJECT));
    if (clear_original)
    {
      const auto clear = [&](const auto& self, auto& feature) -> void {
        feature.getIDMatches().clear();
        feature.getIDQueries().clear();
        if constexpr (std::is_same_v<std::remove_cvref_t<decltype(feature)>, Feature>)
          for (auto& child : feature.getSubordinates())
            self(self, child);
      };
      for (auto& feature : updated)
        clear(clear, feature);
      updated.getIdentificationData().clear();
    }
    map = std::move(updated);
  }
} // namespace
void IdentificationDataConverter::importIDs(ID& data, const std::vector<ProteinIdentification>& proteins, const PeptideIdentificationList& peptides)
{
  const bool generic = std::any_of(peptides.begin(), peptides.end(), [](const auto& query) {
    return std::any_of(query.getHits().begin(), query.getHits().end(), [](const auto& hit) { return hit.metaValueExists("molecule_type"); });
  });
  if (! generic)
  {
    data.merge(Adapter::fromLegacy(proteins, peptides));
    return;
  }
  data.merge(importGeneric(proteins, peptides).data);
}

void IdentificationDataConverter::exportIDs(const ID& data,
                                            std::vector<ProteinIdentification>& proteins,
                                            PeptideIdentificationList& peptides,
                                            bool export_ids_wo_scores)
{
  data.validate();
  if (std::all_of(data.getRuns().begin(), data.getRuns().end(), [](const auto& run) { return run.getMoleculeKind() == ID::MoleculeKind::PEPTIDE; }))
  {
    auto converted = Adapter::toLegacy(data, legacyOptions());
    reportLegacyLosses(converted.losses);
    proteins.insert(proteins.end(), converted.proteins.begin(), converted.proteins.end());
    peptides.insert(peptides.end(), converted.peptides.begin(), converted.peptides.end());
    return;
  }
  // idXML has no RNA/compound sequence type. Preserve the historical label convention
  // explicitly; native persistence retains the full molecular and ion representation.
  std::vector<ProteinIdentification> added_proteins;
  PeptideIdentificationList added_peptides;
  for (const auto& run : data.getRuns())
  {
    auto processing = Adapter::settingsToLegacy(run);
    processing.setIdentifier(run.getIdentifier());
    if (run.getDatabases().size() > 1) reportLegacyLosses({"Legacy export can name only one database per run: " + run.getIdentifier()});
    // Without an inference result, the legacy run takes the primary score as its score type, as search engines write it.
    if (run.getPrimaryScore())
    {
      processing.setScoreType(run.getScoreDefinition(*run.getPrimaryScore()).name);
      processing.setHigherScoreBetter(run.getScoreDefinition(*run.getPrimaryScore()).higher_better);
    }
    std::vector<ProteinHit> database_sequences;
    if (run.getDatabaseSequences())
      for (const auto& sequence : *run.getDatabaseSequences())
      {
        ProteinHit hit;
        static_cast<MetaInfoInterface&>(hit) = sequence;
        hit.setAccession(sequence.accession);
        hit.setSequence(sequence.sequence);
        // ProteinHit keeps its description as metadata, which is absent unless set.
        if (! sequence.description.empty() || hit.metaValueExists("Description")) hit.setDescription(sequence.description);
        if (sequence.metaValueExists("coverage"))
        {
          // the coverage attribute represents it
          hit.setCoverage(static_cast<double>(sequence.getMetaValue("coverage")) * 100.0);
          hit.removeMetaValue("coverage");
        }
        hit.setTargetDecoyType(sequence.target_decoy == ID::TargetDecoy::DECOY    ? ProteinHit::TargetDecoyType::DECOY
                               : sequence.target_decoy == ID::TargetDecoy::TARGET ? ProteinHit::TargetDecoyType::TARGET
                                                                                  : ProteinHit::TargetDecoyType::UNKNOWN);
        database_sequences.push_back(std::move(hit));
      }
    processing.setHits(database_sequences);
    for (const auto& inference : data.getInferenceResults())
      if (std::any_of(inference.inputs.begin(), inference.inputs.end(), [&](const auto& input) { return input.run_uuid == run.getUuid(); }))
      {
        processing.setHits(Adapter::proteinHits(data, inference));
        processing.setScoreType(inference.proteins.getScoreType());
        processing.setHigherScoreBetter(inference.proteins.isHigherScoreBetter());
        processing.getProteinGroups() = inference.proteins.getProteinGroups();
        processing.getIndistinguishableProteins() = inference.proteins.getIndistinguishableProteins();
      }
    const auto files = Adapter::legacyFiles(run);
    if (! files.empty()) processing.setPrimaryMSRunPath(files);
    if (run.getPrimaryScore())
    {
      // Record a producer other than the search engine so the score definition survives the round trip.
      const auto& definition = run.getScoreDefinition(*run.getPrimaryScore());
      if (! definition.name.empty() && ! definition.software.empty()
          && std::make_pair(definition.software, definition.software_version) != processing.getScoreSoftware(definition.name))
        processing.setScoreSoftware(definition.name, definition.software, definition.software_version);
    }
    added_proteins.push_back(std::move(processing));
    // A source with a path is the next file of the legacy file list (see IdentificationDataAdapter::legacyFiles).
    Size file_index = 0;
    for (const auto& source : run.getSources())
    {
      const bool known = ! source.file.path.empty();
      for (const auto& query : source.identifications)
      {
        PeptideIdentification item;
        static_cast<MetaInfoInterface&>(item) = query;
        item.setIdentifier(run.getIdentifier());
        // An empty data ID is no spectrum reference (as in IdentificationDataAdapter).
        if (! query.data_id.empty()) item.setSpectrumReference(query.data_id);
        item.removeMetaValue("id_merge_index");
        if (known && files.size() > 1) item.setMetaValue("id_merge_index", static_cast<Int64>(file_index));
        if (query.rt) item.setRT(*query.rt);
        if (query.mz) item.setMZ(*query.mz);
        if (run.getPrimaryScore())
        {
          const auto& definition = run.getScoreDefinition(*run.getPrimaryScore());
          item.setScoreType(definition.name);
          item.setHigherScoreBetter(definition.higher_better);
        }
        for (const auto& match : query.getMatches())
        {
          const auto scores = run.getScores(match);
          const auto score = run.getPrimaryScore() ? scores.at(run.getPrimaryScore()->value) : std::nullopt;
          if (! score && ! export_ids_wo_scores) continue;
          PeptideHit hit;
          static_cast<MetaInfoInterface&>(hit) = match;
          if (match.encoding == ID::Encoding::AA_SEQUENCE) hit.setSequence(AASequence::fromString(match.representation));
          else
          {
            hit.setMetaValue("label", match.representation);
            hit.setMetaValue("molecule_type", run.getMoleculeKind() == ID::MoleculeKind::OLIGONUCLEOTIDE ? "RNA" : "compound");
          }
          hit.setCharge(match.charge);
          hit.setScore(score.value_or(0));
          hit.setPeakAnnotations(match.peak_annotations);
          exportSequenceEvidence(match.sequence_evidence, hit);
          const auto& details = match.details.value_or_default();
          if (details.adduct)
          {
            hit.setMetaValue("adduct", details.adduct->getName());
            hit.setMetaValue("identification:adduct_formula", details.adduct->getEmpiricalFormula().toString());
            hit.setMetaValue("identification:adduct_multiplier", static_cast<int>(details.adduct->getMolMultiplier()));
          }
          hit.setMetaValue("identification:encoding", static_cast<int>(match.encoding));
          if (! details.identifiers.empty())
          {
            StringList databases, accessions;
            for (const auto& identity : details.identifiers)
            {
              databases.push_back(identity.database);
              accessions.push_back(identity.accession);
            }
            hit.setMetaValue("identification:identifier_databases", databases);
            hit.setMetaValue("identification:identifier_accessions", accessions);
          }
          if (details.formula) hit.setMetaValue("identification:formula", *details.formula);
          if (! details.name.empty()) hit.setMetaValue("identification:name", details.name);
          if (match.calculated_mz) hit.setMetaValue("identification:calculated_mz", *match.calculated_mz);
          if (match.target_decoy != ID::TargetDecoy::UNKNOWN)
            hit.setMetaValue("target_decoy", match.target_decoy == ID::TargetDecoy::DECOY  ? "decoy"
                                             : match.target_decoy == ID::TargetDecoy::BOTH ? "target+decoy"
                                                                                           : "target");
          for (Size i = 0; i < run.getScoreDefinitions().size() && i < scores.size(); ++i)
            if (scores[i] && (! run.getPrimaryScore() || i != run.getPrimaryScore()->value))
              hit.setMetaValue(run.getScoreDefinitions()[i].name, *scores[i]);
          item.insertHit(hit);
        }
        if (! item.getHits().empty() || query.getMatches().empty()) added_peptides.push_back(std::move(item));
      }
      if (known) ++file_index;
    }
  }
  proteins.insert(proteins.end(), added_proteins.begin(), added_proteins.end());
  peptides.insert(peptides.end(), added_peptides.begin(), added_peptides.end());
}

ID::DatabaseId IdentificationDataConverter::importSequences(ID::Run& run, const ID::Database& database, const std::vector<FASTAFile::FASTAEntry>& fasta,
                                                             const std::string& decoy_pattern)
{
  const auto database_id = run.addDatabase(database);
  auto sequences = run.getDatabaseSequences().value_or(std::vector<ID::DatabaseSequence> {});
  for (const auto& entry : fasta)
  {
    ID::DatabaseSequence sequence;
    sequence.database = database_id;
    sequence.accession = entry.identifier;
    sequence.sequence = entry.sequence;
    sequence.description = entry.description;
    sequence.target_decoy
      = ! decoy_pattern.empty() && entry.identifier.find(decoy_pattern) != std::string::npos ? ID::TargetDecoy::DECOY : ID::TargetDecoy::TARGET;
    sequences.push_back(std::move(sequence));
  }
  run.setDatabaseSequences(std::move(sequences));
  return database_id;
}
void IdentificationDataConverter::exportSequenceEvidence(const std::vector<ID::SequenceEvidence>& sequence_evidence, PeptideHit& hit)
{
  std::vector<PeptideEvidence> evidence;
  for (const auto& item : sequence_evidence)
  {
    if ((item.start && *item.start > static_cast<UInt32>(std::numeric_limits<Int>::max()))
        || (item.end && *item.end > static_cast<UInt32>(std::numeric_limits<Int>::max())))
      invalid("Sequence evidence exceeds legacy coordinate limits");
    evidence.emplace_back(item.accession, item.start ? static_cast<Int>(*item.start) : PeptideEvidence::UNKNOWN_POSITION,
                          item.end ? static_cast<Int>(*item.end) : PeptideEvidence::UNKNOWN_POSITION,
                          item.before ? item.before : PeptideEvidence::UNKNOWN_AA, item.after ? item.after : PeptideEvidence::UNKNOWN_AA);
  }
  hit.setPeptideEvidences(evidence);
}
void IdentificationDataConverter::importFeatureIDs(FeatureMap& map, bool clear)
{ importMap(map, clear); }
void IdentificationDataConverter::exportFeatureIDs(FeatureMap& map, bool clear)
{ exportMap(map, clear); }
void IdentificationDataConverter::importConsensusIDs(ConsensusMap& map, bool clear)
{ importMap(map, clear); }
void IdentificationDataConverter::exportConsensusIDs(ConsensusMap& map, bool clear)
{ exportMap(map, clear); }

namespace
{
  template<class Map>
  bool hasLegacyIDs(const Map& map)
  {
    if (! map.getProteinIdentifications().empty() || ! map.getUnassignedPeptideIdentifications().empty()) return true;
    const auto any = [](const auto& self, const auto& feature) -> bool {
      if (! feature.getPeptideIdentifications().empty()) return true;
      if constexpr (std::is_same_v<std::remove_cvref_t<decltype(feature)>, Feature>)
        for (const auto& subordinate : feature.getSubordinates())
          if (self(self, subordinate)) return true;
      return false;
    };
    return std::any_of(map.begin(), map.end(), [&](const auto& feature) { return any(any, feature); });
  }
  template<class Map>
  bool checkedLegacyIDs(const Map& map)
  {
    if (! hasLegacyIDs(map)) return false;
    if (! map.getIdentificationData().empty())
      invalid("The map has peptide identifications and identification data; it can hold its identifications in one of them only");
    return true;
  }
} // namespace

bool IdentificationDataConverter::hasPeptideIdentifications(const FeatureMap& map)
{ return hasLegacyIDs(map); }
bool IdentificationDataConverter::hasPeptideIdentifications(const ConsensusMap& map)
{ return hasLegacyIDs(map); }
const FeatureMap& IdentificationDataConverter::withIdentificationData(const FeatureMap& map, std::optional<FeatureMap>& converted)
{
  if (! checkedLegacyIDs(map)) return map;
  converted = map;
  importMap(*converted, true);
  return *converted;
}
const ConsensusMap& IdentificationDataConverter::withIdentificationData(const ConsensusMap& map, std::optional<ConsensusMap>& converted)
{
  if (! checkedLegacyIDs(map)) return map;
  converted = map;
  importMap(*converted, true);
  return *converted;
}
namespace
{
  template<class Map>
  void editMap(Map& map, const std::function<void(Map&)>& edit)
  {
    if (! checkedLegacyIDs(map))
    {
      edit(map);
      return;
    }
    // The original protein runs (their storage) come back with the exported values.
    std::vector<ProteinIdentification> runs;
    runs.swap(map.getProteinIdentifications());
    map.getProteinIdentifications() = runs;
    try
    {
      importMap(map, true);
    }
    catch (...)
    {
      map.getProteinIdentifications().swap(runs);
      throw;
    }
    const auto restore = [&] {
      exportMap(map, true);
      auto& exported = map.getProteinIdentifications();
      if (exported.size() != runs.size()) return;
      for (Size i = 0; i < runs.size(); ++i)
        runs[i] = std::move(exported[i]);
      exported.swap(runs);
    };
    try
    {
      edit(map);
    }
    catch (...)
    {
      restore();
      throw;
    }
    restore();
  }
} // namespace

void IdentificationDataConverter::editAsIdentificationData(FeatureMap& map, const std::function<void(FeatureMap&)>& edit)
{ editMap(map, edit); }
void IdentificationDataConverter::editAsIdentificationData(ConsensusMap& map, const std::function<void(ConsensusMap&)>& edit)
{ editMap(map, edit); }

namespace
{
  template<class Map>
  bool distinctRuns(Map& map, std::set<std::string>& taken)
  {
    auto& data = map.getIdentificationData();
    std::set<std::string> own;
    for (const auto& run : data.getRuns())
      own.insert(run.getUuid());
    std::map<std::string, std::string> renewed;
    for (const auto& uuid : own)
    {
      if (! taken.contains(uuid)) continue;
      auto fresh = UniqueIdGenerator::getUUID();
      while (taken.contains(fresh) || own.contains(fresh))
        fresh = UniqueIdGenerator::getUUID();
      renewed.emplace(uuid, fresh);
      own.insert(fresh);
    }
    if (renewed.empty())
    {
      taken.insert(own.begin(), own.end());
      return false;
    }
    const auto uuid = [&](const std::string& current) -> const std::string& {
      const auto found = renewed.find(current);
      return found == renewed.end() ? current : found->second;
    };
    ID result;
    for (const auto& run : data.getRuns())
    {
      auto copy = run;
      if (renewed.contains(run.getUuid())) copy.restoreIdentity(uuid(run.getUuid()), copy.getNextQueryId(), copy.getNextMatchId());
      taken.insert(copy.getUuid());
      result.addRun(std::move(copy));
    }
    for (auto inference : data.getInferenceResults())
    {
      for (auto& input : inference.inputs)
        input.run_uuid = uuid(input.run_uuid);
      result.addInferenceResult(std::move(inference));
    }
    data = std::move(result);
    const auto relink = [&](const auto& self, auto& feature) -> void {
      std::set<ID::QueryReference> queries;
      for (const auto& query : feature.getIDQueries())
        queries.insert({uuid(query.run_uuid), query.query});
      feature.getIDQueries() = std::move(queries);
      std::set<ID::MatchReference> matches;
      for (const auto& match : feature.getIDMatches())
        matches.insert({uuid(match.run_uuid), match.match});
      feature.getIDMatches() = std::move(matches);
      if constexpr (std::is_same_v<std::remove_cvref_t<decltype(feature)>, Feature>)
        for (auto& subordinate : feature.getSubordinates())
          self(self, subordinate);
    };
    for (auto& feature : map)
      relink(relink, feature);
    return true;
  }
  template<class Map>
  const std::vector<Map>& distinctMaps(const std::vector<Map>& maps, std::vector<Map>& converted)
  {
    // Decide first, so that maps are only copied if one changes.
    std::set<std::string> taken;
    bool change = false;
    for (const auto& map : maps)
    {
      if (checkedLegacyIDs(map)) change = true;
      for (const auto& run : map.getIdentificationData().getRuns())
        if (! taken.insert(run.getUuid()).second) change = true;
    }
    if (! change) return maps;
    converted.clear();
    converted.reserve(maps.size());
    taken.clear();
    for (const auto& map : maps)
    {
      converted.push_back(map);
      auto& copy = converted.back();
      if (hasLegacyIDs(copy)) importMap(copy, true);
      distinctRuns(copy, taken);
    }
    return converted;
  }
} // namespace

bool IdentificationDataConverter::makeRunsDistinct(FeatureMap& map, std::set<std::string>& taken)
{ return distinctRuns(map, taken); }
bool IdentificationDataConverter::makeRunsDistinct(ConsensusMap& map, std::set<std::string>& taken)
{ return distinctRuns(map, taken); }
const std::vector<FeatureMap>& IdentificationDataConverter::withIdentificationData(const std::vector<FeatureMap>& maps, std::vector<FeatureMap>& converted)
{ return distinctMaps(maps, converted); }
const std::vector<ConsensusMap>& IdentificationDataConverter::withIdentificationData(const std::vector<ConsensusMap>& maps, std::vector<ConsensusMap>& converted)
{ return distinctMaps(maps, converted); }

namespace
{
  /// Meta value of the peptide hits of IdentificationDataConverter::exportWithMatchReferences(): "<run UUID>/<match ID>"
  const std::string MATCH_REFERENCE = "identification:match_reference";
} // namespace

ConsensusMap IdentificationDataConverter::exportWithMatchReferences(const ConsensusMap& map)
{
  if (hasLegacyIDs(map)) invalid("The map has peptide identifications already");
  ConsensusMap copy = map;
  auto& data = copy.getIdentificationData();
  for (const auto& current : data.getRuns())
  {
    std::vector<std::pair<ID::MatchId, ID::MatchData>> named;
    for (const auto& source : current.getSources())
      for (const auto& query : source.identifications)
        for (const auto& match : query.getMatches())
        {
          named.emplace_back(match.getId(), match.getData());
          named.back().second.setMetaValue(MATCH_REFERENCE, current.getUuid() + "/" + std::to_string(match.getId().value));
        }
    auto& run = data.getRun(current.getIdentifier());
    for (const auto& [id, match] : named)
      run.replaceMatch(id, match);
  }
  exportMap(copy, true);
  return copy;
}

std::optional<ID::MatchReference> IdentificationDataConverter::matchReference(const PeptideHit& hit)
{
  if (! hit.metaValueExists(MATCH_REFERENCE)) return std::nullopt;
  const std::string reference = hit.getMetaValue(MATCH_REFERENCE).toString();
  const auto separator = reference.rfind('/');
  if (separator == std::string::npos) invalid("Invalid match reference: " + reference);
  return ID::MatchReference {reference.substr(0, separator), ID::MatchId {std::stoull(reference.substr(separator + 1))}};
}

void IdentificationDataConverter::updateReferencedMatches(ID& data, const PeptideIdentificationList& peptides, bool scores)
{
  struct Update
  {
    std::set<std::string> accessions; // of all hits of the match
    double score = 0.0;               // of its first hit
  };
  std::map<ID::MatchReference, Update> updates;
  for (const auto& peptide : peptides)
  {
    for (const auto& hit : peptide.getHits())
    {
      const auto reference = matchReference(hit);
      if (! reference) continue;
      const auto [update, added] = updates.try_emplace(*reference);
      if (added) update->second.score = hit.getScore();
      const auto referenced = hit.extractProteinAccessionsSet();
      update->second.accessions.insert(referenced.begin(), referenced.end());
    }
  }
  if (updates.empty()) return;
  for (const auto& current : data.getRuns())
  {
    std::vector<std::pair<ID::MatchId, ID::MatchData>> edits;
    std::vector<std::pair<ID::MatchId, double>> scored;
    for (const auto& source : current.getSources())
      for (const auto& query : source.identifications)
        for (const auto& match : query.getMatches())
        {
          const auto found = updates.find({current.getUuid(), match.getId()});
          if (found == updates.end()) continue;
          if (scores) scored.emplace_back(match.getId(), found->second.score);
          ID::MatchData edited = match.getData();
          std::erase_if(edited.sequence_evidence, [&](const ID::SequenceEvidence& evidence) { return ! found->second.accessions.contains(evidence.accession); });
          if (edited.sequence_evidence.size() != match.sequence_evidence.size()) edits.emplace_back(match.getId(), std::move(edited));
        }
    if (edits.empty() && scored.empty()) continue;
    auto& run = data.getRun(current.getIdentifier());
    for (const auto& [id, edited] : edits) run.replaceMatch(id, edited);
    if (scored.empty()) continue;
    if (! run.getPrimaryScore()) invalid("Run '" + run.getIdentifier() + "' has no primary score for the scores of its matches");
    const auto primary = *run.getPrimaryScore();
    for (const auto& [id, score] : scored) run.setScore(id, primary, score);
  }
}

bool IdentificationDataConverter::moveToIdentificationData(FeatureMap& map)
{
  if (! checkedLegacyIDs(map)) return false;
  importMap(map, true);
  return true;
}
bool IdentificationDataConverter::moveToIdentificationData(ConsensusMap& map)
{
  if (! checkedLegacyIDs(map)) return false;
  importMap(map, true);
  return true;
}

MzTab IdentificationDataConverter::exportMzTab(const ID& data)
{
  data.validate();
  if (std::all_of(data.getRuns().begin(), data.getRuns().end(), [](const auto& run) { return run.getMoleculeKind() == ID::MoleculeKind::PEPTIDE; }))
  {
    auto converted = Adapter::toLegacy(data, legacyOptions());
    reportLegacyLosses(converted.losses);
    return MzTab::exportIdentificationsToMzTab(converted.proteins, converted.peptides, "", false, true, true);
  }
  MzTab result;
  MzTabMetaData metadata;
  MzTabNucleicAcidSectionRows nucleic_acids;
  MzTabOligonucleotideSectionRows oligos;
  MzTabOSMSectionRows matches;
  Size file = 0, software = 0;
  std::map<std::string, Size> ms_run_of_file; // a file listed more than once keeps one ms_run
  std::set<std::tuple<std::string, ID::QualifiedAccession, std::optional<UInt32>, std::optional<UInt32>>> seen;
  for (const auto& run : data.getRuns())
  {
    if (run.getMoleculeKind() != ID::MoleculeKind::OLIGONUCLEOTIDE)
      invalid("Use mzTab-M for compounds; mixed peptide/RNA export requires separate outputs");
    const auto& settings = run.getSettings();
    MzTabSoftwareMetaData sw;
    sw.software.setName(settings.software);
    sw.software.setValue(settings.software_version);
    metadata.software[++software] = sw;
    for (Size i = 0; i < run.getScoreDefinitions().size(); ++i)
    {
      const auto& score = run.getScoreDefinitions()[i];
      MzTabParameter definition;
      definition.setName(score.name);
      definition.setAccession(score.accession);
      if (const auto colon = score.accession.find(':'); colon != std::string::npos) definition.setCVLabel(score.accession.substr(0, colon));
      metadata.osm_search_engine_score[i + 1] = definition;
    }
    if (run.getDatabaseSequences())
      for (const auto& sequence : *run.getDatabaseSequences())
      {
        MzTabNucleicAcidSectionRow row;
        row.accession.set(sequence.accession);
        row.description.set(sequence.description);
        MzTabParameter engine;
        engine.setName(settings.software);
        engine.setValue(settings.software_version);
        row.search_engine.set({engine});
        if (sequence.metaValueExists("coverage")) row.coverage.set(static_cast<double>(sequence.getMetaValue("coverage")));
        row.opt_.push_back({"opt_sequence", MzTabString(sequence.sequence)});
        nucleic_acids.push_back(std::move(row));
      }
    for (const auto& source : run.getSources())
    {
      const auto& path = source.file.path;
      auto ms_run = path.empty() ? ms_run_of_file.end() : ms_run_of_file.find(path);
      if (ms_run == ms_run_of_file.end())
      {
        MzTabMSRunMetaData input;
        input.location.set(path);
        metadata.ms_run[++file] = input;
        if (! path.empty()) ms_run = ms_run_of_file.emplace(path, file).first;
      }
      const Size ms_run_index = ms_run == ms_run_of_file.end() ? file : ms_run->second;
      for (const auto& query : source.identifications)
        for (const auto& match : query.getMatches())
        {
          MzTabOSMSectionRow row;
          row.sequence.set(match.representation);
          row.charge.set(match.charge);
          if (query.rt)
          {
            MzTabDouble rt;
            rt.set(*query.rt);
            row.retention_time.set({rt});
          }
          if (query.mz) row.exp_mass_to_charge.set(*query.mz);
          if (match.calculated_mz) row.calc_mass_to_charge.set(*match.calculated_mz);
          else if (match.charge)
            row.calc_mass_to_charge.set(NASequence::fromString(match.representation).getMonoWeight(NASequence::Full, match.charge)
                                        / std::abs(match.charge));
          row.spectra_ref.setMSFile(ms_run_index);
          row.spectra_ref.setSpecRef(query.data_id);
          const auto scores = run.getScores(match);
          for (Size i = 0; i < scores.size(); ++i)
            if (scores[i]) row.search_engine_score[i + 1].set(*scores[i]);
          MzTabParameter engine;
          engine.setName(settings.software);
          engine.setValue(settings.software_version);
          row.search_engine.set({engine});
          if (match.details && match.details->adduct) row.opt_.push_back({"opt_adduct", MzTabString(match.details->adduct->getName())});
          if (match.metaValueExists("isotope_offset"))
            row.opt_.push_back({"opt_isotope_offset", MzTabString(match.getMetaValue("isotope_offset").toString())});
          matches.push_back(std::move(row));
          for (const auto& evidence : match.sequence_evidence)
          {
            if (! seen.emplace(match.representation, run.qualify(evidence.database, evidence.accession), evidence.start, evidence.end).second) continue;
            MzTabOligonucleotideSectionRow oligo;
            oligo.sequence.set(match.representation);
            oligo.accession.set(evidence.accession);
            std::set<std::pair<UInt32, std::string>> sequences;
            for (const auto& item : match.sequence_evidence)
              sequences.emplace(item.database.value, item.accession);
            oligo.unique.set(sequences.size() == 1);
            MzTabParameter engine;
            engine.setName(settings.software);
            engine.setValue(settings.software_version);
            oligo.search_engine.set({engine});
            oligo.pre.set(evidence.before == '[' ? std::string("-") : std::string(evidence.before ? 1 : 0, evidence.before));
            oligo.post.set(evidence.after == ']' ? std::string("-") : std::string(evidence.after ? 1 : 0, evidence.after));
            if (evidence.start) oligo.start.set(*evidence.start + 1);
            if (evidence.end) oligo.end.set(*evidence.end + 1);
            oligos.push_back(std::move(oligo));
          }
          if (match.sequence_evidence.empty() && seen.emplace(match.representation, ID::QualifiedAccession {}, std::nullopt, std::nullopt).second)
          {
            MzTabOligonucleotideSectionRow oligo;
            oligo.sequence.set(match.representation);
            oligos.push_back(std::move(oligo));
          }
        }
    }
  }
  result.setMetaData(metadata);
  result.setNucleicAcidSectionRows(nucleic_acids);
  result.setOligonucleotideSectionRows(oligos);
  result.setOSMSectionRows(matches);
  return result;
}
} // namespace OpenMS
