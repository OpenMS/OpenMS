// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
#include <algorithm>
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
        auto configuration = metadata;
        configuration.getHits().clear();
        configuration.getProteinGroups().clear();
        configuration.getIndistinguishableProteins().clear();
        // The files of the legacy run become the sources of the run.
        StringList files;
        configuration.getPrimaryMSRunPath(files);
        configuration.removeMetaValue("spectra_data");
        run.setProcessingMetadata(configuration);
        Adapter::addLegacySources(run, files);
        std::vector<ID::ParentRecord> parents;
        for (const auto& hit : metadata.getHits())
        {
          ID::ParentRecord parent;
          static_cast<MetaInfoInterface&>(parent) = hit;
          parent.identity = {metadata.getSearchParameters().db, hit.getAccession()};
          parent.sequence = hit.getSequence();
          parent.description = hit.getDescription();
          if (hit.getCoverage() >= 0) parent.setMetaValue("coverage", hit.getCoverage() / 100.0);
          parent.target_decoy = hit.getTargetDecoyType() == ProteinHit::TargetDecoyType::DECOY    ? ID::TargetDecoy::DECOY
                                : hit.getTargetDecoyType() == ProteinHit::TargetDecoyType::TARGET ? ID::TargetDecoy::TARGET
                                                                                                  : ID::TargetDecoy::UNKNOWN;
          parents.push_back(std::move(parent));
        }
        if (! parents.empty()) run.setParents(std::move(parents));
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
            match.identifiers.push_back({databases[i], accessions[i]});
        }
        if (hit.metaValueExists("identification:formula")) match.formula = hit.getMetaValue("identification:formula").toString();
        if (hit.metaValueExists("identification:name")) match.name = hit.getMetaValue("identification:name").toString();
        if (hit.metaValueExists("identification:calculated_mz"))
          match.calculated_mz = static_cast<double>(hit.getMetaValue("identification:calculated_mz"));
        if (hit.metaValueExists("identification:adduct_formula"))
          match.adduct
            = AdductInfo(hit.getMetaValue("adduct").toString(), EmpiricalFormula(hit.getMetaValue("identification:adduct_formula").toString()),
                         match.charge, static_cast<int>(hit.getMetaValue("identification:adduct_multiplier")));
        if (hit.getTargetDecoyType() == PeptideHit::TargetDecoyType::TARGET) match.target_decoy = ID::TargetDecoy::TARGET;
        else if (hit.getTargetDecoyType() == PeptideHit::TargetDecoyType::DECOY)
          match.target_decoy = ID::TargetDecoy::DECOY;
        else if (hit.getTargetDecoyType() == PeptideHit::TargetDecoyType::TARGET_DECOY)
          match.target_decoy = ID::TargetDecoy::BOTH;
        for (const auto& old : hit.getPeptideEvidences())
        {
          ID::ParentEvidence evidence;
          evidence.parent = {run.getProcessingMetadata().getSearchParameters().db, old.getProteinAccession()};
          if (old.getStart() >= 0) evidence.start = old.getStart();
          if (old.getEnd() >= 0) evidence.end = old.getEnd();
          evidence.before = std::string(1, old.getAABefore());
          evidence.after = std::string(1, old.getAAAfter());
          match.parent_evidence.push_back(std::move(evidence));
        }
        run.addMatch(query, match,
                     run.getPrimaryScore() ? std::vector<std::optional<double>> {hit.getScore()} : std::vector<std::optional<double>> {});
      }
    }
    // Reuse the checked protein adapter for inference and protein-only input runs.
    auto parent_data = Adapter::fromLegacy(proteins, {});
    for (const auto& protein : proteins)
      if (std::none_of(runs.begin(), runs.end(), [&](const auto& entry) { return entry.first.first == protein.getIdentifier(); }))
        imported.merge(Adapter::fromLegacy({protein}, {}));
    for (auto inference : parent_data.getInferenceResults())
    {
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
    auto processing = run.getProcessingMetadata();
    processing.setIdentifier(run.getIdentifier());
    std::vector<ProteinHit> parents;
    if (run.getParents())
      for (const auto& parent : *run.getParents())
      {
        ProteinHit hit;
        static_cast<MetaInfoInterface&>(hit) = parent;
        hit.setAccession(parent.identity.accession);
        hit.setSequence(parent.sequence);
        hit.setDescription(parent.description);
        if (parent.metaValueExists("coverage"))
        {
          // the coverage attribute represents it
          hit.setCoverage(static_cast<double>(parent.getMetaValue("coverage")) * 100.0);
          hit.removeMetaValue("coverage");
        }
        hit.setTargetDecoyType(parent.target_decoy == ID::TargetDecoy::DECOY    ? ProteinHit::TargetDecoyType::DECOY
                               : parent.target_decoy == ID::TargetDecoy::TARGET ? ProteinHit::TargetDecoyType::TARGET
                                                                                : ProteinHit::TargetDecoyType::UNKNOWN);
        parents.push_back(std::move(hit));
      }
    processing.setHits(parents);
    for (const auto& inference : data.getInferenceResults())
      if (std::any_of(inference.inputs.begin(), inference.inputs.end(), [&](const auto& input) { return input.run_uuid == run.getUuid(); }))
      {
        processing.setHits(inference.proteins.getHits());
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
        item.setSpectrumReference(query.data_id);
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
          const auto scores = match.getScores();
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
          exportParentMatches(match.parent_evidence, hit);
          if (match.adduct)
          {
            hit.setMetaValue("adduct", match.adduct->getName());
            hit.setMetaValue("identification:adduct_formula", match.adduct->getEmpiricalFormula().toString());
            hit.setMetaValue("identification:adduct_multiplier", static_cast<int>(match.adduct->getMolMultiplier()));
          }
          hit.setMetaValue("identification:encoding", static_cast<int>(match.encoding));
          if (! match.identifiers.empty())
          {
            StringList databases, accessions;
            for (const auto& identity : match.identifiers)
            {
              databases.push_back(identity.database);
              accessions.push_back(identity.accession);
            }
            hit.setMetaValue("identification:identifier_databases", databases);
            hit.setMetaValue("identification:identifier_accessions", accessions);
          }
          if (match.formula) hit.setMetaValue("identification:formula", *match.formula);
          if (! match.name.empty()) hit.setMetaValue("identification:name", match.name);
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

void IdentificationDataConverter::importSequences(ID::Run& run, const std::vector<FASTAFile::FASTAEntry>& fasta, const std::string& decoy_pattern)
{
  std::vector<ID::ParentRecord> parents = run.getParents().value_or(std::vector<ID::ParentRecord> {});
  const auto& database = run.getProcessingMetadata().getSearchParameters().db;
  for (const auto& entry : fasta)
  {
    ID::ParentRecord parent;
    parent.identity = {database, entry.identifier};
    parent.sequence = entry.sequence;
    parent.description = entry.description;
    parent.target_decoy
      = ! decoy_pattern.empty() && entry.identifier.find(decoy_pattern) != std::string::npos ? ID::TargetDecoy::DECOY : ID::TargetDecoy::TARGET;
    parents.push_back(std::move(parent));
  }
  run.setParents(std::move(parents));
}
void IdentificationDataConverter::exportParentMatches(const std::vector<ID::ParentEvidence>& parents, PeptideHit& hit)
{
  std::vector<PeptideEvidence> evidence;
  for (const auto& parent : parents)
  {
    if ((parent.start && *parent.start > std::numeric_limits<Int>::max()) || (parent.end && *parent.end > std::numeric_limits<Int>::max())
        || parent.before.size() > 1 || parent.after.size() > 1)
      invalid("Parent evidence exceeds legacy coordinate or flank limits");
    evidence.emplace_back(parent.parent.accession, parent.start ? static_cast<Int>(*parent.start) : PeptideEvidence::UNKNOWN_POSITION,
                          parent.end ? static_cast<Int>(*parent.end) : PeptideEvidence::UNKNOWN_POSITION,
                          parent.before.empty() ? PeptideEvidence::UNKNOWN_AA : parent.before[0],
                          parent.after.empty() ? PeptideEvidence::UNKNOWN_AA : parent.after[0]);
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
  MzTabNucleicAcidSectionRows parents;
  MzTabOligonucleotideSectionRows oligos;
  MzTabOSMSectionRows matches;
  Size file = 0, software = 0;
  std::map<std::string, Size> ms_run_of_file; // a file listed more than once keeps one ms_run
  std::set<std::tuple<std::string, ID::QualifiedAccession, std::optional<UInt64>, std::optional<UInt64>>> seen;
  for (const auto& run : data.getRuns())
  {
    if (run.getMoleculeKind() != ID::MoleculeKind::OLIGONUCLEOTIDE)
      invalid("Use mzTab-M for compounds; mixed peptide/RNA export requires separate outputs");
    const auto& processing = run.getProcessingMetadata();
    MzTabSoftwareMetaData sw;
    sw.software.setName(processing.getSearchEngine());
    sw.software.setValue(processing.getSearchEngineVersion());
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
    if (run.getParents())
      for (const auto& parent : *run.getParents())
      {
        MzTabNucleicAcidSectionRow row;
        row.accession.set(parent.identity.accession);
        row.description.set(parent.description);
        MzTabParameter engine;
        engine.setName(processing.getSearchEngine());
        engine.setValue(processing.getSearchEngineVersion());
        row.search_engine.set({engine});
        if (parent.metaValueExists("coverage")) row.coverage.set(static_cast<double>(parent.getMetaValue("coverage")));
        row.opt_.push_back({"opt_sequence", MzTabString(parent.sequence)});
        parents.push_back(std::move(row));
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
          for (Size i = 0; i < match.getScoreValues().size(); ++i)
            if (! std::isnan(match.getScoreValues()[i])) row.search_engine_score[i + 1].set(match.getScoreValues()[i]);
          MzTabParameter engine;
          engine.setName(processing.getSearchEngine());
          engine.setValue(processing.getSearchEngineVersion());
          row.search_engine.set({engine});
          if (match.adduct) row.opt_.push_back({"opt_adduct", MzTabString(match.adduct->getName())});
          if (match.metaValueExists("isotope_offset"))
            row.opt_.push_back({"opt_isotope_offset", MzTabString(match.getMetaValue("isotope_offset").toString())});
          matches.push_back(std::move(row));
          for (const auto& evidence : match.parent_evidence)
          {
            if (! seen.emplace(match.representation, evidence.parent, evidence.start, evidence.end).second) continue;
            MzTabOligonucleotideSectionRow oligo;
            oligo.sequence.set(match.representation);
            oligo.accession.set(evidence.parent.accession);
            std::set<ID::QualifiedAccession> parent_ids;
            for (const auto& parent : match.parent_evidence)
              parent_ids.insert(parent.parent);
            oligo.unique.set(parent_ids.size() == 1);
            MzTabParameter engine;
            engine.setName(processing.getSearchEngine());
            engine.setValue(processing.getSearchEngineVersion());
            oligo.search_engine.set({engine});
            oligo.pre.set(evidence.before == "[" ? "-" : evidence.before);
            oligo.post.set(evidence.after == "]" ? "-" : evidence.after);
            if (evidence.start) oligo.start.set(*evidence.start + 1);
            if (evidence.end) oligo.end.set(*evidence.end + 1);
            oligos.push_back(std::move(oligo));
          }
          if (match.parent_evidence.empty() && seen.emplace(match.representation, ID::QualifiedAccession {}, std::nullopt, std::nullopt).second)
          {
            MzTabOligonucleotideSectionRow oligo;
            oligo.sequence.set(match.representation);
            oligos.push_back(std::move(oligo));
          }
        }
    }
  }
  result.setMetaData(metadata);
  result.setNucleicAcidSectionRows(parents);
  result.setOligonucleotideSectionRows(oligos);
  result.setOSMSectionRows(matches);
  return result;
}
} // namespace OpenMS
