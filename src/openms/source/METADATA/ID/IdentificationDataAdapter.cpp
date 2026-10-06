// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/ModificationDefinitionIO.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/ID/IdentificationDataAdapter.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
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

  void loss(Adapter::LegacyResult& result, const Adapter::ExportOptions& options, const std::string& message)
  {
    if (options.loss_policy == Adapter::LossPolicy::STRICT) invalid(message);
    if (std::find(result.losses.begin(), result.losses.end(), message) == result.losses.end()) result.losses.push_back(message);
  }

  ID::TargetDecoy targetDecoy(PeptideHit::TargetDecoyType value)
  {
    switch (value)
    {
      case PeptideHit::TargetDecoyType::TARGET:
        return ID::TargetDecoy::TARGET;
      case PeptideHit::TargetDecoyType::DECOY:
        return ID::TargetDecoy::DECOY;
      case PeptideHit::TargetDecoyType::TARGET_DECOY:
        return ID::TargetDecoy::BOTH;
      default:
        return ID::TargetDecoy::UNKNOWN;
    }
  }

  void registerDefinitions(const ProteinIdentification::SearchParameters& parameters)
  {
    const auto& key = Constants::UserParam::MODIFICATION_DEFINITIONS;
    if (! parameters.metaValueExists(key)) return;
    const auto* db = ModificationsDB::getInstance();
    std::vector<ResidueModification> definitions;
    for (const auto& text : ResidueModification::splitDefinitionRecords(parameters.getMetaValue(key).toString()))
    {
      auto definition = ResidueModification::fromDefinitionString(text);
      if (db->has(definition.getFullId()))
      {
        const auto* old = db->getModification(definition.getFullId());
        if (old->getDiffFormula() != definition.getDiffFormula() || old->getDiffMonoMass() != definition.getDiffMonoMass()
            || old->getDiffAverageMass() != definition.getDiffAverageMass() || old->getOrigin() != definition.getOrigin()
            || old->getTermSpecificity() != definition.getTermSpecificity()
            || old->getNeutralLossDiffFormulas() != definition.getNeutralLossDiffFormulas())
        {
          invalid("Conflicting custom modification definition: " + definition.getFullId());
        }
      }
      definitions.push_back(std::move(definition));
    }
    // Validate every definition before registering any of them. The legacy helper logs
    // malformed definitions and continues; explicit materialization must fail instead.
    for (const auto& definition : definitions)
      db->registerDefinition(definition);
  }

  PeptideHit peptide(const ID::Run& run, const ID::Match& match, ID::ScoreId score)
  {
    if (run.getMoleculeKind() != ID::MoleculeKind::PEPTIDE || match.encoding != ID::Encoding::AA_SEQUENCE)
      invalid("Legacy peptide conversion requires AASequence-encoded peptide matches");
    // Bound access reads the match's dense scores directly instead of building the run's lookup index.
    const auto value = run.bindScore(score)(match);
    if (! value) invalid("Cannot materialize a peptide without the selected score");
    PeptideHit hit;
    static_cast<MetaInfoInterface&>(hit) = match;
    hit.setSequence(AASequence::fromString(match.representation));
    hit.setCharge(match.charge);
    hit.setScore(*value);
    hit.setPeakAnnotations(match.peak_annotations);
    std::vector<PeptideEvidence> evidence;
    for (const auto& parent : match.parent_evidence)
    {
      const auto max_position = static_cast<UInt64>(std::numeric_limits<Int>::max());
      if ((parent.start && *parent.start > max_position) || (parent.end && *parent.end > max_position) || parent.before.size() > 1
          || parent.after.size() > 1)
        invalid("Parent evidence cannot be represented by legacy peptide coordinates/flanking residues");
      evidence.emplace_back(parent.parent.accession, parent.start ? static_cast<Int>(*parent.start) : PeptideEvidence::UNKNOWN_POSITION,
                            parent.end ? static_cast<Int>(*parent.end) : PeptideEvidence::UNKNOWN_POSITION,
                            parent.before.empty() ? PeptideEvidence::UNKNOWN_AA : parent.before.front(),
                            parent.after.empty() ? PeptideEvidence::UNKNOWN_AA : parent.after.front());
    }
    hit.setPeptideEvidences(evidence);
    if (targetDecoy(hit.getTargetDecoyType()) != match.target_decoy)
    {
      switch (match.target_decoy)
      {
        case ID::TargetDecoy::TARGET:
          hit.setTargetDecoyType(PeptideHit::TargetDecoyType::TARGET);
          break;
        case ID::TargetDecoy::DECOY:
          hit.setTargetDecoyType(PeptideHit::TargetDecoyType::DECOY);
          break;
        case ID::TargetDecoy::BOTH:
          hit.setTargetDecoyType(PeptideHit::TargetDecoyType::TARGET_DECOY);
          break;
        default:
          hit.setTargetDecoyType(PeptideHit::TargetDecoyType::UNKNOWN);
          break;
      }
    }
    return hit;
  }

  const ID::InferenceResult* selectInference(const ID& data, const ID::Run& run, const Adapter::ExportOptions& options)
  {
    if (! options.include_inference) return nullptr;
    const ID::InferenceResult* selected = nullptr;
    for (const auto& result : data.getInferenceResults())
    {
      if (options.inference_result && result.identifier != *options.inference_result) continue;
      const bool applies
        = std::any_of(result.inputs.begin(), result.inputs.end(), [&](const auto& input) { return input.run_uuid == run.getUuid(); });
      if (! applies) continue;
      if (selected) invalid("Multiple inference results apply to run " + run.getIdentifier() + "; select a result explicitly");
      selected = &result;
    }
    return selected;
  }

  ProteinIdentification originalParents(const ID::Run& run, Adapter::LegacyResult& result, const Adapter::ExportOptions& options)
  {
    auto proteins = run.getProcessingMetadata();
    if (proteins.getIdentifier().empty()) proteins.setIdentifier(run.getIdentifier());
    StringList paths;
    proteins.getPrimaryMSRunPath(paths);
    if (paths.empty())
    {
      for (const auto& source : run.getSourceBlocks())
      {
        if (! source.source.primary_files.empty())
        {
          paths = source.source.primary_files;
          break;
        }
      }
      if (paths.empty())
        for (const auto& source : run.getSourceBlocks())
          if (! source.source.path.empty() && std::find(paths.begin(), paths.end(), source.source.path) == paths.end())
            paths.push_back(source.source.path);
      if (! paths.empty()) proteins.setPrimaryMSRunPath(paths);
    }
    if (run.getParents())
    {
      std::vector<ProteinHit> hits;
      for (const auto& parent : *run.getParents())
      {
        ProteinHit hit;
        static_cast<MetaInfoInterface&>(hit) = parent;
        hit.setAccession(parent.identity.accession);
        hit.setSequence(parent.sequence);
        hit.setDescription(parent.description);
        if (parent.target_decoy == ID::TargetDecoy::BOTH)
          loss(result, options, "Legacy proteins cannot represent a combined target/decoy parent state");
        else
        {
          const auto state = parent.target_decoy == ID::TargetDecoy::TARGET  ? ProteinHit::TargetDecoyType::TARGET
                             : parent.target_decoy == ID::TargetDecoy::DECOY ? ProteinHit::TargetDecoyType::DECOY
                                                                             : ProteinHit::TargetDecoyType::UNKNOWN;
          if (hit.getTargetDecoyType() != state) hit.setTargetDecoyType(state);
        }
        hits.push_back(std::move(hit));
      }
      proteins.setHits(hits);
    }
    return proteins;
  }

  std::vector<Adapter::FeatureAssociation> makeAssociations(const Adapter::ImportResult& imported, std::vector<Adapter::FeatureAssociation> locations)
  {
    for (Size i = 0; i < locations.size(); ++i)
    {
      locations[i].query = imported.queries[i];
      const auto* run = imported.data.findRunByUuid(locations[i].query.run_uuid);
      for (const auto& match : run->getIdentification(locations[i].query.query).getMatches())
        locations[i].matches.push_back(match.getId());
    }
    return locations;
  }

  void collectFeatures(const std::vector<Feature>& features,
                       std::vector<Size> prefix,
                       PeptideIdentificationList& peptides,
                       std::vector<Adapter::FeatureAssociation>& locations)
  {
    for (Size i = 0; i < features.size(); ++i)
    {
      auto path = prefix;
      path.push_back(i);
      for (const auto& item : features[i].getPeptideIdentifications())
      {
        peptides.push_back(item);
        Adapter::FeatureAssociation association;
        association.feature_path = path;
        locations.push_back(std::move(association));
      }
      collectFeatures(features[i].getSubordinates(), path, peptides, locations);
    }
  }

  void clearFeatures(std::vector<Feature>& features)
  {
    for (auto& feature : features)
    {
      feature.getPeptideIdentifications().clear();
      clearFeatures(feature.getSubordinates());
    }
  }

  PeptideIdentification
  linkedPeptide(const ID& data, const Adapter::LegacyResult& converted, Size index, const Adapter::FeatureAssociation& association)
  {
    auto result = converted.peptides[index];
    const auto* run = data.findRunByUuid(association.query.run_uuid);
    const auto& query = run->getIdentification(association.query.query);
    if (result.getHits().size() != query.getMatches().size()) invalid("Lossy candidate export cannot be used for feature associations");
    std::set<ID::MatchId> retained(association.matches.begin(), association.matches.end());
    std::vector<PeptideHit> hits;
    for (Size i = 0; i < query.getMatches().size(); ++i)
    {
      if (retained.contains(query.getMatches()[i].getId())) hits.push_back(result.getHits().at(i));
    }
    result.setHits(std::move(hits));
    return result;
  }
} // namespace

IdentificationDataAdapter::ImportResult IdentificationDataAdapter::importLegacy(const std::vector<ProteinIdentification>& proteins,
                                                                                const PeptideIdentificationList& peptides)
{
  ImportResult result;
  auto definitions = ModificationDefinitionIO::collect(proteins, peptides);
  std::map<std::string, ProteinIdentification> originals;
  for (auto protein : proteins)
  {
    if (originals.contains(protein.getIdentifier())) invalid("Duplicate legacy protein run identifier: " + protein.getIdentifier());
    auto params = protein.getSearchParameters();
    ModificationDefinitionIO::attach(params, definitions[protein.getIdentifier()]);
    protein.setSearchParameters(params);
    originals.emplace(protein.getIdentifier(), std::move(protein));
  }
  using Contract = std::tuple<std::string, std::string, bool>;
  std::map<Contract, std::string> contracts;
  std::map<std::string, std::vector<std::string>> input_runs;
  std::map<std::string, std::map<SignedSize, UInt32>> sources;
  std::set<std::string> used_names;
  for (const auto& [name, original] : originals)
    used_names.insert(name);

  auto create_run = [&](const Contract& contract) -> ID::Run& {
    const auto& original_id = std::get<0>(contract);
    const auto original = originals.find(original_id);
    if (original == originals.end()) invalid("Peptide identification has no matching protein run: " + original_id);
    auto found = contracts.find(contract);
    if (found != contracts.end()) return result.data.getRun(found->second);
    std::string name = original_id;
    if (! input_runs[original_id].empty())
    {
      Size suffix = 1;
      do
      {
        name = original_id + ":score_" + std::to_string(suffix++);
      } while (used_names.contains(name));
    }
    used_names.insert(name);
    auto& run = result.data.addRun(name);
    auto configuration = original->second;
    configuration.setHits({});
    configuration.getProteinGroups().clear();
    configuration.getIndistinguishableProteins().clear();
    run.setProcessingMetadata(configuration);
    std::vector<ID::ParentRecord> parents;
    for (const auto& hit : original->second.getHits())
    {
      ID::ParentRecord parent;
      static_cast<MetaInfoInterface&>(parent) = hit;
      parent.identity = {configuration.getSearchParameters().db, hit.getAccession()};
      parent.sequence = hit.getSequence();
      parent.description = hit.getDescription();
      if (hit.getTargetDecoyType() == ProteinHit::TargetDecoyType::TARGET) parent.target_decoy = ID::TargetDecoy::TARGET;
      if (hit.getTargetDecoyType() == ProteinHit::TargetDecoyType::DECOY) parent.target_decoy = ID::TargetDecoy::DECOY;
      parents.push_back(std::move(parent));
    }
    run.setParents(std::move(parents));
    ID::ScoreDefinition definition;
    definition.name = std::get<1>(contract);
    definition.higher_better = std::get<2>(contract);
    // A derived score belongs to the tool that recorded itself as its producer, not to the search engine.
    std::tie(definition.software, definition.software_version) = configuration.getScoreSoftware(definition.name);
    if (! definition.name.empty())
    {
      const auto score = run.addScore(definition);
      run.setPrimaryScore(score);
    }
    contracts.emplace(contract, name);
    input_runs[original_id].push_back(name);
    return run;
  };

  for (const auto& item : peptides)
  {
    auto& run = create_run({item.getIdentifier(), item.getScoreType(), item.isHigherScoreBetter()});
    StringList paths;
    originals.at(item.getIdentifier()).getPrimaryMSRunPath(paths);
    SignedSize source_index = -1;
    if (item.metaValueExists("id_merge_index"))
    {
      const auto& value = item.getMetaValue("id_merge_index");
      if (value.valueType() != DataValue::INT_VALUE) invalid("id_merge_index must be an integer");
      source_index = static_cast<SignedSize>(static_cast<Int64>(value));
      if (source_index < 0 || static_cast<Size>(source_index) >= paths.size()) invalid("id_merge_index is outside the primary MS file list");
    }
    else if (paths.size() == 1)
      source_index = 0;
    // Multiple primary files without an explicit merge index remain an unknown
    // source. Neither basename matching nor consensus map indices resolve this.
    auto& by_source = sources[run.getIdentifier()];
    auto source = by_source.find(source_index);
    if (source == by_source.end())
    {
      ID::SourceFile descriptor;
      descriptor.primary_files = paths;
      if (source_index >= 0) descriptor.path = paths[static_cast<Size>(source_index)];
      source = by_source.emplace(source_index, run.addSource(descriptor).value).first;
    }
    ID::Observation observation;
    static_cast<MetaInfoInterface&>(observation) = item;
    observation.data_id = item.getSpectrumReference();
    if (item.hasRT()) observation.rt = item.getRT();
    if (item.hasMZ()) observation.mz = item.getMZ();
    auto query = run.addIdentification(run.getSourceId(source->second), observation);
    result.queries.push_back({run.getUuid(), query});
    for (const auto& hit : item.getHits())
    {
      if (! run.getPrimaryScore()) invalid("A legacy peptide hit requires a nonempty score type");
      ID::MatchData match;
      static_cast<MetaInfoInterface&>(match) = hit;
      match.representation = hit.getSequence().toString();
      match.charge = hit.getCharge();
      match.target_decoy = targetDecoy(hit.getTargetDecoyType());
      match.peak_annotations = hit.getPeakAnnotations();
      for (const auto& item_evidence : hit.getPeptideEvidences())
      {
        ID::ParentEvidence evidence;
        evidence.parent = {run.getProcessingMetadata().getSearchParameters().db, item_evidence.getProteinAccession()};
        if (item_evidence.getStart() != PeptideEvidence::UNKNOWN_POSITION)
        {
          if (item_evidence.getStart() < 0) invalid("Unsupported negative parent evidence start");
          evidence.start = static_cast<UInt64>(item_evidence.getStart());
        }
        if (item_evidence.getEnd() != PeptideEvidence::UNKNOWN_POSITION)
        {
          if (item_evidence.getEnd() < 0) invalid("Unsupported negative parent evidence end");
          evidence.end = static_cast<UInt64>(item_evidence.getEnd());
        }
        evidence.before = std::string(1, item_evidence.getAABefore());
        evidence.after = std::string(1, item_evidence.getAAAfter());
        match.parent_evidence.push_back(std::move(evidence));
      }
      run.addMatch(query, match, {hit.getScore()});
    }
  }
  // Protein-only runs are legitimate and remain representable after an empty search.
  for (const auto& [name, original] : originals)
  {
    if (input_runs[name].empty()) create_run({name, "", true});
    ID::InferenceResult inference;
    inference.identifier = "legacy:" + name;
    inference.proteins = original;
    for (const auto& hit : original.getHits())
      inference.qualified_accessions[hit.getAccession()] = {original.getSearchParameters().db, hit.getAccession()};
    for (const auto& run_name : input_runs[name])
    {
      const auto& run = result.data.getRun(run_name);
      ID::InferenceInput input;
      input.run_identifier = run_name;
      input.run_uuid = run.getUuid();
      input.selection = "Imported legacy run-level provenance";
      inference.inputs.push_back(std::move(input));
    }
    result.data.addInferenceResult(inference);
  }
  result.data.validate();
  return result;
}

IdentificationData IdentificationDataAdapter::fromLegacy(const std::vector<ProteinIdentification>& proteins,
                                                         const PeptideIdentificationList& peptides)
{ return importLegacy(proteins, peptides).data; }

PeptideHit IdentificationDataAdapter::materializePeptide(const ID::Run& run, const ID::Match& match, ID::ScoreId score)
{
  registerDefinitions(run.getProcessingMetadata().getSearchParameters());
  return peptide(run, match, score);
}

IdentificationDataAdapter::LegacyResult IdentificationDataAdapter::toLegacy(const ID& data)
{ return toLegacy(data, ExportOptions {}); }

IdentificationDataAdapter::LegacyResult IdentificationDataAdapter::toLegacy(const ID& data, const ExportOptions& options)
{
  data.validate();
  LegacyResult result;
  if (options.inference_result && ! options.include_inference) invalid("An inference result was selected while inference export is disabled");
  if (options.inference_result && std::none_of(data.getInferenceResults().begin(), data.getInferenceResults().end(), [&](const auto& item) {
        return item.identifier == *options.inference_result;
      }))
    invalid("Selected inference result does not exist");
  std::set<std::string> handled_inference;
  std::map<std::string, Size> protein_indices;
  for (const auto& run : data.getRuns())
  {
    if (run.getMoleculeKind() != ID::MoleculeKind::PEPTIDE)
    {
      loss(result, options, "Legacy peptide/protein export cannot represent non-peptide run " + run.getIdentifier());
      continue;
    }
    const auto* inference = selectInference(data, run, options);
    auto proteins = inference ? inference->proteins : originalParents(run, result, options);
    if (proteins.getIdentifier().empty()) proteins.setIdentifier(run.getIdentifier());
    if (inference && handled_inference.insert(inference->identifier).second)
    {
      if (inference->parent_score || inference->group_score
          || std::any_of(inference->inputs.begin(), inference->inputs.end(), [](const auto& input) { return input.score.has_value(); }))
        loss(result, options, "Legacy export cannot preserve complete inference score definitions: " + inference->identifier);
    }
    const auto primary = run.getPrimaryScore();
    if (primary)
    {
      // Record a producer other than the search engine so the score definition survives the round trip.
      const auto& definition = run.getScoreDefinition(*primary);
      if (! definition.name.empty() && ! definition.software.empty()
          && std::make_pair(definition.software, definition.software_version) != proteins.getScoreSoftware(definition.name))
        proteins.setScoreSoftware(definition.name, definition.software, definition.software_version);
    }
    const auto found = protein_indices.find(proteins.getIdentifier());
    if (found == protein_indices.end())
    {
      protein_indices[proteins.getIdentifier()] = result.proteins.size();
      result.proteins.push_back(proteins);
    }
    else if (result.proteins[found->second] != proteins)
    {
      // Same legacy identifier is only shared when its complete original payload agrees.
      loss(result, options, "Conflicting protein payloads use the same legacy identifier: " + proteins.getIdentifier());
      std::string identifier = run.getIdentifier();
      while (protein_indices.contains(identifier))
        identifier += ":export";
      proteins.setIdentifier(identifier);
      protein_indices[identifier] = result.proteins.size();
      result.proteins.push_back(proteins);
    }
    if (! primary && run.getNumberOfMatches())
    {
      loss(result, options, "Run has matches but no primary score: " + run.getIdentifier());
      continue;
    }
    if (run.getScoreDefinitions().size() > 1)
      loss(result, options, "Legacy export cannot retain complete secondary score definitions: " + run.getIdentifier());
    if (primary)
    {
      const auto& definition = run.getScoreDefinition(*primary);
      if (! definition.accession.empty() || definition.scope != ID::ScoreScope::MATCH || ! definition.calibration.empty()
          || ! definition.aggregation.empty() || ! definition.parameters.isMetaEmpty()
          || std::make_pair(definition.software, definition.software_version) != proteins.getScoreSoftware(definition.name))
        loss(result, options, "Legacy export cannot retain the complete primary score definition: " + run.getIdentifier());
    }
    registerDefinitions(run.getProcessingMetadata().getSearchParameters());
    for (const auto& source : run.getSourceBlocks())
    {
      StringList original_paths;
      proteins.getPrimaryMSRunPath(original_paths);
      if ((! source.source.primary_files.empty() && source.source.primary_files != original_paths)
          || (! source.source.path.empty() && std::find(original_paths.begin(), original_paths.end(), source.source.path) == original_paths.end()))
        loss(result, options, "Source descriptor cannot be represented by the selected legacy protein run: " + run.getIdentifier());
      if (! source.source.isMetaEmpty()) loss(result, options, "Legacy export cannot retain source-level metadata: " + run.getIdentifier());
      if (source.source.path.empty() && original_paths.size() == 1)
        loss(result, options, "An unknown source cannot be represented in a legacy run with exactly one known source: " + run.getIdentifier());
      for (const auto& query : source.identifications)
      {
        PeptideIdentification item;
        static_cast<MetaInfoInterface&>(item) = query;
        item.setIdentifier(proteins.getIdentifier());
        if (primary)
        {
          item.setScoreType(run.getScoreDefinition(*primary).name);
          item.setHigherScoreBetter(run.getScoreDefinition(*primary).higher_better);
        }
        if (query.rt) item.setRT(*query.rt);
        if (query.mz) item.setMZ(*query.mz);
        if (! query.data_id.empty()) item.setSpectrumReference(query.data_id);
        if (item.metaValueExists("id_merge_index"))
        {
          const auto& index_value = item.getMetaValue("id_merge_index");
          bool valid_index = index_value.valueType() == DataValue::INT_VALUE;
          const Int64 index = valid_index ? static_cast<Int64>(index_value) : -1;
          valid_index = valid_index && index >= 0 && static_cast<Size>(index) < original_paths.size() && ! source.source.path.empty()
                        && original_paths[static_cast<Size>(index)] == source.source.path;
          if (! valid_index)
          {
            loss(result, options, "Legacy source index contradicts its source descriptor: " + run.getIdentifier());
            item.removeMetaValue("id_merge_index");
          }
        }
        if (! item.metaValueExists("id_merge_index") && ! source.source.path.empty() && original_paths.size() > 1)
        {
          if (std::count(original_paths.begin(), original_paths.end(), source.source.path) == 1)
            item.setMetaValue("id_merge_index", static_cast<Int64>(std::find(original_paths.begin(), original_paths.end(), source.source.path)
                                                                   - original_paths.begin()));
          else
            loss(result, options, "Duplicate source paths need an explicit legacy file index: " + run.getIdentifier());
        }
        if (query.getSelectedMatch()) loss(result, options, "Legacy export cannot preserve an explicit selected candidate: " + run.getIdentifier());
        for (const auto& match : query.getMatches())
        {
          if (match.calculated_mz || match.adduct || match.formula || ! match.name.empty() || ! match.identifiers.empty())
            loss(result, options, "Legacy export cannot preserve all molecular/ion fields: " + run.getIdentifier());
          for (const auto& evidence : match.parent_evidence)
          {
            if (evidence.parent.database != proteins.getSearchParameters().db)
              loss(result, options, "Legacy export cannot represent the parent database namespace: " + run.getIdentifier());
          }
          if (match.encoding != ID::Encoding::AA_SEQUENCE)
          {
            loss(result, options, "Legacy export cannot materialize the molecular encoding: " + run.getIdentifier());
            continue;
          }
          auto hit = peptide(run, match, *primary);
          const auto scores = match.getScores();
          for (Size score = 0; score < run.getScoreDefinitions().size(); ++score)
          {
            if (score == primary->value || ! scores[score]) continue;
            const auto& name = run.getScoreDefinitions()[score].name;
            if (hit.metaValueExists(name) && hit.getMetaValue(name) != DataValue(*scores[score]))
              loss(result, options, "Secondary score collides with existing metadata: " + name);
            else
              hit.setMetaValue(name, *scores[score]);
          }
          item.insertHit(std::move(hit));
        }
        result.queries.push_back({run.getUuid(), query.getId()});
        result.peptides.push_back(std::move(item));
      }
    }
  }
  if (options.include_inference)
  {
    for (const auto& inference : data.getInferenceResults())
    {
      if (options.inference_result && inference.identifier != *options.inference_result) continue;
      if (! handled_inference.contains(inference.identifier))
        loss(result, options, "Inference result has no representable input run: " + inference.identifier);
    }
  }
  return result;
}

namespace
{
  template<class Map>
  IdentificationDataAdapter::FeatureImportResult nativeMap(const Map& map)
  {
    using ID = IdentificationData;
    using Adapter = IdentificationDataAdapter;
    const auto& data = map.getIdentificationData();
    data.validate();
    std::map<ID::MatchReference, ID::QueryReference> owners;
    for (const auto& run : data.getRuns())
      for (const auto& source : run.getSourceBlocks())
        for (const auto& query : source.identifications)
          for (const auto& match : query.getMatches())
            owners.emplace(ID::MatchReference {run.getUuid(), match.getId()}, ID::QueryReference {run.getUuid(), query.getId()});
    std::vector<Adapter::FeatureAssociation> associations;
    std::map<ID::QueryReference, std::set<ID::MatchId>> assigned;
    std::set<ID::QueryReference> assigned_queries;
    const auto collect = [&](const auto& self, const auto& feature, std::vector<Size> path) -> void {
      std::map<ID::QueryReference, std::set<ID::MatchId>> linked;
      for (const auto& query : feature.getIDQueries())
        linked[query];
      for (const auto& match : feature.getIDMatches())
      {
        auto owner = owners.find(match);
        if (owner == owners.end()) invalid("Feature association refers to a missing match");
        linked[owner->second].insert(match.match);
      }
      for (const auto& [query, matches] : linked)
      {
        const auto* run = data.findRunByUuid(query.run_uuid);
        if (! run || ! run->findIdentification(query.query)) invalid("Feature association refers to a missing query");
        Adapter::FeatureAssociation association;
        association.feature_path = path;
        association.query = query;
        association.matches.assign(matches.begin(), matches.end());
        associations.push_back(std::move(association));
        assigned[query].insert(matches.begin(), matches.end());
        assigned_queries.insert(query);
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
    for (const auto& run : data.getRuns())
      for (const auto& source : run.getSourceBlocks())
        for (const auto& query : source.identifications)
        {
          ID::QueryReference reference {run.getUuid(), query.getId()};
          Adapter::FeatureAssociation unassigned;
          unassigned.unassigned = true;
          unassigned.query = reference;
          for (const auto& match : query.getMatches())
            if (! assigned[reference].contains(match.getId())) unassigned.matches.push_back(match.getId());
          if (! unassigned.matches.empty() || (! assigned_queries.contains(reference) && query.getMatches().empty()))
            associations.push_back(std::move(unassigned));
        }
    return {data, std::move(associations)};
  }
} // namespace

IdentificationDataAdapter::FeatureImportResult IdentificationDataAdapter::fromFeatureMap(const FeatureMap& map)
{
  if (! map.getIdentificationData().empty()) return nativeMap(map);
  PeptideIdentificationList peptides;
  std::vector<FeatureAssociation> locations;
  for (Size i = 0; i < map.size(); ++i)
  {
    for (const auto& item : map[i].getPeptideIdentifications())
    {
      peptides.push_back(item);
      FeatureAssociation location;
      location.feature_path = {i};
      locations.push_back(std::move(location));
    }
    collectFeatures(map[i].getSubordinates(), {i}, peptides, locations);
  }
  for (const auto& item : map.getUnassignedPeptideIdentifications())
  {
    peptides.push_back(item);
    FeatureAssociation location;
    location.unassigned = true;
    locations.push_back(std::move(location));
  }
  auto imported = importLegacy(map.getProteinIdentifications(), peptides);
  auto associations = makeAssociations(imported, std::move(locations));
  return {std::move(imported.data), std::move(associations)};
}

IdentificationDataAdapter::FeatureImportResult IdentificationDataAdapter::fromConsensusMap(const ConsensusMap& map)
{
  if (! map.getIdentificationData().empty()) return nativeMap(map);
  PeptideIdentificationList peptides;
  std::vector<FeatureAssociation> locations;
  for (Size i = 0; i < map.size(); ++i)
  {
    for (const auto& item : map[i].getPeptideIdentifications())
    {
      peptides.push_back(item);
      FeatureAssociation location;
      location.feature_path = {i};
      locations.push_back(std::move(location));
    }
  }
  for (const auto& item : map.getUnassignedPeptideIdentifications())
  {
    peptides.push_back(item);
    FeatureAssociation location;
    location.unassigned = true;
    locations.push_back(std::move(location));
  }
  auto imported = importLegacy(map.getProteinIdentifications(), peptides);
  auto associations = makeAssociations(imported, std::move(locations));
  return {std::move(imported.data), std::move(associations)};
}

Size IdentificationDataAdapter::reconcileAssociations(const ID& data, std::vector<FeatureAssociation>& associations, MissingLinkPolicy policy)
{
  std::vector<FeatureAssociation> reconciled;
  Size removed = 0;
  for (const auto& association : associations)
  {
    const auto* run = data.findRunByUuid(association.query.run_uuid);
    const auto* query = run ? run->findIdentification(association.query.query) : nullptr;
    if (! query)
    {
      if (policy == MissingLinkPolicy::REJECT) invalid("Feature association references a removed query");
      removed += std::max<Size>(1, association.matches.size());
      continue;
    }
    auto retained = association;
    retained.matches.clear();
    for (const auto match : association.matches)
    {
      const bool live = std::any_of(query->getMatches().begin(), query->getMatches().end(), [&](const auto& item) { return item.getId() == match; });
      if (live) retained.matches.push_back(match);
      else
      {
        if (policy == MissingLinkPolicy::REJECT) invalid("Feature association references a removed or foreign candidate");
        ++removed;
      }
    }
    reconciled.push_back(std::move(retained));
  }
  associations.swap(reconciled);
  return removed;
}

std::vector<std::string> IdentificationDataAdapter::applyToFeatureMap(const ID& data,
                                                                      const std::vector<FeatureAssociation>& associations,
                                                                      FeatureMap& map,
                                                                      const ExportOptions& options,
                                                                      MissingLinkPolicy policy)
{
  auto links = associations;
  reconcileAssociations(data, links, policy);
  const auto converted = toLegacy(data, options);
  std::map<QueryReference, Size> query_indices;
  for (Size i = 0; i < converted.queries.size(); ++i)
    query_indices[converted.queries[i]] = i;
  auto updated = map;
  updated.getIdentificationData() = data;
  const auto clear_links = [&](const auto& self, auto& feature) -> void {
    feature.getIDMatches().clear();
    feature.getIDQueries().clear();
    if constexpr (std::is_same_v<std::remove_cvref_t<decltype(feature)>, Feature>)
      for (auto& child : feature.getSubordinates())
        self(self, child);
  };
  for (auto& feature : updated)
    clear_links(clear_links, feature);
  updated.setProteinIdentifications(converted.proteins);
  updated.getUnassignedPeptideIdentifications().clear();
  for (auto& feature : updated)
  {
    feature.getPeptideIdentifications().clear();
    clearFeatures(feature.getSubordinates());
  }
  for (const auto& association : links)
  {
    const auto index = query_indices.find(association.query);
    if (index == query_indices.end()) invalid("Associated query is not representable as a legacy peptide identification");
    auto item = linkedPeptide(data, converted, index->second, association);
    if (association.unassigned) updated.getUnassignedPeptideIdentifications().push_back(std::move(item));
    else
    {
      if (association.feature_path.empty() || association.feature_path.front() >= updated.size()) invalid("Feature association path is out of range");
      auto* feature = &updated[association.feature_path.front()];
      for (Size i = 1; i < association.feature_path.size(); ++i)
      {
        if (association.feature_path[i] >= feature->getSubordinates().size()) invalid("Subordinate feature association path is out of range");
        feature = &feature->getSubordinates()[association.feature_path[i]];
      }
      feature->addIDQuery(association.query);
      for (auto match : association.matches)
        feature->addIDMatch({association.query.run_uuid, match});
      feature->getPeptideIdentifications().push_back(std::move(item));
    }
  }
  map = std::move(updated);
  return converted.losses;
}

std::vector<std::string> IdentificationDataAdapter::applyToConsensusMap(const ID& data,
                                                                        const std::vector<FeatureAssociation>& associations,
                                                                        ConsensusMap& map,
                                                                        const ExportOptions& options,
                                                                        MissingLinkPolicy policy)
{
  auto links = associations;
  reconcileAssociations(data, links, policy);
  const auto converted = toLegacy(data, options);
  std::map<QueryReference, Size> query_indices;
  for (Size i = 0; i < converted.queries.size(); ++i)
    query_indices[converted.queries[i]] = i;
  auto updated = map;
  updated.getIdentificationData() = data;
  const auto clear_links = [&](const auto& self, auto& feature) -> void {
    feature.getIDMatches().clear();
    feature.getIDQueries().clear();
    if constexpr (std::is_same_v<std::remove_cvref_t<decltype(feature)>, Feature>)
      for (auto& child : feature.getSubordinates())
        self(self, child);
  };
  for (auto& feature : updated)
    clear_links(clear_links, feature);
  updated.setProteinIdentifications(converted.proteins);
  updated.getUnassignedPeptideIdentifications().clear();
  for (auto& feature : updated)
    feature.getPeptideIdentifications().clear();
  for (const auto& association : links)
  {
    const auto index = query_indices.find(association.query);
    if (index == query_indices.end()) invalid("Associated query is not representable as a legacy peptide identification");
    auto item = linkedPeptide(data, converted, index->second, association);
    if (association.unassigned) updated.getUnassignedPeptideIdentifications().push_back(std::move(item));
    else
    {
      if (association.feature_path.size() != 1 || association.feature_path.front() >= updated.size())
        invalid("Consensus feature association path is out of range");
      auto& feature = updated[association.feature_path.front()];
      feature.addIDQuery(association.query);
      for (auto match : association.matches)
        feature.addIDMatch({association.query.run_uuid, match});
      feature.getPeptideIdentifications().push_back(std::move(item));
    }
  }
  map = std::move(updated);
  return converted.losses;
}
} // namespace OpenMS
