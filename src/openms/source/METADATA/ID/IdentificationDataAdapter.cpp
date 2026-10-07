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
#include <map>
#include <optional>
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
    const auto files = Adapter::legacyFiles(run);
    if (! files.empty()) proteins.setPrimaryMSRunPath(files);
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

  // Inference score definitions have no legacy field. They travel as prefixed metadata of the
  // legacy protein run, and import removes them again.
  const std::string INFERENCE_SCORE_PREFIX = "identification:inference:";

  void writeScoreDefinition(MetaInfoInterface& target, const std::string& role, const ID::ScoreDefinition& definition)
  {
    const auto prefix = INFERENCE_SCORE_PREFIX + role + ":";
    target.setMetaValue(prefix + "name", definition.name);
    target.setMetaValue(prefix + "higher_better", definition.higher_better ? "true" : "false");
    target.setMetaValue(prefix + "scope", static_cast<Int>(definition.scope));
    for (const auto& [field, value] : {std::pair {"accession", &definition.accession},
                                       std::pair {"software", &definition.software},
                                       std::pair {"software_version", &definition.software_version},
                                       std::pair {"calibration", &definition.calibration},
                                       std::pair {"aggregation", &definition.aggregation}})
      if (! value->empty()) target.setMetaValue(prefix + field, *value);
    std::vector<std::string> keys;
    definition.parameters.getKeys(keys);
    for (const auto& key : keys)
      target.setMetaValue(prefix + "parameter:" + key, definition.parameters.getMetaValue(key));
  }

  std::optional<ID::ScoreDefinition> takeScoreDefinition(MetaInfoInterface& source, const std::string& role)
  {
    const auto prefix = INFERENCE_SCORE_PREFIX + role + ":";
    std::vector<std::string> keys;
    source.getKeys(keys);
    std::erase_if(keys, [&](const auto& key) { return ! key.starts_with(prefix); });
    if (keys.empty()) return std::nullopt;
    if (! source.metaValueExists(prefix + "name")) invalid("Incomplete legacy inference score definition: " + role);
    ID::ScoreDefinition definition;
    for (const auto& key : keys)
    {
      const auto field = key.substr(prefix.size());
      const auto& value = source.getMetaValue(key);
      if (field.starts_with("parameter:"))
        definition.parameters.setMetaValue(field.substr(std::string("parameter:").size()), value);
      else if (field == "higher_better")
      {
        if (value.toString() != "true" && value.toString() != "false") invalid("Invalid legacy inference score direction: " + key);
        definition.higher_better = value.toString() == "true";
      }
      else if (field == "scope")
      {
        const Int scope = value.valueType() == DataValue::INT_VALUE ? static_cast<Int>(value) : -1;
        if (scope < 0 || scope > static_cast<Int>(ID::ScoreScope::OTHER)) invalid("Invalid legacy inference score scope: " + key);
        definition.scope = static_cast<ID::ScoreScope>(scope);
      }
      else if (field == "name") definition.name = value.toString();
      else if (field == "accession") definition.accession = value.toString();
      else if (field == "software") definition.software = value.toString();
      else if (field == "software_version") definition.software_version = value.toString();
      else if (field == "calibration") definition.calibration = value.toString();
      else if (field == "aggregation") definition.aggregation = value.toString();
      else invalid("Unknown legacy inference score field: " + key);
      source.removeMetaValue(key);
    }
    return definition;
  }

  /// The legacy protein run of an inference result. In the legacy model, inference over several runs
  /// works on one merged protein run (as IDMerger or ProteinInference create it), whose PSMs locate
  /// their file through id_merge_index.
  struct LegacyInference
  {
    ProteinIdentification proteins;
    std::map<std::string, Size> file_offsets; ///< Run UUID of each joined run -> index of its first primary MS file.
  };

  LegacyInference legacyInference(const ID& data, const ID::InferenceResult& inference, Adapter::LegacyResult& result,
                                  const Adapter::ExportOptions& options)
  {
    LegacyInference merged {inference.proteins, {}};
    if (merged.proteins.getIdentifier().empty()) merged.proteins.setIdentifier(inference.identifier);
    StringList files;
    std::vector<const ID::Run*> joined;
    const ID::InferenceInput* first_input = nullptr;
    for (const auto& input : inference.inputs)
    {
      const auto* run = data.findRunByUuid(input.run_uuid);
      // Inputs that are no longer part of the dataset contribute no PSMs.
      if (! run || run->getMoleculeKind() != ID::MoleculeKind::PEPTIDE || merged.file_offsets.contains(run->getUuid())) continue;
      if (first_input)
      {
        // A run only joins if its PSMs keep their search settings; otherwise it stays a separate protein run.
        const auto& reference = joined.front()->getProcessingMetadata();
        const auto& processing = run->getProcessingMetadata();
        if (processing.getSearchEngine() != reference.getSearchEngine() || processing.getSearchEngineVersion() != reference.getSearchEngineVersion()
            || ! processing.getSearchParameters().mergeable(reference.getSearchParameters(), "label-free"))
        {
          loss(result, options,
               "Pooled inference input run " + run->getIdentifier() + " cannot share a legacy protein run with " + joined.front()->getIdentifier()
                 + "; it is exported without the inference result " + inference.identifier);
          continue;
        }
        if (processing.getSearchParameters() != reference.getSearchParameters())
          loss(result, options,
               "Pooled inference input run " + run->getIdentifier() + " differs in search settings from " + joined.front()->getIdentifier()
                 + "; the merged legacy protein run keeps those of " + joined.front()->getIdentifier());
        if (input.score != first_input->score)
          loss(result, options, "Legacy export cannot represent different input scores of the pooled inference result " + inference.identifier);
      }
      else
      {
        first_input = &input;
        if (input.score) writeScoreDefinition(merged.proteins, "input_score", *input.score);
      }
      merged.file_offsets[run->getUuid()] = files.size();
      const auto run_files = Adapter::legacyFiles(*run);
      files.insert(files.end(), run_files.begin(), run_files.end());
      joined.push_back(run);
    }
    if (joined.size() > 1)
    {
      for (const auto* run : joined)
        if (Adapter::legacyFiles(*run).empty())
          loss(result, options,
               "Pooled inference input run " + run->getIdentifier() + " has no primary MS file that identifies its PSMs in the merged legacy protein run");
    }
    if (! files.empty()) merged.proteins.setPrimaryMSRunPath(files);
    else merged.proteins.removeMetaValue("spectra_data");
    if (inference.parent_score) writeScoreDefinition(merged.proteins, "parent_score", *inference.parent_score);
    if (inference.group_score) writeScoreDefinition(merged.proteins, "group_score", *inference.group_score);
    // Legacy proteins are identified by accession within the protein run's database.
    for (const auto& [alias, identity] : inference.qualified_accessions)
      if (alias != identity.accession || identity.database != merged.proteins.getSearchParameters().db)
        loss(result, options, "Legacy export cannot represent parents of several databases in the inference result " + inference.identifier);
    return merged;
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
  std::map<std::string, Size> file_counts;
  struct InferenceScores
  {
    std::optional<ID::ScoreDefinition> parent, group, input;
  };
  std::map<std::string, InferenceScores> inference_scores;
  for (auto protein : proteins)
  {
    if (originals.contains(protein.getIdentifier())) invalid("Duplicate legacy protein run identifier: " + protein.getIdentifier());
    // Score definitions of an exported inference result belong to that result, not to the run.
    inference_scores[protein.getIdentifier()] = {takeScoreDefinition(protein, "parent_score"), takeScoreDefinition(protein, "group_score"),
                                                 takeScoreDefinition(protein, "input_score")};
    auto params = protein.getSearchParameters();
    ModificationDefinitionIO::attach(params, definitions[protein.getIdentifier()]);
    protein.setSearchParameters(params);
    file_counts[protein.getIdentifier()] = protein.nrPrimaryMSRunPaths();
    originals.emplace(protein.getIdentifier(), std::move(protein));
  }
  using Contract = std::tuple<std::string, std::string, bool>;
  std::map<Contract, std::string> contracts;
  std::map<std::string, std::vector<std::string>> input_runs;
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
    // The files of the legacy run become the sources of the run.
    StringList files;
    configuration.getPrimaryMSRunPath(files);
    configuration.removeMetaValue("spectra_data");
    run.setProcessingMetadata(configuration);
    addLegacySources(run, files);
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
    // Multiple files without an explicit index leave the file unknown. Neither basename
    // matching nor consensus map indices resolve this.
    const auto source = legacySource(run, file_counts.at(item.getIdentifier()), item);
    ID::Observation observation;
    static_cast<MetaInfoInterface&>(observation) = item;
    // The source is the file, so the index into the legacy file list is not kept.
    observation.removeMetaValue(Constants::UserParam::ID_MERGE_INDEX);
    observation.data_id = item.getSpectrumReference();
    if (item.hasRT()) observation.rt = item.getRT();
    if (item.hasMZ()) observation.mz = item.getMZ();
    auto query = run.addIdentification(source, observation);
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
    // The files of the inference result are those of its input runs.
    inference.proteins.removeMetaValue("spectra_data");
    inference.parent_score = inference_scores[name].parent;
    inference.group_score = inference_scores[name].group;
    for (const auto& hit : original.getHits())
      inference.qualified_accessions[hit.getAccession()] = {original.getSearchParameters().db, hit.getAccession()};
    for (const auto& run_name : input_runs[name])
    {
      const auto& run = result.data.getRun(run_name);
      ID::InferenceInput input;
      input.run_identifier = run_name;
      input.run_uuid = run.getUuid();
      input.selection = "Imported legacy run-level provenance";
      input.score = inference_scores[name].input;
      inference.inputs.push_back(std::move(input));
    }
    result.data.addInferenceResult(inference);
  }
  result.data.validate();
  return result;
}

void IdentificationDataAdapter::addLegacySources(ID::Run& run, const StringList& files)
{
  for (const auto& file : files)
  {
    ID::SourceFile source;
    source.path = file;
    run.addSource(source);
  }
}

ID::SourceId IdentificationDataAdapter::legacySource(ID::Run& run, Size n_files, const PeptideIdentification& item)
{
  if (item.metaValueExists(Constants::UserParam::ID_MERGE_INDEX))
  {
    const auto& value = item.getMetaValue(Constants::UserParam::ID_MERGE_INDEX);
    if (value.valueType() != DataValue::INT_VALUE) invalid("id_merge_index must be an integer");
    const auto index = static_cast<Int64>(value);
    if (index < 0 || static_cast<UInt64>(index) >= n_files) invalid("id_merge_index is outside the file list of its run");
    return run.getSourceId(static_cast<UInt32>(index));
  }
  if (n_files == 1) return run.getSourceId(0);
  const auto& sources = run.getSourceBlocks();
  const auto unknown = std::find_if(sources.begin(), sources.end(), [](const auto& source) { return source.source.path.empty(); });
  if (unknown != sources.end()) return unknown->id;
  return run.addSource(ID::SourceFile {});
}

StringList IdentificationDataAdapter::legacyFiles(const ID::Run& run)
{
  StringList files;
  for (const auto& source : run.getSourceBlocks())
    if (! source.source.path.empty()) files.push_back(source.source.path);
  return files;
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
  std::map<const ID::InferenceResult*, LegacyInference> legacy_inference;
  std::map<std::string, Size> protein_indices;
  for (const auto& run : data.getRuns())
  {
    if (run.getMoleculeKind() != ID::MoleculeKind::PEPTIDE)
    {
      loss(result, options, "Legacy peptide/protein export cannot represent non-peptide run " + run.getIdentifier());
      continue;
    }
    const LegacyInference* shared = nullptr;
    if (const auto* inference = selectInference(data, run, options))
    {
      auto found = legacy_inference.find(inference);
      if (found == legacy_inference.end()) found = legacy_inference.emplace(inference, legacyInference(data, *inference, result, options)).first;
      if (found->second.file_offsets.contains(run.getUuid()))
      {
        shared = &found->second;
        handled_inference.insert(inference->identifier);
      }
    }
    auto proteins = shared ? shared->proteins : originalParents(run, result, options);
    if (proteins.getIdentifier().empty()) proteins.setIdentifier(run.getIdentifier());
    // The files of the legacy protein run. In a merged protein run, those of this run start at file_offset.
    StringList legacy_paths;
    proteins.getPrimaryMSRunPath(legacy_paths);
    const Size file_offset = shared ? shared->file_offsets.at(run.getUuid()) : 0;
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
    // A source with a path is the next file of the run's legacy file list; its identifications point
    // to it with id_merge_index if the legacy run has several files.
    Size file_index = 0;
    for (const auto& source : run.getSourceBlocks())
    {
      const bool known = ! source.source.path.empty();
      if (! source.source.isMetaEmpty()) loss(result, options, "Legacy export cannot retain source-level metadata: " + run.getIdentifier());
      if (! known && legacy_paths.size() == 1)
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
        // The source decides the file; an index in the metadata is not used.
        item.removeMetaValue(Constants::UserParam::ID_MERGE_INDEX);
        if (known && legacy_paths.size() > 1)
          item.setMetaValue(Constants::UserParam::ID_MERGE_INDEX, static_cast<Int64>(file_index + file_offset));
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
      if (known) ++file_index;
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
