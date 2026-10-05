// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
#include <OpenMS/ANALYSIS/ID/BasicProteinInferenceAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/IdentificationDataInference.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/METADATA/ID/IdentificationDataAdapter.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>

namespace OpenMS
{
namespace
{
  using ID = IdentificationData;
  [[noreturn]] void invalidInference(const std::string& message)
  { throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message); }
  void pruneUnusedAliases(ID::InferenceResult& result)
  {
    std::set<std::string> used;
    for (const auto& hit : result.proteins.getHits())
      used.insert(hit.getAccession());
    for (const auto* groups : {&result.proteins.getProteinGroups(), &result.proteins.getIndistinguishableProteins()})
      for (const auto& group : *groups)
        used.insert(group.accessions.begin(), group.accessions.end());
    std::erase_if(result.qualified_accessions, [&](const auto& alias) { return ! used.contains(alias.first); });
  }
} // namespace

ID::InferenceResult IdentificationDataInference::infer(const ID& data, const std::vector<Input>& inputs, const std::string& identifier)
{ return infer(data, inputs, identifier, Param {}); }

ID::InferenceResult
IdentificationDataInference::infer(const ID& data, const std::vector<Input>& inputs, const std::string& identifier, const Param& parameters)
{
  if (identifier.empty() || inputs.empty()) invalidInference("Inference requires a result identifier and at least one input run");
  data.validate();
  ID::InferenceResult result;
  result.identifier = identifier;
  std::map<ID::QualifiedAccession, ProteinHit> parent_hits;
  std::set<ID::QualifiedAccession> referenced_parents;
  std::set<std::string> selected_runs;
  std::optional<std::string> calibration;
  for (const auto& input : inputs)
  {
    const auto* run = data.findRunByUuid(input.run_uuid);
    if (! run) invalidInference("Unknown inference input run UUID: " + input.run_uuid);
    if (! selected_runs.insert(input.run_uuid).second) invalidInference("The same run is selected more than once for inference");
    if (run->getMoleculeKind() != ID::MoleculeKind::PEPTIDE) invalidInference("Basic protein inference only accepts peptide runs");
    const auto& definition = run->getScoreDefinition(input.score);
    if (definition.scope != ID::ScoreScope::MATCH) invalidInference("Protein inference requires match-level probability scores");
    const bool is_pep = input.probability == ProbabilityType::POSTERIOR_ERROR_PROBABILITY;
    if (definition.higher_better == is_pep)
      invalidInference("Probability score direction disagrees with the explicitly selected PEP/PP interpretation");
    if (calibration && *calibration != definition.calibration)
      invalidInference("Input score calibration provenance differs; calibrate comparable scores before pooling");
    calibration = definition.calibration;
    if (result.inputs.empty()) result.proteins = run->getProcessingMetadata();
    ID::InferenceInput provenance;
    provenance.run_identifier = run->getIdentifier();
    provenance.run_uuid = run->getUuid();
    provenance.score = definition;
    provenance.selection = "All candidates in the selected run, in stored order";
    result.inputs.push_back(std::move(provenance));
    if (run->getParents())
    {
      for (const auto& parent : *run->getParents())
      {
        ProteinHit hit;
        static_cast<MetaInfoInterface&>(hit) = parent;
        hit.setAccession(parent.identity.accession);
        hit.setSequence(parent.sequence);
        hit.setDescription(parent.description);
        if (parent.target_decoy == ID::TargetDecoy::BOTH)
          invalidInference("BasicProteinInference cannot represent a combined target/decoy protein state");
        const auto state = parent.target_decoy == ID::TargetDecoy::TARGET  ? ProteinHit::TargetDecoyType::TARGET
                           : parent.target_decoy == ID::TargetDecoy::DECOY ? ProteinHit::TargetDecoyType::DECOY
                                                                           : ProteinHit::TargetDecoyType::UNKNOWN;
        if (hit.getTargetDecoyType() != state) hit.setTargetDecoyType(state);
        const auto existing = parent_hits.find(parent.identity);
        if (existing != parent_hits.end() && existing->second != hit)
          invalidInference("Conflicting parent definitions across inference input runs: " + parent.identity.accession);
        parent_hits[parent.identity] = std::move(hit);
      }
    }
    for (const auto& source : run->getSourceBlocks())
      for (const auto& query : source.identifications)
        for (const auto& match : query.getMatches())
          for (const auto& evidence : match.parent_evidence)
            referenced_parents.insert(evidence.parent);
  }
  // Resolve missing catalogue entries only after every run's catalogue was read.
  // An earlier run with evidence alone must not conflict with a later full record.
  for (const auto& identity : referenced_parents)
  {
    if (! parent_hits.contains(identity))
    {
      ProteinHit hit;
      hit.setAccession(identity.accession);
      parent_hits.emplace(identity, std::move(hit));
    }
  }
  // The legacy algorithm keys proteins by accession. Give cross-database collisions
  // unique aliases while retaining every qualified identity beside the output.
  std::map<std::string, Size> alias_counts;
  for (const auto& [identity, hit] : parent_hits)
    ++alias_counts[identity.accession];
  std::set<std::string> used_aliases;
  for (const auto& [identity, hit] : parent_hits)
    used_aliases.insert(identity.accession);
  std::map<ID::QualifiedAccession, std::string> aliases;
  std::vector<ProteinHit> proteins;
  Size next_alias = 0;
  for (auto& [identity, hit] : parent_hits)
  {
    auto alias = identity.accession;
    if (alias_counts[alias] > 1)
    {
      do
      {
        alias = "__qualified_parent_" + std::to_string(next_alias++);
      } while (used_aliases.contains(alias));
      used_aliases.insert(alias);
    }
    aliases[identity] = alias;
    result.qualified_accessions[alias] = identity;
    hit.setAccession(alias);
    proteins.push_back(std::move(hit));
  }
  result.proteins.setIdentifier(identifier);
  result.proteins.setHits(proteins);
  result.proteins.getProteinGroups().clear();
  result.proteins.getIndistinguishableProteins().clear();
  PeptideIdentificationList peptides;
  const std::string match_key = "__owning_inference_candidate";
  for (Size input_index = 0; input_index < inputs.size(); ++input_index)
  {
    const auto& input = inputs[input_index];
    const auto& run = *data.findRunByUuid(input.run_uuid);
    const auto score = run.bindScore(input.score);
    for (const auto& source : run.getSourceBlocks())
    {
      for (const auto& query : source.identifications)
      {
        PeptideIdentification peptide_id;
        peptide_id.setIdentifier(identifier);
        peptide_id.setScoreType("Posterior Probability");
        peptide_id.setHigherScoreBetter(true);
        for (const auto& match : query.getMatches())
        {
          const auto value = score(match);
          if (! value || ! std::isfinite(*value) || *value < 0.0 || *value > 1.0)
            invalidInference("Inference requires a finite probability in [0, 1] on every candidate");
          auto hit = IdentificationDataAdapter::materializePeptide(run, match, input.score);
          hit.setScore(input.probability == ProbabilityType::POSTERIOR_ERROR_PROBABILITY ? 1.0 - *value : *value);
          auto evidence = hit.getPeptideEvidences();
          for (Size i = 0; i < evidence.size(); ++i)
            evidence[i].setProteinAccession(aliases.at(match.parent_evidence[i].parent));
          hit.setPeptideEvidences(evidence);
          std::set<std::string> unique_parents;
          for (const auto& item : evidence)
            unique_parents.insert(item.getProteinAccession());
          hit.setMetaValue("protein_references", unique_parents.size() == 1 ? "unique" : "non-unique");
          hit.setMetaValue(match_key, static_cast<Int64>(result.assignments.size()));
          peptide_id.insertHit(std::move(hit));
          ID::MatchAssignment assignment;
          assignment.run_identifier = run.getIdentifier();
          assignment.run_uuid = run.getUuid();
          assignment.match = match.getId();
          assignment.input_index = input_index;
          result.assignments.push_back(std::move(assignment));
        }
        if (! peptide_id.getHits().empty()) peptides.push_back(std::move(peptide_id));
      }
    }
  }
  if (peptides.empty()) invalidInference("Protein inference has no candidate matches");
  BasicProteinInferenceAlgorithm algorithm;
  auto configured = parameters;
  const auto defaults = algorithm.getParameters();
  for (auto parameter = configured.begin(); parameter != configured.end(); ++parameter)
    if (! defaults.exists(parameter.getName())) invalidInference("Unknown protein inference parameter: " + parameter.getName());
  configured.setDefaults(defaults);
  // The bridge already selected and normalized the exact score column. A second
  // score switch inside the legacy algorithm would contradict that input contract.
  if (parameters.exists("score_type") && parameters.getValue("score_type").toString() != "")
    invalidInference("Set inference scores through Input; score_type must remain empty");
  configured.setValue("score_type", "");
  algorithm.setParameters(configured);
  std::map<std::pair<std::string, Int>, std::set<std::string>> evidence_contract;
  const bool separate_modifications = configured.getValue("treat_modification_variants_separately").toBool();
  const bool separate_charges = configured.getValue("treat_charge_variants_separately").toBool();
  for (auto& query : peptides)
  {
    query.sort();
    query.getHits().resize(1);
    const auto& hit = query.getHits().front();
    const auto key = std::make_pair(separate_modifications ? hit.getSequence().toString() : hit.getSequence().toUnmodifiedString(),
                                    separate_charges ? hit.getCharge() : 0);
    std::set<std::string> parents;
    for (const auto& evidence : hit.getPeptideEvidences())
      parents.insert(evidence.getProteinAccession());
    const auto [existing, inserted] = evidence_contract.emplace(key, parents);
    if (! inserted && existing->second != parents)
      invalidInference("Conflicting parent mappings for the same inference peptidoform; harmonize search evidence before pooling");
  }
  algorithm.run(peptides, result.proteins);
  for (const auto& peptide_id : peptides)
  {
    for (const auto& hit : peptide_id.getHits())
    {
      auto& assignment = result.assignments.at(static_cast<Size>(static_cast<Int64>(hit.getMetaValue(match_key))));
      for (const auto& evidence : hit.getPeptideEvidences())
        assignment.parents.push_back(result.qualified_accessions.at(evidence.getProteinAccession()));
    }
  }
  ID::ScoreDefinition parent_score;
  parent_score.name = result.proteins.getScoreType();
  parent_score.higher_better = result.proteins.isHigherScoreBetter();
  parent_score.scope = ID::ScoreScope::PROTEIN;
  parent_score.software = "BasicProteinInferenceAlgorithm";
  parent_score.software_version = result.proteins.getInferenceEngineVersion();
  parent_score.calibration = calibration.value_or("");
  parent_score.aggregation = configured.getValue("score_aggregation_method").toString();
  for (auto item = configured.begin(); item != configured.end(); ++item)
    parent_score.parameters.setMetaValue(item.getName(), item->value);
  result.parent_score = parent_score;
  auto group_score = parent_score;
  group_score.scope = ID::ScoreScope::PROTEIN_GROUP;
  result.group_score = group_score;
  pruneUnusedAliases(result);
  return result;
}

void IdentificationDataInference::retainProteins(ID::InferenceResult& result, const std::set<ID::QualifiedAccession>& retained)
{
  auto filtered = result;
  std::set<std::string> removed_aliases;
  for (const auto& hit : filtered.proteins.getHits())
  {
    const auto identity = filtered.qualified_accessions.find(hit.getAccession());
    if (identity == filtered.qualified_accessions.end()) invalidInference("Protein identity is not qualified: " + hit.getAccession());
    if (! retained.contains(identity->second)) removed_aliases.insert(hit.getAccession());
  }
  auto& hits = filtered.proteins.getHits();
  std::erase_if(hits, [&](const auto& hit) { return removed_aliases.contains(hit.getAccession()); });
  auto incomplete = [&](const auto& group) {
    return std::any_of(group.accessions.begin(), group.accessions.end(), [&](const auto& alias) { return removed_aliases.contains(alias); });
  };
  std::erase_if(filtered.proteins.getProteinGroups(), incomplete);
  std::erase_if(filtered.proteins.getIndistinguishableProteins(), incomplete);
  for (auto& assignment : filtered.assignments)
    std::erase_if(assignment.parents, [&](const auto& identity) { return ! retained.contains(identity); });
  // Run-level input provenance and explicit assignment rows remain even when no protein
  // survives. Qualified assignment identities do not depend on the alias dictionary.
  pruneUnusedAliases(filtered);
  result = std::move(filtered);
}
} // namespace OpenMS
