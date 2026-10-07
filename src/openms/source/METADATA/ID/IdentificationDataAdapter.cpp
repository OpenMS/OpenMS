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
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
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

  /// Settings metadata naming the legacy protein run of a run that import split off (e.g. "search:score_1").
  const std::string LEGACY_RUN = "identification:legacy_run";
  /// Settings metadata with the score type and direction ("true"/"false") of a legacy protein run that is kept
  /// without an inference result, if they are not those of the run's primary score (see settingsToLegacy()).
  const std::string LEGACY_PROTEIN_SCORE_TYPE = "identification:legacy_protein_score_type";
  const std::string LEGACY_PROTEIN_HIGHER_BETTER = "identification:legacy_protein_higher_score_better";

  /// Add the database of a legacy search to @p run, if the search names one.
  void addLegacyDatabase(ID::Run& run, const SearchParameters& search)
  {
    const auto database = Adapter::databaseFromLegacy(search);
    if (database != ID::Database {}) run.addDatabase(database);
  }

  /// The database of a legacy run (which has at most one); an unknown database if its search names none.
  ID::DatabaseId legacyDatabase(ID::Run& run)
  { return run.getDatabases().empty() ? run.addDatabase(ID::Database {}) : run.getDatabaseId(0); }

  /// Whether @p metadata has an empty legacy file list ('spectra_data'), which a legacy run without files may carry.
  bool emptyFileList(const MetaInfoInterface& metadata)
  {
    if (! metadata.metaValueExists("spectra_data")) return false;
    const auto& value = metadata.getMetaValue("spectra_data");
    return value.valueType() == DataValue::STRING_LIST && value.toStringList().empty();
  }

  void loss(Adapter::LegacyResult& result, const Adapter::ExportOptions& options, const std::string& message)
  {
    if (options.loss_policy == Adapter::LossPolicy::STRICT) invalid(message);
    if (std::find(result.losses.begin(), result.losses.end(), message) == result.losses.end()) result.losses.push_back(message);
  }

  /// The legacy meta value of a target/decoy state, as PeptideHit and ProteinHit write it ("" if unknown).
  std::string targetDecoyText(ID::TargetDecoy value)
  {
    switch (value)
    {
      case ID::TargetDecoy::TARGET:
        return "target";
      case ID::TargetDecoy::DECOY:
        return "decoy";
      case ID::TargetDecoy::BOTH:
        return "target+decoy";
      default:
        return "";
    }
  }
  /// Remove a legacy meta value that a field of the owning model holds, so it is not stored twice. Export writes the
  /// field back as @p restored; a value it would write differently (other spelling or value type) is kept.
  void dropRestoredMetaValue(MetaInfoInterface& metadata, const std::string& name, const std::string& restored)
  {
    if (restored.empty() || ! metadata.metaValueExists(name)) return;
    const auto& value = metadata.getMetaValue(name);
    if (value.valueType() == DataValue::STRING_VALUE && value.toString() == restored) metadata.removeMetaValue(name);
  }
  /// The target/decoy state of a legacy peptide or protein hit, read from 'target_decoy' as PeptideHit and ProteinHit
  /// read it, but without throwing: a value they do not know (e.g. an empty one) is UNKNOWN and stays metadata.
  ID::TargetDecoy legacyTargetDecoy(const MetaInfoInterface& hit, bool protein)
  {
    if (! hit.metaValueExists("target_decoy")) return ID::TargetDecoy::UNKNOWN;
    const auto text = StringUtils::toLowered(hit.getMetaValue("target_decoy").toString());
    if (text == "target") return ID::TargetDecoy::TARGET;
    if (text == "decoy") return ID::TargetDecoy::DECOY;
    if (! protein && text == "target+decoy") return ID::TargetDecoy::BOTH;
    return ID::TargetDecoy::UNKNOWN;
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
    for (const auto& item : match.sequence_evidence)
    {
      const auto max_position = static_cast<UInt64>(std::numeric_limits<Int>::max());
      if ((item.start && *item.start > max_position) || (item.end && *item.end > max_position) || item.before.size() > 1
          || item.after.size() > 1)
        invalid("Sequence evidence cannot be represented by legacy peptide coordinates/flanking residues");
      evidence.emplace_back(item.accession, item.start ? static_cast<Int>(*item.start) : PeptideEvidence::UNKNOWN_POSITION,
                            item.end ? static_cast<Int>(*item.end) : PeptideEvidence::UNKNOWN_POSITION,
                            item.before.empty() ? PeptideEvidence::UNKNOWN_AA : item.before.front(),
                            item.after.empty() ? PeptideEvidence::UNKNOWN_AA : item.after.front());
    }
    hit.setPeptideEvidences(evidence);
    if (legacyTargetDecoy(hit, false) != match.target_decoy)
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

  ProteinIdentification legacyProteins(const ID::Run& run, Adapter::LegacyResult& result, const Adapter::ExportOptions& options)
  {
    auto proteins = Adapter::settingsToLegacy(run);
    if (run.getDatabases().size() > 1) loss(result, options, "Legacy export can name only one database per run: " + run.getIdentifier());
    proteins.setIdentifier(proteins.metaValueExists(LEGACY_RUN) ? proteins.getMetaValue(LEGACY_RUN).toString() : run.getIdentifier());
    proteins.removeMetaValue(LEGACY_RUN);
    // Without an inference result, the legacy run takes the primary score as its score type, as search engines write it,
    // unless import recorded another one.
    if (run.getPrimaryScore() && ! run.getSettings().metaValueExists(LEGACY_PROTEIN_SCORE_TYPE))
    {
      proteins.setScoreType(run.getScoreDefinition(*run.getPrimaryScore()).name);
      proteins.setHigherScoreBetter(run.getScoreDefinition(*run.getPrimaryScore()).higher_better);
    }
    const auto files = Adapter::legacyFiles(run);
    if (! files.empty()) proteins.setPrimaryMSRunPath(files);
    if (run.getDatabaseSequences())
    {
      std::vector<ProteinHit> hits;
      for (const auto& sequence : *run.getDatabaseSequences())
      {
        ProteinHit hit;
        static_cast<MetaInfoInterface&>(hit) = sequence;
        hit.setAccession(sequence.accession);
        hit.setSequence(sequence.sequence);
        // ProteinHit stores the description as metadata; an empty one is not written.
        if (! sequence.description.empty()) hit.setDescription(sequence.description);
        if (sequence.target_decoy == ID::TargetDecoy::BOTH)
          loss(result, options, "Legacy proteins cannot represent a combined target/decoy database sequence");
        else
        {
          const auto state = sequence.target_decoy == ID::TargetDecoy::TARGET  ? ProteinHit::TargetDecoyType::TARGET
                             : sequence.target_decoy == ID::TargetDecoy::DECOY ? ProteinHit::TargetDecoyType::DECOY
                                                                               : ProteinHit::TargetDecoyType::UNKNOWN;
          if (legacyTargetDecoy(hit, true) != sequence.target_decoy) hit.setTargetDecoyType(state);
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

  /// The identifier of the legacy protein run that @p run was imported from (runs split off by import name it).
  std::string legacyRunName(const ID::Run& run)
  {
    const auto& settings = run.getSettings();
    return settings.metaValueExists(LEGACY_RUN) ? settings.getMetaValue(LEGACY_RUN).toString() : run.getIdentifier();
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
    merged.proteins.setHits(Adapter::proteinHits(data, inference));
    if (merged.proteins.getIdentifier().empty()) merged.proteins.setIdentifier(inference.identifier);
    StringList files;
    Size file_groups = 0; ///< Joined runs with files of their own in the merged list
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
        const auto& reference = joined.front()->getSettings();
        const auto& settings = run->getSettings();
        const auto reference_search = Adapter::settingsToLegacy(*joined.front()).getSearchParameters();
        const auto search = Adapter::settingsToLegacy(*run).getSearchParameters();
        if (settings.software != reference.software || settings.software_version != reference.software_version
            || ! search.mergeable(reference_search, "label-free"))
        {
          loss(result, options,
               "Pooled inference input run " + run->getIdentifier() + " cannot share a legacy protein run with " + joined.front()->getIdentifier()
                 + "; it is exported without the inference result " + inference.identifier);
          continue;
        }
        if (search != reference_search)
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
      // Runs that import split off one legacy run (identifications of another score type) have its files: they
      // share them instead of repeating them.
      const auto run_files = Adapter::legacyFiles(*run);
      const auto split = std::find_if(joined.begin(), joined.end(), [&](const ID::Run* other) {
        return legacyRunName(*other) == legacyRunName(*run) && Adapter::legacyFiles(*other) == run_files;
      });
      if (split != joined.end()) merged.file_offsets[run->getUuid()] = merged.file_offsets.at((*split)->getUuid());
      else
      {
        merged.file_offsets[run->getUuid()] = files.size();
        files.insert(files.end(), run_files.begin(), run_files.end());
        ++file_groups;
      }
      joined.push_back(run);
    }
    if (file_groups > 1)
    {
      for (const auto* run : joined)
        if (Adapter::legacyFiles(*run).empty())
          loss(result, options,
               "Pooled inference input run " + run->getIdentifier() + " has no primary MS file that identifies its PSMs in the merged legacy protein run");
    }
    if (! files.empty()) merged.proteins.setPrimaryMSRunPath(files);
    else if (! joined.empty() && emptyFileList(joined.front()->getSettings())) merged.proteins.setMetaValue("spectra_data", DataValue(StringList()));
    else merged.proteins.removeMetaValue("spectra_data");
    if (inference.protein_score) writeScoreDefinition(merged.proteins, "protein_score", *inference.protein_score);
    if (inference.group_score) writeScoreDefinition(merged.proteins, "group_score", *inference.group_score);
    // Legacy proteins are identified by accession within the protein run's database.
    for (const auto& [alias, identity] : inference.qualified_accessions)
      if (alias != identity.accession || identity.database != merged.proteins.getSearchParameters().db)
        loss(result, options, "Legacy export cannot represent proteins of several databases in the inference result " + inference.identifier);
    return merged;
  }

  /// A query of a dataset at its position: the run and source that own it.
  struct QueryPosition
  {
    Size run;    ///< Index of the run in the dataset
    Size source; ///< Index of the source in the run
    const ID::Identification* query;
  };

  /**
    The queries of @p data in the order in which export writes them: by query ID across all runs if no ID occurs in
    two runs, as after an import, which numbers the queries in their legacy order; otherwise run by run (runs that
    were created independently, e.g. merged), the queries of each run by ID. Either way the queries of a run follow
    their IDs across its sources, i.e. the order in which they were added.
  */
  std::vector<QueryPosition> queryOrder(const ID& data)
  {
    std::vector<QueryPosition> order;
    for (Size r = 0; r < data.getRuns().size(); ++r)
    {
      const auto& sources = data.getRuns()[r].getSources();
      for (Size s = 0; s < sources.size(); ++s)
        for (const auto& query : sources[s].identifications)
          order.push_back({r, s, &query});
    }
    std::stable_sort(order.begin(), order.end(), [](const auto& a, const auto& b) { return a.query->getId() < b.query->getId(); });
    const bool shared_ids = std::adjacent_find(order.begin(), order.end(), [](const auto& a, const auto& b) {
                              return a.query->getId() == b.query->getId();
                            }) != order.end();
    if (shared_ids)
      std::sort(order.begin(), order.end(),
                [](const auto& a, const auto& b) { return std::make_pair(a.run, a.query->getId()) < std::make_pair(b.run, b.query->getId()); });
    return order;
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

  /// Whether export rebuilds the legacy protein run @p original exactly from the settings and database sequences
  /// of @p run, so that an inference result holding @p original would only repeat them.
  bool rebuildsLegacyProteins(const ID::Run& run, const ProteinIdentification& original)
  {
    Adapter::LegacyResult scratch;
    Adapter::ExportOptions options;
    options.loss_policy = Adapter::LossPolicy::ALLOW;
    const auto rebuilt = legacyProteins(run, scratch, options);
    return scratch.losses.empty() && rebuilt == original;
  }

  /**
    Whether @p run represents the legacy protein run @p original without an inference result. The score type
    of a protein run without protein scores is often empty, or still that of the search engine after PSM
    rescoring; if that is all export cannot rebuild, it is recorded in the run settings.
  */
  bool representsLegacyProteins(ID::Run& run, const ProteinIdentification& original)
  {
    if (rebuildsLegacyProteins(run, original)) return true;
    const auto previous = run.getSettings();
    auto settings = previous;
    settings.setMetaValue(LEGACY_PROTEIN_SCORE_TYPE, original.getScoreType());
    settings.setMetaValue(LEGACY_PROTEIN_HIGHER_BETTER, original.isHigherScoreBetter() ? "true" : "false");
    run.setSettings(settings);
    if (rebuildsLegacyProteins(run, original)) return true;
    run.setSettings(previous);
    return false;
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
  std::map<std::string, bool> catalogues;
  struct InferenceScores
  {
    std::optional<ID::ScoreDefinition> protein, group, input;
  };
  std::map<std::string, InferenceScores> inference_scores;
  for (auto protein : proteins)
  {
    if (originals.contains(protein.getIdentifier())) invalid("Duplicate legacy protein run identifier: " + protein.getIdentifier());
    // Score definitions of an exported inference result belong to that result, not to the run.
    inference_scores[protein.getIdentifier()] = {takeScoreDefinition(protein, "protein_score"), takeScoreDefinition(protein, "group_score"),
                                                 takeScoreDefinition(protein, "input_score")};
    auto params = protein.getSearchParameters();
    ModificationDefinitionIO::attach(params, definitions[protein.getIdentifier()]);
    protein.setSearchParameters(params);
    file_counts[protein.getIdentifier()] = protein.nrPrimaryMSRunPaths();
    // Database sequences need distinct accessions. A protein list with empty or repeated accessions (which legacy
    // runs allow) has no catalogue; its inference result keeps the complete protein hits instead.
    std::set<std::string> accessions;
    catalogues[protein.getIdentifier()] = std::all_of(protein.getHits().begin(), protein.getHits().end(), [&](const ProteinHit& hit) {
      return ! hit.getAccession().empty() && accessions.insert(hit.getAccession()).second;
    });
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
    const auto& configuration = original->second;
    auto settings = settingsFromLegacy(configuration);
    if (name != original_id) settings.setMetaValue(LEGACY_RUN, original_id);
    run.setSettings(settings);
    // The files of the legacy run become the sources of the run.
    StringList files;
    configuration.getPrimaryMSRunPath(files);
    addLegacySources(run, files);
    // The database of the legacy run's search becomes the database of the run.
    addLegacyDatabase(run, configuration.getSearchParameters());
    std::vector<ID::DatabaseSequence> sequences;
    for (const auto& hit : catalogues.at(original_id) ? original->second.getHits() : std::vector<ProteinHit> {})
    {
      ID::DatabaseSequence sequence;
      static_cast<MetaInfoInterface&>(sequence) = hit;
      sequence.database = legacyDatabase(run);
      sequence.accession = hit.getAccession();
      sequence.sequence = hit.getSequence();
      sequence.description = hit.getDescription();
      // ProteinHit keeps its description as the meta value "Description"; the field holds it.
      dropRestoredMetaValue(sequence, "Description", sequence.description);
      sequence.target_decoy = legacyTargetDecoy(hit, true);
      dropRestoredMetaValue(sequence, "target_decoy", targetDecoyText(sequence.target_decoy));
      sequences.push_back(std::move(sequence));
    }
    if (catalogues.at(original_id)) run.setDatabaseSequences(std::move(sequences));
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

  // The runs follow the legacy protein runs, so export lists the protein runs in their original order. A run
  // takes the score of the first identification that refers to it; identifications with other scores get runs of their own.
  std::map<std::string, const PeptideIdentification*> first_peptides;
  for (const auto& item : peptides)
    first_peptides.try_emplace(item.getIdentifier(), &item);
  for (const auto& protein : proteins)
  {
    // Protein-only runs are legitimate and remain representable after an empty search.
    const auto first = first_peptides.find(protein.getIdentifier());
    if (first == first_peptides.end()) create_run({protein.getIdentifier(), "", true});
    else create_run({protein.getIdentifier(), first->second->getScoreType(), first->second->isHigherScoreBetter()});
  }

  for (Size index = 0; index < peptides.size(); ++index)
  {
    const auto& item = peptides[index];
    auto& run = create_run({item.getIdentifier(), item.getScoreType(), item.isHigherScoreBetter()});
    // Multiple files without an explicit index leave the file unknown. Neither basename
    // matching nor consensus map indices resolve this.
    const auto source = legacySource(run, file_counts.at(item.getIdentifier()), item);
    ID::Observation observation;
    static_cast<MetaInfoInterface&>(observation) = item;
    // The source is the file, so export writes the index into the legacy file list of a run with several files;
    // in a single-file run, where export writes none, an index (which can only name that file) stays metadata.
    if (file_counts.at(item.getIdentifier()) > 1) observation.removeMetaValue(Constants::UserParam::ID_MERGE_INDEX);
    observation.data_id = item.getSpectrumReference();
    dropRestoredMetaValue(observation, Constants::UserParam::SPECTRUM_REFERENCE, observation.data_id);
    if (item.hasRT()) observation.rt = item.getRT();
    if (item.hasMZ()) observation.mz = item.getMZ();
    // Query IDs follow the legacy order across all runs, so export restores that order (see queryOrder()).
    auto query = run.importIdentification(source, ID::QueryId {index + 1}, std::move(observation));
    result.queries.push_back({run.getUuid(), query});
    for (const auto& hit : item.getHits())
    {
      if (! run.getPrimaryScore()) invalid("A legacy peptide hit requires a nonempty score type");
      ID::MatchData match;
      static_cast<MetaInfoInterface&>(match) = hit;
      match.representation = hit.getSequence().toString();
      match.charge = hit.getCharge();
      match.target_decoy = legacyTargetDecoy(hit, false);
      dropRestoredMetaValue(match, "target_decoy", targetDecoyText(match.target_decoy));
      match.peak_annotations = hit.getPeakAnnotations();
      for (const auto& item_evidence : hit.getPeptideEvidences())
      {
        ID::SequenceEvidence evidence;
        evidence.database = legacyDatabase(run);
        evidence.accession = item_evidence.getProteinAccession();
        if (item_evidence.getStart() != PeptideEvidence::UNKNOWN_POSITION)
        {
          if (item_evidence.getStart() < 0) invalid("Unsupported negative sequence evidence start");
          evidence.start = static_cast<UInt64>(item_evidence.getStart());
        }
        if (item_evidence.getEnd() != PeptideEvidence::UNKNOWN_POSITION)
        {
          if (item_evidence.getEnd() < 0) invalid("Unsupported negative sequence evidence end");
          evidence.end = static_cast<UInt64>(item_evidence.getEnd());
        }
        evidence.before = std::string(1, item_evidence.getAABefore());
        evidence.after = std::string(1, item_evidence.getAAAfter());
        match.sequence_evidence.push_back(std::move(evidence));
      }
      run.addMatch(query, match, {hit.getScore()});
    }
  }
  std::vector<ProteinHit> empty_hits;
  for (const auto& protein : proteins)
  {
    const auto& name = protein.getIdentifier();
    const auto& original = originals.at(name);
    // Search engines list the proteins of their matches, with no inference: if export rebuilds that list from
    // the run's database sequences, the run needs no inference result. (A protein run without peptide
    // identifications keeps one.)
    const auto& scores = inference_scores[name];
    if (first_peptides.contains(name) && input_runs[name].size() == 1 && ! scores.protein && ! scores.group && ! scores.input
        && representsLegacyProteins(result.data.getRun(input_runs[name].front()), original))
      continue;
    ID::InferenceResult inference;
    inference.identifier = "legacy:" + name;
    inference.proteins = original;
    // The files of the inference result are those of its input runs.
    inference.proteins.removeMetaValue("spectra_data");
    // The run's database sequences hold sequence, description and metadata of the proteins; the hits keep their
    // inference values and target/decoy state (protein FDR reads it), and export completes them (proteinHits()).
    for (auto& hit : catalogues.at(name) ? inference.proteins.getHits() : empty_hits)
    {
      hit.setSequence("");
      hit.setDescription("");
      const auto state = hit.metaValueExists("target_decoy") ? std::optional<DataValue>(hit.getMetaValue("target_decoy")) : std::nullopt;
      hit.clearMetaInfo();
      if (state) hit.setMetaValue("target_decoy", *state);
    }
    inference.protein_score = inference_scores[name].protein;
    inference.group_score = inference_scores[name].group;
    for (const auto& hit : original.getHits())
      if (! hit.getAccession().empty()) inference.qualified_accessions[hit.getAccession()] = {original.getSearchParameters().db, hit.getAccession()};
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
  // Queries and candidates were appended one by one.
  std::vector<std::string> names;
  for (const auto& run : result.data.getRuns())
    names.push_back(run.getIdentifier());
  for (const auto& name : names)
    result.data.getRun(name).shrinkToFit();
  result.data.validate();
  return result;
}

ID::RunSettings IdentificationDataAdapter::settingsFromLegacy(const ProteinIdentification& proteins)
{
  ID::RunSettings settings;
  static_cast<MetaInfoInterface&>(settings) = proteins;
  // The files are the sources of the run. An empty list names none and stays, so export can write it again.
  if (! emptyFileList(settings)) settings.removeMetaValue("spectra_data");
  settings.software = proteins.getSearchEngine();
  settings.software_version = proteins.getSearchEngineVersion();
  settings.date = proteins.getDateTime();
  settings.search = proteins.getSearchParameters();
  settings.search.db.clear();
  settings.search.db_version.clear();
  settings.search.taxonomy.clear();
  return settings;
}

ID::Database IdentificationDataAdapter::databaseFromLegacy(const SearchParameters& search)
{
  ID::Database database;
  database.path = search.db;
  database.version = search.db_version;
  database.taxonomy = search.taxonomy;
  return database;
}

ProteinIdentification IdentificationDataAdapter::settingsToLegacy(const ID::Run& run)
{
  const auto& settings = run.getSettings();
  ProteinIdentification proteins;
  static_cast<MetaInfoInterface&>(proteins) = settings;
  proteins.setSearchEngine(settings.software);
  proteins.setSearchEngineVersion(settings.software_version);
  proteins.setDateTime(settings.date);
  auto search = settings.search;
  if (! run.getDatabases().empty())
  {
    const auto& database = run.getDatabases().front();
    search.db = database.path;
    search.db_version = database.version;
    search.taxonomy = database.taxonomy;
  }
  proteins.setSearchParameters(search);
  // The score type of a legacy protein run that import kept without an inference result
  if (proteins.metaValueExists(LEGACY_PROTEIN_SCORE_TYPE))
  {
    proteins.setScoreType(proteins.getMetaValue(LEGACY_PROTEIN_SCORE_TYPE).toString());
    proteins.removeMetaValue(LEGACY_PROTEIN_SCORE_TYPE);
  }
  if (proteins.metaValueExists(LEGACY_PROTEIN_HIGHER_BETTER))
  {
    proteins.setHigherScoreBetter(proteins.getMetaValue(LEGACY_PROTEIN_HIGHER_BETTER).toString() != "false");
    proteins.removeMetaValue(LEGACY_PROTEIN_HIGHER_BETTER);
  }
  return proteins;
}

std::vector<ProteinHit> IdentificationDataAdapter::proteinHits(const ID& data, const ID::InferenceResult& result)
{
  std::map<ID::QualifiedAccession, const ID::DatabaseSequence*> catalog;
  for (const auto& input : result.inputs)
  {
    const auto* run = data.findRunByUuid(input.run_uuid);
    if (! run || ! run->getDatabaseSequences()) continue;
    for (const auto& sequence : *run->getDatabaseSequences())
      catalog.try_emplace(run->qualify(sequence.database, sequence.accession), &sequence);
  }
  auto hits = result.proteins.getHits();
  if (catalog.empty()) return hits;
  for (auto& hit : hits)
  {
    const auto identity = result.qualified_accessions.find(hit.getAccession());
    if (identity == result.qualified_accessions.end()) continue;
    const auto found = catalog.find(identity->second);
    if (found == catalog.end()) continue;
    const auto& sequence = *found->second;
    if (hit.getSequence().empty()) hit.setSequence(sequence.sequence);
    // ProteinHit stores the description as metadata; an empty one is not written.
    if (hit.getDescription().empty() && ! sequence.description.empty()) hit.setDescription(sequence.description);
    std::vector<std::string> keys;
    sequence.getKeys(keys);
    for (const auto& key : keys)
      if (! hit.metaValueExists(key)) hit.setMetaValue(key, sequence.getMetaValue(key));
  }
  return hits;
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
  const auto& sources = run.getSources();
  const auto unknown = std::find_if(sources.begin(), sources.end(), [](const auto& source) { return source.file.path.empty(); });
  if (unknown != sources.end()) return unknown->id;
  return run.addSource(ID::SourceFile {});
}

StringList IdentificationDataAdapter::legacyFiles(const ID::Run& run)
{
  StringList files;
  for (const auto& source : run.getSources())
    if (! source.file.path.empty()) files.push_back(source.file.path);
  return files;
}

IdentificationData IdentificationDataAdapter::fromLegacy(const std::vector<ProteinIdentification>& proteins,
                                                         const PeptideIdentificationList& peptides)
{ return importLegacy(proteins, peptides).data; }

PeptideHit IdentificationDataAdapter::materializePeptide(const ID::Run& run, const ID::Match& match, ID::ScoreId score)
{
  registerDefinitions(run.getSettings().search);
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
  // How the identifications of an exported run refer to its legacy protein run.
  struct RunExport
  {
    std::string identifier;  ///< Of the legacy protein run
    std::string database;    ///< Of the legacy protein run
    std::optional<ID::ScoreId> primary;
    bool merge_index = false; ///< Whether the legacy protein run has several files, which PSMs refer to by id_merge_index
    /// Per source: the index of its file in the legacy file list (in a merged protein run, after the files of
    /// the runs before), or none for a source without a path.
    std::vector<std::optional<Size>> files;
  };
  std::vector<std::optional<RunExport>> exports(data.getRuns().size());
  for (Size run_index = 0; run_index < data.getRuns().size(); ++run_index)
  {
    const auto& run = data.getRuns()[run_index];
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
    auto proteins = shared ? shared->proteins : legacyProteins(run, result, options);
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
      // Same legacy identifier is only shared when its complete original payload agrees. Otherwise (e.g. runs of
      // two bundles merged after a split, which kept the identifier of their original run) the run is exported
      // under its own name: a legacy identifier only links peptide to protein identifications, so nothing is lost.
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
    registerDefinitions(run.getSettings().search);
    // A source with a path is the next file of the run's legacy file list; its identifications point
    // to it with id_merge_index if the legacy run has several files.
    RunExport exported {proteins.getIdentifier(), proteins.getSearchParameters().db, primary, legacy_paths.size() > 1, {}};
    Size file_index = 0;
    for (const auto& source : run.getSources())
    {
      const bool known = ! source.file.path.empty();
      if (! source.file.isMetaEmpty()) loss(result, options, "Legacy export cannot retain source-level metadata: " + run.getIdentifier());
      if (! known && legacy_paths.size() == 1)
        loss(result, options, "An unknown source cannot be represented in a legacy run with exactly one known source: " + run.getIdentifier());
      exported.files.push_back(known ? std::optional<Size>(file_offset + file_index++) : std::nullopt);
    }
    exports[run_index] = std::move(exported);
  }
  // The identifications in the order of their IDs, which is their legacy order after an import (see queryOrder()).
  for (const auto& position : queryOrder(data))
  {
    if (! exports[position.run]) continue;
    const auto& run = data.getRuns()[position.run];
    const auto& exported = *exports[position.run];
    const auto& primary = exported.primary;
    const auto& query = *position.query;
    PeptideIdentification item;
    static_cast<MetaInfoInterface&>(item) = query;
    item.setIdentifier(exported.identifier);
    if (primary)
    {
      item.setScoreType(run.getScoreDefinition(*primary).name);
      item.setHigherScoreBetter(run.getScoreDefinition(*primary).higher_better);
    }
    if (query.rt) item.setRT(*query.rt);
    if (query.mz) item.setMZ(*query.mz);
    if (! query.data_id.empty()) item.setSpectrumReference(query.data_id);
    // The source decides the file: a run with several files writes the index of its file, and an index in the
    // metadata (e.g. of a single-file run) stays only if it names that file.
    const auto& file = exported.files[position.source];
    if (exported.merge_index || ! file
        || (item.metaValueExists(Constants::UserParam::ID_MERGE_INDEX)
            && (item.getMetaValue(Constants::UserParam::ID_MERGE_INDEX).valueType() != DataValue::INT_VALUE
                || static_cast<Int64>(item.getMetaValue(Constants::UserParam::ID_MERGE_INDEX)) != static_cast<Int64>(*file))))
      item.removeMetaValue(Constants::UserParam::ID_MERGE_INDEX);
    if (file && exported.merge_index) item.setMetaValue(Constants::UserParam::ID_MERGE_INDEX, static_cast<Int64>(*file));
    if (query.getSelectedMatch()) loss(result, options, "Legacy export cannot preserve an explicit selected candidate: " + run.getIdentifier());
    for (const auto& match : query.getMatches())
    {
      if (match.calculated_mz || match.adduct || match.formula || ! match.name.empty() || ! match.identifiers.empty())
        loss(result, options, "Legacy export cannot preserve all molecular/ion fields: " + run.getIdentifier());
      for (const auto& evidence : match.sequence_evidence)
      {
        if (run.getDatabase(evidence.database).path != exported.database)
          loss(result, options, "Legacy export cannot represent the database of sequence evidence: " + run.getIdentifier());
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
      for (const auto& source : run.getSources())
        for (const auto& query : source.identifications)
          for (const auto& match : query.getMatches())
            owners.emplace(ID::MatchReference {run.getUuid(), match.getId()}, ID::QueryReference {run.getUuid(), query.getId()});
    std::vector<Adapter::FeatureAssociation> associations;
    std::map<ID::QueryReference, std::set<ID::MatchId>> assigned;
    std::set<ID::QueryReference> assigned_queries;
    std::map<std::string, Size> run_positions;
    for (const auto& run : data.getRuns())
      run_positions.emplace(run.getUuid(), run_positions.size());
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
      // The identifications of a feature in the order of their IDs, then of their runs: the legacy order of an import.
      std::vector<std::pair<ID::QueryReference, std::set<ID::MatchId>>> ordered(linked.begin(), linked.end());
      const auto key = [&](const auto& item) {
        const auto run = run_positions.find(item.first.run_uuid);
        return std::make_pair(item.first.query, run == run_positions.end() ? run_positions.size() : run->second);
      };
      std::stable_sort(ordered.begin(), ordered.end(), [&](const auto& a, const auto& b) { return key(a) < key(b); });
      for (const auto& [query, matches] : ordered)
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
    for (const auto& position : queryOrder(data))
    {
      const auto& query = *position.query;
      ID::QueryReference reference {data.getRuns()[position.run].getUuid(), query.getId()};
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
