// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/ConsensusMapMergerAlgorithm.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/METADATA/ID/IdentificationDataAdapter.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
#include <unordered_map>
#include <unordered_set>

using namespace std;

namespace OpenMS
{
  ConsensusMapMergerAlgorithm::ConsensusMapMergerAlgorithm() :
      ConsensusMapMergerAlgorithm::DefaultParamHandler("ConsensusMapMergerAlgorithm")
    {
      defaults_.setValue("annotate_origin",
                         "true",
                         "If true, adds a map_index MetaValue to the PeptideIDs to annotate the IDRun they came from.");
      defaults_.setValidStrings("annotate_origin", {"true","false"});
      defaultsToParam_();
    }

  //merge proteins across fractions and replicates
  void ConsensusMapMergerAlgorithm::mergeProteinsAcrossFractionsAndReplicates(ConsensusMap& cmap, const ExperimentalDesign& exp_design) const
  {
    const vector<vector<pair<std::string, unsigned>>> toMerge = exp_design.getConditionToPathLabelVector();

    // one of label-free, labeled_MS1, labeled_MS2
    const std::string & experiment_type = cmap.getExperimentType();

    //Not supported because an ID would need to reference multiple protID runs.
    //we could replicate the ID in the future or allow multiple references.
    bool labelfree = true;
    if (experiment_type != "label-free")
    {
      OPENMS_LOG_WARN << "Merging untested for labelled experiments" << endl;
      labelfree = false;
    }

    //out of the path/label combos, construct sets of map indices to be merged
    unsigned lab(0);
    map<unsigned, unsigned> map_idx_2_rep_batch{};
    for (auto& consHeader : cmap.getColumnHeaders())
    {
      bool found = false;
      if (consHeader.second.metaValueExists("channel_id"))
      {
        lab = static_cast<unsigned int>(consHeader.second.getMetaValue("channel_id")) + 1;
      }
      else
      {
        if (!labelfree)
        {
          OPENMS_LOG_WARN << "No channel id annotated in consensusXML. Assuming one channel." << endl;
        }
        lab = 1;
      }
      pair<std::string, unsigned> path_lab{consHeader.second.filename, lab};

      unsigned repBatchIdx(0);
      for (auto& repBatch : toMerge)
      {
        for (const std::pair<std::string, unsigned>& rep : repBatch)
        {
          if (path_lab == rep)
          {
            map_idx_2_rep_batch[consHeader.first] = repBatchIdx;
            found = true;
            break;
          }
        }
        if (found) break;
        repBatchIdx++;
      }
      if (!found)
      {
        throw Exception::MissingInformation(
            __FILE__,
            __LINE__,
            OPENMS_PRETTY_FUNCTION,
            "ConsensusHeader entry ("
            + consHeader.second.filename + ", "
            + consHeader.second.label + ") could not be matched"
            + " to the given experimental design.");
      }
    }

    mergeProteinIDRuns(cmap, map_idx_2_rep_batch);
  }

  namespace
  {
    using ID = IdentificationData;

    /// The identifier of the legacy protein run of @p run: that of the inference result that covers it, else its own
    std::string legacyRun(const ID& data, const ID::Run& run)
    {
      for (const auto& result : data.getInferenceResults())
      {
        for (const auto& input : result.inputs)
        {
          if (input.run_uuid == run.getUuid())
          {
            return result.proteins.getIdentifier().empty() ? result.identifier : result.proteins.getIdentifier();
          }
        }
      }
      return IdentificationDataAdapter::legacyIdentifier(run);
    }

    /// The legacy protein runs of @p data, as export writes them (with their protein hits)
    vector<ProteinIdentification> legacyRuns(const ID& data)
    {
      IdentificationDataAdapter::ExportOptions options;
      options.loss_policy = IdentificationDataAdapter::LossPolicy::ALLOW;
      return IdentificationDataAdapter::toLegacy(data, options).proteins;
    }

    /// Replace the inference results that cover runs that @p results pool with @p results
    void replaceResults(ID& data, vector<ID::InferenceResult> results)
    {
      std::set<std::string> covered;
      for (const auto& result : results)
      {
        for (const auto& item : result.inputs) covered.insert(item.run_uuid);
      }
      vector<ID::InferenceResult> kept;
      for (const auto& result : data.getInferenceResults())
      {
        if (std::none_of(result.inputs.begin(), result.inputs.end(), [&](const auto& item) { return covered.contains(item.run_uuid); }))
        {
          kept.push_back(result);
        }
      }
      data.clearInferenceResults();
      for (auto& result : kept) data.addInferenceResult(std::move(result));
      for (auto& result : results) data.addInferenceResult(std::move(result));
    }

    ID::InferenceInput inferenceInput(const ID::Run& run)
    {
      ID::InferenceInput result;
      result.run_identifier = run.getIdentifier();
      result.run_uuid = run.getUuid();
      return result;
    }
  } // namespace

  void ConsensusMapMergerAlgorithm::mergeProteinIDRuns(ConsensusMap &cmap,
                                             map<unsigned, unsigned> const &mapIdx_to_new_protIDRun) const
  {
    IdentificationDataConverter::editAsIdentificationData(cmap, [&](ConsensusMap& map) { mergeProteinIDRunsNative_(map, mapIdx_to_new_protIDRun); });
  }

  void ConsensusMapMergerAlgorithm::mergeProteinIDRunsNative_(ConsensusMap& cmap, const map<unsigned, unsigned>& mapIdx_to_new_protIDRun) const
  {
    // one of label-free, labeled_MS1, labeled_MS2
    const std::string & experiment_type = cmap.getExperimentType();

    // Not fully supported yet because an ID would need to reference multiple protID runs.
    // we could replicate the ID in the future or allow multiple references.
    if (experiment_type != "label-free")
    {
      OPENMS_LOG_WARN << "Merging untested for labelled experiments" << endl;
    }

    // For the new runs we need newIDRunIdx -> <[file_origins], [mapIdcs]> once to initialize them with metadata
    map<unsigned, pair<set<std::string>,vector<Int>>> new_idcs;
    for (const auto& new_idx : mapIdx_to_new_protIDRun)
    {
      const auto& new_idcs_insert_it = new_idcs.emplace(new_idx.second, make_pair(set<std::string>(), vector<Int>()));
      new_idcs_insert_it.first->second.first.emplace(cmap.getColumnHeaders().at(new_idx.first).filename);
      new_idcs_insert_it.first->second.second.emplace_back(static_cast<Int>(new_idx.first));
    }
    Size new_size = new_idcs.size();

    if (new_size == 1)
    {
      OPENMS_LOG_WARN << "Number of new protein ID runs is one. Consider using mergeAllProteinRuns for some additional speed." << endl;
    }
    else if (new_size >= cmap.getColumnHeaders().size())
      //This also holds for TMT etc. because map_index is a combination of file and label already.
      // even if IDs from the same file are split and replicated, the resulting runs are never more
    {
      throw Exception::InvalidValue(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          "Number of new protein runs after merging"
          " is bigger or equal to the original ones."
          " Aborting. Nothing would be merged.",StringUtils::toStr(new_size));
    }
    else
    {
      OPENMS_LOG_INFO << "Merging into " << new_size << " protein ID runs." << endl;
    }

    // The legacy protein runs (as export writes them: a run, or the runs that an inference result pools) are merged.
    // Mapping from old run ID std::string to new runIDs indices, i.e. calculate from the file/label pairs (=ColumnHeaders),
    // which ProteinIdentifications need to be merged.
    auto& data = cmap.getIdentificationData();
    vector<ProteinIdentification> old_prot_ids = legacyRuns(data);
    map<std::string, set<Size>> run_id_to_new_run_idcs;
    for (const auto& newidx_to_originset_map_idx_pair : new_idcs)
    {
      for (const auto& old_prot_id : old_prot_ids)
      {
        StringList primary_runs;
        old_prot_id.getPrimaryMSRunPath(primary_runs);
        set<std::string> current_content(primary_runs.begin(), primary_runs.end());
        const set<std::string>& merge_request = newidx_to_originset_map_idx_pair.second.first;
        // if this run is fully covered by a requested merged set, use it for it.
        if (std::includes(merge_request.begin(), merge_request.end(), current_content.begin(), current_content.end()))
        {
          run_id_to_new_run_idcs[old_prot_id.getIdentifier()].emplace(newidx_to_originset_map_idx_pair.first);
        }
      }
    }

    vector<ProteinIdentification> new_prot_ids{new_size};
    unordered_map<ProteinHit,set<Size>,hash_type,equal_type> proteins_collected_hits_runs(0, accessionHash_, accessionEqual_);
    // the runs that each new run pools, in the order of their files
    vector<vector<ID::InferenceInput>> new_run_inputs(new_size);
    for (auto& runid2newrunidcs_pair : run_id_to_new_run_idcs)
    {
      // find old run
      auto it = old_prot_ids.begin();
      for (; it != old_prot_ids.end(); ++it)
      {
        if (it->getIdentifier() == runid2newrunidcs_pair.first)
          break;
      }

      for (const auto& newrunid : runid2newrunidcs_pair.second)
      {
        // go through new runs and fill the proteins and update search settings
        // if first time filling this new run:
        if (new_prot_ids.at(newrunid).getIdentifier().empty())
        {
          //initialize new run
          new_prot_ids[newrunid].setSearchEngine(it->getSearchEngine());
          new_prot_ids[newrunid].setSearchEngineVersion(it->getSearchEngineVersion());
          new_prot_ids[newrunid].setSearchParameters(it->getSearchParameters());
          new_prot_ids[newrunid].setIdentifier("condition" + StringUtils::toStr(newrunid));
        }
        // if not, check consistency
        else
        {
          it->peptideIDsMergeable(new_prot_ids[newrunid], experiment_type);
        }
      }
      // the identifications of a run that several new runs use stay with the last one (as merged peptide
      // identifications refer to the last one)
      for (const auto& run : data.getRuns())
      {
        if (run.getMoleculeKind() == ID::MoleculeKind::PEPTIDE && legacyRun(data, run) == runid2newrunidcs_pair.first)
        {
          new_run_inputs[*runid2newrunidcs_pair.second.rbegin()].push_back(inferenceInput(run));
        }
      }

      //Insert hits into collection with empty set (if not present yet) and
      // add destination run indices
      for (auto& hit : it->getHits())
      {
        const auto& foundIt = proteins_collected_hits_runs.emplace(std::move(hit), set<Size>());
        foundIt.first->second.insert(runid2newrunidcs_pair.second.begin(), runid2newrunidcs_pair.second.end());
      }
      it->getHits().clear(); //not needed anymore and moved anyway
    }

    // copy the protein hits into the destination runs
    for (const auto& protToNewRuns : proteins_collected_hits_runs)
    {
      for (Size runID : protToNewRuns.second)
      {
        new_prot_ids.at(runID).getHits().emplace_back(protToNewRuns.first);
      }
    }

    // A new run pools its runs: an inference result (without inference yet), which export writes as the merged
    // protein run, its files in the order of the runs (each identification's 'id_merge_index' points into them).
    for (const auto& run : data.getRuns())
    {
      if (run.getMoleculeKind() == ID::MoleculeKind::PEPTIDE && run.getNumberOfIdentifications() > 0
          && !run_id_to_new_run_idcs.contains(legacyRun(data, run)))
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Identification run '" + legacyRun(data, run) + "' is not part of a merged run.");
      }
    }
    vector<ID::InferenceResult> results;
    for (Size i = 0; i < new_size; ++i)
    {
      ID::InferenceResult result;
      result.identifier = new_prot_ids[i].getIdentifier();
      result.proteins = std::move(new_prot_ids[i]);
      result.inputs = std::move(new_run_inputs[i]);
      results.push_back(std::move(result));
    }
    replaceResults(data, std::move(results));
  }

  //merge proteins across fractions and replicates
  void ConsensusMapMergerAlgorithm::mergeAllIDRuns(ConsensusMap& cmap) const
  {
    IdentificationDataConverter::editAsIdentificationData(cmap, [&](ConsensusMap& map) { mergeAllIDRunsNative_(map); });
  }

  void ConsensusMapMergerAlgorithm::mergeAllIDRunsNative_(ConsensusMap& cmap) const
  {
    auto& data = cmap.getIdentificationData();
    // the legacy protein runs: runs, or the runs that an inference result pools
    vector<ProteinIdentification> old_prot_runs = legacyRuns(data);
    if (old_prot_runs.size() <= 1)
      return;

    // Everything needs to agree
    checkOldRunConsistency_(old_prot_runs, cmap.getExperimentType());

    ProteinIdentification new_prot_id_run;
    //TODO create better ID
    new_prot_id_run.setIdentifier("merged");
    //TODO merge SearchParams e.g. in case of SILAC
    new_prot_id_run.setSearchEngine(old_prot_runs[0].getSearchEngine());
    new_prot_id_run.setSearchEngineVersion(old_prot_runs[0].getSearchEngineVersion());
    new_prot_id_run.setSearchParameters(old_prot_runs[0].getSearchParameters());
    std::string old_inference_engine = old_prot_runs[0].getInferenceEngine();
    if (!old_inference_engine.empty())
    {
      OPENMS_LOG_WARN << "Inference was already performed on the runs in this ConsensusXML."
                         " Merging their proteins, will invalidate correctness of the inference."
                         " You should redo it.\n";
      // deliberately do not take over old inference settings.
    }

    unordered_set<ProteinHit,hash_type,equal_type> proteins_collected_hits(0, accessionHash_, accessionEqual_);
    typedef std::vector<ProteinHit>::iterator iter_t;
    for (auto& prot_run : old_prot_runs)
    {
      auto& hits = prot_run.getHits();
      proteins_collected_hits.insert(
          std::move_iterator<iter_t>(hits.begin()),
          std::move_iterator<iter_t>(hits.end())
      );
      hits.clear();
    }
    auto& hits = new_prot_id_run.getHits();
    for (auto& prot : proteins_collected_hits)
    {
      hits.emplace_back(std::move(const_cast<ProteinHit&>(prot))); //careful this completely invalidates the set
    }
    proteins_collected_hits.clear();

    // The merged run pools all runs: an inference result (without inference yet), which export writes as the merged
    // protein run, its files in the order of the runs (each identification's 'id_merge_index' points into them).
    ID::InferenceResult result;
    result.identifier = new_prot_id_run.getIdentifier();
    result.proteins = std::move(new_prot_id_run);
    for (const auto& run : data.getRuns())
    {
      if (run.getMoleculeKind() == ID::MoleculeKind::PEPTIDE) result.inputs.push_back(inferenceInput(run));
    }
    replaceResults(data, {std::move(result)});
  }

  bool ConsensusMapMergerAlgorithm::checkOldRunConsistency_(const vector<ProteinIdentification>& protRuns, const std::string& experiment_type) const
  {
    return checkOldRunConsistency_(protRuns, protRuns[0], experiment_type);
  }

  //TODO refactor the next two functions
  bool ConsensusMapMergerAlgorithm::checkOldRunConsistency_(const vector<ProteinIdentification>& protRuns, const ProteinIdentification& ref, const std::string& experiment_type) const
  {
    bool ok = true;
    for (const auto& idRun : protRuns)
    {
      // collect warnings and throw at the end if at least one failed
      ok = ok && ref.peptideIDsMergeable(idRun, experiment_type);
    }
    if (!ok /*&& TODO and no force flag*/)
    {
      throw Exception::MissingInformation(__FILE__,
                                     __LINE__,
                                     OPENMS_PRETTY_FUNCTION,
                                     "Search settings are not matching across IdentificationRuns. "
                                     "See warnings. Aborting..");
    }
    return ok;
  }
} // namespace OpenMS
