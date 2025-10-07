// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/ConsensusMapMergerAlgorithm.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <unordered_map>

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
    const vector<vector<pair<String, unsigned>>> toMerge = exp_design.getConditionToPathLabelVector();

    // one of label-free, labeled_MS1, labeled_MS2
    const String & experiment_type = cmap.getExperimentType();

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
      pair<String, unsigned> path_lab{consHeader.second.filename, lab};

      unsigned repBatchIdx(0);
      for (auto& repBatch : toMerge)
      {
        for (const std::pair<String, unsigned>& rep : repBatch)
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

  void ConsensusMapMergerAlgorithm::mergeProteinIDRuns(ConsensusMap &cmap,
                                             map<unsigned, unsigned> const &mapIdx_to_new_protIDRun) const
  {
    // one of label-free, labeled_MS1, labeled_MS2
    const String & experiment_type = cmap.getExperimentType();

    // Not fully supported yet because an ID would need to reference multiple protID runs.
    // we could replicate the ID in the future or allow multiple references.
    if (experiment_type != "label-free")
    {
      OPENMS_LOG_WARN << "Merging untested for labelled experiments" << endl;
    }

    // Unfortunately we need a kind of bimap here.
    // For the features we need oldMapIdx -> newIDRunIdx
    // For the new runs we need newIDRunIdx -> <[file_origins], [mapIdcs]> once to initialize them with metadata
    // TODO should we instead better try to collect the primaryMSRuns from the old Runs?
    // TODO I just saw that in the columnHeaders might have the featureXMLs as origins but we should enforce that
    //  this will be changed to the mzML by all tools
    //  Therefore we somehow need to check consistency of ColumnHeaders and ProteinIdentification (file_origins).
    map<unsigned, pair<set<String>,vector<Int>>> new_idcs;
    for (const auto& new_idx : mapIdx_to_new_protIDRun)
    {
      const auto& new_idcs_insert_it = new_idcs.emplace(new_idx.second, make_pair(set<String>(), vector<Int>()));
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
          " Aborting. Nothing would be merged.", String(new_size));
    }
    else
    {
      OPENMS_LOG_INFO << "Merging into " << new_size << " protein ID runs." << endl;
    }

    // Mapping from old run ID String to new runIDs indices, i.e. calculate from the file/label pairs (=ColumnHeaders),
    // which ProteinIdentifications need to be merged.
    map<String, set<Size>> run_id_to_new_run_idcs;
    // this is to check how many old runs contribute to the new runs
    // this can help save time and we can double check
    vector<Size> nr_inputs_for_new_run_ids(new_size, 0);
    for (const auto& newidx_to_originset_map_idx_pair : new_idcs)
    {
      for (auto& old_prot_id : cmap.getProteinIdentifications())
      {
        StringList primary_runs;
        old_prot_id.getPrimaryMSRunPath(primary_runs);
        set<String> current_content(primary_runs.begin(), primary_runs.end());
        const set<String>& merge_request = newidx_to_originset_map_idx_pair.second.first;
        // if this run is fully covered by a requested merged set, use it for it.
        Size count = 1;
        if (std::includes(merge_request.begin(), merge_request.end(), current_content.begin(), current_content.end()))
        {
          auto it = run_id_to_new_run_idcs.emplace(old_prot_id.getIdentifier(), set<Size>());
          if (!it.second)
          {
            OPENMS_LOG_WARN << "Duplicate protein run ID found. Uniquifying it." << endl;
            old_prot_id.setIdentifier(old_prot_id.getIdentifier() + "_" + count);
            it = run_id_to_new_run_idcs.emplace(old_prot_id.getIdentifier(), set<Size>());
          }
          it.first->second.emplace(newidx_to_originset_map_idx_pair.first);
          nr_inputs_for_new_run_ids.at(newidx_to_originset_map_idx_pair.first)++;
        }
      }
    }

    vector<ProteinIdentification> new_prot_ids{new_size};
    //TODO preallocate with sum of proteins?
    unordered_map<ProteinHit,set<Size>,hash_type,equal_type> proteins_collected_hits_runs(0, accessionHash_, accessionEqual_);

    // we only need to store an offset if we append the primaryRunPaths
    //(oldRunID, newRunIdx) -> newMergeIdxOffset
    map<pair<String,Size>, Size> oldrunid_newrunidx_pair2newmergeidx_offset;

    for (auto& runid2newrunidcs_pair : run_id_to_new_run_idcs)
    {
      // find old run
      auto it = cmap.getProteinIdentifications().begin();
      for (; it != cmap.getProteinIdentifications().end(); ++it)
      {
        if (it->getIdentifier() == runid2newrunidcs_pair.first)
          break;
      }

      for (const auto& newrunid : runid2newrunidcs_pair.second)
      {
        // go through new runs and fill the proteins and update search settings
        // if first time filling this new run:
        //TODO safe to check for empty identifier?

        if (new_prot_ids.at(newrunid).getIdentifier().empty())
        {
          //initialize new run
          new_prot_ids[newrunid].setSearchEngine(it->getSearchEngine());
          new_prot_ids[newrunid].setSearchEngineVersion(it->getSearchEngineVersion());
          new_prot_ids[newrunid].setSearchParameters(it->getSearchParameters());
          StringList toFill;
          it->getPrimaryMSRunPath(toFill);
          new_prot_ids[newrunid].setPrimaryMSRunPath(toFill);
          new_prot_ids[newrunid].setIdentifier("condition" + String(newrunid));
          oldrunid_newrunidx_pair2newmergeidx_offset.emplace(std::piecewise_construct,
                                                             std::forward_as_tuple(runid2newrunidcs_pair.first, newrunid),
                                                             std::forward_as_tuple(0));
        }
        // if not, merge settings or check consistency
        else
        {
          //check consistency and add origins
          it->peptideIDsMergeable(new_prot_ids[newrunid], experiment_type);
          Size offset = new_prot_ids[newrunid].nrPrimaryMSRunPaths();
          StringList toFill; it->getPrimaryMSRunPath(toFill); // new ones
          new_prot_ids[newrunid].addPrimaryMSRunPath(toFill); //add to previous
          oldrunid_newrunidx_pair2newmergeidx_offset.emplace(std::piecewise_construct,
                                                             std::forward_as_tuple(runid2newrunidcs_pair.first, newrunid),
                                                             std::forward_as_tuple(offset));
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

    //Now update the references in the PeptideHits
    //TODO double check the PrimaryRunPaths with the initial requested merge

    function<void(PeptideIdentification&)> fun = [&run_id_to_new_run_idcs, &oldrunid_newrunidx_pair2newmergeidx_offset, &new_prot_ids](PeptideIdentification& pid)
    {
      const set<Size>& runs_to_put = run_id_to_new_run_idcs.at(pid.getIdentifier());
      //TODO check that in the beginning until we support it
      if (runs_to_put.size() > 1)
      {
       OPENMS_LOG_WARN << "Warning: Merging parts of IDRuns currently untested. If it is not a TMT/iTraq sample,"
                    "something is wrong anyway.";
        // in this case you would need to copy the PeptideID
        // should only happen in TMT/itraq
      }

      Size old_merge_idx = 0;
      //TODO we could lookup the old protein ID and see if there were multiple MSruns. If so, we should fail if not
      // exist
      if (pid.metaValueExists(Constants::UserParam::ID_MERGE_INDEX))
      {
        old_merge_idx = pid.getMetaValue(Constants::UserParam::ID_MERGE_INDEX);
      }

      for (const auto& run_to_put : runs_to_put)
      {
        const ProteinIdentification& new_prot_id_run = new_prot_ids[run_to_put];
        pid.setIdentifier(new_prot_id_run.getIdentifier());
        pid.setMetaValue(Constants::UserParam::ID_MERGE_INDEX,
            old_merge_idx + oldrunid_newrunidx_pair2newmergeidx_offset[{pid.getIdentifier(), run_to_put}]);
      }
    };

    cmap.applyFunctionOnPeptideIDs(fun);
    cmap.setProteinIdentifications(std::move(new_prot_ids));
  }

  //merge proteins across fractions and replicates
  void ConsensusMapMergerAlgorithm::mergeAllIDRuns(ConsensusMap& cmap) const
  {
    if (cmap.getProteinIdentifications().size() == 1)
      return;

    // Everything needs to agree
    checkOldRunConsistency_(cmap.getProteinIdentifications(), cmap.getExperimentType());

    ProteinIdentification new_prot_id_run;
    //TODO create better ID
    new_prot_id_run.setIdentifier("merged");
    //TODO merge SearchParams e.g. in case of SILAC
    new_prot_id_run.setSearchEngine(cmap.getProteinIdentifications()[0].getSearchEngine());
    new_prot_id_run.setSearchEngineVersion(cmap.getProteinIdentifications()[0].getSearchEngineVersion());
    new_prot_id_run.setSearchParameters(cmap.getProteinIdentifications()[0].getSearchParameters());
    String old_inference_engine = cmap.getProteinIdentifications()[0].getInferenceEngine();
    if (!old_inference_engine.empty())
    {
      OPENMS_LOG_WARN << "Inference was already performed on the runs in this ConsensusXML."
                         " Merging their proteins, will invalidate correctness of the inference."
                         " You should redo it.\n";
      // deliberately do not take over old inference settings.
    }

    //we do it based on the IDRuns since ID Runs maybe different from quantification in e.g. TMT
    vector<String> merged_origin_files{};
    map<String,pair<Size,bool>> oldrunid2offset_multi_pair;
    for (const auto& pid : cmap.getProteinIdentifications())
    {
      vector<String> out;
      pid.getPrimaryMSRunPath(out);
      Size offset = merged_origin_files.size();
      merged_origin_files.insert(merged_origin_files.end(), out.begin(), out.end());
      oldrunid2offset_multi_pair.emplace(std::piecewise_construct,
                                         std::forward_as_tuple(pid.getIdentifier()),
                                         std::forward_as_tuple(offset, out.size() > 1));
    }
    new_prot_id_run.setPrimaryMSRunPath(merged_origin_files);

    unordered_set<ProteinHit,hash_type,equal_type> proteins_collected_hits(0, accessionHash_, accessionEqual_);

    std::vector<ProteinIdentification>& old_prot_runs = cmap.getProteinIdentifications();
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

    std::map<String, Size> run_id_to_run_idx;
    for (Size old_prot_run_idx = 0; old_prot_run_idx < old_prot_runs.size(); ++old_prot_run_idx)
    {
      ProteinIdentification& protIDRun = old_prot_runs[old_prot_run_idx];
      run_id_to_run_idx[protIDRun.getIdentifier()] = old_prot_run_idx;
    }

    const String& new_prot_id_run_string = new_prot_id_run.getIdentifier();

    function<void(PeptideIdentification &)> fun =
    [&new_prot_id_run_string, &oldrunid2offset_multi_pair](PeptideIdentification& pid) -> void
    {
      const auto& p = oldrunid2offset_multi_pair[pid.getIdentifier()];
      pid.setIdentifier(new_prot_id_run_string);
      Size old = 0;
      if (pid.metaValueExists(Constants::UserParam::ID_MERGE_INDEX))
      {
        old = pid.getMetaValue(Constants::UserParam::ID_MERGE_INDEX);
      }
      else
      {
        if (p.second)
        {
          throw Exception::MissingInformation(
              __FILE__,
              __LINE__,
              OPENMS_PRETTY_FUNCTION,
              "No id_merge_index value in a merged ID run."); //TODO add more info about where.
        }
      }
      pid.setMetaValue(Constants::UserParam::ID_MERGE_INDEX, old + p.first);
    };

    cmap.applyFunctionOnPeptideIDs(fun);

    auto& hits = new_prot_id_run.getHits();
    for (auto& prot : proteins_collected_hits)
    {
      hits.emplace_back(std::move(const_cast<ProteinHit&>(prot))); //careful this completely invalidates the set
    }
    proteins_collected_hits.clear();
    cmap.getProteinIdentifications().resize(1);
    swap(cmap.getProteinIdentifications()[0], new_prot_id_run);
    //TODO remove unreferenced proteins? Can this happen when merging all? I think not.
  }

  bool ConsensusMapMergerAlgorithm::checkOldRunConsistency_(const vector<ProteinIdentification>& protRuns, const String& experiment_type) const
  {
    return checkOldRunConsistency_(protRuns, protRuns[0], experiment_type);
  }

  //TODO refactor the next two functions
  bool ConsensusMapMergerAlgorithm::checkOldRunConsistency_(const vector<ProteinIdentification>& protRuns, const ProteinIdentification& ref, const String& experiment_type) const
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
