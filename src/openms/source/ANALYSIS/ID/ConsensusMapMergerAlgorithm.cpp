// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/IDMergerAlgorithm.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <unordered_map>
#include <include/OpenMS/ANALYSIS/ID/ConsensusMapMergerAlgorithm.h>


using namespace std;
namespace OpenMS
{
  ConsensusMapMergerAlgorithm::ConsensusMapMergerAlgorithm() :
      ConsensusMapMergerAlgorithm::DefaultParamHandler("ConsensusMapMergerAlgorithm")
    {
      defaults_.setValue("annotate_origin",
                         "true",
                         "If true, adds a map_index MetaValue to the PeptideIDs to annotate the IDRun they came from.");
      defaults_.setValidStrings("annotate_origin", ListUtils::create<String>("true,false"));
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
    map<unsigned, unsigned> mapIdx2repBatch{};
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
      pair<String, unsigned> pathLab{consHeader.second.filename, lab};

      unsigned repBatchIdx(0);
      for (auto& repBatch : toMerge)
      {
        for (const std::pair<String, unsigned>& rep : repBatch)
        {
          if (pathLab == rep)
          {
            mapIdx2repBatch[consHeader.first] = repBatchIdx;
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

    mergeProteinIDRuns(cmap, mapIdx2repBatch);
  }

  void ConsensusMapMergerAlgorithm::mergeProteinIDRuns(ConsensusMap &cmap,
                                             map<unsigned, unsigned> const &mapIdx_to_new_protIDRun) const
  {
    // one of label-free, labeled_MS1, labeled_MS2
    const String & experiment_type = cmap.getExperimentType();

    // Not fully supported yet because an ID would need to reference multiple protID runs.
    // we could replicate the ID in the future or allow multiple references.
    bool labelfree = true;
    if (experiment_type != "label-free")
    {
      OPENMS_LOG_WARN << "Merging untested for labelled experiments" << endl;
      labelfree = false;
    }

    // Unfortunately we need a kind of bimap here.
    // For the features we need oldMapIdx -> newIDRunIdx
    // For the new runs we need newIDRunIdx -> <[file_origins], [mapIdcs]> once to initialize them with metadata
    // TODO should we instead better try to collect the primaryMSRuns from the old Runs?
    // TODO I just saw that in the columnHeaders might have the featureXMLs as origins but we should enforce that
    //  this will be changed to the mzML by all tools
    //  Therefore we somehow need to check consistency of ColumnHeaders and ProteinIdentification (file_origins).
    map<unsigned, pair<set<String>,vector<Int>>> newIdcs;
    for (const auto& newIdx : mapIdx_to_new_protIDRun)
    {
      const auto& newIdcsInsertIt = newIdcs.emplace(newIdx.second, make_pair(set<String>(),vector<Int>()));
      newIdcsInsertIt.first->second.first.emplace(cmap.getColumnHeaders().at(newIdx.first).filename);
      newIdcsInsertIt.first->second.second.emplace_back(static_cast<Int>(newIdx.first));
    }
    Size new_size = newIdcs.size();

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
    map<String,set<Size>> runIDToNewRunIdcs;
    // this is to check how many old runs contribute to the new runs
    // this can help save time and we can double check
    vector<Size> nrInputsForNewRunIDs(new_size, 0);
    for (const auto& newidxToOriginsetMapIdxPair : newIdcs)
    {
      for (auto& oldProtID : cmap.getProteinIdentifications())
      {
        set<String> currentContent;
        oldProtID.getPrimaryMSRunPath(currentContent);
        const set<String>& mergeRequest = newidxToOriginsetMapIdxPair.second.first;
        // if this run is fully covered by a requested merged set, use it for it.
        Size count = 1;
        if (std::includes(mergeRequest.begin(), mergeRequest.end(), currentContent.begin(), currentContent.end()))
        {
          auto it = runIDToNewRunIdcs.emplace(oldProtID.getIdentifier(), set<Size>());
          if (!it.second)
          {
            OPENMS_LOG_WARN << "Duplicate protein run ID found. Uniquifying it." << endl;
            oldProtID.setIdentifier(oldProtID.getIdentifier()+"_"+count);
            it = runIDToNewRunIdcs.emplace(oldProtID.getIdentifier(), set<Size>());
          }
          it.first->second.emplace(newidxToOriginsetMapIdxPair.first);
          nrInputsForNewRunIDs.at(newidxToOriginsetMapIdxPair.first)++;
        }
      }
    }

    vector<ProteinIdentification> newProtIDs{new_size};
    //TODO preallocate with sum of proteins?
    unordered_map<ProteinHit,set<Size>,hash_type,equal_type> proteinsCollectedHitsRuns(0,accessionHash,accessionEqual);

    for (auto& runIDToNewRunIdcsPair : runIDToNewRunIdcs)
    {
      // find old run
      auto it = cmap.getProteinIdentifications().begin();
      for (; it != cmap.getProteinIdentifications().end(); ++it)
      {
        if (it->getIdentifier() == runIDToNewRunIdcsPair.first)
          break;
      }

      for (const auto& newrunid : runIDToNewRunIdcsPair.second)
      {
        // go through new runs and fill the proteins and update search settings
        // if first time filling this new run:
        //TODO safe to check for empty identifier?
        if (newProtIDs.at(newrunid).getIdentifier().empty())
        {
          //initialize new run
          newProtIDs[newrunid].setSearchEngine(it->getSearchEngine());
          newProtIDs[newrunid].setSearchEngineVersion(it->getSearchEngineVersion());
          newProtIDs[newrunid].setSearchParameters(it->getSearchParameters());
          StringList toFill;
          it->getPrimaryMSRunPath(toFill);
          newProtIDs[newrunid].setPrimaryMSRunPath(toFill);
          newProtIDs[newrunid].setIdentifier("condition" + String(newrunid));
        }
        // if not, merge settings or check consistency
        else
        {
          //check consistency and add origins
          checkRunSettings_(*it, newProtIDs[newrunid], experiment_type);
          StringList toFill; it->getPrimaryMSRunPath(toFill); // new ones
          newProtIDs[newrunid].addPrimaryMSRunPath(toFill); //add to previous
        }
      }

      //Insert hits into collection with empty set (if not present yet) and
      // add destination run indices
      for (auto& hit : it->getHits())
      {
        const auto& foundIt = proteinsCollectedHitsRuns.emplace(std::move(hit), set<Size>());
        foundIt.first->second.insert(runIDToNewRunIdcsPair.second.begin(), runIDToNewRunIdcsPair.second.end());
      }
      it->getHits().clear(); //not needed anymore and moved anyway
    }

    // copy the protein hits into the destination runs
    for (const auto& protToNewRuns : proteinsCollectedHitsRuns)
    {
      for (Size runID : protToNewRuns.second)
      {
        newProtIDs.at(runID).getHits().emplace_back(protToNewRuns.first);
      }
    }

    //Now update the references in the PeptideHits
    //TODO double check the PrimaryRunPaths with the initial requested merge

    function<void(PeptideIdentification&)> fun = [&runIDToNewRunIdcs, &newProtIDs](PeptideIdentification& pid)
    {
      const set<Size>& runsToPut = runIDToNewRunIdcs.at(pid.getIdentifier());
      //TODO check that in the beginning until we support it
      if (runsToPut.size() > 1)
      {
       OPENMS_LOG_WARN << "Warning: Merging parts of IDRuns currently untested. If it is not a TMT/iTraq sample,"
                    "something is wrong anyway.";
        // in this case you would need to copy the PeptideID
        // should only happen in TMT/itraq
      }

      for (const auto& runToPut : runsToPut)
      {
        const ProteinIdentification& newProtIDRun = newProtIDs[runToPut];
        pid.setIdentifier(newProtIDRun.getIdentifier());
      }
    };

    cmap.applyFunctionOnPeptideIDs(fun);
    cmap.setProteinIdentifications(std::move(newProtIDs));
  }

  //merge proteins across fractions and replicates
  void ConsensusMapMergerAlgorithm::mergeAllIDRuns(ConsensusMap& cmap) const
  {
    // Everything needs to agree
    checkOldRunConsistency_(cmap.getProteinIdentifications(), cmap.getExperimentType());

    ProteinIdentification newProtIDRun;
    //TODO create better ID
    newProtIDRun.setIdentifier("merged");
    //TODO merge SearchParams e.g. in case of SILAC
    newProtIDRun.setSearchEngine(cmap.getProteinIdentifications()[0].getSearchEngine());
    newProtIDRun.setSearchEngineVersion(cmap.getProteinIdentifications()[0].getSearchEngineVersion());
    newProtIDRun.setSearchParameters(cmap.getProteinIdentifications()[0].getSearchParameters());
    //TODO based on old IDRuns or based on consensusHeaders?
    //vector<String> mergedOriginFiles = cmap.getProteinIdentifications()...
    vector<String> mergedOriginFiles{};
    for (const auto& chead : cmap.getColumnHeaders())
    {
      mergedOriginFiles.push_back(chead.second.filename);
    }
    newProtIDRun.setPrimaryMSRunPath(mergedOriginFiles);

    unordered_set<ProteinHit,hash_type,equal_type> proteinsCollectedHits(0,accessionHash,accessionEqual);

    std::vector<ProteinIdentification>& oldProtRuns = cmap.getProteinIdentifications();
    typedef std::vector<ProteinHit>::iterator iter_t;
    for (auto& protRun : oldProtRuns)
    {
      auto& hits = protRun.getHits();
      proteinsCollectedHits.insert(
          std::move_iterator<iter_t>(hits.begin()),
          std::move_iterator<iter_t>(hits.end())
      );
      hits.clear();
    }

    std::map<String, Size> runIDToRunIdx;
    for (Size oldProtRunIdx = 0; oldProtRunIdx < oldProtRuns.size(); ++oldProtRunIdx)
    {
      ProteinIdentification& protIDRun = oldProtRuns[oldProtRunIdx];
      runIDToRunIdx[protIDRun.getIdentifier()] = oldProtRunIdx;
    }

    const String& newProtIDRunString = newProtIDRun.getIdentifier();

    function<void(PeptideIdentification &)> fun =
    [&newProtIDRunString](PeptideIdentification& pid) -> void
    {
      pid.setIdentifier(newProtIDRunString);
    };

    cmap.applyFunctionOnPeptideIDs(fun);

    auto& hits = newProtIDRun.getHits();
    for (auto& prot : proteinsCollectedHits)
    {
      hits.emplace_back(std::move(const_cast<ProteinHit&>(prot))); //careful this completely invalidates the set
    }
    proteinsCollectedHits.clear();
    cmap.getProteinIdentifications().resize(1);
    swap(cmap.getProteinIdentifications()[0],newProtIDRun);
    //TODO remove unreferenced proteins? Can this happen when merging all? I think not.
  }

  //merge proteins across fractions and replicates only copying over referenced proteins
/*
  void ConsensusMapMergerAlgorithm::mergeAllIDRunsOld(ConsensusMap& cmap) const
  {
    checkOldRunConsistency_(cmap.getProteinIdentifications(), cmap.getExperimentType());

    ProteinIdentification newProtIDRun;
    newProtIDRun.setIdentifier("merged");
    //TODO merge SearchParams e.g. in case of SILAC
    newProtIDRun.setSearchEngine(cmap.getProteinIdentifications()[0].getSearchEngine());
    newProtIDRun.setSearchEngineVersion(cmap.getProteinIdentifications()[0].getSearchEngineVersion());
    newProtIDRun.setSearchParameters(cmap.getProteinIdentifications()[0].getSearchParameters());
    //TODO based on old IDRuns or based on consensusHeaders?
    //vector<String> mergedOriginFiles = cmap.getProteinIdentifications()...
    vector<String> mergedOriginFiles{};
    for (const auto& chead : cmap.getColumnHeaders())
    {
      mergedOriginFiles.push_back(chead.second.filename);
    }
    newProtIDRun.setPrimaryMSRunPath(mergedOriginFiles);

    unordered_set<string> proteinsCollected{};
    for (auto& cf : cmap)
    {
      for (auto& pid : cf.getPeptideIdentifications())
      {
        ProteinIdentification* oldProtIDRun(nullptr);
        const String& runID = pid.getIdentifier();
        for (auto& protIDRun : cmap.getProteinIdentifications())
        {
          if (protIDRun.getIdentifier() == runID)
          {
            oldProtIDRun = &protIDRun;
          }
        }
        if (!oldProtIDRun)
        {
          throw Exception::MissingInformation(
              __FILE__,
              __LINE__,
              OPENMS_PRETTY_FUNCTION,
              "Old IdentificationRun not found for PeptideIdentification "
              "(" + String(pid.getMZ()) + ", " + String(pid.getRT()) + ") in feature " + String(cf.getUniqueId()) +  ".");

        }

        pid.setIdentifier(newProtIDRun.getIdentifier());

        for (auto& phit : pid.getHits())
        {
          for (auto& pev : phit.getPeptideEvidences())
          {
            auto acc = proteinsCollected.find(pev.getProteinAccession());
            if (acc == proteinsCollected.end())
            {
              newProtIDRun.getHits().emplace_back(*oldProtIDRun->findHit(pev.getProteinAccession()));
              proteinsCollected.insert(pev.getProteinAccession());
            }
          }
        }
      }
    }

    for (auto& upid : cmap.getUnassignedPeptideIdentifications())
    {
      ProteinIdentification* oldProtIDRun(nullptr);
      const String& runID = upid.getIdentifier();
      for (auto& protIDRun : cmap.getProteinIdentifications())
      {
        if (protIDRun.getIdentifier() == runID)
        {
          oldProtIDRun = &protIDRun;
        }
      }
      if (!oldProtIDRun)
      {
        throw Exception::MissingInformation(
            __FILE__,
            __LINE__,
            OPENMS_PRETTY_FUNCTION,
            "Old IdentificationRun not found for UnassignedPeptideIdentification "
            "(" + String(upid.getMZ()) + ", " + String(upid.getRT()) + ").");

      }

      upid.setIdentifier(newProtIDRun.getIdentifier());

      for (auto& phit : upid.getHits())
      {
        for (auto& pev : phit.getPeptideEvidences())
        {
          auto acc = proteinsCollected.find(pev.getProteinAccession());
          if (acc == proteinsCollected.end())
          {
            newProtIDRun.getHits().emplace_back(*oldProtIDRun->findHit(pev.getProteinAccession()));
            proteinsCollected.insert(pev.getProteinAccession());
          }
        }
      }
    }
    cmap.setProteinIdentifications(vector<ProteinIdentification>{newProtIDRun});
  }
*/

  //merge proteins across fractions and replicates. Sort peptides into the same amount of vec<PepIDs>
  //I think this should only ever work for labelfree since we do not store the labels in the idXML structures
  //TODO we could need a method to "not split" the peptides (or you just append the results)
  //TODO we could need a method to merge by origin filenames. Then we could just construct a mapping here,
  // and assume that all the origin filename in one run map to the same new run and in the end call this method.
  //TODO we could need a method to split and merge based on the origin filenames. But this requires major
  // rewriting and should be a new method.

  //TODO docu that oldrunToNewRun should be continuous in the mapped-to indices
  void ConsensusMapMergerAlgorithm::mergeIDRunsAndSplitPeptides(
      vector<ProteinIdentification>& oldProtRuns,
      vector<PeptideIdentification>& pepIDs,
      const map<Size, Size>& oldrunToNewrun,
      vector<vector<PeptideIdentification>>& splitPepIDs) const
  {
    //Without any exp. design we assume label-free for checking mods
    checkOldRunConsistency_(oldProtRuns, "label-free");

    vector<map<Size, Size>> oldToNewFileIdx;
    vector<ProteinIdentification> newProtIDRuns;
    initNewRunsAndFileMappings_(oldProtRuns, oldrunToNewrun, oldToNewFileIdx, newProtIDRuns);

    // for each new run collect the set of proteins to quickly check for existence
    vector<unordered_set<string>> proteinsCollected{newProtIDRuns.size()};

    movePepIDsAndRefProteinsToResult_(
      pepIDs,
      oldProtRuns,
      newProtIDRuns,
      splitPepIDs,
      oldrunToNewrun,
      oldToNewFileIdx,
      proteinsCollected
    );

    oldProtRuns.swap(newProtIDRuns);
    newProtIDRuns.clear();
  }

  //same as above only for vecs of vecs of pepIDs (e.g. from multiple files)
  void ConsensusMapMergerAlgorithm::mergeIDRunsAndSplitPeptides(
      vector<ProteinIdentification>& oldProtRuns,
      vector<vector<PeptideIdentification>>& pepIDs,
      const map<Size, Size>& oldrunToNewrun,
      vector<vector<PeptideIdentification>>& splitPepIDs) const
  {
    //Without any exp. design we assume label-free for checking mods
    checkOldRunConsistency_(oldProtRuns, "label-free");

    vector<map<Size, Size>> oldToNewFileIdx;
    vector<ProteinIdentification> newProtIDRuns;
    initNewRunsAndFileMappings_(oldProtRuns, oldrunToNewrun, oldToNewFileIdx, newProtIDRuns);

    // for each new run collect the set of proteins to quickly check for existence
    vector<unordered_set<string>> proteinsCollected{newProtIDRuns.size()};

    for (auto& pepID : pepIDs)
    {
      movePepIDsAndRefProteinsToResult_(
          pepID,
          oldProtRuns,
          newProtIDRuns,
          splitPepIDs,
          oldrunToNewrun,
          oldToNewFileIdx,
          proteinsCollected
      );
    }

    oldProtRuns.swap(newProtIDRuns);
    newProtIDRuns.clear();
  }

  void ConsensusMapMergerAlgorithm::movePepIDsAndRefProteinsToResult_(
      vector<PeptideIdentification>& pepIDs,
      vector<ProteinIdentification>& oldProtRuns,
      vector<ProteinIdentification>& newProtIDRuns,
      vector<vector<PeptideIdentification>>& splitPepIDs,
      const map<Size,Size>& oldrunToNewrun,
      const vector<map<Size,Size>>& oldToNewFileIdx,
      vector<unordered_set<string>> proteinsCollected
      ) const
  {
    bool annotate_origin(param_.getValue("annotate_origin").toBool());

    for (auto &pid : pepIDs)
    {
      const String &runID = pid.getIdentifier();

      //TODO maybe create lookup table in the beginning
      Size oldProtRunIdx = 0;
      for (; oldProtRunIdx < oldProtRuns.size(); ++oldProtRunIdx)
      {
        ProteinIdentification &protIDRun{oldProtRuns[oldProtRunIdx]};
        if (protIDRun.getIdentifier() == runID)
        {
          break;
        }
      }
      if (oldProtRunIdx == oldProtRuns.size())
      {
        throw Exception::MissingInformation(
            __FILE__,
            __LINE__,
            OPENMS_PRETTY_FUNCTION,
            "Old IdentificationRun not found for PeptideIdentification "
            "(" + String(pid.getMZ()) + ", " + String(pid.getRT()) + ").");

      }

      Size newRunIdx = oldrunToNewrun.at(oldProtRunIdx);
      bool annotated = pid.metaValueExists("map_index");
      if (annotate_origin || annotated)
      {
        Size oldFileIdx(0);
        if (annotated)
        {
          oldFileIdx = pid.getMetaValue("map_index");
        }
          // If there is more than one possible file it might be from
          // and it is not annotated -> fail
        else if (oldToNewFileIdx[oldProtRunIdx].size() > 1)
        {
          throw Exception::MissingInformation(
              __FILE__,
              __LINE__,
              OPENMS_PRETTY_FUNCTION,
              "Trying to annotate new map_index for PeptideIdentification "
              "(" + String(pid.getMZ()) + ", " + String(pid.getRT()) + ") but"
                                                                       "no old map_index present");
        }
        pid.setMetaValue("map_index", oldToNewFileIdx[oldProtRunIdx].at(oldFileIdx));
      }
      pid.setIdentifier(newProtIDRuns[newRunIdx].getIdentifier());
      for (auto &phit : pid.getHits())
      {
        //TODO think about getting the set first and then look for each acc
        for (auto &acc : phit.extractProteinAccessionsSet())
        {
          const auto &it = proteinsCollected[newRunIdx].emplace(acc);
          if (it.second) // was newly inserted
          {
            newProtIDRuns[newRunIdx].getHits().emplace_back(std::move(*oldProtRuns[oldProtRunIdx].findHit(acc)));
          }
        }
      }
      //move peptides into right vector
      splitPepIDs[newRunIdx].emplace_back(std::move(pid));
    }
  }

  void ConsensusMapMergerAlgorithm::initNewRunsAndFileMappings_(
      const vector<ProteinIdentification>& oldProtRuns,
      const map<Size,Size>& oldrunToNewrun,
      vector<map<Size, Size>>& oldToNewFileIdx,
      vector<ProteinIdentification> & newProtIDRuns) const
  {
    //TODO do some precondition checks for the map etc.
    Size nrNewRuns = 1;
    for (const auto& runmap : oldrunToNewrun)
    {
      if (runmap.second > nrNewRuns)
      {
        nrNewRuns = runmap.second;
      }
    }
    nrNewRuns++;


    // For each of the files in the old runs, we construct a mapping to the future index
    // in the new runs
    //TODO actually could just be a vector instead of map
    oldToNewFileIdx.clear();
    oldToNewFileIdx.resize(oldProtRuns.size());
    // holds the proteins collected for each of the new runs
    newProtIDRuns.clear();
    newProtIDRuns.resize(nrNewRuns);

    // For each new set of files of the new runs, we build a map with its future index
    // If a file is not yet present, it gets inserted with a new index
    vector<map<String, Size>> newFileToIdx{nrNewRuns};
    // collect origin files for all new runs
    for (auto &runmap : oldrunToNewrun)
    {
      StringList toFill;
      oldProtRuns[runmap.first].getPrimaryMSRunPath(toFill);
      Size c = 0;
      for (auto &f : toFill)
      {
        auto it_inserted = newFileToIdx[runmap.second].emplace(std::move(f), newFileToIdx[runmap.second].size());
        oldToNewFileIdx[runmap.first].emplace(c, it_inserted.first->second);
        c++;
      }
      toFill.clear();
    }

    vector<StringList> newRunOriginFiles{nrNewRuns};
    //transform the map to vectors that will be inserted into the new protein runs
    // at the index that is stored in the maps that we created
    for (Size s = 0; s < nrNewRuns; ++s)
    {
      StringList &newFiles = newRunOriginFiles[s];
      newFiles.resize(newFileToIdx[s].size());
      for (auto &newRunFileIdxMapEntry : newFileToIdx.at(s))
      {
        newFiles.at(newRunFileIdxMapEntry.second) = std::move(newRunFileIdxMapEntry.first);
      }
    }
    newFileToIdx.clear();


    Size idcount = 0;
    for (auto &newProtIDRun : newProtIDRuns)
    {
      //TODO better ID?
      newProtIDRun.setIdentifier("merged" + String(idcount));
      //TODO merge SearchParams e.g. in case of SILAC mods
      newProtIDRun.setSearchEngine(oldProtRuns[0].getSearchEngine());
      newProtIDRun.setSearchEngineVersion(oldProtRuns[0].getSearchEngineVersion());
      newProtIDRun.setSearchParameters(oldProtRuns[0].getSearchParameters());
      newProtIDRun.setPrimaryMSRunPath(std::move(newRunOriginFiles[idcount]));
      idcount++;
    }
    newRunOriginFiles.clear();
  }

  bool ConsensusMapMergerAlgorithm::checkOldRunConsistency_(const vector<ProteinIdentification>& protRuns, const String& experiment_type) const
  {
    return checkOldRunConsistency_(protRuns, protRuns[0], experiment_type);
  }

  //TODO refactor the next two functions
  bool ConsensusMapMergerAlgorithm::checkOldRunConsistency_(const vector<ProteinIdentification>& protRuns, const ProteinIdentification& ref, const String& experiment_type) const
  {
    const String& engine = ref.getSearchEngine();
    const String& version = ref.getSearchEngineVersion();
    ProteinIdentification::SearchParameters params = ref.getSearchParameters();
    set<String> fixed_mods(params.fixed_modifications.begin(), params.fixed_modifications.end());
    set<String> var_mods(params.variable_modifications.begin(), params.variable_modifications.end());
    bool ok = false;
    unsigned runID = 0;
    for (const auto& idRun : protRuns)
    {
      ok = true;
      if (idRun.getSearchEngine() != engine || idRun.getSearchEngineVersion() != version)
      {
        ok = false;
       OPENMS_LOG_WARN << "Search engine " + idRun.getSearchEngine() + "from IDRun " + String(runID) + " does not match "
                                                                                                 "with the others. You probably do not want to merge the results with this tool.";
        break;
      }
      const ProteinIdentification::SearchParameters& sp = idRun.getSearchParameters();
      if (params.precursor_mass_tolerance != sp.precursor_mass_tolerance ||
          params.precursor_mass_tolerance_ppm != sp.precursor_mass_tolerance_ppm ||
          params.db != sp.db ||
          params.db_version != sp.db_version ||
          params.fragment_mass_tolerance != sp.fragment_mass_tolerance ||
          params.fragment_mass_tolerance_ppm != sp.fragment_mass_tolerance_ppm ||
          params.charges != sp.charges ||
          params.digestion_enzyme != sp.digestion_enzyme ||
          params.taxonomy != sp.taxonomy)
      {
        ok = false;
       OPENMS_LOG_WARN << "Searchengine settings from IDRun " + String(runID) + " does not match with the others."
                                                                          " You probably do not want to merge the results with this tool if they differ significantly.";
        break;
      }

      set<String> curr_fixed_mods(sp.fixed_modifications.begin(), sp.fixed_modifications.end());
      set<String> curr_var_mods(sp.variable_modifications.begin(), sp.variable_modifications.end());
      if (fixed_mods != curr_fixed_mods ||
          var_mods != curr_var_mods)
      {
        if (experiment_type != "labeled_MS1")
        {
          ok = false;
         OPENMS_LOG_WARN << "Used modification settings from IDRun " + String(runID) + " does not match with the others."
                                                                                 " Since the experiment is not annotated as MS1-labeled you probably do not want to merge the results with this tool.";
          break;
        }
        else
        {
          //TODO actually introduce a flag for labelling modifications in the Mod datastructures?
          //OR put a unique ID for the used mod as a UserParam to the mapList entries (consensusHeaders)
          //TODO actually you would probably need an experimental design here, because
          //settings have to agree exactly in a FractionGroup but can slightly differ across runs.
         OPENMS_LOG_WARN << "Used modification settings from IDRun " + String(runID) + " does not match with the others."
                                                                                 " Although it seems to be an MS1-labeled experiment, check carefully that only non-labelling mods differ.";
        }
      }
    }
    if (!ok /*&& TODO and no force flag*/)
    {
      throw Exception::BaseException(__FILE__,
                                     __LINE__,
                                     OPENMS_PRETTY_FUNCTION,
                                     "InvalidData",
                                     "Search settings are not matching across IdentificationRuns. "
                                     "See warnings. Aborting..");
    }
    return ok;
  }

  bool ConsensusMapMergerAlgorithm::checkRunSettings_(const ProteinIdentification& idRun, const ProteinIdentification& ref, const String& experiment_type) const
  {
    const String& warn = " You probably do not want to merge the results with this tool."
                         " For merging searches with different engines/settings please use ConsensusID or PercolatorAdapter"
                         " to create a comparable score.";
    const String& engine = ref.getSearchEngine();
    const String& version = ref.getSearchEngineVersion();
    ProteinIdentification::SearchParameters params = ref.getSearchParameters();
    set<String> fixed_mods(params.fixed_modifications.begin(), params.fixed_modifications.end());
    set<String> var_mods(params.variable_modifications.begin(), params.variable_modifications.end());
    bool ok = true;

    if (idRun.getSearchEngine() != engine || idRun.getSearchEngineVersion() != version)
    {
      ok = false;
     OPENMS_LOG_WARN << "Search engine " + idRun.getSearchEngine() + "from IDRun " + idRun.getIdentifier()
      + " does not match with the others." + warn;
    }
    const ProteinIdentification::SearchParameters& sp = idRun.getSearchParameters();
    if (params.precursor_mass_tolerance != sp.precursor_mass_tolerance ||
        params.precursor_mass_tolerance_ppm != sp.precursor_mass_tolerance_ppm ||
        params.db != sp.db ||
        params.db_version != sp.db_version ||
        params.fragment_mass_tolerance != sp.fragment_mass_tolerance ||
        params.fragment_mass_tolerance_ppm != sp.fragment_mass_tolerance_ppm ||
        params.charges != sp.charges ||
        params.digestion_enzyme != sp.digestion_enzyme ||
        params.taxonomy != sp.taxonomy)
    {
      ok = false;
     OPENMS_LOG_WARN << "Searchengine settings from IDRun " + idRun.getIdentifier() + " does not match with the others." + warn;
    }

    set<String> curr_fixed_mods(sp.fixed_modifications.begin(), sp.fixed_modifications.end());
    set<String> curr_var_mods(sp.variable_modifications.begin(), sp.variable_modifications.end());
    if (fixed_mods != curr_fixed_mods ||
        var_mods != curr_var_mods)
    {
      if (experiment_type != "labeled_MS1")
      {
        ok = false;
       OPENMS_LOG_WARN << "Used modification settings from IDRun " + idRun.getIdentifier() + " does not match with the others." + warn;
      }
      else
      {
        //TODO actually introduce a flag for labelling modifications in the Mod datastructures?
        //OR put a unique ID for the used mod as a UserParam to the mapList entries (consensusHeaders)
        //TODO actually you would probably need an experimental design here, because
        //settings have to agree exactly in a FractionGroup but can slightly differ across runs.
        //Or just ignore labelling mods during the check
       OPENMS_LOG_WARN << "Used modification settings from IDRun " + idRun.getIdentifier() + " does not match with the others."
        " Although it seems to be an MS1-labeled experiment, check carefully that only non-labelling mods differ.";
      }
    }

    if (!ok /*&& TODO and no force flag*/)
    {
      throw Exception::BaseException(__FILE__,
                                     __LINE__,
                                     OPENMS_PRETTY_FUNCTION,
                                     "InvalidData",
                                     "Search settings are not matching across IdentificationRuns. "
                                     "See warnings. Aborting..");
    }
    // TODO else merge as far as possible (mainly mods I guess)
    return ok;
  }
} // namespace OpenMS
