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
#include <unordered_set>

using namespace std;
namespace OpenMS
{
  IDMergerAlgorithm::IDMergerAlgorithm() :
    IDMergerAlgorithm::DefaultParamHandler("IDMergerAlgorithm")
    {
      defaults_.setValue("annotate_origin",
                         "true",
                         "If true, adds a map_index MetaValue to the PeptideIDs to annotate the IDRun they came from.");
      defaults_.setValidStrings("annotate_origin", ListUtils::create<String>("true,false"));
      defaultsToParam_();
    }

  //merge proteins across fractions and replicates
  void IDMergerAlgorithm::mergeProteinsAcrossFractionsAndReplicates(ConsensusMap& cmap, const ExperimentalDesign& exp_design) const
  {
    const vector<vector<pair<String, unsigned>>> toMerge = exp_design.getSampleWOReplicatesToMSFilesMapping();

    // one of label-free, labeled_MS1, labeled_MS2
    const String & experiment_type = cmap.getExperimentType();

    //Not supported because an ID would need to reference multiple protID runs.
    //we could replicate the ID in the future or allow multiple references.
    bool labelfree = true;
    if (experiment_type != "label-free")
    {
      LOG_WARN << "Merging untested for labelled experiments" << endl;
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
        if (labelfree)
        {
          LOG_WARN << "No channel id annotated in consensusXML. Assuming one channel." << endl;
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

  void IDMergerAlgorithm::mergeProteinIDRuns(ConsensusMap &cmap,
                                             map<unsigned, unsigned> const &mapIdx_to_new_protIDRun) const
  {
    // one of label-free, labeled_MS1, labeled_MS2
    const String & experiment_type = cmap.getExperimentType();

    checkOldRunConsistency_(cmap.getProteinIdentifications(), experiment_type);


    //Not supported because an ID would need to reference multiple protID runs.
    //we could replicate the ID in the future or allow multiple references.
    bool labelfree = true;
    if (experiment_type != "label-free")
    {
      LOG_WARN << "Merging untested for labelled experiments" << endl;
      labelfree = false;
    }

    // Unfortunately we need a kind of bimap here.
    // For the features we need map -> IDRun
    // For the new runs we need IDRun -> <maps> once to initialize them with metadata
    // TODO should we instead better try to collect the primaryMSRuns from the old Runs?
    // Or are the filenames in the consensusHeaders always equivalent?
    map<unsigned, pair<vector<String>,vector<Int>>> newIdcs;
    for (const auto& newIdx : mapIdx_to_new_protIDRun)
    {
      const auto& newIdcsIt = newIdcs.find(newIdx.second);
      if (newIdcsIt == newIdcs.end())
      {
        newIdcs[newIdx.second] = make_pair(vector<String>{cmap.getColumnHeaders()[newIdx.first].filename},vector<Int>{static_cast<Int>(newIdx.first)});
      }
      else
      {
        (*newIdcsIt).second.first.emplace_back(cmap.getColumnHeaders()[newIdx.first].filename);
        (*newIdcsIt).second.second.emplace_back(static_cast<Int>(newIdx.first));
      }
    }
    Size new_size = newIdcs.size();

    if (new_size == 1)
    {
      LOG_WARN << "Number of new protein ID runs is one. Consider using mergeAllProteinRuns for some additional speed." << endl;
    }
    else if (new_size >= cmap.getColumnHeaders().size())
    {
      throw Exception::InvalidValue(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          "Number of new protein runs after merging"
          " is smaller or equal to the original ones."
          " Aborting. Nothing would be merged.", String(new_size));
    }

    vector<ProteinIdentification> newProtIDs{new_size};
    unsigned j = 0;
    for (auto& pid : newProtIDs)
    {
      pid.setIdentifier("condition" + String(j));
      //TODO merge SearchParams e.g. in case of SILAC
      pid.setSearchEngine(cmap.getProteinIdentifications()[0].getSearchEngine());
      pid.setSearchEngineVersion(cmap.getProteinIdentifications()[0].getSearchEngineVersion());
      pid.setSearchParameters(cmap.getProteinIdentifications()[0].getSearchParameters());
      pid.setPrimaryMSRunPath(newIdcs.at(j).first);
      //TODO is setting a map_index or a label necessary?
      pid.setMetaValue("map_indices", newIdcs.at(j).second);
      j++;
    }

    //stores the index in the newProtIDs ProteinHits + 1. So an entry of 0 means not yet in there.
    unordered_map<string, vector<Size>> acc2ProtHitIdxPerRun{};
    for (auto& cf : cmap)
    {
      for (auto& pid : cf.getPeptideIdentifications())
      {
        int mapIdx(-1);
        //set<int> mapIdcs; //would be needed for itraq
        if (pid.metaValueExists("map_index"))
        {
          //mapIdcs.insert(static_cast<int>(pid.getMetaValue("map_index")));
          mapIdx = static_cast<int>(pid.getMetaValue("map_index"));
        }
        else
        {
          if (labelfree)
          {
            throw Exception::MissingInformation(
                __FILE__,
                __LINE__,
                OPENMS_PRETTY_FUNCTION,
                "UserParam 'map_index' not found for PeptideIdentification"
                "(" + String(pid.getMZ()) + ", " + String(pid.getRT()) + ") in feature " + String(cf.getUniqueId()) +  "."
                " Cannot associate it to the experimental design");
          }
          else //TODO how to check between ITraq and SILAC? Does it make sense?
          {
            //if SILAC, add to an extra prot id run or ignore
            // should not really happen if searched with fixed mod.

            //if itraq/tmt, infer the map_indices by looking at non-zero quants
            // currently not supported
          }
        }

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

        //for (const auto& mapIdx : mapIdcs)
        //{
        unsigned runToPut;
        const auto& foundIt = mapIdx_to_new_protIDRun.find(mapIdx);
        if (foundIt == mapIdx_to_new_protIDRun.end())
        {
          throw Exception::BaseException(
              __FILE__,
              __LINE__,
              OPENMS_PRETTY_FUNCTION,
              "ElementNotFound",
              "Mapping for map_index " + String(mapIdx) + " not found. Check experimental design.");
        }
        else
        {
          runToPut = foundIt->second;
        }
        ProteinIdentification& newProtIDRun = newProtIDs[runToPut];
        pid.setIdentifier(newProtIDRun.getIdentifier());

        for (auto& phit : pid.getHits())
        {
          for (auto& pev : phit.getPeptideEvidences())
          {
            auto acc2ProtHitIdxPerRunIt = acc2ProtHitIdxPerRun.find(pev.getProteinAccession());
            if (acc2ProtHitIdxPerRunIt != acc2ProtHitIdxPerRun.end())
            {
              if (acc2ProtHitIdxPerRunIt->second[runToPut] != 0) //already there
              {
                //TODO update target decoy info? Can this happen? Should not!
                //Remember to access the proteinID with
                //newProtIDs[runToPut][acc2ProtHitIdxPerRunIt->second[runToPut] - 1]
              }
              else
              {
                newProtIDRun.getHits().emplace_back(*oldProtIDRun->findHit(pev.getProteinAccession()));
                acc2ProtHitIdxPerRunIt->second[runToPut] = newProtIDRun.getHits().size();
              }
            }
            else
            {
              newProtIDRun.getHits().emplace_back(*oldProtIDRun->findHit(pev.getProteinAccession()));
              vector<Size> initVec{new_size, 0u};
              initVec[runToPut] = newProtIDRun.getHits().size();
              acc2ProtHitIdxPerRun[pev.getProteinAccession()] = move(initVec);
            }
          }
        }
        //}
      }
    }

    for (auto& upid : cmap.getUnassignedPeptideIdentifications())
    {
      int mapIdx(-1);
      //set<int> mapIdcs; //would be needed for itraq
      if (upid.metaValueExists("map_index"))
      {
        //mapIdcs.insert(static_cast<int>(pid.getMetaValue("map_index")));
        mapIdx = static_cast<int>(upid.getMetaValue("map_index"));
      }
      else
      {
        if (labelfree)
        {
          throw Exception::MissingInformation(
              __FILE__,
              __LINE__,
              OPENMS_PRETTY_FUNCTION,
              "UserParam 'map_index' not found for UnassignedPeptideIdentification "
              "(" + String(upid.getMZ()) + ", " + String(upid.getRT()) + ")."
              " Cannot associate it to the experimental design");
        }
        else //TODO how to check between ITraq and SILAC? Does it make sense?
        {
          //if SILAC, add to an extra prot id run or ignore unlabeled
          // should not really happen if searched with fixed mod.

          //if itraq/tmt, infer the map_indices by looking at non-zero quants
          // currently not supported as we probably would have to replicate IDs
          // for every map
        }
      }

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
            "Old IdentificationRun for UnassignedPeptideIdentification"
            "(" + String(upid.getMZ()) + ", " + String(upid.getRT()) + ") not found.");
      }

      //for (const auto& mapIdx : mapIdcs)
      //{
      unsigned runToPut;
      const auto& foundIt = mapIdx_to_new_protIDRun.find(mapIdx);
      if (foundIt == mapIdx_to_new_protIDRun.end())
      {
        throw Exception::BaseException(
            __FILE__,
            __LINE__,
            OPENMS_PRETTY_FUNCTION,
            "ElementNotFound",
            "Mapping for map_index " + String(mapIdx) + " not found. Check experimental design.");
      }
      else
      {
        runToPut = foundIt->second;
      }
      ProteinIdentification& newProtIDRun = newProtIDs[runToPut];
      upid.setIdentifier(newProtIDRun.getIdentifier());

      for (auto& phit : upid.getHits())
      {
        for (auto& pev : phit.getPeptideEvidences())
        {
          auto acc2ProtHitIdxPerRunIt = acc2ProtHitIdxPerRun.find(pev.getProteinAccession());
          if (acc2ProtHitIdxPerRunIt != acc2ProtHitIdxPerRun.end())
          {
            if (acc2ProtHitIdxPerRunIt->second[runToPut] != 0) //already there
            {
              //TODO update target decoy info? Can this happen? Should not!
              //Remember to access the proteinID with
              //newProtIDs[repBatchToPut][acc2ProtHitIdxPerRunIt->second[repBatchToPut] - 1]
            }
            else
            {
              newProtIDRun.getHits().emplace_back(*oldProtIDRun->findHit(pev.getProteinAccession()));
              acc2ProtHitIdxPerRunIt->second[runToPut] = newProtIDRun.getHits().size();
            }
          }
          else
          {
            newProtIDRun.getHits().emplace_back(*oldProtIDRun->findHit(pev.getProteinAccession()));
            vector<Size> initVec{new_size, 0u};
            initVec[runToPut] = newProtIDRun.getHits().size();
            acc2ProtHitIdxPerRun[pev.getProteinAccession()] = move(initVec);
          }
        }
      }
      //}
    }
    cmap.setProteinIdentifications(newProtIDs);
  }

  //merge proteins across fractions and replicates
  void IDMergerAlgorithm::mergeAllIDRuns(ConsensusMap& cmap) const
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

  //merge proteins across fractions and replicates
  void IDMergerAlgorithm::mergeAllIDRuns(vector<ProteinIdentification>& protRuns, vector<PeptideIdentification>& pepIDs) const
  {
    //Without any exp. design we assume label-free for checking mods
    checkOldRunConsistency_(protRuns, "label-free");

    bool annotate_origin(param_.getValue("annotate_origin").toBool());

    ProteinIdentification newProtIDRun;
    //TODO better ID?
    newProtIDRun.setIdentifier("merged");
    //TODO merge SearchParams e.g. in case of SILAC mods
    newProtIDRun.setSearchEngine(protRuns[0].getSearchEngine());
    newProtIDRun.setSearchEngineVersion(protRuns[0].getSearchEngineVersion());
    newProtIDRun.setSearchParameters(protRuns[0].getSearchParameters());

    //TODO construct mapping from RunPath to index in set
    set<String> mergedOriginFiles{};
    for (const auto& protIDRun: protRuns)
    {
      StringList toFill;
      protIDRun.getPrimaryMSRunPath(toFill);
      //TODO I think we can remove that.  We have sensible default now.
      if (toFill.size() > 1 && annotate_origin)
      {
        throw Exception::IllegalArgument(
            __FILE__,
            __LINE__,
            OPENMS_PRETTY_FUNCTION,
            "Annotating origin not yet supported when merging already merged ID runs.");
      }
      for (String& s : toFill)
      {
        mergedOriginFiles.insert(s);
      }
      toFill.clear();
    }
    vector<String> mergedOriginFilesVec(mergedOriginFiles.begin(), mergedOriginFiles.end());
    newProtIDRun.setPrimaryMSRunPath(mergedOriginFilesVec);

    unordered_set<string> proteinsCollected{};
    for (auto& pid : pepIDs)
    {
        ProteinIdentification* oldProtIDRun(nullptr);
        int oldProtRunIdx = -1;
        const String& runID = pid.getIdentifier();
        for (auto& protIDRun : protRuns)
        {
          if (protIDRun.getIdentifier() == runID)
          {
            oldProtIDRun = &protIDRun;

            if (annotate_origin)
            {
              StringList runs;
              oldProtIDRun->getPrimaryMSRunPath(runs);
              String run = runs[0];

              if (runs.size() > 1 && pid.metaValueExists("map_index"))
              {
                run = runs[pid.getMetaValue("map_index")];
              }

              //TODO maybe create lookup table in the beginning
              for (Size u = 0; u < mergedOriginFilesVec.size(); ++u)
              {
                if (mergedOriginFilesVec[u] == run)
                {
                  oldProtRunIdx = u;
                  break;
                }
              }
            }
          }
        }
        if (!oldProtIDRun)
        {
          throw Exception::MissingInformation(
              __FILE__,
              __LINE__,
              OPENMS_PRETTY_FUNCTION,
              "Old IdentificationRun not found for PeptideIdentification "
              "(" + String(pid.getMZ()) + ", " + String(pid.getRT()) + ").");

        }

        pid.setIdentifier(newProtIDRun.getIdentifier());

        //TODO think about a better way to annotate the originating run.
        // this is borrowed from consensusXML where map_index is used.
        if (annotate_origin) pid.setMetaValue("map_index", oldProtRunIdx);

        for (auto& phit : pid.getHits())
        {
          //TODO think about getting the set first and then look for each acc
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

    protRuns = vector<ProteinIdentification>{newProtIDRun};
  }

  bool IDMergerAlgorithm::checkOldRunConsistency_(const vector<ProteinIdentification> protRuns, String experiment_type) const
  {
    String engine = protRuns[0].getSearchEngine();
    String version = protRuns[0].getSearchEngineVersion();
    ProteinIdentification::SearchParameters params = protRuns[0].getSearchParameters();
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
        LOG_WARN << "Search engine " + idRun.getSearchEngine() + "from IDRun " + String(runID) + " does not match "
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
        LOG_WARN << "Searchengine settings from IDRun " + String(runID) + " does not match with the others."
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
          LOG_WARN << "Used modification settings from IDRun " + String(runID) + " does not match with the others."
          " Since the experiment is not annotated as MS1-labeled you probably do not want to merge the results with this tool.";
          break;
        }
        else
        {
          //TODO actually introduce a flag for labelling modifications in the Mod datastructures?
          //OR put a unique ID for the used mod as a UserParam to the mapList entries (consensusHeaders)
          //TODO actually you would probably need an experimental design here, because
          //settings have to agree exactly in a FractionGroup but can slightly differ across runs.
          LOG_WARN << "Used modification settings from IDRun " + String(runID) + " does not match with the others."
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
} // namespace OpenMS
