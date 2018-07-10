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

using namespace std;
namespace OpenMS
{
  IDMergerAlgorithm::IDMergerAlgorithm() :
    IDMergerAlgorithm::DefaultParamHandler("ConsensusIDAlgorithm")
    {

    }
  //merge proteins across fractions and replicates
  void IDMergerAlgorithm::mergeProteinsAcrossFractionsAndReplicates(ConsensusMap& cmap, const ExperimentalDesign& exp_design)
  {
    const auto toMerge = exp_design.getSampleWOReplicatesToMSFilesMapping();

    //Not supported because an ID would need to reference multiple protID runs.
    //we could replicate the ID in the future or allow multiple references.
    if (cmap.getExperimentType() == "labelled") return; //TODO throw error.

    //out of the path/label combos, construct sets of map indices to be merged
    unsigned int repBatchIdx(0);
    map<unsigned, unsigned> mapIdx2repBatch{};
    for (auto& repBatch : toMerge)
    {
      std::cout << "Merging: (";
      for (const std::pair<String, unsigned>& rep : repBatch)
      {
        std::cout << rep.first << ", ";
        for (auto& consHeader : cmap.getColumnHeaders())
        {
          pair<String, unsigned> pathLab{consHeader.second.filename, consHeader.second.unique_id};
          if (pathLab == rep)
          {
            mapIdx2repBatch[consHeader.first] = repBatchIdx;
          }
        }
      }
      std::cout << ")" << std::endl;
      repBatchIdx++;
    }

    vector<ProteinIdentification> newProtIDs{toMerge.size()};
    int j = 0;
    for (auto& pid : newProtIDs)
    {
      pid.setIdentifier("condition" + String(j));
    }
    //stores the index in the newProtIDs ProteinHits + 1. So an entry of 0 means not yet in there.
    unordered_map<string, vector<Size>> acc2ProtHitIdxPerBatch{};
    for (auto& cf : cmap)
    {
      for (auto& pid : cf.getPeptideIdentifications())
      {
        int mapIdx(0);
        //set<int> mapIdcs; //would be needed for itraq
        if (pid.metaValueExists("map_index"))
        {
          //mapIdcs.insert(static_cast<int>(pid.getMetaValue("map_index")));
          mapIdx = static_cast<int>(pid.getMetaValue("map_index"));
        }
        else
        {
          if (cmap.getExperimentType() == "label-free")
          {
            //error, should not happen
          }
          else if (cmap.getExperimentType() == "labelled") //TODO how to check between ITraq and SILAC?
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
          std::cout << "Skipped PID because old ProtID Run ID could not be found."
          continue; //TODO Throw warning or exception that the old ProtIDrun could not be found.
        }

        //for (const auto& mapIdx : mapIdcs)
        //{
        unsigned repBatchToPut = mapIdx2repBatch[mapIdx];
        ProteinIdentification& newProtIDRun = newProtIDs[repBatchToPut];
        pid.setIdentifier(newProtIDRun.getIdentifier());

        for (auto& phit : pid.getHits())
        {
          for (auto& pev : phit.getPeptideEvidences())
          {
            auto acc2ProtHitIdxPerBatchIt = acc2ProtHitIdxPerBatch.find(pev.getProteinAccession());
            if (acc2ProtHitIdxPerBatchIt != acc2ProtHitIdxPerBatch.end())
            {
              if (acc2ProtHitIdxPerBatchIt->second[repBatchToPut] != 0) //already there
              {
                //TODO update target decoy info? Can this happen? Should not!
                //Remember to access the proteinID with
                //newProtIDs[repBatchToPut][acc2ProtHitIdxPerBatchIt->second[repBatchToPut] - 1]
              }
              else
              {
                newProtIDRun.getHits().emplace_back(*oldProtIDRun->findHit(pev.getProteinAccession()));
              }
            }
            else
            {
              newProtIDRun.getHits().emplace_back(*oldProtIDRun->findHit(pev.getProteinAccession()));
              vector<Size> initVec{toMerge.size(), 0u};
              initVec[repBatchToPut] = newProtIDRun.getHits().size();
              acc2ProtHitIdxPerBatch[pev.getProteinAccession()] = move(initVec);
            }
          }
        }
        //}
      }
    }
    cmap.setProteinIdentifications(newProtIDs);
  }
} // namespace OpenMS
