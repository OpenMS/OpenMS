// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h>

// Derived classes are included here
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmLabeled.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmUnlabeled.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmQT.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmKD.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Factory.h>

#include <unordered_set>
#include <regex>

using namespace std;

namespace OpenMS
{
  //register products here
  void FeatureGroupingAlgorithm::registerChildren()
  {
    Factory<FeatureGroupingAlgorithm>::registerProduct(FeatureGroupingAlgorithmLabeled::getProductName(), &FeatureGroupingAlgorithmLabeled::create);
    Factory<FeatureGroupingAlgorithm>::registerProduct(FeatureGroupingAlgorithmUnlabeled::getProductName(), &FeatureGroupingAlgorithmUnlabeled::create);
    Factory<FeatureGroupingAlgorithm>::registerProduct(FeatureGroupingAlgorithmQT::getProductName(), &FeatureGroupingAlgorithmQT::create);
    Factory<FeatureGroupingAlgorithm>::registerProduct(FeatureGroupingAlgorithmKD::getProductName(), &FeatureGroupingAlgorithmKD::create);

  }

  FeatureGroupingAlgorithm::FeatureGroupingAlgorithm() :
    DefaultParamHandler("FeatureGroupingAlgorithm")
  {
  }

  void FeatureGroupingAlgorithm::group(const vector<ConsensusMap>& maps, ConsensusMap& out)
  {
    OPENMS_LOG_WARN << "FeatureGroupingAlgorithm::group() does not support ConsensusMaps directly. Converting to FeatureMaps." << endl;

    vector<FeatureMap> maps_f;
    for (Size i = 0; i < maps.size(); ++i)
    {
      FeatureMap fm;
      MapConversion::convert(maps[i], true, fm);
      maps_f.push_back(fm);
    }
    // call FeatureMap version of group()
    group(maps_f, out);
  }

  void FeatureGroupingAlgorithm::transferSubelements(const vector<ConsensusMap>& maps, ConsensusMap& out) const
  {
    // accumulate file descriptions from the input maps:
    // cout << "Updating file descriptions..." << endl;
    out.getColumnHeaders().clear();
    // mapping: (input file index / map index assigned by the linkers, old map index) -> new map index
    map<pair<Size, UInt64>, Size> mapid_table;
    for (Size i = 0; i < maps.size(); ++i)
    {
      const ConsensusMap& consensus = maps[i];
      for (ConsensusMap::ColumnHeaders::const_iterator desc_it = consensus.getColumnHeaders().begin(); desc_it != consensus.getColumnHeaders().end(); ++desc_it)
      {
        Size counter = mapid_table.size();
        mapid_table[make_pair(i, desc_it->first)] = counter;
        out.getColumnHeaders()[counter] = desc_it->second;
      }
    }

    // look-up table: input map -> unique ID -> consensus feature
    // cout << "Creating look-up table..." << endl;
    vector<map<UInt64, ConsensusMap::ConstIterator> > feat_lookup(maps.size());
    for (Size i = 0; i < maps.size(); ++i)
    {
      const ConsensusMap& consensus = maps[i];
      for (ConsensusMap::ConstIterator feat_it = consensus.begin();
           feat_it != consensus.end(); ++feat_it)
      {
        // do NOT use "id_lookup[i][feat_it->getUniqueId()] = feat_it;" here as
        // you will get "attempt to copy-construct an iterator from a singular
        // iterator" in STL debug mode:
        feat_lookup[i].insert(make_pair(feat_it->getUniqueId(), feat_it));
      }
    }
    // adjust the consensus features:
    // cout << "Adjusting consensus features..." << endl;
    for (ConsensusMap::iterator cons_it = out.begin(); cons_it != out.end(); ++cons_it)
    {
      ConsensusFeature adjusted = ConsensusFeature(
        static_cast<BaseFeature>(*cons_it)); // remove sub-features
      for (ConsensusFeature::HandleSetType::const_iterator sub_it = cons_it->getFeatures().begin(); sub_it != cons_it->getFeatures().end(); ++sub_it)
      {
        UInt64 id = sub_it->getUniqueId();
        Size map_index = sub_it->getMapIndex();
        ConsensusMap::ConstIterator origin = feat_lookup[map_index][id];
        for (ConsensusFeature::HandleSetType::const_iterator handle_it = origin->getFeatures().begin(); handle_it != origin->getFeatures().end(); ++handle_it)
        {
          FeatureHandle handle = *handle_it;
          Size new_id = mapid_table[make_pair(map_index, handle.getMapIndex())];
          handle.setMapIndex(new_id);
          adjusted.insert(handle);
        }
      }
      *cons_it = adjusted;

      for (auto& id : cons_it->getPeptideIdentifications())
      {
        // if old_map_index is not present, there was no map_index in the beginning,
        // therefore the newly assigned map_index cannot be "corrected"
        // -> remove the MetaValue to be consistent.
        if (id.metaValueExists("old_map_index"))
        {
          Size old_map_index = id.getMetaValue("old_map_index");
          Size file_index = id.getMetaValue("map_index");
          Size new_idx = mapid_table[make_pair(file_index, old_map_index)];
          id.setMetaValue("map_index", new_idx);
          id.removeMetaValue("old_map_index");
        }
        else
        {
          id.removeMetaValue("map_index");
        }
      }
    }
    for (auto& id : out.getUnassignedPeptideIdentifications())
    {
      // if old_map_index is not present, there was no map_index in the beginning,
      // therefore the newly assigned map_index cannot be "corrected"
      // -> remove the MetaValue to be consistent.
      if (id.metaValueExists("old_map_index"))
      {
        Size old_map_index = id.getMetaValue("old_map_index");
        Size file_index = id.getMetaValue("map_index");
        Size new_idx = mapid_table[make_pair(file_index, old_map_index)];
        id.setMetaValue("map_index", new_idx);
        id.removeMetaValue("old_map_index");
      }
      else
      {
        id.removeMetaValue("map_index");
      }
    }
  }

  void FeatureGroupingAlgorithm::annotateIonIdentityNetworks(ConsensusMap& out) const
  {
    // set row ID for each feature and reformat best ion
    set<size_t> all;
    for (size_t i = 0; i < out.size(); i++)
    {
      all.insert(i);
      out[i].setMetaValue(Constants::UserParam::ROW_ID, i+1);
      if (out[i].metaValueExists(Constants::UserParam::BEST_ION))
      {
        String best_ion = out[i].getMetaValue(Constants::UserParam::BEST_ION);
        String formatted_best_ion = "[M";
        std::smatch matches;
        std::regex_search(best_ion, matches, std::regex(".*[^\\d]+"));
        if (matches.ready())
        {
          if (matches.str(0).rfind("-", 0) == 0)
          {
            formatted_best_ion += matches.str(0);
          } else
          {
            formatted_best_ion += "+";
            formatted_best_ion += matches.str(0);
          }
        }
        formatted_best_ion += "]";
        std::regex_search(best_ion, matches, std::regex("\\d$"));
        if (matches.ready())
        {
          if (matches.str(0) == "1")
          {
            formatted_best_ion += "+";
          } else
          {
            formatted_best_ion += matches.str(0);
            formatted_best_ion += "+";
          }
        }
        out[i].setMetaValue(Constants::UserParam::BEST_ION, formatted_best_ion);
      }
    }

    Int annotation_network_number = 1;

    while (!all.empty())
    {
      size_t n = *all.begin();
      all.erase(n);
      if (!out[n].metaValueExists(Constants::UserParam::LINKED_GROUPS)) continue;
      set<size_t> connected;
      set<size_t> connected_queue = {n};
      while (!connected_queue.empty())
      {
        size_t i = *connected_queue.begin();
        connected_queue.erase(i);
        connected.insert(i);
        vector<String> groups_i = out[i].getMetaValue(Constants::UserParam::LINKED_GROUPS).toStringList();
        set<size_t> to_erase;
        for (const auto& j: all)
        {
          if (!out[j].metaValueExists(Constants::UserParam::LINKED_GROUPS)) continue;
          vector<String> groups_j = out[j].getMetaValue(Constants::UserParam::LINKED_GROUPS).toStringList();
          vector<String> intersection;
          set_intersection(groups_i.begin(), groups_i.end(), groups_j.begin(), groups_j.end(), back_inserter(intersection));
          if (!intersection.empty())
          {
            connected_queue.insert(j);
            to_erase.insert(j);
          }
        }
        for (const auto& e: to_erase) all.erase(e);
      }

      cout << "connected components: ";
      for (const auto& xx: connected) cout << xx << " ";
      cout << endl;

      for (const auto& i : connected)
      {
        out[i].setMetaValue(Constants::UserParam::ANNOTATION_NETWORK_NUMBER, annotation_network_number);
        const vector<String>& groups_i = out[i].getMetaValue(Constants::UserParam::LINKED_GROUPS).toStringList();
        for (const auto& j : connected)
        {
          const vector<String>& groups_j = out[j].getMetaValue(Constants::UserParam::LINKED_GROUPS).toStringList();
          vector<String> intersection;
          set_intersection(groups_i.begin(), groups_i.end(), groups_j.begin(), groups_j.end(), back_inserter(intersection));
          
          if (i == j || intersection.empty()) continue;
          if (out[i].metaValueExists(Constants::UserParam::ADDUCT_PARTNERS))
          {
            out[i].setMetaValue(Constants::UserParam::ADDUCT_PARTNERS, out[i].getMetaValue(Constants::UserParam::ADDUCT_PARTNERS).toString()
                                +";"
                                +out[j].getMetaValue(Constants::UserParam::ROW_ID).toString());
          } else
          {
            out[i].setMetaValue(Constants::UserParam::ADDUCT_PARTNERS, out[j].getMetaValue(Constants::UserParam::ROW_ID));
          }
        }
      }
      annotation_network_number++;
    }
    // remove Constants::UserParam::LINKED_GROUPS mv
    // for (size_t i = 0; i < out.size(); i++)
    // {
    //   out[i].removeMetaValue(Constants::UserParam::LINKED_GROUPS);
    // }
 }

  FeatureGroupingAlgorithm::~FeatureGroupingAlgorithm()
  {
  }

} //namespace OpenMS
