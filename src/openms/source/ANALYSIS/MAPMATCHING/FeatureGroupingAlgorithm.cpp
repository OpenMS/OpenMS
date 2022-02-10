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
    // set row ID for each feature
    for (size_t i = 0; i < out.size(); i++)
    {
      out[i].setMetaValue("row ID", i+1);
    }
    // set partners for each feature (row IDs of other features with intersecting adduct groups, separated by semicolon)
    // store all related groups in a map with key index (equal to i) and value a set of all groups from other features with
    // intersecting groups, wil be used later to create the annotation network number.
    map<int, set<String>> related_groups;
    // loop over each feature (i) ...
    for (size_t i = 0; i < out.size(); ++i)
    {
      // if a feature has no groups mv, skip it
      if (!out[i].metaValueExists("groups"))
      {
        continue;
      }
      // get the groups of the feature as a const & to list of strings
      const vector<String>& groups = out[i].getMetaValue("groups").toStringList();
      // add each group in feature groups to all related groups with key i
      for (String group: groups)
      {
        related_groups[i].insert(group);
      }
      // ...  and compare to all other features (j)
      for (size_t j = 0; j < out.size(); ++j)
      {
        // same here, if other feature has no groups mv, skip it, also if feature and other feature are the same (i == j)
        if (!out[j].metaValueExists("groups") || (i == j))
        {
          continue;
        }
        // get the groups of the feature as a const & to list of strings
        const vector<String>& other_groups = out[j].getMetaValue("groups").toStringList();
        // create a string of lists: intersections, to store intersecting groups between feature and other feature
        vector<String> intersection;
        std::set_intersection(groups.begin(), groups.end(), other_groups.begin(), other_groups.end(), std::back_inserter(intersection));
        // if intersections not empty the features are related, and row IDs added to partners mv of feature (i)
        if (!intersection.empty())
        {
          // add each group in other_groups to all related groups with key i
          for (String group: other_groups)
          {
            related_groups[i].insert(group);
          }
          // check if feature has partners, if so add  other feature row ID with semi colon, else create new partners mv with row ID
          if (out[i].metaValueExists("partners"))
          {
            out[i].setMetaValue("partners", out[i].getMetaValue("partners").toString()
                                +";"
                                +out[j].getMetaValue("row ID").toString());
          } else
          {
            out[i].setMetaValue("partners", out[j].getMetaValue("row ID"));
          }
        }
      }
    }

    // all via partners and groups related features will get the same network annotation number mv
    Int annotation_network_number = 1;
    // if feature has been annotated already (processed), it's index (i or j) gets added to this list
    set<size_t> already_processed;
    // loop over every feature (i) ...
    for (size_t i=0; i<out.size(); i++)
    {
      // only proceed if feature index (i) in related_groups map feature has not been processed
      if (related_groups.count(i) > 0 && !already_processed.count(i))
      {
        // annotate feature
        out[i].setMetaValue("annotation network number", annotation_network_number);
        // ... and compare to all following features (j) (not to every feature!)
        for (size_t j=i+1; j<out.size(); j++)
        {
          if (related_groups.count(j) > 0)
          {
            // check for intersection between the related groups of feature (i) and other feature (j)
            vector<String> intersection;
            std::set_intersection(related_groups[i].begin(), related_groups[i].end(),
                                  related_groups[j].begin(), related_groups[j].end(),
                                  std::back_inserter(intersection));
            if (!intersection.empty())
            {
              // annotate other feature and add j to already_processed
              out[j].setMetaValue("annotation network number", annotation_network_number);
              already_processed.insert(j);
            }
          }
        }
        // increment annotation network number
        annotation_network_number++;
      }
    }

    // remove "groups" mv
    for (size_t i = 0; i < out.size(); i++)
    {
      out[i].removeMetaValue("groups");
    }
 }

  FeatureGroupingAlgorithm::~FeatureGroupingAlgorithm()
  {
  }

} //namespace OpenMS
