// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
      for (const auto& desc : consensus.getColumnHeaders())
      {
        Size counter = mapid_table.size();
        mapid_table[make_pair(i, desc.first)] = counter;
        out.getColumnHeaders()[counter] = desc.second;
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
    for (auto& cons : out)
    {
      ConsensusFeature adjusted = ConsensusFeature(
        static_cast<BaseFeature>(cons)); // remove sub-features
      for (const auto& sub : cons.getFeatures())
      {
        UInt64 id = sub.getUniqueId();
        Size map_index = sub.getMapIndex();
        ConsensusMap::ConstIterator origin = feat_lookup[map_index][id];
        for (const auto& hdl : origin->getFeatures())
        {
          FeatureHandle handle = hdl;
          Size new_id = mapid_table[make_pair(map_index, handle.getMapIndex())];
          handle.setMapIndex(new_id);
          adjusted.insert(handle);
        }
      }
      cons = adjusted;

      for (auto& id : cons.getPeptideIdentifications())
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

  FeatureGroupingAlgorithm::~FeatureGroupingAlgorithm() = default;

} //namespace OpenMS
