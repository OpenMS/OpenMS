// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Stephan Aiche $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/SIMULATION/LABELING/BaseLabeler.h>
#include <OpenMS/SIMULATION/LABELING/BaseLabeler_impl.h>

#include <map>

using std::vector;
using std::pair;
using std::set;

namespace OpenMS
{

  BaseLabeler::BaseLabeler()
    : DefaultParamHandler("BaseLabeler"),
      rng_(0)
  {
  }

  String BaseLabeler::getChannelIntensityName(const Size channel_index) const
  {
    return String("channel_")+ String(channel_index) + "_intensity";
  }

  FeatureMapSim BaseLabeler::mergeProteinIdentificationsMaps_(const FeatureMapSimVector &maps)
  {
    // we do not have any features yet (or at least we ignore them), so simply iterate over the protein
    // identifications
    std::map<String, ProteinHit> prot_hits;
    Size channel_index = 1;
    for(FeatureMapSimVector::const_iterator maps_iterator = maps.begin() ; maps_iterator != maps.end() ; ++maps_iterator)
    {
      for(std::vector<ProteinHit>::const_iterator protein_hit = (*maps_iterator).getProteinIdentifications()[0].getHits().begin();
        protein_hit != (*maps_iterator).getProteinIdentifications()[0].getHits().end();
        ++protein_hit)
      {
        if(prot_hits.count((*protein_hit).getSequence()))
        { // we already know this protein -- sum up abundances
          SimIntensityType new_intensity = prot_hits[(*protein_hit).getSequence()].getMetaValue("intensity");

          // remember channel intensity
          prot_hits[(*protein_hit).getSequence()].setMetaValue("intensity_" + String(channel_index),new_intensity);

          new_intensity += static_cast<SimIntensityType>((*protein_hit).getMetaValue("intensity"));
          prot_hits[(*protein_hit).getSequence()].setMetaValue("intensity",new_intensity);
        }
        else
        { // new protein hit .. remember
          ProteinHit protHit(*protein_hit);
          protHit.setMetaValue("intensity_" + String(channel_index), protHit.getMetaValue("intensity"));
          prot_hits.insert(std::pair<String, ProteinHit>((*protein_hit).getSequence(), protHit));
        }
      }
      ++channel_index;
    }

    FeatureMapSim final_map;
    ProteinIdentification protIdent;

    for(std::map<String, ProteinHit>::iterator prot_hit_iter = prot_hits.begin() ; prot_hit_iter != prot_hits.end() ; ++prot_hit_iter)
    {
      protIdent.insertHit(prot_hit_iter->second);
    }
    std::vector<ProteinIdentification> protIdents;
    protIdents.push_back(protIdent);
    final_map.setProteinIdentifications(protIdents);

    return final_map;
  }

  void BaseLabeler::mergeProteinAccessions_(Feature& target, const Feature& source) const
  {
    std::vector<String> target_acc (target.getPeptideIdentifications()[0].getHits()[0].getProteinAccessions());
    std::vector<String> source_acc (source.getPeptideIdentifications()[0].getHits()[0].getProteinAccessions());

    std::set<String> unique_acc;
    std::pair<std::set<String>::iterator, bool> result;

    for(vector<String>::iterator target_acc_iterator = target_acc.begin() ; target_acc_iterator != target_acc.end() ; ++target_acc_iterator)
    {
      unique_acc.insert(*target_acc_iterator);
    }

    for(vector<String>::iterator source_acc_iterator = source_acc.begin() ; source_acc_iterator != source_acc.end() ; ++source_acc_iterator)
    {
      result = unique_acc.insert(*source_acc_iterator);

      if(result.second)
      {
        target_acc.push_back(*source_acc_iterator);
      }
    }

    PeptideHit pepHit(target.getPeptideIdentifications()[0].getHits()[0]);
    pepHit.setProteinAccessions(target_acc);

    std::vector<PeptideHit> pepHits;
    pepHits.push_back(pepHit);

    target.getPeptideIdentifications()[0].setHits(pepHits);
  }

  const ConsensusMap& BaseLabeler::getConsensus() const
  {
    return consensus_;
  }
  
} // namespace OpenMS
