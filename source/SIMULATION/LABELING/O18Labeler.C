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

#include <OpenMS/SIMULATION/LABELING/O18Labeler.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>

using std::vector;

// TODO: change postRawHook to FeatureMapSim?
// TODO: implement correct consensus in postRaw


namespace OpenMS
{
  O18Labeler::O18Labeler()
    : BaseLabeler()
  {
    setName("O18Labeler");
    defaults_.setValue("labeling_efficiency", 1.0, "Describes the distribution of the labeled peptide over the different states (unlabeled, mono- and dilabeled)");
    defaults_.setMinFloat("labeling_efficiency", 0.0);
    defaults_.setMaxFloat("labeling_efficiency", 1.0);
    defaultsToParam_();
  }

  O18Labeler::~O18Labeler()
  {
  }

  void O18Labeler::preCheck(Param & param) const
  {
    // check for trypsin
    // TODO: add other enzymes
    if(param.getValue("Digestion:enzyme") != "Trypsin")
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "18 O Labeling requires digestion with Trypsin");
    }
  }

  void O18Labeler::setUpHook(FeatureMapSimVector & features)
  {
    // no action here .. just check for 2 channels
    if(features.size() != 2)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "18 O Labeling only works with 2 channels.");
    }
  }

  /// Labeling between digestion and rt simulation
  void O18Labeler::postDigestHook(FeatureMapSimVector & features_to_simulate )
  {
    SimIntensityType labeling_efficiency = param_.getValue("labeling_efficiency");

    // index unlabeled map
    // merge channel one and two into a single feature map
    FeatureMapSim final_feature_map = mergeProteinIdentificationsMaps_(features_to_simulate);
    FeatureMapSim& unlabeled_features = features_to_simulate[0];

    std::map<AASequence, Feature> unlabeled_features_index;

    for(FeatureMapSim::iterator unlabeled_features_iter = unlabeled_features.begin() ;
        unlabeled_features_iter != unlabeled_features.end() ;
        ++unlabeled_features_iter)
    {
      (*unlabeled_features_iter).ensureUniqueId();
      unlabeled_features_index.insert(std::make_pair(
          (*unlabeled_features_iter).getPeptideIdentifications()[0].getHits()[0].getSequence()
          ,
          *unlabeled_features_iter
          ));

    }

    // iterate over second map
    FeatureMapSim& labeled_features = features_to_simulate[1];

    for(FeatureMapSim::iterator lf_iter = labeled_features.begin() ; lf_iter != labeled_features.end() ; ++lf_iter)
    {
      AASequence unmodified_sequence = (*lf_iter).getPeptideIdentifications()[0].getHits()[0].getSequence();

      // check if feature has tryptic c-terminus
      PeptideHit ph = (*lf_iter).getPeptideIdentifications()[0].getHits()[0];
      if(ph.getSequence().getResidue(ph.getSequence().size() - 1) == 'R'
         ||
         ph.getSequence().getResidue(ph.getSequence().size() - 1) == 'K')
      {
        // this one will be modified since it shows a Trypsin-C-Term
        // relevant unimod modifications are
        // Label:18O(1) -- 258
        // Label:18O(2) -- 193

        if(labeling_efficiency != 1.0)
        {                  
          Feature b1(*lf_iter);
          b1.ensureUniqueId();
          Feature b2(*lf_iter);
          b2.ensureUniqueId();

          SimIntensityType total_intensity = (*lf_iter).getIntensity();

          // dilabled
          addModificationToPeptideHit_(b2,"UniMod:193");
          b2.setIntensity(total_intensity * labeling_efficiency  * labeling_efficiency);

          final_feature_map.push_back(b2);

          // mono labeled
          addModificationToPeptideHit_(b1, "UniMod:258");
          b1.setIntensity(total_intensity * 2.0 * (1 - labeling_efficiency) * labeling_efficiency);

          final_feature_map.push_back(b1);

          // merge unlabeled with possible labeled feature
          // modify unlabeled intensity
          (*lf_iter).setIntensity(total_intensity * (1 - labeling_efficiency) * (1 - labeling_efficiency));

          // all three partial intensities from above should add up to 1 now

          // generate consensus feature
          ConsensusFeature cf;
          // add mono and &dilabeled variant to ConsensusFeature
          cf.insert(0, b1);
          cf.insert(0, b2);

          // merge unlabeled with unlabeled from other channel (if it exists)
          Feature final_unlabeled_feature = mergeFeatures_(*lf_iter, unmodified_sequence, unlabeled_features_index);
          final_unlabeled_feature.ensureUniqueId();
          cf.insert(0,final_unlabeled_feature);

          consensus_.push_back(cf);
          final_feature_map.push_back(final_unlabeled_feature);

          // remove unlabeled feature
          unlabeled_features_index.erase(unmodified_sequence);
        }
        else
        {
          // generate labeled feature
          // labeling_efficiency is 100% so we transform the complete
          // feature in a dilabeled feature
          addModificationToPeptideHit_(*lf_iter, "UniMod:193");
          (*lf_iter).ensureUniqueId();
          final_feature_map.push_back(*lf_iter);

          // add corresponding feature if it exists
          // and generate consensus feature for the unlabeled/labeled pair
          if(unlabeled_features_index.count(unmodified_sequence) != 0)
          {
            ConsensusFeature cf;
            final_feature_map.push_back(unlabeled_features_index[unmodified_sequence]);
            cf.insert(0, *lf_iter);
            cf.insert(0, unlabeled_features_index[unmodified_sequence]);

            // remove unlabeled feature
            unlabeled_features_index.erase(unmodified_sequence);

            consensus_.push_back(cf);
          }
        }
      }
      else
      {
        Feature final_feature = mergeFeatures_(*lf_iter, unmodified_sequence, unlabeled_features_index);
        final_feature_map.push_back(final_feature);
      }
    }

    // add remaining feature from first channel
    for(std::map<AASequence, Feature>::iterator remaining_features_iter = unlabeled_features_index.begin() ; remaining_features_iter != unlabeled_features_index.end() ; ++remaining_features_iter)
    {
      final_feature_map.push_back(remaining_features_iter->second);
    }

    features_to_simulate.clear();
    features_to_simulate.push_back(final_feature_map);

    consensus_.setProteinIdentifications(final_feature_map.getProteinIdentifications());
    ConsensusMap::FileDescription map_description;
    map_description.label = "Simulation (Labeling Consensus)";
    map_description.size = features_to_simulate.size();
    consensus_.getFileDescriptions()[0] = map_description;
  }


  Feature O18Labeler::mergeFeatures_(Feature& labeled_channel_feature, const AASequence& unmodified_sequence, std::map<AASequence, Feature>& unlabeled_features_index) const
  {
    // merge with feature from first map (if it exists)
    if(unlabeled_features_index.count(unmodified_sequence) != 0)
    {
      // we only merge abundance and use feature from first map
      Feature new_f = unlabeled_features_index[unmodified_sequence];

      new_f.setMetaValue(getChannelIntensityName(1), new_f.getIntensity());
      new_f.setMetaValue(getChannelIntensityName(2), labeled_channel_feature.getIntensity());

      new_f.setIntensity(new_f.getIntensity() + labeled_channel_feature.getIntensity());

      mergeProteinAccessions_(new_f, labeled_channel_feature);

      // remove feature from index
      unlabeled_features_index.erase(unmodified_sequence);

      return new_f;
    }
    else
    {
      // simply add feature from labeled channel, since we
      // have no corresponding feature in the unlabeled channel
      return labeled_channel_feature;
    }
  }

  void O18Labeler::addModificationToPeptideHit_(Feature& feature, const String& modification) const
  {
    vector<PeptideHit> pep_hits(feature.getPeptideIdentifications()[0].getHits());
    AASequence modified_sequence(pep_hits[0].getSequence());
    modified_sequence.setCTerminalModification(modification);
    pep_hits[0].setSequence(modified_sequence);
    feature.getPeptideIdentifications()[0].setHits(pep_hits);
  }

  /// Labeling between RT and Detectability
  void O18Labeler::postRTHook(FeatureMapSimVector & /* features_to_simulate */)
  {
  }

  /// Labeling between Detectability and Ionization
  void O18Labeler::postDetectabilityHook(FeatureMapSimVector & /* features_to_simulate */)
  {
  }

  /// Labeling between Ionization and RawMS
  void O18Labeler::postIonizationHook(FeatureMapSimVector & /* features_to_simulate */)
  {
  }

  /// Labeling after RawMS
  void O18Labeler::postRawMSHook(FeatureMapSimVector & features_to_simulate )
  {
    recomputeConsensus_(features_to_simulate[0]);
  }  

  void O18Labeler::postRawTandemMSHook(FeatureMapSimVector &, MSSimExperiment &)
  {
  }
} // namespace OpenMS
