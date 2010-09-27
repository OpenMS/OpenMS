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

#include <OpenMS/SIMULATION/LABELING/SILACLabeler.h>

#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

#include <vector>

using std::vector;

namespace OpenMS
{

  SILACLabeler::SILACLabeler()
  {
    setName("SILACLabeler");
    defaults_.setValue("SILAC_modifications", "UniMod:944", "Comma separated list of SILAC modifications that will be added to the labeled channel.");
    defaultsToParam_();
  }

  SILACLabeler::~SILACLabeler()
  {

  }

  void SILACLabeler::preCheck(Param & /* param */) const
  {
    // we do not have any prerequisite
  }

  void SILACLabeler::setUpHook(FeatureMapSimVector & features)
  {
    // check for 2 channels
    if(features.size() != 2)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "We currently support only 2-channel SILAC");
    }

    // parse modifications
    String mod_list = param_.getValue("SILAC_modifications");
    StringList modifications = StringList::create(mod_list);

    Map<String, const ResidueModification*> site_resmod_map;

    for(StringList::iterator it = modifications.begin() ; it != modifications.end() ; ++it)
    {
      try
      {
        const ResidueModification* mod = &ModificationsDB::getInstance()->getModification(*it);
        if(site_resmod_map.has(mod->getOrigin()))
        {
          throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "MSSimulator can only handle one modification specific for the residue \"" + mod->getOrigin() + "\"");
        }
        else
        {
          site_resmod_map.insert(std::make_pair<String, const ResidueModification*>(mod->getOrigin(), mod));
        }
      }
      catch (Exception::ElementNotFound ex)
      {
        // nothing to clean up here
        ex.setMessage("The modification \"" + *it + "\" could not be found in the local UniMod DB! Please check if you used the correct format (e.g. UniMod:Accession#)");
        throw ex;
      }
    }

    // label second channel
    for(std::vector<ProteinHit>::iterator protein_hit = features[1].getProteinIdentifications()[0].getHits().begin();
      protein_hit != features[1].getProteinIdentifications()[0].getHits().end();
      ++protein_hit)
    {
      AASequence aa(protein_hit->getSequence());

      for(Size r = 0 ; r < aa.size() ; ++r)
      {
        if(site_resmod_map.has(aa[r].getOneLetterCode()))
        {
          aa.setModification(r, site_resmod_map[aa[r].getOneLetterCode()]->getUniModAccession());
        }
      }

      protein_hit->setSequence(aa.toString());
    }

  }

  void SILACLabeler::postDigestHook(FeatureMapSimVector & features_to_simulate)
  {
    // merge the generated feature maps and create consensus
    FeatureMapSim final_feature_map = mergeProteinIdentificationsMaps_(features_to_simulate);
    FeatureMapSim& unlabeled_features = features_to_simulate[0];

    Map<String, Feature> unlabeled_features_index;

    for(FeatureMapSim::iterator unlabeled_features_iter = unlabeled_features.begin() ;
        unlabeled_features_iter != unlabeled_features.end() ;
        ++unlabeled_features_iter)
    {
      (*unlabeled_features_iter).ensureUniqueId();
      unlabeled_features_index.insert(std::make_pair(
          (*unlabeled_features_iter).getPeptideIdentifications()[0].getHits()[0].getSequence().toString()
          ,
          *unlabeled_features_iter
          ));

    }

    // iterate over second map
    FeatureMapSim& labeled_features = features_to_simulate[1];

    for(FeatureMapSim::iterator lf_iter = labeled_features.begin() ; lf_iter != labeled_features.end() ; ++lf_iter)
    {
      AASequence unmodified_sequence = (*lf_iter).getPeptideIdentifications()[0].getHits()[0].getSequence();

      // guarantee uniqueness
      (*lf_iter).ensureUniqueId();

      // check if we have a pair ..
      if(unlabeled_features_index.has(unmodified_sequence.toUnmodifiedString()))
      {
        // guarantee uniqueness
        unlabeled_features_index[unmodified_sequence.toUnmodifiedString()].ensureUniqueId();

        // add features to final map
        final_feature_map.push_back(*lf_iter);
        final_feature_map.push_back(unlabeled_features_index[unmodified_sequence.toUnmodifiedString()]);

        // create consensus feature
        ConsensusFeature cf;
        cf.insert(0, *lf_iter);
        cf.insert(0, unlabeled_features_index[unmodified_sequence.toUnmodifiedString()]);

        consensus_.push_back(cf);

        // remove unlabeled feature
        unlabeled_features_index.erase(unmodified_sequence.toUnmodifiedString());
      }
      else // no SILAC pair, just add the labeled one
      {
        final_feature_map.push_back(*lf_iter);
      }
    }

    // add remaining feature from first channel
    for(Map<String, Feature>::iterator remaining_features_iter = unlabeled_features_index.begin() ; remaining_features_iter != unlabeled_features_index.end() ; ++remaining_features_iter)
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

  void SILACLabeler::postRTHook(FeatureMapSimVector & /* features_to_simulate */)
  {
    // no changes to the features .. nothing todo here
  }

  void SILACLabeler::postDetectabilityHook(FeatureMapSimVector & /* features_to_simulate */)
  {
    // no changes to the features .. nothing todo here
  }

  void SILACLabeler::postIonizationHook(FeatureMapSimVector & /* features_to_simulate */)
  {
    // no changes to the features .. nothing todo here
  }

  void SILACLabeler::postRawMSHook(FeatureMapSimVector & features_to_simulate)
  {
    recomputeConsensus_(features_to_simulate[0]);
  }

  void SILACLabeler::postRawTandemMSHook(FeatureMapSimVector & /* features_to_simulate */, MSSimExperiment & /* simulated map */)
  {
    // no changes to the features .. nothing todo here
  }
} // namespace OpenMS
