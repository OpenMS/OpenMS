// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Stephan Aiche, Frederic Lehnert, Fabian Kriegel $
// --------------------------------------------------------------------------

#include <OpenMS/SIMULATION/LABELING/ICPLLabeler.h>
//#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <vector>

using std::vector;

namespace OpenMS
{

	ICPLLabeler::ICPLLabeler()
    : BaseLabeler()
	{
		setName("ICPLLabeler");

		//defaults for RT-Fix and protein-labeling
		defaults_.setValue("ICPL_fixed_rtshift", 0.0, "Fixed retention time shift between labeled pairs. If set to 0.0 only the retention times, computed by the RT model step are used.");
		//defaults for protein-labeling
    defaults_.setValue("label_proteins","true","Enables protein-labeling. (select 'false' if you only need peptide-labeling)");
    defaults_.setValidStrings("label_proteins", StringList::create("true,false"));

    // labels
    defaults_.setValue("ICPL_light_channel_label","UniMod:365","UniMod Id of the light channel ICPL label.",StringList::create("advanced"));
    defaults_.setValue("ICPL_medium_channel_label","UniMod:687","UniMod Id of the medium channel ICPL label.",StringList::create("advanced"));
    defaults_.setValue("ICPL_heavy_channel_label","UniMod:364","UniMod Id of the heavy channel ICPL label.",StringList::create("advanced"));

		defaultsToParam_();
	}

	ICPLLabeler::~ICPLLabeler()
	{
	}

	void ICPLLabeler::preCheck(Param & /* param */) const
	{
		// we do not have any prerequisite
	}

  void ICPLLabeler::addLabelToProteinHits_(FeatureMapSim &features, const String &label) const
  {
    // check if proteinIdentification exists before accessing it
    if(features.getProteinIdentifications().empty()) return;
    for(std::vector<ProteinHit>::iterator protein_hit = features.getProteinIdentifications()[0].getHits().begin();
      protein_hit != features.getProteinIdentifications()[0].getHits().end();
      ++protein_hit)
    {
      AASequence aa(protein_hit->getSequence());
      // modify only if the term is accessible
      if(aa.getNTerminalModification() == "")
      {
        aa.setNTerminalModification(label);
        protein_hit->setSequence(aa.toString());
      }
    }
  }

	void ICPLLabeler::setUpHook(FeatureMapSimVector & features)
	{
		// channel check
		if(features.size() < 2 || features.size() > 3)
		{
			throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "We currently support only 2- or 3-channel ICPL");
		}

    if (param_.getValue("label_proteins") == "true")	// loop for protein-labeling (pre-digest-labeling)
		{
			// first channel labeling (light)
      addLabelToProteinHits_(features[0], light_channel_label_);

			// second channel labeling (medium)
      addLabelToProteinHits_(features[1], medium_channel_label_);

      // third channel labeling (heavy) .. if exists
      if(features.size() == 3)
			{
        addLabelToProteinHits_(features[2], heavy_channel_label_);
      }
		}
	}

	void ICPLLabeler::postDigestHook(FeatureMapSimVector & features_to_simulate)
	{
		FeatureMapSim& light_labeled_features = features_to_simulate[0];
		FeatureMapSim& medium_labeled_features = features_to_simulate[1];

    if (param_.getValue("label_proteins") == "false")	// loop for peptide-labeling (post-digest-labeling)
		{
			// iterate over first map for light labeling
			for(FeatureMapSim::iterator lf_iter = light_labeled_features.begin() ; lf_iter != light_labeled_features.end() ; ++lf_iter)
			{
				lf_iter->ensureUniqueId();
        addModificationToPeptideHit_(*lf_iter, light_channel_label_);
			}

			// iterate over second map for medium labeling
			for(FeatureMapSim::iterator lf_iter = medium_labeled_features.begin() ; lf_iter != medium_labeled_features.end() ; ++lf_iter)
			{
				lf_iter->ensureUniqueId();
        addModificationToPeptideHit_(*lf_iter,medium_channel_label_);
			}

			if(features_to_simulate.size() == 3) //third channel labeling can only be done, if a third channel exist
			{
				FeatureMapSim& heavy_labeled_features = features_to_simulate[2];

				// iterate over third map
				for(FeatureMapSim::iterator lf_iter = heavy_labeled_features.begin() ; lf_iter != heavy_labeled_features.end() ; ++lf_iter)
				{
					lf_iter->ensureUniqueId();
          addModificationToPeptideHit_(*lf_iter, heavy_channel_label_);
				}
			}
		}

		// merge the generated feature maps and create consensus
		FeatureMapSim final_feature_map = mergeProteinIdentificationsMaps_(features_to_simulate);

    if (features_to_simulate.size() == 2)	// merge_modus for two FeatureMaps
		{
      // create index of light channel features for easy mapping of medium-to-light channel
			Map<String, Feature> light_labeled_features_index;
			for(FeatureMapSim::iterator light_labeled_features_iter = light_labeled_features.begin();
			light_labeled_features_iter != light_labeled_features.end();
			++light_labeled_features_iter)
			{
				(*light_labeled_features_iter).ensureUniqueId();
				light_labeled_features_index.insert(std::make_pair(
                                              getUnmodifiedAASequence_((*light_labeled_features_iter), light_channel_label_),
                                              *light_labeled_features_iter
                                              ));
			}

			// iterate over second map
			for(FeatureMapSim::iterator medium_labeled_feature_iter = medium_labeled_features.begin() ; medium_labeled_feature_iter != medium_labeled_features.end() ; ++medium_labeled_feature_iter)
			{
				AASequence medium_labeled_feature_sequence = (*medium_labeled_feature_iter).getPeptideIdentifications()[0].getHits()[0].getSequence();

				// guarantee uniqueness
				(*medium_labeled_feature_iter).ensureUniqueId();

				// check if we have a pair
        if(light_labeled_features_index.has(getUnmodifiedAASequence_((*medium_labeled_feature_iter), medium_channel_label_)))
				{
					// own scope as we don't know what happens to 'f_modified' once we call erase() below
          Feature& light_labeled_feature = light_labeled_features_index[getUnmodifiedAASequence_((*medium_labeled_feature_iter), medium_channel_label_)];
					// guarantee uniquenes
					light_labeled_feature.ensureUniqueId();

					if(medium_labeled_feature_sequence.isModified()) // feature has a medium ICPL-Label and is not equal to light-labeled
					{
						// add features to final map
						final_feature_map.push_back(*medium_labeled_feature_iter);
						final_feature_map.push_back(light_labeled_feature);

						// create consensus feature
						ConsensusFeature cf;
						cf.insert(0, *medium_labeled_feature_iter);
						cf.insert(0, light_labeled_feature);

						consensus_.push_back(cf);
						
						// remove light-labeled feature
            light_labeled_features_index.erase(getUnmodifiedAASequence_((*medium_labeled_feature_iter), medium_channel_label_));
					}
					else
					{
						// merge features since they are equal
						Feature final_feature = mergeFeatures_(*medium_labeled_feature_iter, medium_labeled_feature_sequence, light_labeled_features_index);
						final_feature_map.push_back(final_feature);
					}
				}
				else // no ICPL pair, just add the medium-labeled one
				{
					final_feature_map.push_back(*medium_labeled_feature_iter);
				}
			}

			// add singletons from light-labeled channel
			// clean up light-labeled_index
			for(Map<String, Feature>::iterator light_labeled_index_iter = light_labeled_features_index.begin() ; light_labeled_index_iter != light_labeled_features_index.end() ; ++light_labeled_index_iter)
			{
				// the single ones from c0
				final_feature_map.push_back(light_labeled_index_iter->second);
			}
		}
    else if (features_to_simulate.size() == 3) // merge_modus for three Channels
		{
      // create index of light channel features for easy mapping of heavy-to-medium-to-light channel
			Map<String, Feature> light_labeled_features_index;
			for(FeatureMapSim::iterator light_labeled_features_iter = light_labeled_features.begin();
			light_labeled_features_iter != light_labeled_features.end();
			++light_labeled_features_iter)
			{
				(*light_labeled_features_iter).ensureUniqueId();
        light_labeled_features_index.insert(std::make_pair(
                                              getUnmodifiedAASequence_(*light_labeled_features_iter,light_channel_label_),
                                              *light_labeled_features_iter
                                              ));
			}

      // create index of medium channel features for easy mapping of heavy-to-medium-to-light channel
			Map<String, Feature> medium_labeled_features_index;
			for(FeatureMapSim::iterator medium_labeled_features_iter = medium_labeled_features.begin();
			medium_labeled_features_iter != medium_labeled_features.end();
			++medium_labeled_features_iter)
			{
				(*medium_labeled_features_iter).ensureUniqueId();
        medium_labeled_features_index.insert(std::make_pair(
                                               getUnmodifiedAASequence_((*medium_labeled_features_iter), medium_channel_label_),
                                               *medium_labeled_features_iter
                                               ));
			}

      for(FeatureMapSim::iterator heavy_labeled_feature_iter = features_to_simulate[2].begin() ; heavy_labeled_feature_iter != features_to_simulate[2].end() ; ++heavy_labeled_feature_iter)
			{
				Feature& heavy_feature = *heavy_labeled_feature_iter;
        String heavy_feature_unmodified_sequence = getUnmodifiedAASequence_(heavy_feature, heavy_channel_label_);
				heavy_feature.ensureUniqueId();

				if(light_labeled_features_index.has(heavy_feature_unmodified_sequence) && medium_labeled_features_index.has(heavy_feature_unmodified_sequence))
				{
          // 1st case .. it is a triplet
          if(heavy_feature.getPeptideIdentifications()[0].getHits()[0].getSequence().isModified())
					{
            // if heavy feature is modified, all light and medium channel are also

						// add features to final map
						final_feature_map.push_back(heavy_feature);
						final_feature_map.push_back(medium_labeled_features_index[heavy_feature_unmodified_sequence]);
						final_feature_map.push_back(light_labeled_features_index[heavy_feature_unmodified_sequence]);

            // create triplet consensus feature
						ConsensusFeature c_triplet;
						c_triplet.insert(0,heavy_feature);
						c_triplet.insert(0,light_labeled_features_index[heavy_feature_unmodified_sequence]);
						c_triplet.insert(0,medium_labeled_features_index[heavy_feature_unmodified_sequence]);

						consensus_.push_back(c_triplet);
					}
					else
					{
						// merge all three channels
						Feature c2c1 = mergeFeatures_(heavy_feature, heavy_feature_unmodified_sequence, medium_labeled_features_index);
						Feature completeMerge = mergeFeatures_(c2c1, heavy_feature_unmodified_sequence, light_labeled_features_index);

						final_feature_map.push_back(completeMerge);
					}
					// remove features from indices
					light_labeled_features_index.erase(heavy_feature_unmodified_sequence);
					medium_labeled_features_index.erase(heavy_feature_unmodified_sequence);
				}
				else if(light_labeled_features_index.has(heavy_feature_unmodified_sequence))
				{
					// 2.Fall -> c0 - c2
					if(heavy_feature.getPeptideIdentifications()[0].getHits()[0].getSequence().isModified())
					{
						// add features to final map
						final_feature_map.push_back(heavy_feature);
						final_feature_map.push_back(light_labeled_features_index[heavy_feature_unmodified_sequence]);

						ConsensusFeature c_triplet;
						c_triplet.insert(0,heavy_feature);
						c_triplet.insert(0,light_labeled_features_index[heavy_feature_unmodified_sequence]);

						consensus_.push_back(c_triplet);
					}
					else
					{
						// merge all three channels
						Feature completeMerge = mergeFeatures_(heavy_feature, heavy_feature_unmodified_sequence, light_labeled_features_index);
						final_feature_map.push_back(completeMerge);
					}
					// remove features from indices
					light_labeled_features_index.erase(heavy_feature_unmodified_sequence);
				}
				else if(medium_labeled_features_index.has(heavy_feature_unmodified_sequence))
				{
					// 3.Fall -> c1 - c2
					if(heavy_feature.getPeptideIdentifications()[0].getHits()[0].getSequence().isModified())
					{
						// add features to final map
						final_feature_map.push_back(heavy_feature);
						final_feature_map.push_back(medium_labeled_features_index[heavy_feature_unmodified_sequence]);

						ConsensusFeature c_triplet;
						c_triplet.insert(0,heavy_feature);
						c_triplet.insert(0,medium_labeled_features_index[heavy_feature_unmodified_sequence]);

						consensus_.push_back(c_triplet);
					}
					else
					{
						// merge all
						Feature completeMerge = mergeFeatures_(heavy_feature, heavy_feature_unmodified_sequence, medium_labeled_features_index);
						final_feature_map.push_back(completeMerge);
					}
					// remove features from indices
					medium_labeled_features_index.erase(heavy_feature_unmodified_sequence);
				}
				else
				{
					// 4.Fall -> alleine
					final_feature_map.push_back(heavy_feature);
				}
			}

			// clean up medium-labeled_index
			for(Map<String, Feature>::iterator medium_labeled_index_iter = medium_labeled_features_index.begin() ; medium_labeled_index_iter != medium_labeled_features_index.end() ; ++medium_labeled_index_iter)
			{
				Feature& medium_labeled_feature = medium_labeled_index_iter->second;
				medium_labeled_feature.ensureUniqueId();

				String medium_labeled_feature_unmodified_sequence = medium_labeled_feature.getPeptideIdentifications()[0].getHits()[0].getSequence().toUnmodifiedString();

				if(light_labeled_features_index.has(medium_labeled_feature_unmodified_sequence))
				{
					// 1.Fall paar zwischen c0 und c1
					if(medium_labeled_feature.getPeptideIdentifications()[0].getHits()[0].getSequence().isModified())
					{
						// add features to final map
						final_feature_map.push_back(medium_labeled_feature);
						final_feature_map.push_back(light_labeled_features_index[medium_labeled_feature_unmodified_sequence]);

						ConsensusFeature c_triplet;
						c_triplet.insert(0,medium_labeled_feature);
						c_triplet.insert(0,light_labeled_features_index[medium_labeled_feature_unmodified_sequence]);

						consensus_.push_back(c_triplet);
					}
					else
					{
						// merge
						Feature completeMerge = mergeFeatures_(medium_labeled_feature, medium_labeled_feature_unmodified_sequence, light_labeled_features_index);
						final_feature_map.push_back(completeMerge);
					}
					// remove features from indices
					light_labeled_features_index.erase(medium_labeled_feature_unmodified_sequence);
				}
				else
				{
					// c1 ist alleine
					final_feature_map.push_back(medium_labeled_feature);
				}
			}

			// clean up light-labeled_index
			for(Map<String, Feature>::iterator light_labeled_index_iter = light_labeled_features_index.begin() ; light_labeled_index_iter != light_labeled_features_index.end() ; ++light_labeled_index_iter)
			{
				// the single ones from c0
				final_feature_map.push_back(light_labeled_index_iter->second);
			}
		}

		features_to_simulate.clear();
		features_to_simulate.push_back(final_feature_map);

		consensus_.setProteinIdentifications(final_feature_map.getProteinIdentifications());
		ConsensusMap::FileDescription map_description;
		map_description.label = "Simulation (Labeling Consensus)";
		map_description.size = features_to_simulate.size();
		consensus_.getFileDescriptions()[0] = map_description;
	}


  Feature ICPLLabeler::mergeFeatures_(Feature& feature_to_merge, const AASequence& labeled_feature_sequence, Map<String, Feature>& feature_index) const
	{
		// merge with feature from first map (if it exists)
    if(feature_index.count(labeled_feature_sequence.toString()) != 0)
		{
			// we only merge abundance and use feature from first map
      Feature new_f = feature_index[labeled_feature_sequence.toString()];

			new_f.setMetaValue(getChannelIntensityName(1), new_f.getIntensity());
      new_f.setMetaValue(getChannelIntensityName(2), feature_to_merge.getIntensity());
      new_f.setIntensity(new_f.getIntensity() + feature_to_merge.getIntensity());

      mergeProteinAccessions_(new_f, feature_to_merge);

			// remove feature from index
      feature_index.erase(labeled_feature_sequence.toString());

			return new_f;
		}
		else
		{
			// simply add feature from second channel, since we have no corresponding feature in the first channel
      return feature_to_merge;
		}
	}

	void ICPLLabeler::addModificationToPeptideHit_(Feature& feature, const String & modification) const
	{
		vector<PeptideHit> pep_hits(feature.getPeptideIdentifications()[0].getHits());
		AASequence modified_sequence(pep_hits[0].getSequence());
    if(modified_sequence.getNTerminalModification() == "")
    {
      // attach label only if the nterm is accessible
      modified_sequence.setNTerminalModification(modification);
      pep_hits[0].setSequence(modified_sequence);
      feature.getPeptideIdentifications()[0].setHits(pep_hits);
    }
	}


	void ICPLLabeler::postRTHook(FeatureMapSimVector & features_to_simulate)
	{
		DoubleReal rt_shift = param_.getValue("ICPL_fixed_rtshift");

		if(rt_shift != 0.0)
		{
			Map<UInt64, Feature*> id_map;
			FeatureMapSim& feature_map = features_to_simulate[0];
			for(FeatureMapSim::Iterator it = feature_map.begin() ; it != feature_map.end() ; ++it)
			{
				id_map.insert(std::make_pair<UInt64,Feature*>(it->getUniqueId(), &(*it)));
			}

			// recompute RT of pairs
			for(ConsensusMap::Iterator consensus_it = consensus_.begin() ; consensus_it != consensus_.end() ; ++consensus_it)
			{
				ConsensusFeature& cf = *consensus_it;

				// check if these features are still available and were not removed during RT sim
				bool complete = true;

				for(ConsensusFeature::iterator cfit = cf.begin() ; cfit != cf.end() ; ++cfit)
				{
					complete &= id_map.has(cfit->getUniqueId());
				}

				if(complete)
				{
					Feature* f1 = id_map[(cf.begin())->getUniqueId()];
					Feature* f2 = id_map[(++cf.begin())->getUniqueId()];

					// the lighter one should be the unmodified one
					EmpiricalFormula ef1 = (f1->getPeptideIdentifications())[0].getHits()[0].getSequence().getFormula();
					EmpiricalFormula ef2 = (f2->getPeptideIdentifications())[0].getHits()[0].getSequence().getFormula();

					if(ef1.getMonoWeight() < ef2.getMonoWeight())
					{
						f2->setRT(f1->getRT() + rt_shift);
					}
					else
					{
						f1->setRT(f2->getRT() + rt_shift);
					}
				}
			}
		}
	}

	void ICPLLabeler::postDetectabilityHook(FeatureMapSimVector & /* features_to_simulate */)
	{
		// no changes to the features .. nothing todo here
	}

	void ICPLLabeler::postIonizationHook(FeatureMapSimVector & /* features_to_simulate */)
	{
		// no changes to the features .. nothing todo here
	}

	void ICPLLabeler::postRawMSHook(FeatureMapSimVector & features_to_simulate)
	{
		recomputeConsensus_(features_to_simulate[0]);
	}

	void ICPLLabeler::postRawTandemMSHook(FeatureMapSimVector & /* features_to_simulate */, MSSimExperiment & /* simulated map */)
	{
		// no changes to the features .. nothing todo here
	}

  void ICPLLabeler::updateMembers_()
  {
    light_channel_label_ = param_.getValue("ICPL_light_channel_label");
    medium_channel_label_ = param_.getValue("ICPL_medium_channel_label");
    heavy_channel_label_ = param_.getValue("ICPL_heavy_channel_label");
  }

  String ICPLLabeler::getUnmodifiedAASequence_(const Feature &feature, const String &label) const
  {
    AASequence unmodified = feature.getPeptideIdentifications()[0].getHits()[0].getSequence();
    if(unmodified.getNTerminalModification() == label)
    {
      unmodified.setNTerminalModification(""); // remove terminal modification, if it is the channel specific one
    }
    return unmodified.toString();
  }
} // namespace OpenMS
