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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/SIMULATION/LABELING/ITRAQLabeler.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>

using std::vector;
using std::pair;
using std::set;

namespace OpenMS
{
  ITRAQLabeler::ITRAQLabeler()
    : BaseLabeler(),
		  itraq_type_(),
		  channel_map_(),
		  isotope_corrections_()
  {
		// this needs to come first!
		isotope_corrections_.resize(2);
		isotope_corrections_[0].setMatrix<4,4>(ItraqConstants::ISOTOPECORRECTIONS_FOURPLEX);
		isotope_corrections_[1].setMatrix<8,4>(ItraqConstants::ISOTOPECORRECTIONS_EIGHTPLEX);

		// iTRAQ
		defaults_.setValue("iTRAQ", "4plex", "4plex or 8plex iTRAQ?");
		defaults_.setValidStrings("iTRAQ", StringList::create("4plex,8plex"));

		defaults_.setValue("reporter_mass_shift", 0.1, "Allowed shift (uniformly distributed - left to right) in Da from the expected position (of e.g. 114.1, 115.1)"); 
		defaults_.setMinFloat ("reporter_mass_shift", 0);
		defaults_.setMaxFloat ("reporter_mass_shift", 0.5);

		defaults_.setValue("channel_active_4plex", StringList::create("114:myReference"), "Four-plex only: Each channel that was used in the experiment and its description (114-117) in format <channel>:<name>, e.g. \"114:myref\",\"115:liver\"."); 
		defaults_.setValue("channel_active_8plex", StringList::create("113:myReference"), "Eight-plex only: Each channel that was used in the experiment and its description (113-121) in format <channel>:<name>, e.g. \"113:myref\",\"115:liver\",\"118:lung\"."); 

		StringList isotopes = ItraqConstants::getIsotopeMatrixAsStringList(itraq_type_, isotope_corrections_);
		defaults_.setValue("isotope_correction_values", isotopes, "override default values (see Documentation); use the following format: <channel>:<-2Da>/<-1Da>/<+1Da>/<+2Da> ; e.g. '114:0/0.3/4/0' , '116:0.1/0.3/3/0.2' ", StringList::create("advanced"));
		
    defaultsToParam_();
  }

  ITRAQLabeler::~ITRAQLabeler()
  {
  }

  void ITRAQLabeler::updateMembers_()
  {
		StringList channels_active;

		if (param_.getValue("iTRAQ") == "4plex") 
		{
			itraq_type_ = ItraqConstants::FOURPLEX;
			channels_active = param_.getValue("channel_active_4plex");
		}
		else if (param_.getValue("iTRAQ") == "8plex") 
		{
			itraq_type_ = ItraqConstants::EIGHTPLEX;
			channels_active = param_.getValue("channel_active_8plex");
		}

		ItraqConstants::initChannelMap(itraq_type_, channel_map_);
		ItraqConstants::updateChannelMap(channels_active, channel_map_);


		// update isotope_corrections_ Matrix with custom values
		StringList channels = param_.getValue("isotope_correction_values");
		if (channels.size()>0)
		{
			ItraqConstants::updateIsotopeMatrixFromStringList(itraq_type_, channels, isotope_corrections_);
		}

	}
  void ITRAQLabeler::preCheck(Param & param) const
  {
    // check for valid MS/MS method
    if(!StringList::create("disabled,precursor").contains( param.getValue("RawTandemSignal:status")))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "iTRAQ Labeling does not work with the chosen MS/MS type");
    }
  }

  void ITRAQLabeler::setUpHook(FeatureMapSimVector & features)
  {
    // no action here .. just check for correct # of channels
    Size active_channel_count=0;
    for (ChannelMapType::ConstIterator it=channel_map_.begin();it!=channel_map_.end();++it) 
    {
      if (it->second.active) ++active_channel_count;
    }
    if(features.size() != active_channel_count)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("iTRAQ Labeling received wrong number of channels: ") + String(active_channel_count) + " defined, but " + String(features.size()) + " given as FASTA files.");
    }
  }

  /// Labeling between digestion and rt simulation
  void ITRAQLabeler::postDigestHook(FeatureMapSimVector & features_to_simulate )
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
          b1.setIntensity(total_intensity * 2.0 * (1 - labeling_efficiency));

          final_feature_map.push_back(b1);

          // merge unlabeled with possible labeled feature
          // modify unlabeled intensity
          (*lf_iter).setIntensity(total_intensity * (1 - labeling_efficiency) * (1 - labeling_efficiency));

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

        }
        else
        {
          // generate labeled feature
          // labeling_efficiency is 100% so we transform the complete
          // feature in a dilabeled feature
          addModificationToPeptideHit_(*lf_iter, "UniMod:193");
          final_feature_map.push_back(*lf_iter);

          // add corresponding feature if it exists
          // and generate consensus feature for the unlabeled/labeled pair
          if(unlabeled_features_index.count(unmodified_sequence) != 0)
          {
            ConsensusFeature cf;
            final_feature_map.push_back(unlabeled_features_index[unmodified_sequence]);
            (*lf_iter).ensureUniqueId();
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
  }

  void ITRAQLabeler::mergeProteinAccessions_(Feature& target, const Feature& source) const
  {
    std::vector<String> target_acc (target.getPeptideIdentifications()[0].getHits()[0].getProteinAccessions());
    std::vector<String> source_acc (source.getPeptideIdentifications()[0].getHits()[0].getProteinAccessions());

    std::set<String> unique_acc;
    std::pair<set<String>::iterator, bool> result;

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

  Feature ITRAQLabeler::mergeFeatures_(Feature& labeled_channel_feature, const AASequence& unmodified_sequence, std::map<AASequence, Feature>& unlabeled_features_index) const
  {
    // merge with feature from first map (if it exists)
    if(unlabeled_features_index.count(unmodified_sequence) != 0)
    {
      // we only merge abundance and use feature from first map
      Feature new_f = unlabeled_features_index[unmodified_sequence];

      new_f.setMetaValue("channel_1_intensity", new_f.getIntensity());
      new_f.setMetaValue("channel_2_intensity", labeled_channel_feature.getIntensity());

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

  void ITRAQLabeler::addModificationToPeptideHit_(Feature& feature, const String& modification) const
  {
    vector<PeptideHit> pepHits(feature.getPeptideIdentifications()[0].getHits());
    AASequence modified_sequence(pepHits[0].getSequence());
    //modified_sequence.setModification(modified_sequence.size() - 1, modification);
    modified_sequence.setCTerminalModification(modification);
    pepHits[0].setSequence(modified_sequence);
    feature.getPeptideIdentifications()[0].setHits(pepHits);
  }


  Matrix<SimIntensityType> ITRAQLabeler::getItraqIntensity_(const Feature & f) const
  {
		StringList keys;
		f.getKeys(keys);

		// prepare map
		Map <Int, SimIntensityType> channel_intensities;
		std::vector< Matrix<Int> > channel_names(2);
		channel_names[0].setMatrix<4,1>(ItraqConstants::CHANNELS_FOURPLEX);
		channel_names[1].setMatrix<8,1>(ItraqConstants::CHANNELS_EIGHTPLEX);
		for (Int i=0; i<ItraqConstants::CHANNEL_COUNT[itraq_type_]; ++i)
		{
			channel_intensities[channel_names[itraq_type_].getValue(i,0)] = 0;
		}		
		
		// fill map with values present (all missing ones remain 0)
		for (StringList::const_iterator it_key = keys.begin(); it_key != keys.end(); it_key++)
		{
			if (!it_key->hasPrefix("intensity_itraq")) continue;
			Int ch = it_key->substr(String("intensity_itraq").length()).toInt();
			channel_intensities[ch] = f.getMetaValue(*it_key);
			std::cout << "raw itraq intensity: " << ch << "->" << f.getMetaValue(*it_key) << "\n";
		}	
		
		// fill the matrix
		Matrix<SimIntensityType> m(ItraqConstants::CHANNEL_COUNT[itraq_type_], 1, 0);
		Size index=0;
		for (Map <Int, SimIntensityType>::const_iterator it=channel_intensities.begin(); it!=channel_intensities.end(); ++it)
		{
			m.setValue(index ,0, it->second);
			++index;
		}

		return m;
		
  }

  /// Labeling between RT and Detectability
  void ITRAQLabeler::postRTHook(FeatureMapSimVector & /* features_to_simulate */)
  {
  }

  /// Labeling between Detectability and Ionization
  void ITRAQLabeler::postDetectabilityHook(FeatureMapSimVector & /* features_to_simulate */)
  {
  }

  /// Labeling between Ionization and RawMS
  void ITRAQLabeler::postIonizationHook(FeatureMapSimVector & /* features_to_simulate */)
  {
  }

  /// Labeling after RawMS
  void ITRAQLabeler::postRawMSHook(FeatureMapSimVector & /* features_to_simulate */)
  {
  }

  void ITRAQLabeler::postRawTandemMSHook(FeatureMapSimVector &, MSSimExperiment &)
  {
/*
			std::cout << "Matrix used: \n" << ItraqConstants::translateIsotopeMatrix(itraq_type_, isotope_corrections_) << "\n\n";
				
			gsl_matrix* channel_frequency = ItraqConstants::translateIsotopeMatrix(itraq_type_, isotope_corrections_).toGslMatrix();
			gsl_matrix* itraq_intensity_observed = Matrix<SimIntensityType>(ItraqConstants::CHANNEL_COUNT[itraq_type_],1).toGslMatrix();
			gsl_matrix* itraq_intensity_sum = Matrix<SimIntensityType>(ItraqConstants::CHANNEL_COUNT[itraq_type_],1).toGslMatrix();
			
			std::vector< Matrix<Int> > channel_names(2);
			channel_names[0].setMatrix<4,1>(ItraqConstants::CHANNELS_FOURPLEX);
			channel_names[1].setMatrix<8,1>(ItraqConstants::CHANNELS_EIGHTPLEX);

			// add signal...
			for (MSSimExperiment::iterator it=ms2.begin(); it!=ms2.end(); ++it)
			{
				// reset sum matrix to 0
				gsl_matrix_scale (itraq_intensity_sum, 0);
				
				// add up signal of all features
				// TODO: take care of actual position of feature relative to precursor!
				IntList parent_fs = (IntList) it->getMetaValue("parent_feature_ids");
				for (Size i_f=0; i_f < parent_fs.size(); ++i_f)
				{
					// apply isotope matrix to active channels
					gsl_matrix* row = getItraqIntensity_(features[i_f]).toGslMatrix();
					// row * channel_frequency = observed iTRAQ intensities
					gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
													1.0, channel_frequency, row,
													0.0, itraq_intensity_observed);
          // add result to sum
          gsl_matrix_add (itraq_intensity_sum, itraq_intensity_observed);
					gsl_matrix_free (row);
				}
				
				// add signal to MS2 spectrum
				for (Int i_channel=0; i_channel< ItraqConstants::CHANNEL_COUNT[itraq_type_]; ++i_channel)
				{
					MSSimExperiment::SpectrumType::PeakType p;
					// dummy
					p.setMZ(channel_names[itraq_type_].getValue(i_channel,0) + 0.1);
					p.setIntensity(gsl_matrix_get(itraq_intensity_sum, i_channel, 0));
					std::cout << "inserted iTRAQ peak: " << p << "\n";
					it->push_back(p);
				}
			}
			
			gsl_matrix_free (channel_frequency);
			gsl_matrix_free (itraq_intensity_observed);
			gsl_matrix_free (itraq_intensity_sum);

      */
  }
} // namespace OpenMS
