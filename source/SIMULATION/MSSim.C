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
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche, Chris Bielow$
// --------------------------------------------------------------------------

#include<OpenMS/SIMULATION/MSSim.h>

#include <OpenMS/SIMULATION/DigestSimulation.h>
#include <OpenMS/SIMULATION/DetectabilitySimulation.h>
#include <OpenMS/SIMULATION/RawMSSignalSimulation.h>
#include <OpenMS/SIMULATION/RawTandemMSSignalSimulation.h>
#include <OpenMS/SIMULATION/IonizationSimulation.h>
#include <OpenMS/SIMULATION/RTSimulation.h>

#include <OpenMS/SIMULATION/LABELING/BaseLabeler.h>

#define _DEBUG

namespace OpenMS {

  void verbosePrintFeatureMap(FeatureMapSimVector feature_maps, String stage)
  {
#ifdef _DEBUG
    std::cout << "############## DEBUG (" << stage << ") -- FEATURE MAPS ##############" << std::endl;

    Size map_count = 1;
    for(FeatureMapSimVector::iterator map_iter = feature_maps.begin() ; map_iter != feature_maps.end() ; ++map_iter)
    {
      std::cout << "FEATURE MAP #" << map_count << std::endl;

      std::cout << "contained proteins" << std::endl;
      ProteinIdentification protIdent = (*map_iter).getProteinIdentifications()[0];
      for(std::vector<ProteinHit>::iterator proteinHit = protIdent.getHits().begin();
          proteinHit != protIdent.getHits().end();
          ++proteinHit)
      {
        std::cout << "- " << proteinHit->getAccession() << std::endl;
      }
      std::cout << "----------------------------------------------" << std::endl;
      for(FeatureMapSim::const_iterator feat = (*map_iter).begin();
          feat != (*map_iter).end();
          ++feat)
      {
        std::cout << " RT: " << (*feat).getRT()
                  << " MZ: " << (*feat).getMZ()
                  << " INT: " << (*feat).getIntensity()
                  << " CHARGE: " << (*feat).getCharge()
                  << " Det: " << (*feat).getMetaValue("detectibility")
                  << " Pep: " << feat->getPeptideIdentifications()[0].getHits()[0].getSequence () .toString() << ::std::endl;
        std::cout << "derived from protein(s): ";
        for(std::vector<String>::const_iterator it = (*feat).getPeptideIdentifications()[0].getHits()[0].getProteinAccessions().begin();
            it != (*feat).getPeptideIdentifications()[0].getHits()[0].getProteinAccessions().end();
            ++it)
        {
          std::cout << (*it) << " ";
        }
        std::cout << std::endl << "----------------------------------------------" << std::endl;
      }
      std::cout << std::endl;
      ++map_count;
    }
    std::cout << "############## END DEBUG -- FEATURE MAPS ##############" << std::endl;
#else
		if (feature_map.size()==0) std::cout << stage; // just to avoid warnings of unused parameters
#endif
  }

  MSSim::MSSim()
    : DefaultParamHandler("MSSim"),
			experiment_(),
      feature_maps_(),
			consensus_map_()
  {
		// section params
    defaults_.insert("Digestion:", DigestSimulation().getDefaults());
    defaults_.insert("RTSimulation:",RTSimulation(NULL).getDefaults());
    defaults_.insert("PeptideDetectabilitySimulation:",DetectabilitySimulation().getDefaults());
    defaults_.insert("Ionization:",IonizationSimulation(NULL).getDefaults());
    defaults_.insert("RawSignal:",RawMSSignalSimulation(NULL).getDefaults());
		defaults_.insert("RawTandemSignal:",RawTandemMSSignalSimulation(NULL).getDefaults());

    subsections_.push_back("Labeling");

		//sync params (remove duplicates from modules and put them in a global module)
		syncParams_(defaults_, true);
    defaultsToParam_();
  }

  MSSim::~MSSim()
  {}

  Param MSSim::getParameters(const String &labeling_name) const
  {
    Param tmp;
    tmp.insert("", this->param_); // get non-labeling options

    if(labeling_name != "")
    {
      BaseLabeler* labeler = Factory<BaseLabeler>::create(labeling_name);
      tmp.insert("Labeling:", labeler->getDefaultParameters());
      delete(labeler);
    }
    return tmp;
  }

  void MSSim::simulate(gsl_rng* const rnd_gen, SampleChannels& channels, const String &labeling_name)
  {
    // TODO: add method to read contaminants
    // TODO: add method to select contaminants

    /*
      General progress should be
        1. Digest Proteins
        2. Predict retention times
        3. predict detectibility
        4. simulate ionization
        5. simulate the ms signal
        6. select features for MS2
        7. generate MS2 signals for selected features
     */

    BaseLabeler* labeler = Factory<BaseLabeler>::create(labeling_name);
    Param labeling_parameters = param_.copy("Labeling:",true);
    labeler->setParameters(labeling_parameters);
    labeler->setRnd(rnd_gen);

    // check parameters ..
    labeler->preCheck(param_);

		// re-distribute synced parameters:
		syncParams_(param_, false);
		//param_.store("c:/mssim_param.ini"); // test reconstruction

    // convert sample proteins into an empty FeatureMap with ProteinHits
    for(SampleChannels::const_iterator channel_iterator = channels.begin() ; channel_iterator != channels.end() ; ++channel_iterator)
    {
      FeatureMapSim map;
      createFeatureMap_(*channel_iterator, map);
      feature_maps_.push_back(map);
    }

    // Call setUpHook
    labeler->setUpHook(feature_maps_);

		// digest
    DigestSimulation digest_sim;
    digest_sim.setParameters(param_.copy("Digestion:",true));
    for(FeatureMapSimVector::iterator map_iterator = feature_maps_.begin() ; map_iterator != feature_maps_.end() ; ++map_iterator)
    {
      digest_sim.digest(*map_iterator);
    }

    // post digest labeling
    labeler->postDigestHook(feature_maps_);

    // debug
    verbosePrintFeatureMap(feature_maps_, "digested");

		// RT prediction
		RTSimulation rt_sim(rnd_gen);
		rt_sim.setParameters(param_.copy("RTSimulation:",true));
    for(FeatureMapSimVector::iterator map_iterator = feature_maps_.begin() ; map_iterator != feature_maps_.end() ; ++map_iterator)
    {
      rt_sim.predictRT(*map_iterator);
    }
    rt_sim.createExperiment(experiment_);

    // post rt sim labeling
    labeler->postRTHook(feature_maps_);

    // debug
    verbosePrintFeatureMap(feature_maps_, "RT sim done");

		// Detectability prediction
		DetectabilitySimulation dt_sim;
		dt_sim.setParameters(param_.copy("PeptideDetectabilitySimulation:",true));
    for(FeatureMapSimVector::iterator map_iterator = feature_maps_.begin() ; map_iterator != feature_maps_.end() ; ++map_iterator)
    {
      dt_sim.filterDetectability(*map_iterator);
    }

    // post detectability labeling
    labeler->postDetectabilityHook(feature_maps_);

    // debug
    verbosePrintFeatureMap(feature_maps_, "DT sim done");

    // at this point all feature maps should be combined to one?
    IonizationSimulation ion_sim(rnd_gen);
    ion_sim.setParameters(param_.copy("Ionization:", true));
    ion_sim.ionize(feature_maps_.front(), consensus_map_, experiment_);

    // post ionization labeling
    labeler->postIonizationHook(feature_maps_);

    // debug
    verbosePrintFeatureMap(feature_maps_, "ION sim done");

    RawMSSignalSimulation raw_sim(rnd_gen);
    raw_sim.setParameters(param_.copy("RawSignal:", true));
    raw_sim.generateRawSignals(feature_maps_.front(), experiment_);

    // post raw sim labeling
    labeler->postRawMSHook(feature_maps_);

    // debug
    verbosePrintFeatureMap(feature_maps_, "RawSignal sim done");

    RawTandemMSSignalSimulation raw_tandemsim(rnd_gen);
    raw_tandemsim.setParameters(param_.copy("RawTandemSignal:", true));
    raw_tandemsim.generateRawTandemSignals(feature_maps_.front(), experiment_);

    labeler->postRawTandemMSHook(feature_maps_,experiment_);

    LOG_INFO << "Final number of simulated features: " << feature_maps_[0].size() << "\n";

  }

	void MSSim::createFeatureMap_(const SampleProteins& proteins, FeatureMapSim& feature_map)
	{
    // clear feature map
    feature_map.clear(true);
    ProteinIdentification protIdent;

		for (SampleProteins::const_iterator it=proteins.begin(); it!=proteins.end(); ++it)
		{
      std::cout << (it->first).identifier << " " << (it->first).sequence << " " << (it->second) << ::std::endl;
      // add new ProteinHit to ProteinIdentification
      ProteinHit protHit(0.0, 1, (it->first).identifier, (it->first).sequence);
      protHit.setMetaValue("description", it->first.description);
      // add intensity
      protHit.setMetaValue("intensity", it->second);
      protIdent.insertHit(protHit);

		}
    std::vector<ProteinIdentification> vec_protIdent;
    vec_protIdent.push_back(protIdent);
    feature_map.setProteinIdentifications(vec_protIdent);
	}

  void MSSim::syncParams_(Param& p, bool to_outer)
  {
		std::vector<StringList> globals;
		// here the globals params are listed that require to be in sync across several modules
		// - first the global param name and following that the module names where this param occurs
		// - Warning: the module params must have unchanged names and restrictions! (descriptions can differ though)
		globals.push_back(StringList::create("ionization_type,Ionization,RawTandemSignal"));
		
		String global_prefix = "Global";
		// remove or add local params
		if (to_outer)
		{	// remove local params and merge to global
			for (Size i = 0; i < globals.size(); ++i)
			{
				// set the global param:
				OPENMS_PRECONDITION(globals[i].size()>=2, "Param synchronisation aborting due to missing local parameters!");
				p.insert(global_prefix+":"+globals[i][0], p.copy(globals[i][1]+":"+globals[i][0],true));
				// remove local params
				for (Size i_local=1;i_local<globals[i].size(); ++i_local)
				{
					p.remove(globals[i][i_local]+":"+globals[i][0]);
				}
			}
		}
		else // restore local params from global one
		{
			for (Size i = 0; i < globals.size(); ++i)
			{
				// get the global param:
				OPENMS_PRECONDITION(globals[i].size()>=2, "Param synchronisation aborting due to missing local parameters!");

				Param p_global = p.copy(global_prefix + ":" + globals[i][0],true);
				// insert into local params
				for (Size i_local=1;i_local<globals[i].size(); ++i_local)
				{
					p.insert(globals[i][i_local]+":"+globals[i][0], p_global);
				}
			}		
		}
		
  }

  void MSSim::updateMembers_()
  {
  }

  MSSimExperiment const & MSSim::getExperiment() const
  {
    return experiment_;
  }

  FeatureMapSim const & MSSim::getSimulatedFeatures() const
  {
    OPENMS_PRECONDITION(feature_maps_.size()==1, "More than one feature map remains after simulation. The channels should however be merged by now. Check!")
    return feature_maps_[0];
  }

  ConsensusMap const & MSSim::getSimulatedConsensus() const
  {
		return consensus_map_;
  }
  
 

}
