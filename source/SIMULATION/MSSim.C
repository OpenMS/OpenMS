// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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

//#define OPENMS_DEBUG_SIM_

namespace OpenMS {

  void verbosePrintFeatureMap(FeatureMapSimVector feature_maps, String stage)
  {
#ifdef OPENMS_DEBUG_SIM_
    std::cout << "############## DEBUG (" << stage << ") -- FEATURE MAPS ##############" << std::endl;

    Size map_count = 1;
    StringList keys;
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
        (*feat).getKeys (keys);
        std::cout << " RT: " << (*feat).getRT()
                  << " MZ: " << (*feat).getMZ()
                  << " INT: " << (*feat).getIntensity()
                  << " CHARGE: " << (*feat).getCharge()
                  << " Det: " << (*feat).getMetaValue("detectibility") << std::endl
                  << " Pep: " << feat->getPeptideIdentifications()[0].getHits()[0].getSequence () .toString() << std::endl
                  << " ID: " << feat->getUniqueId() << std::endl
                  << " Meta: " << keys.concatenate(",") << std::endl;
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
		if (feature_maps.empty()) std::cout << stage; // just to avoid warnings of unused parameters
#endif
  }

  MSSim::MSSim()
    : DefaultParamHandler("MSSim"),
			experiment_(),
      feature_maps_(),
      consensus_map_(),
      labeler_(0)
  {
		// section params
    defaults_.insert("Digestion:", DigestSimulation().getDefaults());
    defaults_.insert("RT:",RTSimulation(SimRandomNumberGenerator()).getDefaults());
    defaults_.insert("Detectability:",DetectabilitySimulation().getDefaults());
    defaults_.insert("Ionization:",IonizationSimulation(SimRandomNumberGenerator()).getDefaults());
    defaults_.insert("RawSignal:",RawMSSignalSimulation(SimRandomNumberGenerator()).getDefaults());
    defaults_.insert("RawTandemSignal:",RawTandemMSSignalSimulation(SimRandomNumberGenerator()).getDefaults());

    subsections_.push_back("Labeling");

		//sync params (remove duplicates from modules and put them in a global module)
		syncParams_(defaults_, true);
    defaultsToParam_();
  }

  MSSim::~MSSim()
  {
    delete labeler_;
  }

  Param MSSim::getParameters() const
  {
    Param tmp;
    tmp.insert("", this->param_); // get non-labeling options

    std::vector<String> products = Factory<BaseLabeler>::registeredProducts();

    tmp.setValue("Labeling:type", "labelfree", "Select the labeling type you want for your experiment");
    tmp.setValidStrings("Labeling:type", products);

    for(std::vector<String>::iterator product_name = products.begin() ; product_name != products.end() ; ++product_name)
    {
      BaseLabeler* labeler = Factory<BaseLabeler>::create(*product_name);
      tmp.insert("Labeling:" + *product_name + ":", labeler->getDefaultParameters());
      delete(labeler);
    }

    return tmp;
  }

  void MSSim::simulate(const SimRandomNumberGenerator & rnd_gen, SampleChannels& channels)
  {
    /*todo: move to a global config file or into INI file */
    Log_fatal.setPrefix("%S: ");
    Log_error.setPrefix("%S: ");
    Log_warn.setPrefix("%S: ");
    Log_info.setPrefix("%S: ");
    Log_debug.setPrefix("%S: ");

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

    // re-distribute synced parameters:
    //param_.store("c:/mssim_param.ini"); // test reconstruction
    syncParams_(param_, false);
    
    // instanciate and pass params before doing any actual work
    // ... this way, each module can throw an Exception when the parameters
    // ... do not fit, and the users gets an immediate feedback
    DigestSimulation digest_sim;
    digest_sim.setParameters(param_.copy("Digestion:",true));
		RTSimulation rt_sim(rnd_gen);
		rt_sim.setParameters(param_.copy("RT:",true));
		DetectabilitySimulation dt_sim;
		dt_sim.setParameters(param_.copy("Detectability:",true));
    IonizationSimulation ion_sim(rnd_gen);
    ion_sim.setParameters(param_.copy("Ionization:", true));
    ion_sim.setLogType(this->getLogType());
    RawMSSignalSimulation raw_sim(rnd_gen);
    raw_sim.setParameters(param_.copy("RawSignal:", true));
    raw_sim.setLogType(this->getLogType());
    raw_sim.loadContaminants(); // check if the file is valid (if not, an error is raised here instead of half-way through simulation)


    String labeling = param_.getValue("Labeling:type");
    labeler_ = Factory<BaseLabeler>::create(labeling);
    Param labeling_parameters = param_.copy("Labeling:" + labeling + ":",true);
    labeler_->setParameters(labeling_parameters);
    labeler_->setRnd(rnd_gen);

    // check parameters ..
    labeler_->preCheck(param_);

    // convert sample proteins into an empty FeatureMap with ProteinHits
    for(SampleChannels::const_iterator channel_iterator = channels.begin() ; channel_iterator != channels.end() ; ++channel_iterator)
    {
      FeatureMapSim map;
      createFeatureMap_(*channel_iterator, map, feature_maps_.size());
      feature_maps_.push_back(map);
    }

    // Call setUpHook
    labeler_->setUpHook(feature_maps_);

		// digest
    for(FeatureMapSimVector::iterator map_iterator = feature_maps_.begin() ; map_iterator != feature_maps_.end() ; ++map_iterator)
    {
      digest_sim.digest(*map_iterator);
    }

    // post digest labeling
    labeler_->postDigestHook(feature_maps_);

    // debug
    verbosePrintFeatureMap(feature_maps_, "digested");

		// RT prediction
    for(FeatureMapSimVector::iterator map_iterator = feature_maps_.begin() ; map_iterator != feature_maps_.end() ; ++map_iterator)
    {
      rt_sim.predictRT(*map_iterator);
    }
    rt_sim.createExperiment(experiment_);

    peak_map_ = experiment_; // initial Ground Truth for peak map is the same as for raw data

    // post rt sim labeling
    labeler_->postRTHook(feature_maps_);

    // debug
    verbosePrintFeatureMap(feature_maps_, "RT sim done");

		// Detectability prediction
    for(FeatureMapSimVector::iterator map_iterator = feature_maps_.begin() ; map_iterator != feature_maps_.end() ; ++map_iterator)
    {
      dt_sim.filterDetectability(*map_iterator);
    }

    // post detectability labeling
    labeler_->postDetectabilityHook(feature_maps_);

    // debug
    verbosePrintFeatureMap(feature_maps_, "DT sim done");

    // at this point all feature maps should be combined to one?
    ion_sim.ionize(feature_maps_.front(), consensus_map_, experiment_);

    // post ionization labeling
    labeler_->postIonizationHook(feature_maps_);

    // debug
    verbosePrintFeatureMap(feature_maps_, "ION sim done");

    raw_sim.generateRawSignals(feature_maps_.front(), experiment_, peak_map_, contaminants_map_);

    // post raw sim labeling
    labeler_->postRawMSHook(feature_maps_);

    // debug
    verbosePrintFeatureMap(feature_maps_, "RawSignal sim done");

    RawTandemMSSignalSimulation raw_tandemsim(rnd_gen);
    raw_tandemsim.setParameters(param_.copy("RawTandemSignal:", true));
    raw_tandemsim.generateRawTandemSignals(feature_maps_.front(), experiment_);

    labeler_->postRawTandemMSHook(feature_maps_, experiment_);

    
    // some last fixing of meta-values (this is impossible to do before as we do not know the final number of scans)
    for (Size i = 0 ; i < feature_maps_[0].size() ; ++i)
    {
      Feature& f = feature_maps_[0][i];
      PeptideIdentification& pi = f.getPeptideIdentifications()[0];
      // search for closest scan index:
      MSSimExperiment::ConstIterator it_rt = experiment_.RTBegin(f.getRT());
      SignedSize scan_index = std::distance<MSSimExperiment::ConstIterator> (experiment_.begin(), it_rt);
      pi.setMetaValue("RT_index", scan_index);
      pi.setMetaValue("RT", f.getRT());
    }

    LOG_INFO << "Final number of simulated features: " << feature_maps_[0].size() << "\n";

    // reindex spectra to avoid naming conflicts
    Size id = 1;
    experiment_.sortSpectra();
    for(MSSimExperiment::Iterator spectrum_iterator = experiment_.begin() ; spectrum_iterator != experiment_.end() ; ++spectrum_iterator)
    {
      MSSimExperiment::SpectrumType& spectrum = *spectrum_iterator;
      String spec_id = String("scan=") + id++;
      spectrum.setNativeID(spec_id);
    }
  }

  void MSSim::createFeatureMap_(const SampleProteins& proteins, FeatureMapSim& feature_map, Size map_index)
	{
    // clear feature map
    feature_map.clear(true);
    ProteinIdentification protIdent;

		for (SampleProteins::const_iterator it=proteins.begin(); it!=proteins.end(); ++it)
		{
      // add new ProteinHit to ProteinIdentification
      ProteinHit protHit(0.0, 1, (it->first).identifier, (it->first).sequence);
      // copy all meta values from FASTA file parsing
      protHit=(it->second);
      // additional meta values:
      protHit.setMetaValue("description", it->first.description);
      protHit.setMetaValue("map_index", map_index);
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
		globals.push_back(StringList::create("ionization_type,Ionization,RawSignal,RawTandemSignal"));
		
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

  ConsensusMap & MSSim::getChargeConsensus()
  {
		return consensus_map_;
  }

  ConsensusMap & MSSim::getLabelingConsensus()
  {
    return labeler_->getConsensus();
  }
  
  FeatureMapSim const & MSSim::getContaminants() const
  {
		return contaminants_map_;
  }
 
  MSSimExperiment const & MSSim::getPeakMap() const
  {
    return peak_map_;
  }


}
