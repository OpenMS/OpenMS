// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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

#include<OpenMS/SIMULATION/DigestSimulation.h>
#include<OpenMS/SIMULATION/DetectabilitySimulation.h>
#include <OpenMS/SIMULATION/RawMSSignalSimulation.h>
#include <OpenMS/SIMULATION/RawTandemMSSignalSimulation.h>
#include <OpenMS/SIMULATION/IonizationSimulation.h>
#include <OpenMS/SIMULATION/PTMSimulation.h>
#include <OpenMS/SIMULATION/RTSimulation.h>

namespace OpenMS {

  void verbosePrintFeatureMap(FeatureMapSim feature_map)
  {
#ifdef _DEBUG
    std::cout << "############## DEBUG -- FEATURE MAP ##############" << std::endl;

    std::cout << "contained proteins" << std::endl;
    ProteinIdentification protIdent = feature_map.getProteinIdentifications()[0];
    for(std::vector<ProteinHit>::iterator proteinHit = protIdent.getHits().begin();
        proteinHit != protIdent.getHits().end();
        ++proteinHit)
    {
      std::cout << "- " << proteinHit->getAccession() << std::endl;
    }
    std::cout << "----------------------------------------------" << std::endl;
    for(FeatureMapSim::const_iterator feat = feature_map.begin();
        feat != feature_map.end();
        ++feat)
    {
      std::cout << " RT: " << (*feat).getRT() << " MZ: " << (*feat).getMZ() << " INT: " << (*feat).getIntensity() << " CHARGE: " << (*feat).getCharge() << " Det: " << (*feat).getMetaValue("detectibility") << ::std::endl;
      std::cout << "derived from protein(s): ";
      for(std::vector<String>::const_iterator it = (*feat).getPeptideIdentifications()[0].getHits()[0].getProteinAccessions().begin();
          it != (*feat).getPeptideIdentifications()[0].getHits()[0].getProteinAccessions().end();
          ++it)
      {
        std::cout << (*it) << " ";
      }
      std::cout << std::endl << "----------------------------------------------" << std::endl;
    }

    std::cout << "############## END DEBUG -- FEATURE MAP ##############" << std::endl;
#endif
  }

  MSSim::MSSim()
    : DefaultParamHandler("MSSim"),
			experiment_(),
			features_(),
			consensus_map_()
  {
    setDefaultParams_();
  }

  MSSim::MSSim(const MSSim& source)
    : DefaultParamHandler(source),
			experiment_(source.experiment_),
			features_(source.features_),
			consensus_map_(source.consensus_map_)
  {
    setParameters( source.getParameters() );

    updateMembers_();
  }

  MSSim& MSSim::operator = (const MSSim& source)
  {
    setParameters( source.getParameters() );
    experiment_ = source.experiment_;
    features_ = source.features_;
		consensus_map_ = source.consensus_map_;
    updateMembers_();
    return *this;
  }

  MSSim::~MSSim()
  {}

  void MSSim::simulate(const gsl_rng* rnd_gen, const SampleProteins& proteins)
  {
    // TODO: add method to read contaminants
    // TODO: add method to select contaminants

    /*
      General progress should be
        1. Digest Proteins
        2. add Post Translational modifications
        3. Predict retention times
        4. predict detectibility
        5. simulate ionization
        6. simulate the (lc)ms signal -> TODO: integrate parameter for signal in lc direction
        7. select features for MS2
        8. generate MS2 signals for selected features
     */

		// re-distribute synced parameters:
		syncParams_(param_, false);
		//param_.store("c:/mssim_param.ini"); // test reconstruction

    // convert sample proteins into an empty FeatureMap with ProteinHits
    createFeatureMap_(proteins, features_);

		// digest
    DigestSimulation digest_sim;
    digest_sim.setParameters(param_.copy("Digestion:",true));
    digest_sim.digest(features_);

    // debug
    std::cout << "digested" << std::endl;
    verbosePrintFeatureMap(features_);

		// add PTM's
		PTMSimulation ptm_sim(rnd_gen);
    
    std::cout << param_.copy("PostTranslationalModifications:",true) << std::endl;
		
    ptm_sim.setParameters(param_.copy("PostTranslationalModifications:",true));
		ptm_sim.predictPTMs(features_);

    // debug
    std::cout << "ptms added" << std::endl;
    verbosePrintFeatureMap(features_);

		// RT prediction
		RTSimulation rt_sim(rnd_gen);
		rt_sim.setParameters(param_.copy("RTSimulation:",true));
		rt_sim.predictRT(features_, experiment_);

    // debug
    std::cout << "rt simulated" << std::endl;
    verbosePrintFeatureMap(features_);

		// Detectability prediction
		DetectabilitySimulation dt_sim;
		dt_sim.setParameters(param_.copy("PeptideDetectibilitySimulation:",true));
		dt_sim.filterDetectability(features_);

    // debug
    std::cout << "pd filtered" << std::endl;
    verbosePrintFeatureMap(features_);

    IonizationSimulation ion_sim(rnd_gen);
    ion_sim.setParameters(param_.copy("Ionization:", true));
    ion_sim.ionize(features_, consensus_map_);

    // debug
    std::cout << "ionized" << std::endl;
    verbosePrintFeatureMap(features_);

    RawMSSignalSimulation raw_sim(rnd_gen);
    raw_sim.setParameters(param_.copy("RawSignal:", true));
    raw_sim.generateRawSignals(features_, experiment_);

    // debug
    std::cout << "ionized" << std::endl;
    verbosePrintFeatureMap(features_);

    RawTandemMSSignalSimulation raw_tandemsim(rnd_gen);
    raw_tandemsim.setParameters(param_.copy("RawTandemSignal:", true));
    raw_tandemsim.generateRawTandemSignals(features_, experiment_);


  }

	void MSSim::createFeatureMap_(const SampleProteins& proteins, FeatureMapSim& feature_map)
	{
    // clear feature map
    feature_map.clear();
    ProteinIdentification protIdent;

		for (SampleProteins::const_iterator it=proteins.begin(); it!=proteins.end(); ++it)
		{
      std::cout << (it->first).identifier << " " << (it->first).sequence << " " << (it->second["intensity"]) << ::std::endl;
      // add new ProteinHit to ProteinIdentification
      ProteinHit protHit(0.0, 1, (it->first).identifier, (it->first).sequence);
      protHit.setMetaValue("description", it->first.description);
      // add intensity (global, iTRAQ,...) to Protein
			for (FASTAEntryEnhanced::const_iterator it_q = it->second.begin(); it_q!=it->second.end(); ++it_q)
			{
				protHit.setMetaValue(it_q->first, it_q->second);
			}
      protIdent.insertHit(protHit);

		}
    std::vector<ProteinIdentification> vec_protIdent;
    vec_protIdent.push_back(protIdent);
    feature_map.setProteinIdentifications(vec_protIdent);
	}

  void MSSim::setDefaultParams_()
  {
		// section params
    defaults_.insert("Digestion:", DigestSimulation().getDefaults());
    defaults_.insert("PostTranslationalModifications:",PTMSimulation(NULL).getDefaults());
    defaults_.insert("RTSimulation:",RTSimulation(NULL).getDefaults());
    defaults_.insert("PeptideDetectibilitySimulation:",DetectabilitySimulation().getDefaults());
    defaults_.insert("Ionization:",IonizationSimulation(NULL).getDefaults());
    defaults_.insert("RawSignal:",RawMSSignalSimulation(NULL).getDefaults());
		defaults_.insert("RawTandemSignal:",RawTandemMSSignalSimulation(NULL).getDefaults());

		//sync params (remove duplicates from modules and put them in a global module)
		syncParams_(defaults_, true);

    defaultsToParam_();
  }
  
  void MSSim::syncParams_(Param& p, bool to_outer)
  {
		std::vector<StringList> globals;
		// here the globals params are listed that require to be in sync across several modules
		// - first the global param name and following that the module names where this param occurs
		// - Warning: the module params must have unchanged names and restrictions! (descriptions can differ though)
		globals.push_back(StringList::create("iTRAQ,PostTranslationalModifications,RawTandemSignal:iTRAQ"));
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
    return features_;
  }

  ConsensusMap const & MSSim::getSimulatedConsensus() const
  {
		return consensus_map_;
  }
  
 

}
