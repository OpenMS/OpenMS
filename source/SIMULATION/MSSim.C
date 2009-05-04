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
#include <OpenMS/SIMULATION/RawSignalSimulation.h>
#include <OpenMS/SIMULATION/IonizationSimulation.h>
#include <OpenMS/SIMULATION/PTMSimulation.h>
#include <OpenMS/SIMULATION/RTSimulation.h>

#define _VFeature(feat) std::cout << __LINE__ << " RT: " << (feat).getRT() << " MZ: " << (feat).getMZ() << " INT: " << (feat).getIntensity() << " CHARGE: " << (feat).getCharge() << " Det: " << (feat).getMetaValue("detectibility") << ::std::endl; 

namespace OpenMS {

  MSSim::MSSim()
    : DefaultParamHandler("MSSim"), experiment_(), features_()
  {
    setDefaultParams_();
  }

  MSSim::MSSim(const MSSim& source)
    : DefaultParamHandler(source)
  {
    setParameters( source.getParameters() );
    experiment_ = source.experiment_;
    features_ = source.features_;
    updateMembers_();  
  }

  MSSim& MSSim::operator = (const MSSim& source)
  {
    setParameters( source.getParameters() );
    experiment_ = source.experiment_;
    features_ = source.features_;
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
    
		// digest
    DigestSimulation digest_sim;
    digest_sim.setParameters(param_.copy("Digestion:",true)); 
    
    // read proteins from protein file  
		SamplePeptides peptides;
		digest_sim.digest(proteins, peptides);
		
    // debug
    for (SamplePeptides::const_iterator peptide = peptides.begin();
         peptide != peptides.end();
         ++peptide)
    {
      std::cout << "Generated " << (*peptide).first << " with abundance " << (*peptide).second << " via tryptic digestion" << std::endl; 
    }
      
    
		// convert
    createFeatureMap_(peptides);

    // debug
		for(FeatureMapSim::const_iterator feature = features_.begin();
        feature != features_.end();
        ++feature)
    {
      _VFeature(*feature)
    }
    
		// add PTM's
		PTMSimulation ptm_sim(rnd_gen);
		ptm_sim.setParameters(param_.copy("PostTranslationalModifications:",true));
		ptm_sim.predict_ptms(features_);

		// RT prediction
		RTSimulation rt_sim(rnd_gen);
		rt_sim.setParameters(param_.copy("RTSimulation:",true));
		rt_sim.predict_rt(features_);
       
		// Detectability prediction
		DetectabilitySimulation dt_sim;
		dt_sim.setParameters(param_.copy("PeptideDetectibilitySimulation:",true));
		dt_sim.filterDetectability(features_);
    
    IonizationSimulation ion_sim(rnd_gen);
    ion_sim.setParameters(param_.copy("Ionization:", true));
    ion_sim.ionize(features_);
    
    // debug
		for(FeatureMapSim::const_iterator feature = features_.begin();
        feature != features_.end();
        ++feature)
    {
      _VFeature(*feature)
    }
    
    RawSignalSimulation raw_sim(rnd_gen);
    raw_sim.setParameters(param_.copy("RawSignal:", true));
    createExperiment_(rt_sim.getGradientTime(), raw_sim.getRTSamplingRate());
    
    raw_sim.generateRawSignals(features_, experiment_);

/**
			...
        7. select features for MS2
        8. generate MS2 signals for selected features
**/

    
    
  }

	void MSSim::createFeatureMap_(const SamplePeptides& peptides)
	{
    features_.clear();
		features_.reserve(peptides.size());

		for (SamplePeptides::const_iterator it=peptides.begin(); it!=peptides.end(); ++it)
		{
			Feature f;
			PeptideIdentification pep_id;
			pep_id.insertHit(PeptideHit(1.0, 1, 1, it->first));
			f.getPeptideIdentifications().push_back(pep_id);
			f.setIntensity(it->second);
			features_.push_back(f);
		}
	}

  void MSSim::createExperiment_(const DoubleReal gradient_time, const DoubleReal rt_sampling_rate)
  {
    experiment_.clear();
    Size number_of_scans = static_cast<Size>(gradient_time / rt_sampling_rate);
    experiment_.resize(number_of_scans);
    
    DoubleReal current_scan_rt = rt_sampling_rate;
    for(MSSimExperiment::iterator exp_it = experiment_.begin();
        exp_it != experiment_.end();
        ++exp_it)
    {
      // TODO: maybe we should also apply an error here like Ole did it in the original MapSimulator
      // double n = gsl_ran_gaussian(rand_gen_, 0.05);
      (*exp_it).setRT(current_scan_rt);
      current_scan_rt += rt_sampling_rate;
    }
  }
  
  void MSSim::setDefaultParams_()
  {   
    defaults_.insert("Digestion:", DigestSimulation().getDefaults());  
    defaults_.insert("PostTranslationalModifications:",PTMSimulation(NULL).getDefaults());
    defaults_.insert("RTSimulation:",RTSimulation(NULL).getDefaults());
    defaults_.insert("PeptideDetectibilitySimulation:",DetectabilitySimulation().getDefaults());
    defaults_.insert("Ionization:",IonizationSimulation(NULL).getDefaults());
    defaults_.insert("RawSignal:",RawSignalSimulation(NULL).getDefaults());
    
    defaultsToParam_();  
  }
  
  void MSSim::updateMembers_()
  {}

  MSSimExperiment const & MSSim::getExperiment() const
  {
    return experiment_;
  }
  
  FeatureMapSim const & MSSim::getSimulatedFeatures() const
  {
    return features_;
  }
  
}
