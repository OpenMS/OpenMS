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

#include <OpenMS/SIMULATION/IonizationSimulation.h>

namespace OpenMS {

  IonizationSimulation::IonizationSimulation()
    : DefaultParamHandler("IonizationSimulation")
  {
    setDefaultParams_();
    updateMembers_();    
  }

  IonizationSimulation::IonizationSimulation(const IonizationSimulation& source)
    : DefaultParamHandler(source)
  {
    setParameters( source.getParameters() );
    rnd_gen_ = source.rnd_gen_;
    updateMembers_();  
  }

  IonizationSimulation::IonizationSimulation(const gsl_rng * random_generator)
    : DefaultParamHandler("IonizationSimulation"), rnd_gen_(random_generator)
  {
    setDefaultParams_();
    updateMembers_();
  }
  
  IonizationSimulation& IonizationSimulation::operator = (const IonizationSimulation& source)
  {
    setParameters( source.getParameters() );
    rnd_gen_ = source.rnd_gen_;
    updateMembers_();
    return *this;
  }
  
  IonizationSimulation::~IonizationSimulation()
  {}

  void IonizationSimulation::ionize(FeatureMapSim & features)
  {
    switch (ionization_type) {
      case MALDI:
        ionize_maldi(features);
        break;
      case ESI:
        ionize_esi(features);
        break;
    }
  }
  
  void IonizationSimulation::setDefaultParams_()
  {
    defaults_.setValue("ionization_type", "ESI", "Type of Ionization (MALDI or ESI)");
    StringList valid_ionization_types = StringList::create("MALDI,ESI");
    defaults_.setValidStrings("ionization_type", valid_ionization_types);

    defaults_.setValue("esi:ionized_residues", StringList::create("Arg,Lys,His"), "List of residues (as three letter code) that will be considered during ESI ionization. This parameter will be ignores during MALDI ionization.");
    StringList valid_ionized_residues = StringList::create("Ala,Cys,Asp,Glu,Phe,Gly,His,Ile,Lys,Leu,Met,Asn,Pro,Gln,Arg,Sec,Ser,Thr,Val,Trp,Tyr");
    defaults_.setValidStrings("esi:ionized_residues", valid_ionized_residues);
    
    // ionization probabilities
    defaults_.setValue("esi:ionization_probability",0.8, "Probability for the binomial distribution of the ESI charge states");
    defaults_.setValue("maldi:ionization_probabilities", DoubleList::create("0.9,0.1") , "List of probabilities for the different charge states during MALDI ionization (the list must sum up to 1.0)");
    defaultsToParam_();
  }
  
  void IonizationSimulation::updateMembers_()
  {
    String type = param_.getValue("ionization_type");
    if(type == "ESI")
    {
      ionization_type = ESI;    
    }
    else if(type == "MALDI")
    {
      ionization_type = MALDI;
    }
    else
    {
      /// unsupported ionization model
      throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, "IonizationSimulation got invalid Ionization type" + type);
    }
    
    // get basic residues from params
    basic_residues_.clear();
    StringList basic_residues = (StringList) defaults_.getValue("esi:ionized_residues");
    for(StringList::const_iterator it = basic_residues.begin(); it != basic_residues.end(); ++it)
    {
      basic_residues_.insert(*it);
    }
    
    // get probs
    maldi_probabilities = param_.getValue("maldi:ionization_probabilities");
    // cumulate probabilities in list
    for(Size i = 1 ; i < maldi_probabilities.size() ; ++i )
    {
      maldi_probabilities[i] += maldi_probabilities[i-1];
    }
    
    esi_probability = param_.getValue("esi:ionization_probability");
  }
  
  void IonizationSimulation::ionize_esi(FeatureMapSim & features)
  {
    FeatureMap< > copyMap;

    // iterate over all features
    for(FeatureMapSim::iterator feature_it = features.begin();
        feature_it != features.end();
        ++feature_it)
    {
      // iterate on abundance
      Int abundance = ceil( (*feature_it).getIntensity() );
      UInt basic_residues_c = countIonizedResidues_((*feature_it).getPeptideIdentifications()[0].getHits()[0].getSequence());
      
      std::vector<UInt> charge_states(basic_residues_c + 1, 0);
      
      // sample different charge states
      for(Int j = 0; j < abundance ; ++j)
      {
        // TODO: currently we could also loose some molecules here
        // sample charge state from binomial
        UInt charge = gsl_ran_binomial(rnd_gen_,esi_probability,basic_residues_c);
        // add 1 to abundance of sampled charge state
        ++charge_states[ charge ];
      }

      // only consider charged (charge >= 1) ions
      for(UInt c = 1 ; c < charge_states.size() ; ++c)
      {
        // empty charge states won't be generated
        if(charge_states[c] == 0) { continue; }
        else
        {
          Feature chargedFeature((*feature_it));
          chargedFeature.setCharge(c);
          chargedFeature.setIntensity(charge_states[c]);
          copyMap.push_back(chargedFeature);
        }

      }
    }
    
    // swap feature maps
    features.swap(copyMap);
  }
  
  
  UInt IonizationSimulation::countIonizedResidues_(const AASequence& seq) const
  {
    UInt count = 0;
    for (Size i = 0; i<seq.size(); ++i)
    {
      // check for basic residues
      if (basic_residues_.count(seq[i].getShortName()) == 1)
      {
        ++count;
      }
    }
    
    return count;
  }
  
  void IonizationSimulation::ionize_maldi(FeatureMapSim & features)
  {
    FeatureMap< > copyMap;
    for(FeatureMap< >::iterator feature_it = features.begin();
        feature_it != features.end();
        ++feature_it)
    {
      Int abundance = ceil( (*feature_it).getIntensity() );
      std::vector<UInt> charge_states(maldi_probabilities.size() + 1);
      // sample different charge states
      for(Int j = 0; j < abundance ; ++j)
      {
        // sample charge from discrete distribution
        // TODO: maybe we should switch to gsl_ran_discrete .. but this needs preprocessing
        Real pr = gsl_rng_uniform(rnd_gen_);
        Size charge = 1;
        for(Size pi = 0 ; pi < maldi_probabilities.size() ; ++pi)
        {
          if(pr < maldi_probabilities[pi])
          {
            charge = pi + 1;
            break;
          }
        }

        // add 1 to abundance of sampled charge state
        ++charge_states[ charge ];
      }
      
      // only consider charged (charge >= 1) ions
      for(UInt c = 1 ; c < charge_states.size() ; ++c)
      {
        // empty charge states won't be generated
        if(charge_states[c] == 0) { continue; }
        else
        {
          Feature chargedFeature((*feature_it));
          chargedFeature.setCharge(c);
          chargedFeature.setIntensity(charge_states[c]);
          copyMap.push_back(chargedFeature);
        }
        
      }
    }
    
    // swap feature maps
    features.swap(copyMap);
  }
}

