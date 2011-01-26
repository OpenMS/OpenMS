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
// $Maintainer: Chris Bielow$
// $Authors: Stephan Aiche, Chris Bielow$
// --------------------------------------------------------------------------

#include <OpenMS/SIMULATION/IonizationSimulation.h>

#include <OpenMS/DATASTRUCTURES/Compomer.h>
#include <OpenMS/CONCEPT/Constants.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS {


  IonizationSimulation::IonizationSimulation(const SimRandomNumberGenerator& random_generator)
    : DefaultParamHandler("IonizationSimulation"),
      ProgressLogger(),
      ionization_type_(),
      basic_residues_(),
      esi_probability_(),
      esi_impurity_probabilities_(),
      esi_adducts_(),
      max_adduct_charge_(),
      maldi_probabilities_(),
      rnd_gen_(&random_generator)
  {
    setDefaultParams_();
    updateMembers_();
  }

  IonizationSimulation::IonizationSimulation(const IonizationSimulation& source)
    : DefaultParamHandler(source),
      ProgressLogger(source),
      ionization_type_(source.ionization_type_ ),
      basic_residues_(source.basic_residues_ ),
      esi_probability_(source.esi_probability_ ),
      esi_impurity_probabilities_(source.esi_impurity_probabilities_),
      esi_adducts_(source.esi_adducts_ ),
      max_adduct_charge_(source.max_adduct_charge_ ),
      maldi_probabilities_(source.maldi_probabilities_),
      rnd_gen_(source.rnd_gen_)
  {
    //updateMembers_();  
  }

  
  IonizationSimulation& IonizationSimulation::operator = (const IonizationSimulation& source)
  {
		DefaultParamHandler::operator=(source);
		ionization_type_ = source.ionization_type_;
		basic_residues_ = source.basic_residues_;
		esi_probability_ = source.esi_probability_;
		esi_impurity_probabilities_ = source.esi_impurity_probabilities_;
		esi_adducts_ = source.esi_adducts_;
		max_adduct_charge_ = source.max_adduct_charge_;
		maldi_probabilities_ = source.maldi_probabilities_;
		rnd_gen_ = source.rnd_gen_;
    //updateMembers_();
    return *this;
  }
  
  IonizationSimulation::~IonizationSimulation()
  {}

  void IonizationSimulation::ionize(FeatureMapSim & features, ConsensusMap & charge_consensus, MSSimExperiment & experiment)
  {
    LOG_INFO << "Ionization Simulation ... started" << std::endl;

		// clear the consensus map
		charge_consensus = ConsensusMap();
		charge_consensus.setProteinIdentifications(features.getProteinIdentifications());

    switch (ionization_type_) 
		{
      case MALDI:
        ionizeMaldi_(features, charge_consensus);
        break;
      case ESI:
        ionizeEsi_(features, charge_consensus);
        break;
    }

		// add params for subsequent modules
		ScanWindow sw;
		sw.begin = minimal_mz_measurement_limit_;
		sw.end = maximal_mz_measurement_limit_;
		for (Size i=0;i<experiment.size();++i)
		{
			experiment[i].getInstrumentSettings().getScanWindows().push_back(sw);
		}

		ConsensusMap::FileDescription map_description;
    map_description.label = "Simulation (Charge Consensus)";
		map_description.size = features.size();
		charge_consensus.getFileDescriptions()[0] = map_description;
  }
  
  void IonizationSimulation::setDefaultParams_()
  {
    defaults_.setValue("ionization_type", "ESI", "Type of Ionization (MALDI or ESI)");
    defaults_.setValidStrings("ionization_type", StringList::create("MALDI,ESI"));

    defaults_.setValue("esi:ionized_residues", StringList::create("Arg,Lys,His"), "List of residues (as three letter code) that will be considered during ESI ionization. This parameter will be ignored during MALDI ionization.");
    StringList valid_ionized_residues = StringList::create("Ala,Cys,Asp,Glu,Phe,Gly,His,Ile,Lys,Leu,Met,Asn,Pro,Gln,Arg,Sec,Ser,Thr,Val,Trp,Tyr");
    defaults_.setValidStrings("esi:ionized_residues", valid_ionized_residues);
		defaults_.setValue("esi:charge_impurity", StringList::create("H+:1"), "List of charged ions that contribute to charge with weight of occurence (which must not sum to 1), e.g. ['H:1'] or ['H:0.7' 'Na:0.3']");

		defaults_.setValue("esi:max_impurity_set_size", 3, "Maximal #combinations of charge impurities allowed (each generating one feature) per charge state. E.g. assuming charge=3 and this parameter is 2, then we could choose to allow '3H+, 2H+Na+' features (given certain 'charge_impurity' constaints), but no '3H+, 2H+Na+, 3Na+'", StringList::create("advanced"));

    // ionization probabilities
    defaults_.setValue("esi:ionization_probability",0.8, "Probability for the binomial distribution of the ESI charge states");
    defaults_.setValue("maldi:ionization_probabilities", DoubleList::create("0.9,0.1") , "List of probabilities for the different charge states during MALDI ionization (the list must sum up to 1.0)");
    
    // maximal size of map in mz dimension
    defaults_.setValue("mz:upper_measurement_limit",2500.0,"Upper m/z detecter limit.");
    defaults_.setValue("mz:lower_measurement_limit",200.0,"Lower m/z detecter limit.");
    
    defaultsToParam_();
  }
  
  void IonizationSimulation::updateMembers_()
  {
    String type = param_.getValue("ionization_type");
    if(type == "ESI")
    {
      ionization_type_ = ESI;    
    }
    else if(type == "MALDI")
    {
      ionization_type_ = MALDI;
    }
    else
    {
      /// unsupported ionization model
      throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, "IonizationSimulation got invalid Ionization type '" + type + "'");
    }
    
    // get basic residues from params
    basic_residues_.clear();
    StringList basic_residues = (StringList) param_.getValue("esi:ionized_residues");
    for(StringList::const_iterator it = basic_residues.begin(); it != basic_residues.end(); ++it)
    {
      basic_residues_.insert(*it);
    }

		// parse possible ESI adducts
    StringList esi_charge_impurity = param_.getValue("esi:charge_impurity");
		StringList components;
		max_adduct_charge_ = 0;
		// reset internal state:
		esi_impurity_probabilities_.clear();
		esi_adducts_.clear();
    // cumulate probabilities in list
    for(Size i = 0 ; i < esi_charge_impurity.size() ; ++i )
    {
			esi_charge_impurity[i].split(':',components);
			if (components.size() != 2) throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, String("IonizationSimulation got invalid esi:charge_impurity (") + esi_charge_impurity[i] + ")with " + components.size() + " components instead of 2.");
			// determine charge of adduct (by # of '+')
			Size l_charge = components[0].size();
			l_charge -= components[0].remove('+').size();
			EmpiricalFormula ef(components[0].remove('+'));
			// effectively substract electrones
			ef.setCharge(l_charge); ef -= String("H")+String(l_charge);
			// create adduct
			Adduct a((Int)l_charge, 1, ef.getMonoWeight(), components[0].remove('+'), log(components[1].toDouble()),0);
			esi_adducts_.push_back(a);
			esi_impurity_probabilities_.push_back(components[1].toDouble());

			max_adduct_charge_ = std::max(max_adduct_charge_, l_charge);
    }

    // MALDI charge distribution
		maldi_probabilities_ = param_.getValue("maldi:ionization_probabilities");
    
    esi_probability_ = param_.getValue("esi:ionization_probability");
    
    // detector ranges
    maximal_mz_measurement_limit_ = param_.getValue("mz:upper_measurement_limit");
    minimal_mz_measurement_limit_ = param_.getValue("mz:lower_measurement_limit");
    
  }
  
  void IonizationSimulation::ionizeEsi_(FeatureMapSim & features, ConsensusMap & charge_consensus)
  {
		// we need to do this locally to avoid memory leaks (copying this stuff in C'tors is not wise)
		gsl_ran_discrete_t * gsl_ran_lookup_esi_charge_impurity = gsl_ran_discrete_preproc (esi_impurity_probabilities_.size(), &esi_impurity_probabilities_[0]);

		try
		{
			// map for charged features
			FeatureMapSim copy_map = features;
			// but leave meta information & other stuff intact
			copy_map.clear(false);

			// features which are not ionized
			Size uncharged_feature_count = 0;
			// features discarded - out of mz detection range
			Size undetected_features_count = 0;
			
      LOG_INFO << "Simulating " << features.size() << " features" << std::endl;

      this->startProgress(0,features.size(),"Ionization");
      Size progress(0);

			// iterate over all features
      #pragma omp parallel for reduction(+: uncharged_feature_count, undetected_features_count)
			for(SignedSize index = 0; index < (SignedSize)features.size(); ++index)
			{
        // no barrier here .. only an atomic update of progress value
        #pragma omp atomic
        ++progress;

#ifdef _OPENMP
        // progresslogger, only master thread sets progress (no barrier here)
        if (omp_get_thread_num() == 0) this->setProgress(progress);
#else
        this->setProgress(progress);
#endif

				ConsensusFeature cf;

				// iterate on abundance
				Int abundance = (Int) ceil( features[index].getIntensity() );
				UInt basic_residues_c = countIonizedResidues_(features[index].getPeptideIdentifications()[0].getHits()[0].getSequence());
	      
        if (basic_residues_c==0)
        {
					++uncharged_feature_count; // OMP
					//std::cout << "  not ionized: " << feature_it -> getPeptideIdentifications()[0].getHits()[0].getSequence().toUnmodifiedString() << "\n";
					continue;
        }

        // precompute random numbers:
        std::vector<UInt> prec_rndbin(abundance);
        #pragma omp critical (OPENMS_gsl)
        {
          for(Int j = 0; j < abundance ; ++j)
				  {
            prec_rndbin[j] = gsl_ran_binomial(rnd_gen_->technical_rng,esi_probability_,basic_residues_c);
          }
        }

        std::vector<Size> prec_rnduni(50); // uniform numbers container
        Size prec_rnduni_remaining(0);

				// assumption: each basic residue can hold one charged adduct
				Map<Compomer, UInt> charge_states;
        Size adduct_index;
        UInt charge;

				// sample different charge states (dice for each peptide molecule separately)
				for(Int j = 0; j < abundance ; ++j)
				{
					// currently we might also loose some molecules here (which is ok?)
					// sample charge state from binomial

          charge = prec_rndbin[j]; // get precomputed rnd

					if (charge==0)
					{
						continue;
					}

					Compomer cmp;
					// distribute charges across adduct types
					for (UInt charge_site=0;charge_site<charge;++charge_site)
					{
            if (prec_rnduni_remaining == 0)
            {
              #pragma omp critical (OPENMS_gsl)
              {
                for (Size i_rnd=0;i_rnd<prec_rnduni.size();++i_rnd)
                {
                  prec_rnduni[i_rnd] = gsl_ran_discrete (rnd_gen_->technical_rng, gsl_ran_lookup_esi_charge_impurity);
                }
                prec_rnduni_remaining = prec_rnduni.size();
              }
            }
            adduct_index = prec_rnduni[--prec_rnduni_remaining];
            cmp.add(esi_adducts_[adduct_index],Compomer::RIGHT);
					}

					// add 1 to abundance of sampled charge state
					++charge_states[ cmp ];
				}

        // no charges > 0 selected (this should be really rare)
				if (charge_states.size()==0) 
				{
					++uncharged_feature_count; // OMP!
					//std::cout << "  not ionized: " << feature_it -> getPeptideIdentifications()[0].getHits()[0].getSequence().toUnmodifiedString() << "\n";
					continue;
				}

				// transform into a set (for sorting by abundance)
				Int max_observed_charge(0);
				std::set< std::pair<UInt, Compomer > > charge_states_sorted;
				for (Map<Compomer, UInt>::const_iterator it_m=charge_states.begin(); it_m!=charge_states.end();++it_m)
				{ // create set of pair(value, key)
					charge_states_sorted.insert( std::make_pair(it_m->second,it_m->first) );
					// update maximal observed charge
					max_observed_charge = std::max(max_observed_charge, it_m->first.getNetCharge());
				}

				Int max_compomer_types = param_.getValue("esi:max_impurity_set_size");
				std::vector<Int> allowed_entities_of_charge(max_observed_charge+1, max_compomer_types);
				// start at highest abundant ions
				for(std::set< std::pair<UInt, Compomer > >::reverse_iterator it_s=charge_states_sorted.rbegin();
						it_s!=charge_states_sorted.rend();
						++it_s)
				{
					Int charge = it_s->second.getNetCharge();
					if (allowed_entities_of_charge[charge]>0)
					{
						Feature charged_feature(features[index]);            

            setFeatureProperties_(charged_feature, it_s->second.getMass(), it_s->second.getAdductsAsString(1), charge, it_s->first, index);
	
            // remember the original feature as parent feature (needed for labeling consensus)
            charged_feature.setMetaValue("parent_feature", String( features[index].getUniqueId() ));

						if (!isFeatureValid_(charged_feature))
						{
							++undetected_features_count; // OMP!
							continue;
						}

            #pragma omp critical (OPENMS_copy_map)
            {
              copy_map.push_back(charged_feature);
            }
						// add to consensus
						cf.insert(0, charged_feature);

						// decrease # of allowed compomers of current compomer's charge
						--allowed_entities_of_charge[charge];
					}
				}

				// add consensus element containing all charge variants just created
        #pragma omp critical (OPENMS_charge_consensus)
        {
          charge_consensus.push_back(cf);
        }

			} // ! for feature  (parallel)

      this->endProgress();	    

      for (Size i=0;i<charge_consensus.size(); ++i)
      { // this cannot be done inside the parallel-for as the copy_map might be populated meanwhile, which changes the internal uniqueid-map (used in below function)
        charge_consensus[i].computeDechargeConsensus(copy_map);
      }

			// swap feature maps
			features.swap(copy_map);

	    LOG_INFO << "#Peptides not ionized: " << uncharged_feature_count << "\n";
	    LOG_INFO << "#Peptides outside mz range: " << undetected_features_count << "\n";
		}
		catch (std::exception& e)
		{
			// before leaving: free
			gsl_ran_discrete_free (gsl_ran_lookup_esi_charge_impurity);
			throw e;
		}

		// all ok: free
		gsl_ran_discrete_free (gsl_ran_lookup_esi_charge_impurity);

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
  
  void IonizationSimulation::ionizeMaldi_(FeatureMapSim & features, ConsensusMap & charge_consensus)
  {
		gsl_ran_discrete_t * gsl_ran_lookup_maldi = gsl_ran_discrete_preproc (maldi_probabilities_.size(), &maldi_probabilities_[0]);

		try
		{
			// features discarded - out of mz detection range
			Size undetected_features_count = 0;
			Size feature_index = 0;
					
			FeatureMapSim copy_map(features);
      copy_map.clear(false);
      DoubleReal h_mono_weight = Constants::PROTON_MASS_U;

      this->startProgress(0,features.size(),"Ionization");
      Size progress=0;

			for(SignedSize index=0; index<(SignedSize)features.size();++index)
			{
				Int abundance = (Int) ceil( features[index].getIntensity() );
				std::vector<UInt> charge_states(((DoubleList) param_.getValue("maldi:ionization_probabilities")).size() + 1);
				// sample different charge states
				for(Int j = 0; j < abundance ; ++j)
				{
					// sample charge from discrete distribution
          Size charge = gsl_ran_discrete (rnd_gen_->technical_rng, gsl_ran_lookup_maldi) + 1;

					// add 1 to abundance of sampled charge state
					++charge_states[ charge ];
				}
	      
				ConsensusFeature cf;
				// only consider charged (charge >= 1) ions
				for(UInt c = 1 ; c < charge_states.size() ; ++c)
				{
					// empty charge states won't be generated
					if (charge_states[c] == 0) { continue; }
					else
					{
						Feature charged_feature(features[index]);
						
						setFeatureProperties_(charged_feature, h_mono_weight*c, String("H")+String(c), c, charge_states[c], feature_index);

            // remember the original feature as parent feature (needed for labeling consensus)
            charged_feature.setMetaValue("parent_feature", String( features[index].getUniqueId() ));

						if (!isFeatureValid_(charged_feature))
						{
							++undetected_features_count;
							continue;
						}
						
						copy_map.push_back(charged_feature);
						
						cf.insert(0, charged_feature);
					}
				}
				// add consensus element containing all charge variants just created
				cf.computeDechargeConsensus(copy_map);
				charge_consensus.push_back(cf);
				
        this->setProgress(progress);
				++feature_index;
			} // ! feature loop (parallel)
      
      this->endProgress();

			// swap feature maps
			features.swap(copy_map);
			
	    std::cout << "#Peptides outside mz range: " << undetected_features_count << "\n";
		}
		catch (std::exception& e)
		{
			// before leaving: free
			gsl_ran_discrete_free (gsl_ran_lookup_maldi);
			throw e;
		}
		// all ok: free
		gsl_ran_discrete_free (gsl_ran_lookup_maldi);
  }
  
	void IonizationSimulation::setFeatureProperties_(Feature & f, 
																									 const DoubleReal & adduct_mass, 
																									 const String & adduct_formula, 
																									 const SimChargeType charge,
																									 const SimIntensityType new_intensity,
																									 const Size parent_index)
	{
		EmpiricalFormula feature_ef = f.getPeptideIdentifications()[0].getHits()[0].getSequence().getFormula();

		f.setMZ( (feature_ef.getMonoWeight() + adduct_mass ) / charge);
		f.setCharge(charge);
    std::vector<PeptideHit> hits = f.getPeptideIdentifications()[0].getHits();
    hits[0].setCharge(charge);
    f.getPeptideIdentifications()[0].setHits(hits);
	  // set "main" intensity
    SimIntensityType old_intensity = f.getIntensity();
	  f.setIntensity(new_intensity);
	  DoubleReal factor = new_intensity / old_intensity;

    #pragma omp critical (OPENMS_setfeatureprop)
    {
      // ensure uniquenes
      f.setUniqueId();
		  // add meta information on compomer (mass)
		  f.setMetaValue("charge_adduct_mass", adduct_mass );
		  f.setMetaValue("charge_adducts", adduct_formula );
		  f.setMetaValue("parent_feature_number", parent_index );

		  // adapt "other" intensities (iTRAQ...) by the factor we just decreased real abundance
		  StringList keys;
		  f.getKeys(keys);
		  for (StringList::const_iterator it_key = keys.begin(); it_key != keys.end(); ++it_key)
		  {
			  if (it_key->hasPrefix("intensity"))
			  {
				  f.setMetaValue(*it_key, SimIntensityType(f.getMetaValue(*it_key)) * factor);
			  }
		  }
    } // ! pragma
	}

  bool IonizationSimulation::isFeatureValid_(const Feature & feature)
	{
		if (feature.getMZ() > maximal_mz_measurement_limit_ || feature.getMZ() < minimal_mz_measurement_limit_)
		{ // remove feature
			return false;
		}
		else
		{
			return true;
		}
	}

}

