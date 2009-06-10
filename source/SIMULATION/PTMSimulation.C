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
// $Maintainer: Chris Bielow$
// $Authors: Stephan Aiche, Chris Bielow$
// --------------------------------------------------------------------------

#include<OpenMS/SIMULATION/PTMSimulation.h>

#include <OpenMS/CHEMISTRY/ResidueDB.h>

namespace OpenMS {

  PTMSimulation::PTMSimulation(const gsl_rng* rnd_gen)
    : DefaultParamHandler("PTMSimulation"), rnd_gen_(rnd_gen)
  {
    setDefaultParams_();
	}

  PTMSimulation::PTMSimulation(const PTMSimulation& source)
    : DefaultParamHandler(source)
  {
    rnd_gen_ = source.rnd_gen_;
		//ptms_ will be done by updateMembers_
		updateMembers_();
  }

  PTMSimulation& PTMSimulation::operator = (const PTMSimulation& source)
  {
		if (this != &source)
		{
			DefaultParamHandler::operator=(source);
	    rnd_gen_ = source.rnd_gen_;
			//ptms_ will be done by updateMembers_
			updateMembers_();
		}
    return *this;
  }
	
	PTMSimulation::~PTMSimulation()
  {}
  
  void PTMSimulation::setDefaultParams_()
  {
		defaults_.setValue("modification_bound", 3, "no more modifications are added to a peptide when this number is reached. (set to 0 to disable this module)");
		defaults_.setMinInt("modification_bound",0);
		defaults_.setValue("potential_modifications", StringList::create("MOD:00071|0.13,MOD:00076|0.12,MOD:00078|0.31,MOD:00130|0.21"), "List of PSI::MOD modifications (each with probability of occurence)");
		
		defaults_.setValue("iTRAQ","off","add iTRAQ modifications?");
		defaults_.setValidStrings("iTRAQ",StringList::create("off,4plex,8plex"));
				
		defaultsToParam_();	  
  }
  
	void PTMSimulation::updateMembers_()
	{
		ptms_.clear();

		// test if modification exists; will throw exception otherwise
		StringList mods = (StringList) param_.getValue("potential_modifications");
		StringList mod_info;
		ProbabilityType p;

		// used to make sure identifier is unique
		std::set<String> unique_ids;

		for (StringList::const_iterator it=mods.begin();it!=mods.end();++it)
		{
			// two elements required
			if ( !(it->split('|',mod_info) && mod_info.size()==2) )
			{
				throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, "PTMSimulation got invalid parameter 'allowed_modifications' for entry: " + *it + " (invalid identifier)");
			}

			// just test if modification exists (will throw exception)
			const Residue* residue = ResidueDB::getInstance()->getModifiedResidue(mod_info[0].trim());

			// check if we've seen the identifier before
			if (!unique_ids.insert(mod_info[0]).second)
			{
				throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, "PTMSimulation got invalid parameter 'allowed_modifications' for entry: " + *it + " (double occurence)");
			}

			// try to convert the probability
			p = mod_info[1].toFloat();
			if (!(0<=p && p<=1))
			{
				throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, "PTMSimulation got invalid parameter 'allowed_modifications' for entry: " + *it + " (invalid probability)");
			}
			// store & update table
			ptms_[residue->getOneLetterCode()].add( PTM(*residue,p) );

		}

	}
 
  void PTMSimulation::predict_ptms(FeatureMapSim & map)
  {
		// TODO: handle fixed modifications (they need to show up in every feature)
  
		if (ptms_.size()==0) return;

		Size max_mod_count = (UInt) param_.getValue("modification_bound");
		if (max_mod_count == 0) return;

		FeatureMapSim map_ptm;

		// each peptide has 0-1 modified versions currently
		for (FeatureMapSim::iterator it_f=map.begin(); it_f!=map.end(); ++it_f)
		{
			// assume there is *one* peptide attached to current feature
			AASequence seq_aa = it_f->getPeptideIdentifications()[0].getHits()[0].getSequence();
			// get AA frequencies
			Map<String, Size> aa_table;
			seq_aa.getAAFrequencies(aa_table);
			String seq = seq_aa.toUnmodifiedString();

			
			// structures to randomly choose AA (only required due to 'modification_bound' param)
			Map<UInt,String> index2AA;
			UInt index=0;
			std::vector<UInt> AAcandidates;
			
			// draw from binomial for each type of AA 
			for (Map<String, Size>::const_iterator it_aa=aa_table.begin(); it_aa!=aa_table.end();++it_aa)
			{
				// can this type of AA be modified?
				if (!ptms_.has(it_aa->first)) continue;

				// dice how many are modified
				ptms_[it_aa->first].rnd_amount = gsl_ran_binomial (rnd_gen_, 1-ptms_[it_aa->first].probability_none, (UInt)aa_table[it_aa->first]);

				// map AA index to AA name
				index2AA[index] = it_aa->first;
				// add AA indizes for later shuffling
				AAcandidates.insert(AAcandidates.end(), ptms_[it_aa->first].rnd_amount, index);
				++index;
			}

			// pick candidates
			Size bound = std::min(max_mod_count, AAcandidates.size());

			if (bound==0) continue; // not even a single modification selected

			std::vector<UInt> AAcandidates_picked(bound);
			// choose #bound candidates (no replacement) and put them into AAcandidates_picked
			gsl_ran_choose (rnd_gen_, &(AAcandidates_picked[0]), bound, &(AAcandidates[0]), AAcandidates.size(), sizeof (UInt));

			// this is now the list of AA types to modify, e.g. [0,0,3,6]
			std::sort(AAcandidates_picked.begin(),AAcandidates_picked.end());
			
			// now add the modifications
			UInt last_aa_type=AAcandidates_picked[0];
			UInt current_aa_type;
			Size search_offset = 0;
			for (Size i_mods=0;i_mods<AAcandidates_picked.size();++i_mods)
			{
				current_aa_type = AAcandidates_picked[i_mods];
				String aa_name = index2AA[current_aa_type];
				// new AA --> start searching from start of peptide
				if (current_aa_type!=last_aa_type) search_offset=0;
				Size hit_pos = seq.find(aa_name, search_offset);
				
				// we know the AA is there because we counted it
				OPENMS_POSTCONDITION(hit_pos!=std::string::npos, "AA selected for modification could not be found! check code!");

				// modify AA
				seq_aa.setModification (hit_pos, ptms_[aa_name].draw(rnd_gen_).getModification() );

				search_offset = hit_pos+1;
			}

			// create modified peptide
			FeatureMapSim::FeatureType f = *it_f;
			std::vector< PeptideIdentification > pep_ids = f.getPeptideIdentifications();
			std::vector<PeptideHit> p_hits = pep_ids[0].getHits();
			p_hits[0].setSequence(seq_aa);
			pep_ids[0].setHits(p_hits);
			f.setPeptideIdentifications(pep_ids);
			map_ptm.push_back(f);

		} //! for each peptide

		// append PTM's to feature map
		map.insert(map.end(), map_ptm.begin(), map_ptm.end());

  }
}
