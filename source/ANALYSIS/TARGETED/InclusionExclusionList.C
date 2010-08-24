// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Alexandra Zerck $
// $Authors: Alexandra Zerck $
// --------------------------------------------------------------------------
//

#include <OpenMS/ANALYSIS/TARGETED/InclusionExclusionList.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/SIMULATION/RTSimulation.h>


#include <fstream>

namespace OpenMS
{
  InclusionExclusionList::InclusionExclusionList()
  {

  }
    

//   void InclusionExclusionList::loadTargets(FeatureMap<>& map, std::vector<IncludeExcludeTarget>& targets,TargetedExperiment& exp)
//   {

//   }

//   void InclusionExclusionList::loadTargets(std::vector<FASTAFile::FASTAEntry>& fasta_entries, std::vector<IncludeExcludeTarget>& targets,TargetedExperiment& exp, Size missed_cleavages)
//   {

//   }

	void InclusionExclusionList::writeTargets(std::vector<FASTAFile::FASTAEntry>& fasta_entries, String& out_path,IntList& charges,String rt_model_path,
																						DoubleReal rel_rt_window_size,bool rt_in_seconds,Size missed_cleavages)
	{
		EnzymaticDigestion digest;
		digest.setMissedCleavages(missed_cleavages);
		std::vector<FASTAFile::FASTAEntry>::iterator entry_iter = fasta_entries.begin();
		std::ofstream outs(out_path.c_str());
		outs.precision(8);
		if(!outs)
		{
			throw Exception::UnableToCreateFile(__FILE__, __LINE__,__PRETTY_FUNCTION__,"Cannot open output file.");
		}

    SimRandomNumberGenerator rnd_gen;
    RTSimulation rt_sim(rnd_gen);
		Param rt_param;
		rt_param.setValue("HPLC:model_file",rt_model_path);
		rt_sim.setParameters(rt_param);
		std::vector<AASequence> pep_seqs;
		for(;entry_iter != fasta_entries.end();++entry_iter)
		{
			// digest sequence
			AASequence aa_seq(entry_iter->sequence);
			std::vector<AASequence> vec;
			digest.digest(aa_seq,vec);

			// copy 
			pep_seqs.insert(pep_seqs.begin(),vec.begin(),vec.end());
		
			// TODO: enter modifications

					// // enter mod
					// if(fixed_mods_)
					// 	{
// 								// go through peptide sequence and check if AA is modified
// 								for(Size aa = 0; aa < vec_iter->size();++aa)
// 									{
// 										if(fixed_modifications_.find((vec_iter->toUnmodifiedString())[aa])!= fixed_modifications_.end())
// 											{
// #ifdef DEBUG_PISP
// 												std::cout << "w/o Mod "<<*vec_iter<<" "
// 																	<<vec_iter->getMonoWeight(Residue::Full,1)<<std::endl;
// #endif
// 												std::vector<String> & mods = fixed_modifications_[(vec_iter->toUnmodifiedString())[aa]];
// 												for(Size m = 0; m < mods.size();++m)
// 													{
// 														vec_iter->setModification(aa,mods[m]);
// 													}
// #ifdef DEBUG_PISP														
// 												std::cout << "set Mods "<<*vec_iter<<" "
// 																	<<vec_iter->getMonoWeight(Residue::Full,1)<<std::endl;
// #endif
// 											}
// 									}
//							}
					
		}
		std::vector<DoubleReal> rts;
		rt_sim.wrapSVM(pep_seqs,rts);
		for(Size i = 0; i < pep_seqs.size();++i)
		{
			for(Size c = 0; c < charges.size();++c)
			{
				// calculate m/z
				DoubleReal mz = pep_seqs[i].getMonoWeight(Residue::Full,charges[c])/(DoubleReal)charges[c];
				DoubleReal rt_start,rt_stop;
				if(rt_in_seconds)
				{
					rt_start = (rts[i] - rel_rt_window_size * rts[i]); // RT in minutes
					if(rt_start < 0.) rt_start = 0.;
					rt_stop = (rts[i] + rel_rt_window_size * rts[i]) ; // RT in minutes
				}
				else
				{
					rt_start = (rts[i] - rel_rt_window_size * rts[i]) / 60.; // RT in minutes
					if(rt_start < 0.) rt_start = 0.;
					rt_stop = (rts[i] + rel_rt_window_size * rts[i]) / 60.; // RT in minutes
				}
				outs << mz << "\t"<< rt_start << "\t" << rt_stop << "\n";
			}
		}

		outs.close();
	}
    

	void InclusionExclusionList::writeTargets(FeatureMap<>& map,String& out_path,DoubleReal rel_rt_window_size,bool rt_in_seconds)
	{
		std::ofstream outs(out_path.c_str());
		if(!outs)
		{
			throw Exception::UnableToCreateFile(__FILE__, __LINE__,__PRETTY_FUNCTION__,"Cannot open output file.");
		}

    DoubleReal min_to_s_factor = rt_in_seconds ? 1.0 : (1.0/60.0) ;
		for(Size f = 0; f < map.size(); ++f)
		{
			DoubleReal rt_start =  map[f].getRT() - map[f].getRT() * rel_rt_window_size;
			if(rt_start < 0.) rt_start = 0.;
			DoubleReal rt_end =  map[f].getRT() + map[f].getRT() * rel_rt_window_size;
			
			outs << map[f].getMZ() << "\t" << (rt_start*min_to_s_factor) <<"\t" << (rt_end*min_to_s_factor) << "\n";
		}
		outs.close();
	}
	
	void InclusionExclusionList::writeTargets(std::vector<PeptideIdentification>& pep_ids,String& out_path,
																						DoubleReal rel_rt_window_size,IntList& charges,bool rt_in_seconds)
	{
		std::ofstream outs(out_path.c_str());
		outs.precision(8);
		if(!outs)
		{
			throw Exception::UnableToCreateFile(__FILE__, __LINE__,__PRETTY_FUNCTION__,"Cannot open output file.");
		}

    Size charge_invalid_count(0);

		std::vector<PeptideIdentification>::const_iterator pep_id_iter = pep_ids.begin();
		for(;pep_id_iter != pep_ids.end();++pep_id_iter)
		{
			if(pep_id_iter->getHits().size() > 1)
			{
				Exception::InvalidSize(__FILE__, __LINE__,__PRETTY_FUNCTION__,pep_id_iter->getHits().size());
			}
			if(!pep_id_iter->metaValueExists("RT"))
			{
				Exception::MissingInformation(__FILE__, __LINE__,__PRETTY_FUNCTION__,"Peptide identification contains no RT information.");
			}
			DoubleReal rt = pep_id_iter->getMetaValue("RT");
			DoubleReal rt_start,rt_stop;
			if(rt_in_seconds)
			{
				rt_start = (rt - rel_rt_window_size * rt);
				if(rt_start < 0.) rt_start = 0.;
				rt_stop = (rt + rel_rt_window_size * rt); 
			}
			else
			{
				rt_start = (rt - rel_rt_window_size * rt) / 60.; // RT in minutes
				if(rt_start < 0.) rt_start = 0.;
				rt_stop = (rt + rel_rt_window_size * rt) / 60.; // RT in minutes
			}
			std::vector<PeptideHit>::const_iterator pep_hit_iter = pep_id_iter->getHits().begin();
			for(;pep_hit_iter != pep_id_iter->getHits().end();++pep_hit_iter)
			{
				Int charge = pep_hit_iter->getCharge();
        if (charge == 0)
        {
          ++charge_invalid_count;
          //fix charge
          charge = 2;
        }

        bool charge_found = false;
				for(Size c = 0; c < charges.size();++c)
				{
					DoubleReal mz = pep_hit_iter->getSequence().getMonoWeight(Residue::Full,charges[c])/(DoubleReal)charges[c];
					outs << mz <<"\t"<<rt_start<<"\t"<<rt_stop<<"\n";
					if(charges[c] == charge)
					{
						charge_found = true;
					}
				}
        if(!charge_found) // if not already done, consider annotated charge of peptide (unless its 0)
				{
					DoubleReal mz = pep_hit_iter->getSequence().getMonoWeight(Residue::Full,charge)/(DoubleReal)charge;
					outs << mz <<"\t"<<rt_start<<"\t"<<rt_stop<<"\n";
				}
			}
		}
    
    if (charge_invalid_count>0) LOG_WARN << "Warning: " << charge_invalid_count << " peptides with charge=0 were found, and assumed to have charge=2.\n";

		outs.close();
						
	}

} // namespace
