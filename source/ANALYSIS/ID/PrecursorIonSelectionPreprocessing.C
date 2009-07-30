// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Alexandra Zerck $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/ANALYSIS/ID/PrecursorIonSelectionPreprocessing.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/SIMULATION/DetectabilitySimulation.h>
#include <OpenMS/SIMULATION/RTSimulation.h>

using namespace std;
//#define PISP_DEBUG
#undef PISP_DEBUG
namespace OpenMS
{
	PrecursorIonSelectionPreprocessing::PrecursorIonSelectionPreprocessing()
		: DefaultParamHandler("PrecursorIonSelectionPreprocessing"),
			f_max_(0)
	{
		defaults_.setValue("precursor_mass_tolerance", 10., "Precursor mass tolerance which is used to query the peptide database for peptides");
		defaults_.setMinFloat("precursor_mass_tolerance",0.);
		defaults_.setValue("rt_weighting:total_gradient_time",7640.,"the total gradient time in seconds, needed for normalization.");
		defaults_.setValue("rt_weighting:gradient_offset",600,"the total gradient time in seconds, needed for normalization.");

		defaults_.setValue("rt_weighting:gauss_amplitude",100.,"amplitude at the gauss_mean");
		defaults_.setValue("rt_weighting:gauss_mean",0.0,"mean of the gauss curve");
		defaults_.setValue("rt_weighting:gauss_std",0.01,"std of the gauss curve");
		//		defaults_.setMinFloat("rt_tolerance",0.);
		defaults_.setValue("precursor_mass_tolerance_unit", "ppm", "Precursor mass tolerance unit.");
		defaults_.setValidStrings("precursor_mass_tolerance_unit",StringList::create("ppm,Da"));
		defaults_.setValue("preprocessing:preprocessed_db_path","","Path where the preprocessed database should be stored");
		defaults_.setValue("preprocessing:preprocessed_db_pred_rt_path","","Path where the predicted rts of the preprocessed database should be stored");
		defaults_.setValue("preprocessing:preprocessed_db_pred_dt_path","","Path where the predicted rts of the preprocessed database should be stored");
		defaults_.setValue("preprocessing:max_peptides_per_run",100000,"Number of peptides for that the pt and rt are parallely predicted.");
		defaults_.setMinInt("preprocessing:max_peptides_per_run",1);
		defaults_.setValue("missed_cleavages",1,"Number of allowed missed cleavages.");
		defaults_.setMinInt("missed_cleavages",0);
		defaults_.setValue("preprocessing:taxonomy","","Taxonomy");
		defaults_.setValue("tmp_dir","","Absolute path to tmp data directory used to store files needed for rt and dt prediction.");
		defaultsToParam_();
	}

	PrecursorIonSelectionPreprocessing::PrecursorIonSelectionPreprocessing(const PrecursorIonSelectionPreprocessing& source)
		: DefaultParamHandler(source),		
			sequences_(source.sequences_),
			prot_masses_(source.prot_masses_),
			bin_masses_(source.bin_masses_),
			f_max_(source.f_max_)
	{
		
	}

	PrecursorIonSelectionPreprocessing::~PrecursorIonSelectionPreprocessing()
	{
	  //????
	}

	PrecursorIonSelectionPreprocessing& PrecursorIonSelectionPreprocessing::operator = (const PrecursorIonSelectionPreprocessing& source)
	{
	  if (&source != this)
 	 {
		 DefaultParamHandler::operator=(source);
 	   sequences_ = source.sequences_;
 	   prot_masses_ = source.prot_masses_;
 	   bin_masses_ = source.bin_masses_;
 	   f_max_ = source.f_max_;
	  }
	  return *this;
	}
	
	
// 	const std::set<AASequence>& PrecursorIonSelectionPreprocessing::getSequences() const 
// 	{
// 	  return sequences_; 
// 	}
	
	
	const std::map<String,std::vector<DoubleReal> >& PrecursorIonSelectionPreprocessing::getProtMasses() const
	{
	  return prot_masses_; 
	}

	const std::vector<DoubleReal> & PrecursorIonSelectionPreprocessing::getMasses(String acc)
	{
	  return prot_masses_[acc]; 
	}
	
	DoubleReal PrecursorIonSelectionPreprocessing::getRT(String peptide)
	{
		if(rt_map_.size()>0)
			{
				if(rt_map_.find(peptide)!=rt_map_.end()) return rt_map_[peptide];
				else return -1;
			}
		std::cout << "rt_map is empty, no rts predicted!"<<std::endl;
		return -1;
	}

	
	DoubleReal PrecursorIonSelectionPreprocessing::getPT(String peptide)
	{
		if(pt_map_.size()>0)
			{
				if(pt_map_.find(peptide)!=pt_map_.end()) return pt_map_[peptide];
				else return 1;
			}
		std::cout << "pt_map is empty, no detectabilities predicted!"<<std::endl;
		return 1;
	}

	DoubleReal PrecursorIonSelectionPreprocessing::getRT(String prot_id,Size peptide_index)
	{
		if(rt_prot_map_.size()>0)
			{
// 				if(rt_prot_map_.find(prot_id)!=rt_prot_map_.end())
// 					{
						if(rt_prot_map_[prot_id].size() >peptide_index)	return rt_prot_map_[prot_id][peptide_index];
						else return -1;
// 					}
// 				else return -1;
			}
		std::cout << "rt_map is empty, no rts predicted!"<<std::endl;
		return -1;
	}

	DoubleReal PrecursorIonSelectionPreprocessing::getPT(String prot_id,Size peptide_index)
	{
		if(pt_prot_map_.size()>0)
			{
// 				if(pt_prot_map_.find(prot_id)!=pt_prot_map_.end())
// 					{
						if(pt_prot_map_[prot_id].size() >peptide_index)	return pt_prot_map_[prot_id][peptide_index];
						else return 1;
// 					}
// 				else return 1;
			}
		std::cout << "pt_map is empty, no detectabilities predicted!"<<std::endl;
		return 1;
	}

	
	DoubleReal PrecursorIonSelectionPreprocessing::getRTWeight(String prot_id, Size peptide_index,DoubleReal meas_rt)
	{
		DoubleReal pred_rt = getRT(prot_id,peptide_index);
		//		std::cout << "pred rt: "<<pred_rt << std::endl;
		// TODO: what to return if no rt was predicted for this peptide?
		if(pred_rt == -1) return 1.;
		// determine difference of measured and predicted rt and normalize by total gradient time
		DoubleReal diff = (meas_rt - pred_rt)/(DoubleReal)param_.getValue("rt_weighting:total_gradient_time");
		// get parameters for gauss curve representing the distribution of the rt differences
		DoubleReal a = param_.getValue("rt_weighting:gauss_amplitude");
		DoubleReal m = param_.getValue("rt_weighting:gauss_mean");
		DoubleReal s = param_.getValue("rt_weighting:gauss_std");
		// get gauss value for the specific rt difference
		DoubleReal gauss_diff = a*exp(-1.0 *pow(diff-m,2)/(2*pow(s,2)));
		return gauss_diff;
	}


	
	DoubleReal PrecursorIonSelectionPreprocessing::getRTWeight(String peptide,DoubleReal meas_rt)
	{
		DoubleReal pred_rt = getRT(peptide);
		std::cout << "pep: "<<peptide << "\t rt: "<<pred_rt << std::endl;
		// TODO: what to return if no rt was predicted for this peptide?
		if(pred_rt == -1) return 1.;
		// determine difference of measured and predicted rt and normalize by total gradient time
		DoubleReal diff = (meas_rt - pred_rt)/(DoubleReal)param_.getValue("rt_weighting:total_gradient_time");
		// get parameters for gauss curve representing the distribution of the rt differences
		DoubleReal a = param_.getValue("rt_weighting:gauss_amplitude");
		DoubleReal m = param_.getValue("rt_weighting:gauss_mean");
		DoubleReal s = param_.getValue("rt_weighting:gauss_std");
		// get gauss value for the specific rt difference
		DoubleReal gauss_diff = a*exp(-1.0 *pow(diff-m,2)/(2*pow(s,2)));
		return gauss_diff;
	}


	
	DoubleReal PrecursorIonSelectionPreprocessing::getWeight(DoubleReal mass)
	{
		if(param_.getValue("precursor_mass_tolerance_unit") == "Da")
			{
				return (DoubleReal)counter_[(Size) floor((mass - masses_[0])/(DoubleReal)param_.getValue("precursor_mass_tolerance") +0.5)]/(DoubleReal)f_max_;
			}
		else // 
			{
#ifdef PISP_DEBUG
				std::cout << bin_masses_.size() << " "<< mass << std::endl;
				std::cout << *(bin_masses_.begin()) << " "<< *(bin_masses_.end()-1) << std::endl;
#endif
				std::vector<DoubleReal>::iterator tmp_iter = bin_masses_.begin();
				while(tmp_iter!=bin_masses_.end() && *tmp_iter<mass)
					{
						++tmp_iter;
					}
				if(tmp_iter!=bin_masses_.begin()) --tmp_iter;
#ifdef PISP_DEBUG
				if(tmp_iter!=bin_masses_.end())	std::cout << "tmp_iter "<<*tmp_iter << std::endl;
				else std::cout << "tmp_iter am ende"<< std::endl;
				std::cout << "fmax "<<f_max_ << std::endl;
#endif
				if((tmp_iter+1)==bin_masses_.end()
					 || fabs(*tmp_iter - mass) < fabs(*(tmp_iter+1) - mass))
					{
						return (DoubleReal)counter_[distance(bin_masses_.begin(),tmp_iter)]/(DoubleReal)f_max_;
					}
				else return (DoubleReal)counter_[distance(bin_masses_.begin(),tmp_iter+1)]/(DoubleReal)f_max_;
			}
	}

	void PrecursorIonSelectionPreprocessing::loadPreprocessing()
	{
		// first check if preprocessed db already exists
		String path = param_.getValue("preprocessing:preprocessed_db_path");
		
		// check if file exists
		std::ifstream test(path.c_str());
		if(test)
			{
				loadPreprocessedDB_(path);
			}
		else Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, path);
	}

	
	void PrecursorIonSelectionPreprocessing::dbPreprocessing(String db_path,String rt_model_path,
																													 String dt_model_path,bool save)
	{
		//#ifdef PISP_DEBUG
		std::cout << "Parameters: "<< param_.getValue("preprocessing:preprocessed_db_path")
							<< "\t" << param_.getValue("precursor_mass_tolerance")
							<< " " << param_.getValue("precursor_mass_tolerance_unit")
			//<< "\t"<<param_.getValue("rt_tolerance")
							<< "\t"<<param_.getValue("rt_weighting:total_gradient_time")
							<< "\t"<<param_.getValue("rt_weighting:gauss_amplitude")
							<< "\t"<<param_.getValue("rt_weighting:gauss_mean")
							<< "\t"<<param_.getValue("rt_weighting:gauss_std")
							<< "\t" << param_.getValue("missed_cleavages")
							<< "\t" << param_.getValue("preprocessing:taxonomy")
							<< "\t" << param_.getValue("tmp_dir") << "---"
							<< std::endl;
		//#endif

		FASTAFile fasta_file;
		std::vector<FASTAFile::FASTAEntry> entries;
		fasta_file.load(db_path,entries);
		EnzymaticDigestion digest;
		digest.setMissedCleavages((UInt)param_.getValue("missed_cleavages"));
		String tmp_filename = (String)param_.getValue("tmp_dir") + "/sequences";

		// first get all protein sequences and calculate digest
		for(UInt e=0;e<entries.size();++e)
			{
				
				// filter for taxonomy
				if(entries[e].description.toUpper().hasSubstring(((String)param_.getValue("preprocessing:taxonomy")).toUpper())) 
					{
#ifdef PISP_DEBUG
						std::cout << entries[e].identifier << std::endl;
#endif
						if(entries[e].identifier.hasPrefix("sp|") || entries[e].identifier.hasPrefix("tr|"))
							{
								entries[e].identifier = entries[e].identifier.suffix(entries[e].identifier.size()-3);
							}
						entries[e].identifier = entries[e].identifier.prefix('|');
						String& seq = entries[e].sequence;
						// check for unallowed characters
						if(seq.hasSubstring("X") || seq.hasSubstring("B") ||  seq.hasSubstring("Z") )
							{
								continue;
							}
						std::vector<DoubleReal> prot_masses;
						// digest sequence
						AASequence aa_seq(seq);
						std::vector<AASequence> vec;
						digest.digest(aa_seq,vec);

						// enter peptide sequences in map
						std::vector<AASequence>::iterator vec_iter = vec.begin();
						for(;vec_iter != vec.end();++vec_iter)
							{
								// write peptide seq in temporary file, for rt prediction
								//seq_file << *vec_iter << "\n";
								DoubleReal mass = vec_iter->getMonoWeight(Residue::Full,1);
								prot_masses.push_back(mass);
								if(sequences_.count(*vec_iter)==0) // peptide sequences are considered only once
									{
										sequences_.insert(*vec_iter);
										masses_.push_back(mass);
									}
							}
						prot_masses_.insert(make_pair(entries[e].identifier,prot_masses));
					}

			}
		entries.clear();
		std::cout << "now make the rt and pt predictions"<<std::endl;
		std::set<AASequence>::iterator seq_it = sequences_.begin();
		std::vector<String> peptide_sequences;
		UInt index = 0;
		DoubleReal total_gradient_time = param_.getValue("rt_weighting:total_gradient_time");
		DoubleReal gradient_offset = param_.getValue("rt_weighting:gradient_offset");
		Param rt_param;
		rt_param.setValue("rt_model_file",rt_model_path);
		Param dt_param;
		dt_param.setValue("dt_simulation_on","true");
		dt_param.setValue("dt_model_file",dt_model_path);
		gsl_rng* random_generator = gsl_rng_alloc(gsl_rng_mt19937);// not needed for this rt prediction, but needed for RTSimulation
		// this is needed, as too many sequences require too much memory for the rt and dt prediction
		Size max_peptides_per_run = (Int)param_.getValue("preprocessing:max_peptides_per_run");
		peptide_sequences.resize(max_peptides_per_run);
		std::cout << sequences_.size()<<" peptides for predictions."<<std::endl;
		std::set<AASequence>::iterator seq_it_end = sequences_.end();
		if(seq_it != seq_it_end)  --seq_it_end;
		for(;seq_it!=sequences_.end();++seq_it)
			{
				peptide_sequences[index] = seq_it->toUnmodifiedString();
				++index;
				if(index == max_peptides_per_run || seq_it == seq_it_end)
					{
						// now make RTPrediction using the RTSimulation class of the simulator
						RTSimulation rt_sim(random_generator);
						rt_sim.setParameters(rt_param);
						std::vector<DoubleReal> rts;
						rt_sim.wrapSVM(peptide_sequences,rts);
						
						for(Size index2 = 0; index2 < rts.size();++index2)
							{
								rt_map_.insert(std::make_pair(peptide_sequences[index2],rts[index2]*total_gradient_time+gradient_offset));
							}
						rts.clear();
						
						// now make DTPrediction using the DetectabilitySimulation class of the simulator
						DetectabilitySimulation dt_sim;
						dt_sim.setParameters(dt_param);
						std::vector<DoubleReal> labels;
						std::vector<DoubleReal> detectabilities;
						dt_sim.predictDetectabilities(peptide_sequences,labels,detectabilities);
						
						for(Size index2 = 0; index2 < detectabilities.size();++index2)
							{
								pt_map_.insert(make_pair(peptide_sequences[index2],detectabilities[index2]));
							}
						peptide_sequences.clear();
						Int size = std::min((int)distance(seq_it,sequences_.end())-1,(int)max_peptides_per_run);
						std::cout << "peptide_sequences.resize(size) "<<size<<std::endl;
 						peptide_sequences.resize(size);
						index = 0;
					}
			}
		peptide_sequences.clear();
		sequences_.clear();
		std::cout << "Finished predictions!"<<std::endl;
		if(masses_.size() == 0)
			{
				std::cout << "no masses entered" << std::endl;
				return;
			}
		std::sort(masses_.begin(),masses_.end());
		// now get minimal and maximal mass and create counter_-vectors
		// count mass occurences using bins
#ifdef PISP_DEBUG
		std::cout << "min\tmax "<<masses_[0] << "\t"<<*(masses_.end()-1)<<std::endl;
		std::cout << "prot_masses.size() "<<prot_masses_.size()<<std::endl;
#endif
		// if the precursor mass tolerance is given in Da
		// we have equidistant bins
		if(param_.getValue("precursor_mass_tolerance_unit") == "Da")
			{
				counter_.resize((UInt)(ceil((*(masses_.end()-1))-masses_[0]) / (DoubleReal)param_.getValue("precursor_mass_tolerance")) +2,0);
				for(UInt i=0;i<masses_.size();++i)
					{
						// get bin index
						DoubleReal tmp = (masses_[i] - masses_[0] ) / (DoubleReal)param_.getValue("precursor_mass_tolerance");
						++counter_[(Size) ceil(tmp)];
					}
				UInt max = 0;
				for(UInt i=0;i<counter_.size();++i)
					{
						if(counter_[i] >  max )  max = counter_[i];
					}
				// store maximal frequency
				f_max_ = max;
			}
		else // with ppm, we have bins of increasing size and need to save the bin limits
			{
				bin_masses_.clear();
				counter_.clear();
				// so first we calculate the boundings of the bins
				DoubleReal curr_mass = masses_[0];
#ifdef PISP_DEBUG 
				std::cout << "min_max_curr "<<masses_[0] << " "<<*(masses_.end()-1) << " "<< curr_mass << std::endl;
#endif
				UInt size = 0;
				while(curr_mass < *(masses_.end()-1))
					{
						/// store the lower bound of the current bin
						bin_masses_.push_back(curr_mass); 
						++size;
						/// calculate lower bound for next bin
						curr_mass = curr_mass + curr_mass*(DoubleReal)param_.getValue("precursor_mass_tolerance")/1e06;
					}
#ifdef PISP_DEBUG
				std::cout <<"bin_masses_.size() " <<  bin_masses_.size() << " "<<size<< std::endl;
#endif
				counter_.resize(size,0);
#ifdef PISP_DEBUG
				std::cout << "masses_.size() "<< masses_.size()
									<< " counter_.size() "<<counter_.size()
									<< " bin_masses_.size() "<<bin_masses_.size()<<std::endl;

#endif
				std::vector<DoubleReal>::iterator old_begin = bin_masses_.begin();
				// then we put the peptide masses into the right bins
				for(UInt i=0;i<masses_.size();++i)
					{
#ifdef PISP_DEBUG
						std::cout << "i "<<i << std::endl;
						StopWatch timer;
						timer.start();
#endif
						std::vector<DoubleReal>::iterator tmp_iter = old_begin;

#ifdef PISP_DEBUG
						timer.start();
						std::cout << "old_begin "<<*old_begin << " masses_[i] "<<masses_[i]<<std::endl;
						if(old_begin != bin_masses_.begin()) std::cout << "old_begin-1 "<<*(old_begin-1)<<std::endl;
#endif
						while(tmp_iter != bin_masses_.end() &&  *tmp_iter < masses_[i])
							{
								++tmp_iter;
							}
#ifdef PISP_DEBUG
						timer.stop();
						std::cout << "while schleife\t"; 
						std::cout << timer.getCPUTime()<<"\t";
						timer.reset();
#endif
						
						
						old_begin = tmp_iter;
						if(tmp_iter== bin_masses_.end())
							{
								++counter_[distance(bin_masses_.begin(),tmp_iter-1)];
							}
						else if((tmp_iter+1)==bin_masses_.end() // last entry
										|| fabs(*tmp_iter - masses_[i]) < fabs(*(tmp_iter+1) - masses_[i])) // or closer to left entry
							{
								// increase counter for bin to the left
								++counter_[distance(bin_masses_.begin(),tmp_iter)];
							}
						else ++counter_[distance(bin_masses_.begin(),tmp_iter+1)]; // increase right counter
#ifdef PISP_DEBUG
						timer.stop();
						std::cout << timer.getCPUTime ()<<std::endl;
#endif
					}
				// determine maximal frequency for normalization
				UInt max = 0;
				for(UInt i=0;i<counter_.size();++i)
					{
						if(counter_[i] >  max )  max = counter_[i];
					}
				f_max_ = max;
#ifdef PISP_DEBUG
				std::cout << "f_max_ "<<f_max_ <<std::endl;
#endif
					
			}
		if(save)
			{
				savePreprocessedDBWithRT_(db_path,(String)param_.getValue("preprocessing:preprocessed_db_path"));
			}
		
	}
	
	void PrecursorIonSelectionPreprocessing::dbPreprocessing(String db_path,bool save)
	{

#ifdef PISP_DEBUG
		std::cout << "Parameters: "<< param_.getValue("preprocessing:preprocessed_db_path")
							<< "\t" << param_.getValue("precursor_mass_tolerance")
							<< " " << param_.getValue("precursor_mass_tolerance_unit")
							<< "\t" << param_.getValue("missed_cleavages")
							<< "\t" << param_.getValue("preprocessing:taxonomy") <<std::endl;
#endif

		FASTAFile fasta_file;
		std::vector<FASTAFile::FASTAEntry> entries;
		fasta_file.load(db_path,entries);
		EnzymaticDigestion digest;
		digest.setMissedCleavages((UInt)param_.getValue("missed_cleavages"));

		// first get all protein sequences and calculate digest
		for(UInt e=0;e<entries.size();++e)
			{
				// filter for taxonomy
				if(entries[e].description.toUpper().hasSubstring(((String)param_.getValue("preprocessing:taxonomy")).toUpper())) 
					{
#ifdef PISP_DEBUG
						std::cout << entries[e].identifier << std::endl;
#endif
						entries[e].identifier = entries[e].identifier.prefix('|');
						String& seq = entries[e].sequence;
						// check for unallowed characters
						if(seq.hasSubstring("X") || seq.hasSubstring("B") ||  seq.hasSubstring("Z") )
							{
								continue;
							}
						std::vector<DoubleReal> prot_masses;
						// digest sequence
						AASequence aa_seq(seq);
						std::vector<AASequence> vec;
						digest.digest(aa_seq,vec);

						// enter peptide sequences in map
						std::vector<AASequence>::iterator vec_iter = vec.begin();
						for(;vec_iter != vec.end();++vec_iter)
							{
								DoubleReal mass = vec_iter->getMonoWeight(Residue::Full,1);
								prot_masses.push_back(mass);
								if(sequences_.count(*vec_iter)==0) // peptide sequences are considered only once
									{
										sequences_.insert(*vec_iter);
										masses_.push_back(mass);
									}
							}
						prot_masses_.insert(make_pair(entries[e].identifier,prot_masses));
					}

			}

		if(masses_.size() == 0)
			{
				std::cout << "no masses entered" << std::endl;
				return;
			}
		std::sort(masses_.begin(),masses_.end());
		// now get minimal and maximal mass and create counter_-vectors
		// count mass occurences using bins
#ifdef PISP_DEBUG
		std::cout << "min\tmax "<<masses_[0] << "\t"<<*(masses_.end()-1)<<std::endl;
		std::cout << "prot_masses.size() "<<prot_masses_.size()<<std::endl;
#endif
		// if the precursor mass tolerance is given in Da
		// we have equidistant bins
		if(param_.getValue("precursor_mass_tolerance_unit") == "Da")
			{
				counter_.resize((UInt)(ceil((*(masses_.end()-1))-masses_[0]) / (DoubleReal)param_.getValue("precursor_mass_tolerance")) +2,0);
				for(UInt i=0;i<masses_.size();++i)
					{
						// get bin index
						DoubleReal tmp = (masses_[i] - masses_[0] ) / (DoubleReal)param_.getValue("precursor_mass_tolerance");
						++counter_[(Size) ceil(tmp)];
					}
				UInt max = 0;
				for(UInt i=0;i<counter_.size();++i)
					{
						if(counter_[i] >  max )  max = counter_[i];
					}
				// store maximal frequency
				f_max_ = max;
			}
		else // with ppm, we have bins of increasing size and need to save the bin limits
			{
				bin_masses_.clear();
				counter_.clear();
				// so first we calculate the boundings of the bins
				DoubleReal curr_mass = masses_[0];
#ifdef PISP_DEBUG 
				std::cout << "min_max_curr "<<masses_[0] << " "<<*(masses_.end()-1) << " "<< curr_mass << std::endl;
#endif
				UInt size = 0;
				while(curr_mass < *(masses_.end()-1))
					{
						/// store the lower bound of the current bin
						bin_masses_.push_back(curr_mass); 
						++size;
						/// calculate lower bound for next bin
						curr_mass = curr_mass + curr_mass*(DoubleReal)param_.getValue("precursor_mass_tolerance")/1e06;
					}
#ifdef PISP_DEBUG
				std::cout <<"bin_masses_.size() " <<  bin_masses_.size() << " "<<size<< std::endl;
#endif
				counter_.resize(size,0);
#ifdef PISP_DEBUG
				std::cout << "masses_.size() "<< masses_.size()
									<< " counter_.size() "<<counter_.size()
									<< " bin_masses_.size() "<<bin_masses_.size()<<std::endl;

#endif
				std::vector<DoubleReal>::iterator old_begin = bin_masses_.begin();
				// then we put the peptide masses into the right bins
				for(UInt i=0;i<masses_.size();++i)
					{
#ifdef PISP_DEBUG
						std::cout << "i "<<i << std::endl;
						StopWatch timer;
						timer.start();
#endif
						std::vector<DoubleReal>::iterator tmp_iter = old_begin;

#ifdef PISP_DEBUG
						timer.start();
						std::cout << "old_begin "<<*old_begin << " masses_[i] "<<masses_[i]<<std::endl;
						if(old_begin != bin_masses_.begin()) std::cout << "old_begin-1 "<<*(old_begin-1)<<std::endl;
#endif
						while(tmp_iter != bin_masses_.end() &&  *tmp_iter < masses_[i])
							{
								++tmp_iter;
							}
#ifdef PISP_DEBUG
						timer.stop();
						std::cout << "while schleife\t"; 
						std::cout << timer.getCPUTime()<<"\t";
						timer.reset();
#endif
						
						
						old_begin = tmp_iter;
						if(tmp_iter== bin_masses_.end())
							{
								++counter_[distance(bin_masses_.begin(),tmp_iter-1)];
							}
						else if((tmp_iter+1)==bin_masses_.end() // last entry
										|| fabs(*tmp_iter - masses_[i]) < fabs(*(tmp_iter+1) - masses_[i])) // or closer to left entry
							{
								// increase counter for bin to the left
								++counter_[distance(bin_masses_.begin(),tmp_iter)];
							}
						else ++counter_[distance(bin_masses_.begin(),tmp_iter+1)]; // increase right counter
#ifdef PISP_DEBUG
						timer.stop();
						std::cout << timer.getCPUTime ()<<std::endl;
#endif
					}
				// determine maximal frequency for normalization
				UInt max = 0;
				for(UInt i=0;i<counter_.size();++i)
					{
						if(counter_[i] >  max )  max = counter_[i];
					}
				f_max_ = max;
#ifdef PISP_DEBUG
				std::cout << "f_max_ "<<f_max_ <<std::endl;
#endif
					
			}
		if(save)
			{
				savePreprocessedDB_(db_path,(String)param_.getValue("preprocessing:preprocessed_db_path"));
			}

	}

	void PrecursorIonSelectionPreprocessing::savePreprocessedDB_(String db_path,String path)
	{
		std::ofstream out(path.c_str());
		out.precision(10);
		if (!out)
			{
				throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, path);
			}
		// header: db_name  precursor_mass_tolerance precursor_mass_tolerance_unit  taxonomy);
		Size pos1 = db_path.rfind("/") +1;
		// get db-name
		Size pos2 = db_path.rfind(".");
		String db_name = db_path.substr(pos1,pos2-pos1);
		out << db_name <<"\t" <<param_.getValue("precursor_mass_tolerance")  << "\t"
				<< param_.getValue("precursor_mass_tolerance_unit")
				<< "\t"<< (String)param_.getValue("preprocessing:taxonomy");
		// first save protein_masses_map
		out << prot_masses_.size() <<std::endl;
#ifdef PISP_DEBUG
		std::cout << prot_masses_.size() << " "<<counter_.size() << " "<< bin_masses_.size()<< std::endl;
#endif
		std::map<String,std::vector<DoubleReal> >::iterator pm_iter = prot_masses_.begin();
		for(;pm_iter!=prot_masses_.end();++pm_iter)
			{
					out << pm_iter->second.size() << "\t"<< pm_iter->first;
					for(UInt i = 0; i<pm_iter->second.size();++i)
					{
							out << "\t" << pm_iter->second[i];
					}
					out << "\n";
			}
		//	out.close();
#ifdef PISP_DEBUG
		std::cout << (path + "_counter") << std::endl;
		std::cout  <<"counter: "  << counter_.size() << "\t" << masses_[0] << "\t"<<*(masses_.end()-1) <<"\n";
#endif
// 		// now save counter
// 		std::ofstream out2((path + "_counter").c_str());
// 		if (!out2)
// 			{
// 				throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, path+ "_counter");
// 			}
// 		out2.precision(10);
		
		// header: number of entries, min, max
		out << "###\n";
		out << counter_.size() << "\t" << masses_[0] << "\t"<<*(masses_.end()-1) <<"\n";
		for(UInt i = 0; i < counter_.size(); ++i)
			{
					out << counter_[i] << "\t";
			}
		out << "\n";
		if(path.hasSubstring("ppm"))
			{
#ifdef PISP_DEBUG
				std::cout << (path + "_bin_masses") << std::endl;
#endif
// 				// now save bin_masses
// 				std::ofstream out3((path + "_bin_masses").c_str());
// 				if (!out3)
// 					{
// 						throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, path+ "_bin_masses");
// 					}
				
// 				out3.precision(10);				
				// header: number of entries
				out << "###\n";
				out << bin_masses_.size() << "\n";
				for(UInt i = 0; i < bin_masses_.size(); ++i)
					{
						out<< bin_masses_[i] << "\n";
					}
			}

	}


	void PrecursorIonSelectionPreprocessing::savePreprocessedDBWithRT_(String db_path,String path)
	{
		std::ofstream out(path.c_str());
		out.precision(10);
		if (!out)
			{
				throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, path);
			}
		// header: db_name  precursor_mass_tolerance precursor_mass_tolerance_unit  taxonomy);
		Size pos1 = db_path.rfind("/") +1;
		// get db-name
		Size pos2 = db_path.rfind(".");
		String db_name = db_path.substr(pos1,pos2-pos1);
		out << db_name <<"\t" <<param_.getValue("precursor_mass_tolerance")  << "\t"
				<< param_.getValue("precursor_mass_tolerance_unit")
				<< "\t"<< (String)param_.getValue("preprocessing:taxonomy");
		// first save protein_masses_map
		out << prot_masses_.size() <<std::endl;
#ifdef PISP_DEBUG
		std::cout << prot_masses_.size() << " "<<counter_.size() << " "<< bin_masses_.size()<< std::endl;
#endif
		FASTAFile fasta_file;
		std::vector<FASTAFile::FASTAEntry> entries;
		fasta_file.load(db_path,entries);
		EnzymaticDigestion digest;
		digest.setMissedCleavages((UInt)param_.getValue("missed_cleavages"));
		
		// first get all protein sequences and calculate digest
		for(UInt e=0;e<entries.size();++e)
			{
				// filter for taxonomy
				if(entries[e].description.toUpper().hasSubstring(((String)param_.getValue("preprocessing:taxonomy")).toUpper())) 
					{
#ifdef PISP_DEBUG
						std::cout << entries[e].identifier << std::endl;
#endif
						if(entries[e].identifier.hasPrefix("sp|") || entries[e].identifier.hasPrefix("tr|"))
							{
								entries[e].identifier = entries[e].identifier.suffix(entries[e].identifier.size()-3);
							}
						entries[e].identifier = entries[e].identifier.prefix('|');
						String& seq = entries[e].sequence;
						// check for unallowed characters
						if(seq.hasSubstring("X") || seq.hasSubstring("B") ||  seq.hasSubstring("Z") )
							{
								continue;
							}
						std::vector<DoubleReal> prot_masses;
						// digest sequence
						AASequence aa_seq(seq);
						std::vector<AASequence> vec;
						digest.digest(aa_seq,vec);
						out << vec.size() << "\t"<<entries[e].identifier;
						// enter peptide sequences in map
						std::vector<AASequence>::iterator vec_iter = vec.begin();
						for(;vec_iter != vec.end();++vec_iter)
							{
								// write peptide seq in temporary file, for rt prediction
								//seq_file << *vec_iter << "\n";
								DoubleReal mass = vec_iter->getMonoWeight(Residue::Full,1);
								// out : masse, rt, pt
								out << "\t"<<mass << ","<<getRT(vec_iter->toString())<<","<<getPT(vec_iter->toString());
							}
						out << "\n";
					}

			}
		
#ifdef PISP_DEBUG
		std::cout << (path + "_counter") << std::endl;
		std::cout  <<"counter: "  << counter_.size() << "\t" << masses_[0] << "\t"<<*(masses_.end()-1) <<"\n";
#endif
// 		// now save counter
// 		std::ofstream out2((path + "_counter").c_str());
// 		if (!out2)
// 			{
// 				throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, path+ "_counter");
// 			}
// 		out2.precision(10);
		
		// header: number of entries, min, max
		out << "###\n";
		out << counter_.size() << "\t" << masses_[0] << "\t"<<*(masses_.end()-1) <<"\n";
		for(UInt i = 0; i < counter_.size(); ++i)
			{
					out << counter_[i] << "\t";
			}
		out << "\n";
		if(path.hasSubstring("ppm"))
			{
#ifdef PISP_DEBUG
				std::cout << (path + "_bin_masses") << std::endl;
#endif
// 				// now save bin_masses
// 				std::ofstream out3((path + "_bin_masses").c_str());
// 				if (!out3)
// 					{
// 						throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, path+ "_bin_masses");
// 					}
				
// 				out3.precision(10);				
				// header: number of entries
				out << "###\n";
				out << bin_masses_.size() << "\n";
				for(UInt i = 0; i < bin_masses_.size(); ++i)
					{
						out<< bin_masses_[i] << "\n";
					}
			}

	}


	
	void PrecursorIonSelectionPreprocessing::loadPreprocessedDB_(String path)
	{
		// first get protein_masses_map
		TextFile file;
		file.load(path,true);
		//#ifdef PISP_DEBUG
		std::cout << "load " << path << std::endl;
		//#endif
		TextFile::Iterator iter = file.begin();
		++iter;
		for(; iter != file.end() && !iter->hasPrefix("###"); ++iter)
	    {
 				std::vector<String> parts;
				iter->split('\t',parts);
				std::vector<DoubleReal> masses;
				masses.reserve(parts[0].toInt());
				std::vector<String> line_parts;
				std::vector<DoubleReal> rts;
				std::vector<DoubleReal> pts;
				for(UInt i = 2; i < parts.size();++i)
					{
						// check if rts are stored also
						if(parts[i].hasSubstring(","))
							{
								parts[i].split(',',line_parts);
								masses.push_back(line_parts[0].toDouble());
								if(line_parts.size()>1)
									{
										rts.push_back(line_parts[1].toDouble());
										if(line_parts.size()==3)	pts.push_back(line_parts[2].toDouble());
									}
							}
						else masses.push_back(parts[i].toDouble());
					}
				if(parts[1].hasSubstring(".")) parts[1] = parts[1].prefix(11);
				prot_masses_.insert(make_pair(parts[1],masses));
				if(rts.size()>0) rt_prot_map_.insert(make_pair(parts[1],rts));
				if(pts.size()>0)	pt_prot_map_.insert(make_pair(parts[1],pts));

#ifdef PISP_DEBUG
				std::cout << parts[1] << " "<< masses.size()<< std::endl;
#endif
			}
#ifdef PISP_DEBUG
		std::cout << "loaded"<<std::endl;
#endif
		//		if(pt_prot_map)
// 		file.clear();
// 		// now get the counter_ vector
// 		TextFile file2;
// 		file2.load(path+"_counter",true);
// 		iter = file2.begin();
		++iter;
		std::vector<String> header;
		iter->split('\t',header);
		// store minimal mass
		masses_.push_back(header[1].toFloat());
		
		f_max_ = 0;
		
		++iter;
		std::vector<String> freqs;
		iter->split('\t',freqs);

		for(std::vector<String>::const_iterator f_iter=freqs.begin(); f_iter != freqs.end(); ++f_iter)
	    {
				counter_.push_back(f_iter->toInt());
				if((UInt)f_iter->toInt() > f_max_) f_max_ = f_iter->toInt();
			}
		++iter;
		
		if(param_.getValue("precursor_mass_tolerance_unit")=="ppm")
			{
				if(iter == file.end())  	throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__,
																																		"ppm is used as precursor_mass_tolerance_unit, which requires the file "
																																		+ path+"_bin_masses"+ ", that could not be found." );
				else if(iter->hasPrefix("###"))
					{
						++iter;
#ifdef PISP_DEBUG
						std::cout << "load "<<path<<"_bin_masses"<<std::endl;
#endif
						bin_masses_.reserve(iter->toInt());
						++iter;
						for(; iter != file.end(); ++iter)
							{
								bin_masses_.push_back(iter->toDouble());
							}
						
					}
				else throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__,
																							 "ppm is used as precursor_mass_tolerance_unit, which requires the file "
																							 + path+"_bin_masses"+ ", that could not be found." );
			}


// 		// load rt predictions
// 		if(param_.getValue("preprocessing:preprocessed_db_pred_rt_path")!= ""
// 			 &&File().exists(param_.getValue("preprocessing:preprocessed_db_pred_rt_path")))
// 			{
// 				rt_map_.clear();
// 				TextFile rt_pred_file;
// 				rt_pred_file.load(param_.getValue("preprocessing:preprocessed_db_pred_rt_path"),true);
// 				TextFile::Iterator text_it = rt_pred_file.begin();
// 				for(;text_it != rt_pred_file.end();++text_it)
// 					{
// 						std::vector<String> parts;
// 						text_it->split(' ',parts);
// 						rt_map_.insert(make_pair(parts[0],parts[1].toDouble()));
// 					}
// 			}
// 		std::cout <<"rt_map.size: "<<rt_map_.size()<<std::endl;
		
	}


	
} //namespace

