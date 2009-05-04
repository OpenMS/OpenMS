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
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>

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
		defaults_.setValue("precursor_mass_tolerance_unit", "ppm", "Precursor mass tolerance unit.");
		defaults_.setValidStrings("precursor_mass_tolerance_unit",StringList::create("ppm,Da"));
		defaults_.setValue("preprocessing:preprocessed_db_path","","Path where the preprocessed database should be stored");
		defaults_.setValue("missed_cleavages",1,"Number of allowed missed cleavages.");
		defaults_.setMinInt("missed_cleavages",0);
		defaults_.setValue("preprocessing:taxonomy","","Taxonomy");
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
	
	void PrecursorIonSelectionPreprocessing::loadPreprocessedDB_(String path)
	{
		// first get protein_masses_map
		TextFile file;
		file.load(path,true);
#ifdef PISP_DEBUG
		std::cout << "load " << path << std::endl;
#endif
		TextFile::Iterator iter = file.begin();
		++iter;
		for(; iter != file.end() && !iter->hasPrefix("###"); ++iter)
	    {
				std::vector<String> parts;
				iter->split('\t',parts);
				std::vector<DoubleReal> masses;
				masses.reserve(parts[0].toInt());

				for(UInt i = 2; i < parts.size();++i)
					{
						masses.push_back(parts[i].toDouble());
					}
				if(parts[1].hasSubstring(".")) parts[1] = parts[1].prefix(11);
				prot_masses_.insert(make_pair(parts[1],masses));
#ifdef PISP_DEBUG
				std::cout << parts[1] << " "<< masses.size()<< std::endl;
#endif
			}
#ifdef PISP_DEBUG
		std::cout << "loaded"<<std::endl;
#endif
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
		
	}


	
} //namespace

