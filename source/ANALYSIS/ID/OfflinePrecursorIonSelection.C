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
// $Maintainer: Alexandra Zerck $
// $Authors: $
// --------------------------------------------------------------------------
#include <OpenMS/ANALYSIS/ID/OfflinePrecursorIonSelection.h>

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <OpenMS/KERNEL/ComparatorUtils.h>

using namespace OpenMS;

OfflinePrecursorIonSelection::OfflinePrecursorIonSelection() : DefaultParamHandler("OfflinePrecursorIonSelection")
{
	defaults_.setValue("peptides_per_protein",2,"Minimal number of peptides selected for each protein.");
	defaults_.setMinInt("peptides_per_protein",1);
	defaults_.setValue("ms2_spectra_per_rt_bin",5,"Number of allowed MS/MS spectra in a retention time bin.");
	defaults_.setMinInt("ms2_spectra_per_rt_bin",1);
	defaultsToParam_();
}

OfflinePrecursorIonSelection::~OfflinePrecursorIonSelection()
{

}


std::vector<PeptideIdentification> OfflinePrecursorIonSelection::filterPeptideIds_(std::vector<PeptideIdentification>& pep_ids)
	{
		std::vector<PeptideIdentification> filtered_pep_ids;
		
		for(UInt id_c = 0; id_c < pep_ids.size();++id_c)
			{
				std::vector<PeptideHit> tmp_hits;
				if(pep_ids[id_c].getHits()[0].metaValueExists("Rank"))
					{
						for(UInt hit_c=0; hit_c < pep_ids[id_c].getHits().size();++hit_c)
							{
								if(pep_ids[id_c].getHits()[hit_c].getScore() >= pep_ids[id_c].getSignificanceThreshold() && 
									 (Int)pep_ids[id_c].getHits()[hit_c].getMetaValue("Rank") == 1 )
									{
										tmp_hits.push_back(pep_ids[id_c].getHits()[hit_c]);
									}
							}
					}
				else // if meta value rank doesn't exist, take highest scoring peptide hit
					{
						if(pep_ids[id_c].getHits().size() == 1 &&
							 pep_ids[id_c].getHits()[0].getScore() >= pep_ids[id_c].getSignificanceThreshold())
							{
								tmp_hits.push_back(pep_ids[id_c].getHits()[0]);
							}
						else if(pep_ids[id_c].getHits().size()>1)
							{
								UInt max_score_idx=0;
								for(UInt hit_c=1; hit_c < pep_ids[id_c].getHits().size();++hit_c)
									{
										if(pep_ids[id_c].getHits()[hit_c].getScore() >
											 pep_ids[id_c].getHits()[max_score_idx].getScore())
											{
												max_score_idx = hit_c;
											}
									}
								// check if highest scoring peptide hit is significant
								if(pep_ids[id_c].getHits()[max_score_idx].getScore() >= pep_ids[id_c].getSignificanceThreshold())
									{
										tmp_hits.push_back(pep_ids[id_c].getHits()[max_score_idx]);
									}
							}
					}
				
				if(!tmp_hits.empty()) // if there were significant hits save them
					{
						PeptideIdentification tmp_id = pep_ids[id_c];
						tmp_id.setHits(tmp_hits);
						filtered_pep_ids.push_back(tmp_id);
					}  
			}
	
	return filtered_pep_ids;
}


void OfflinePrecursorIonSelection::computeOptimalSolution(std::vector<ProteinIdentification>& prot_ids,
																													std::vector<PeptideIdentification>& pep_ids,
																													FeatureMap<>& features,
																													FeatureMap<>& /*optimal_set*/,bool /*filter*/)
{
	//TODO: filter peptides???
	// TODO: ensure that prot_id contains one proteinid and several prot_hits
	
  // first map the ids onto the features
	IDMapper mapper;
	mapper.annotate(features,pep_ids,prot_ids);
	
	features.sortByRT();
	DoubleReal min_rt = features[0].getRT();
	DoubleReal max_rt = (features.end()-1)->getRT();
	DoubleReal rt_bin_span = ((String) param_.getValue("rt_prediction:rt_bin_size")).toFloat();
	// create bin_map
	rt_bins_.clear();
	rt_bins_.resize((Size) ceil((max_rt-min_rt)/rt_bin_span));

	// create protein acc map
	std::map<String,Size> accessions;
	Size index = 0;
	for(Size p_id = 0; p_id < prot_ids.size();++p_id)
		{
			for(Size p_h = 0; p_h < prot_ids[p_id].getHits().size();++p_h)
				{
					if(accessions.find(prot_ids[p_id].getHits()[p_h].getAccession()) == accessions.end())
						{
							accessions.insert(make_pair(prot_ids[p_id].getHits()[p_h].getAccession(),index));
							++index;
						}
				}
		}
	
	// create protein_peptide map
	protein_precursor_map_.clear();
	// usually we one protein id with many hits
	//	protein_precursor_map_.resize(accessions.size());
	
	FeatureMap<> filtered_features;
	for(Size i=0; i<features.size();++i) // not really elegant
		{
			// filter for features with an id
			if(features[i].getPeptideIdentifications().size()>0)
				{
					filtered_features.push_back(features[i]);

					// enter index in filtered_features in rt_bins_
					rt_bins_[(Size) floor((features[i].getRT()-min_rt) / rt_bin_span)].push_back((Int)filtered_features.size()-1);

					// find protein id for peptide id
					for(Size id = 0; id < features[i].getPeptideIdentifications().size();++id)
						{
							for(Size h = 0; h < features[i].getPeptideIdentifications()[id].getHits().size(); ++h)
								{
									const std::vector<String>& accs = features[i].getPeptideIdentifications()[id].getHits()[h].getProteinAccessions();

									// enter feature index in the corresponding protein_precursor vector
									for(Size a = 0; a < accs.size();++a)
										{
											if(accessions.find(accs[a]) == accessions.end())
												{
													std::cout << "huch accession war nicht in den prot_ids"<<std::endl;
												}
											else
												{
													protein_precursor_map_[accessions[accs[a]]].push_back(filtered_features.size()-1);
												}
										}
								}
						}

				}
		}




// 	// build model
// 	ILPWrapper ilp_wrapper;
//   ilp_wrapper.encodeModelForOptimalSolution(filtered_features,rt_bins_,protein_precursor_map_);


// 	// compute solution
// 	std::vector<int> solution_indices;
//   ilp_wrapper.solve(solution_indices);
// 	for(Size i = 0; i < solution_indices.size();++i)
// 		{
// 			std::cout << filtered_features[solution_indices[i]].getRT() <<  " "
// 								<< filtered_features[solution_indices[i]].getMZ() <<  "\n";
// 			optimal_set.push_back(filtered_features[solution_indices[i]]);
// 		}
	
}



