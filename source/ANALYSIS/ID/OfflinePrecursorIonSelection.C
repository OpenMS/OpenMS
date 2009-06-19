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
#include <OpenMS/ANALYSIS/ID/ILPWrapper.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
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

void OfflinePrecursorIonSelection::getMassRanges(FeatureMap<>& features, MSExperiment<>& experiment,
																								 std::vector<std::vector<std::pair<Size,Size> > > & indices)
{
	for(Size f = 0; f < features.size();++f)
		{
			std::vector<std::pair<Size,Size> > vec;

			for(Size rt = 0; rt < experiment.size();++rt)
				{
					// is scan relevant?
					if(!features[f].encloses(experiment[rt].getRT(),features[f].getMZ()))
						{
							continue;
						}
					std::pair<Size,Size> start;
					std::pair<Size,Size> end;
					bool start_found = false;
					bool end_found = false;
					for(Size mz = 0; mz < experiment[rt].size();++mz)
						{
							if(!start_found && features[f].encloses(experiment[rt].getRT(),experiment[rt][mz].getMZ()))
								{
									start_found = true;
									start.first = rt;
									start.second = mz;
								}
							if(start_found && features[f].encloses(experiment[rt].getRT(),experiment[rt][mz].getMZ()))
								{
									end_found = true;
									end.first = rt;
									end.second = mz;
								}
						}
					if(start_found && end_found)
						{
							vec.push_back(start);
							vec.push_back(end);
						}
#ifdef DEBUG_OPS
					else
						{
							std::cout << "start "<<start_found<<" end "<<end_found<<std::endl;
							std::cout << "feature: "<<f << " rt: "<<rt<<std::endl;
						}
#endif
				}
#ifdef DEBUG_OPS
			if(vec.size()>0)
				{
					std::cout << "Feature "<< f<< " RT von: "<<vec[0].first << " "<<vec[vec.size()-1].first
										<< " MZ : "<<experiment[vec[0].first][vec[0].second].getMZ() << " "
										<< experiment[vec[1].first][vec[1].second].getMZ() << std::endl;
				}
#endif
			indices.push_back(vec);
		}
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
																													MSExperiment<>& experiment,
																													FeatureMap<>& features,
																													FeatureMap<>& optimal_set,
																													bool /*filter*/)
{
	std::vector<PeptideIdentification> filtered_pep_ids = filterPeptideIds_(pep_ids);

	
  // first map the ids onto the features
	IDMapper mapper;
	mapper.annotate(features,filtered_pep_ids,prot_ids);

	// get the mass ranges for each features for each scan it occurs in
	std::vector<std::vector<std::pair<Size,Size> > >  indices;
	getMassRanges(features,experiment,indices);
	

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
							std::vector<Size> vec;
							protein_precursor_map_.insert(make_pair(prot_ids[p_id].getHits()[p_h].getAccession(),vec));
							++index;
						}
				}
		}
	
	// create protein_peptide map
	protein_precursor_map_.clear();
	// usually we one protein id with many hits
	//	protein_precursor_map_.resize(accessions.size());
	

	for(Size i=0; i<features.size();++i) // not really elegant
		{
			// filter for features with an id
			if(features[i].getPeptideIdentifications().size()>0)
				{
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
													// store feature index
													protein_precursor_map_[accs[a]].push_back(i);
												}
										}
								}
						}
				}
		}

	// build model
	ILPWrapper ilp_wrapper;
  ilp_wrapper.encodeModelForOptimalSolution(features,experiment,indices,protein_precursor_map_,
																						param_.getValue("ms2_spectra_per_rt_bin"));


	// compute solution
	std::vector<int> solution_indices;
  ilp_wrapper.solve(solution_indices);
	for(Size i = 0; i < solution_indices.size();++i)
		{
			std::cout << features[solution_indices[i]].getRT() <<  " "
								<< features[solution_indices[i]].getMZ() <<  "\n";
			optimal_set.push_back(features[solution_indices[i]]);
		}
	
}


void OfflinePrecursorIonSelection::calculateXICs_(FeatureMap<> &features,
																									std::vector<std::vector<std::pair<Size,Size> > >& mass_ranges,
																									std::vector<std::vector<std::pair<Size,DoubleReal> > >& xics,
																									MSExperiment<>& experiment,
																									std::set<Int>& charges_set)
{
	xics.clear();
	xics.resize(experiment.size());
	// for each feature
	for(Size f = 0; f < mass_ranges.size();++f)
		{
			// is charge valid
			if(charges_set.count(features[f].getCharge()) < 1)
				{
					continue;
				}
			// go through all scans where the feature occurs
			for(Size s = 0; s < mass_ranges[f].size();s+=2)
				{
					// sum intensity over all raw datapoints belonging to the feature in the current scan
					DoubleReal weight = 0.;
					for(Size j = mass_ranges[f][s].second;j <= mass_ranges[f][s+1].second;++j)
						{
							weight += experiment[mass_ranges[f][s].first][j].getIntensity();
						}
					// enter xic in the vector for scan s
					xics[s].push_back(std::make_pair(f,weight));
				}
		}

	for(Size s = 0; s < xics.size(); ++s)
		{
			sort(xics[s].begin(),xics[s].end(),PairComparatorSecondElement<std::pair<Size,DoubleReal> >());
		}
}


void OfflinePrecursorIonSelection::makePrecursorSelectionForKnownLCMSMap(FeatureMap<>& features,
																																				 MSExperiment< Peak1D > & experiment,
																																				 MSExperiment< Peak1D > & ms2,
																																				 std::set<Int>& charges_set,
																																				 bool feature_based)
{
	
	// get the mass ranges for each features for each scan it occurs in
	std::vector<std::vector<std::pair<Size,Size> > >  indices;
	getMassRanges(features,experiment,indices);

	// feature based selection (e.g. with LC-MALDI)
	if(feature_based)
		{
			// create ILP
			ILPWrapper ilp_wrapper;
			
			std::vector<IndexTriple> variable_indices;
			ilp_wrapper.encodeModelForKnownLCMSMapFeatureBased(features, experiment,variable_indices,
																												 indices,charges_set,
																												 param_.getValue("ms2_spectra_per_rt_bin"));
			sort(variable_indices.begin(),variable_indices.end(),OfflinePrecursorIonSelection::IndexLess());
			
			// solve it
			std::vector<int> solution_indices;
			ilp_wrapper.solve(solution_indices);
			
			std::cout << "best_solution "<<std::endl;
			// print best solution
			// create inclusion list
			for(Size i = 0; i < solution_indices.size();++i)
				{
					Size feature_index = variable_indices[solution_indices[i]].feature;
					Size feature_scan_idx = variable_indices[solution_indices[i]].scan;
					MSExperiment<Peak1D>::iterator scan = experiment.begin()+feature_scan_idx;
					MSExperiment<Peak1D>::SpectrumType ms2_spec;
					Precursor p;
					std::vector< Precursor > pcs;
					p.setIntensity(features[feature_index].getIntensity());
					p.setMZ(features[feature_index].getMZ());
					p.setCharge(features[feature_index].getCharge());
					pcs.push_back(p);
					ms2_spec.setPrecursors(pcs);
					ms2_spec.setRT(scan->getRT());
					ms2.push_back(ms2_spec);
					// link ms2 spectrum with features overlapping its precursor
					// Warning: this depends on the current order of features in the map
					// Attention: make sure to name ALL features that overlap, not only one!
					ms2.setMetaValue("parent_feature_ids", IntList::create(String(feature_index)));
					std::cout << " MS2 spectra generated at: " << scan->getRT() << " x " << p.getMZ() << "\n";
					
				}
			std::cout << solution_indices.size() << " out of " << features.size()
								<< " precursors are in best solution.\n";
			
		}
	else // scan based selection (take the x highest signals for each spectrum)
		{
#ifdef DEBUG_OPS
			std::cout << "scan based precursor selection"<<std::endl;
#endif
			// if the highest signals for each scan shall be selected we don't need an ILP formulation
			std::vector<std::vector<std::pair<Size,DoubleReal> > > xics;			
			calculateXICs_(features,indices,xics,experiment,charges_set);

			Size max_spec = (Int)param_.getValue("ms2_spectra_per_rt_bin");
			// get best x signals for each scan
			for(Size i = 0; i < experiment.size();++i)
				{
#ifdef DEBUG_OPS
					std::cout << "scan "<<experiment[i].getRT() << ":";
#endif
					for(Size j = 0; j < xics[i].size() && j <= max_spec; ++j)
						{
							MSExperiment<Peak1D>::SpectrumType ms2_spec;
							Precursor p;
							std::vector< Precursor > pcs;
							p.setIntensity(features[(xics[i].end()-1-j)->first].getIntensity());
							p.setMZ(features[(xics[i].end()-1-j)->first].getMZ());
							p.setCharge(features[(xics[i].end()-1-j)->first].getCharge());
							pcs.push_back(p);
							ms2_spec.setPrecursors(pcs);
							ms2_spec.setRT(experiment[i].getRT());
							ms2.push_back(ms2_spec);
							// link ms2 spectrum with features overlapping its precursor
							// Warning: this depends on the current order of features in the map
							// Attention: make sure to name ALL features that overlap, not only one!
							ms2.setMetaValue("parent_feature_ids", IntList::create(String((xics[i].end()-1-j)->first)));
#ifdef DEBUG_OPS
							std::cout << " MS2 spectra generated at: " << experiment[i].getRT() << " x " << p.getMZ()
												<< " int: "<<(xics[i].end()-1-j)->second<< "\n";
#endif
						}
				}
		}

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
	rt_bins_.resize(ceil((max_rt-min_rt)/rt_bin_span));

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
					rt_bins_[floor((features[i].getRT()-min_rt) / rt_bin_span)].push_back((Int)filtered_features.size()-1);

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



