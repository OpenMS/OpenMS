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

#ifndef OPENMS_ANALYSIS_ID_OFFLINEPRECURSORIONSELECTOR_H
#define OPENMS_ANALYSIS_ID_OFFLINEPRECURSORIONSELECTOR_H


#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/ANALYSIS/ID/ILPWrapper.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>

namespace OpenMS
{
	class PeptideIdentification;
	class ProteinIdentification;
	class String;


	/**
		 @brief Implements different algorithms for precursor ion selection

		 Implements different algorithms for precursor ion selection,
		 either based on a whole FeatureMap (e.g. like with LC-MALDI MS data)
		  or based on single scans (e.g. with LC-ESI MS data).
			
		 @htmlinclude OpenMS_OfflinePrecursorIonSelection.parameters
  */
  class OPENMS_DLLAPI OfflinePrecursorIonSelection: public DefaultParamHandler
  {
  public:
		typedef ILPWrapper::IndexTriple IndexTriple;

    OfflinePrecursorIonSelection();
    virtual ~OfflinePrecursorIonSelection();

		/**
			 @brief Determines the minimal set of features needed to obtain all
			 protein identifications.

			 

		 */
    void computeOptimalSolution(std::vector<ProteinIdentification>& prot_ids,
																std::vector<PeptideIdentification>& pep_ids,
																FeatureMap<>& features,FeatureMap<>& optimal_set,bool filter);

		/**
			 @brief Makes the precursor selection for a given feature map, either feature or scan based.

			 

		 */
		template <typename InputPeakType>
		void makePrecursorSelectionForKnownLCMSMap(FeatureMap<>& features,MSExperiment< InputPeakType > & experiment,
																							 MSExperiment< InputPeakType > & ms2,std::set<Int>& charges_set,
																							 bool feature_based);

		/**
			 @brief Calculates the mass ranges for each feature and stores them as indices of the raw data.
			 
		*/
		template <typename InputPeakType>
		void getMassRanges(FeatureMap<>& features, MSExperiment<InputPeakType>& experiment,
											 std::vector<std::vector<std::pair<Size,Size> > > & indices);

		template <typename InputPeakType>
		void computeOptimalSolution(std::vector<ProteinIdentification>& prot_ids,
																std::vector<PeptideIdentification>& pep_ids,
																MSExperiment<InputPeakType>& experiment,
																FeatureMap<>& features,
																FeatureMap<>& optimal_set,
																bool filter);
	
		
	private:
		/**
			 @brief Calculate the sum of intensities of relevant features for each scan separately.

		 */
		template <typename InputPeakType>
		void calculateXICs_(FeatureMap<> &features,std::vector<std::vector<std::pair<Size,Size> > >& mass_ranges,
												std::vector<std::vector<std::pair<Size,DoubleReal> > >& xics,MSExperiment<InputPeakType>& experiment,
												std::set<Int>& charges_set);

		std::vector<PeptideIdentification> filterPeptideIds_(std::vector<PeptideIdentification>& pep_ids);
		std::map<String,std::vector<Size> > protein_precursor_map_;
		std::vector<std::vector<Int> > rt_bins_;

		
		
  };

	template <typename InputPeakType>
	void OfflinePrecursorIonSelection::getMassRanges(FeatureMap<>& features, MSExperiment<InputPeakType>& experiment,
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
						typename MSSpectrum<InputPeakType>::Iterator mz_iter = experiment[rt].MZBegin(features[f].getMZ());
						typename MSSpectrum<InputPeakType>::Iterator mz_end = mz_iter;
						if(mz_iter == experiment[rt].end()) continue;
						// check to the left
						while(features[f].encloses(experiment[rt].getRT(),mz_iter->getMZ()))
							{
								start_found = true;
								start.first = rt;
								start.second = distance(experiment[rt].begin(),mz_iter);
								if(mz_iter == experiment[rt].begin()) break;
								--mz_iter;
							}
						// and now to the right
						while(mz_end != experiment[rt].end() && features[f].encloses(experiment[rt].getRT(),mz_end->getMZ()))
							{
								end_found = true;
								end.first = rt;
								end.second = distance(experiment[rt].begin(),mz_end);
								++mz_end;
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
						std::cout << vec.size() << " / 2 scans"<<std::endl;
						for(Size i = 0; i < vec.size(); i+=2)
							{
								std::cout << "Feature "<< f<< " RT : "<<vec[i].first 
													<< " MZ : "<<experiment[vec[i].first][vec[i].second].getMZ() << " "
													<< experiment[vec[i+1].first][vec[i+1].second].getMZ() << std::endl;
							}
					}
#endif
				indices.push_back(vec);
			}
	}

	template <typename InputPeakType>
	void OfflinePrecursorIonSelection::computeOptimalSolution(std::vector<ProteinIdentification>& prot_ids,
																														std::vector<PeptideIdentification>& pep_ids,
																														MSExperiment<InputPeakType>& experiment,
																														FeatureMap<>& features,
																														FeatureMap<>& optimal_set,
																														bool /*filter*/)
	{
		std::vector<PeptideIdentification> filtered_pep_ids = filterPeptideIds_(pep_ids);

		std::cout << "filtered"<<std::endl;
		// first map the ids onto the features
		IDMapper mapper;
		mapper.annotate(features,filtered_pep_ids,prot_ids);
		std::cout << "mapped"<<std::endl;
		// get the mass ranges for each features for each scan it occurs in
		std::vector<std::vector<std::pair<Size,Size> > >  indices;
		getMassRanges(features,experiment,indices);
		std::cout << "got mass ranges"<<std::endl;

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
		//protein_precursor_map_.clear();
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
		std::cout << "created protein_precursor_map"<<std::endl;
		std::vector<IndexTriple> variable_indices;
		// build model
		ILPWrapper ilp_wrapper;
		ilp_wrapper.encodeModelForOptimalSolution(features,experiment,indices,protein_precursor_map_,variable_indices,
																							param_.getValue("ms2_spectra_per_rt_bin"));

		std::cout << "encoded problem"<<std::endl;
		// compute solution
		std::vector<int> solution_indices;
		ilp_wrapper.solve(solution_indices);
		sort(variable_indices.begin(),variable_indices.end(),ILPWrapper::VariableIndexLess());
		std::cout << solution_indices.size() <<" features in solution"<<std::endl;
		std::map<String,std::vector<Size> > prot_feat_map;
		FeatureMap<> map_out;
		map_out.setProteinIdentifications(features.getProteinIdentifications());
		for(Size i = 0; i < solution_indices.size();++i)
			{	
				Int index = variable_indices[solution_indices[i]].variable;
				std::cout << features[index].getRT() <<  " " << variable_indices[solution_indices[i]].scan<< " " 
									<< features[index].getMZ() <<	"\n";
				optimal_set.push_back(features[index]);
				map_out.push_back(features[index]);
				if(features[index].getPeptideIdentifications().size()>0)
					{
						if(features[index].getPeptideIdentifications()[0].getHits().size() > 0)
							{
								std::vector<String> accs = features[index].getPeptideIdentifications()[0].getHits()[0].getProteinAccessions();
								if(accs.size() >0) std::cout << "more than 1 protein hit for this peptide" << std::endl;
								if(features[index].getPeptideIdentifications()[0].getHits().size()>1) std::cout << "more than 1 peptide hit"<<std::endl;
							}
						else std::cout << "ATTENTION: feature without pep id in optimal solution"<<std::endl;
					}
			}
		FeatureXMLFile().store("/home/zerck/data/presentations/poster/ismb09/bruker_optimal_features_5msmsperspot.featureXML",map_out);
	}


	template <typename InputPeakType>
	void OfflinePrecursorIonSelection::calculateXICs_(FeatureMap<> &features,
																										std::vector<std::vector<std::pair<Size,Size> > >& mass_ranges,
																										std::vector<std::vector<std::pair<Size,DoubleReal> > >& xics,
																										MSExperiment<InputPeakType>& experiment,
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
								// 							std::cout <<"exp["<< mass_ranges[f][s].first << " "<<j<<"]="<<std::endl;
								// 							std::cout << experiment[mass_ranges[f][s].first][j].getIntensity()<<std::endl;
								weight += experiment[mass_ranges[f][s].first][j].getIntensity();
							}
						// enter xic in the vector for scan s/2
						xics[s/2].push_back(std::make_pair(f,weight));
					}
			}

		for(Size s = 0; s < xics.size(); ++s)
			{
				sort(xics[s].begin(),xics[s].end(),PairComparatorSecondElement<std::pair<Size,DoubleReal> >());
			}
	}


	template <typename InputPeakType>
	void OfflinePrecursorIonSelection::makePrecursorSelectionForKnownLCMSMap(FeatureMap<>& features,
																																					 MSExperiment< InputPeakType > & experiment,
																																					 MSExperiment< InputPeakType > & ms2,
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
				sort(variable_indices.begin(),variable_indices.end(),ILPWrapper::IndexLess());
			
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
						typename MSExperiment<InputPeakType>::iterator scan = experiment.begin()+feature_scan_idx;
						typename MSExperiment<InputPeakType>::SpectrumType ms2_spec;
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
						for(Size j = 0; j < xics[i].size() && j < max_spec; ++j)
							{
								typename MSExperiment<InputPeakType>::SpectrumType ms2_spec;
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

	
}

#endif //  OPENMS_ANALYSIS_ID_OFFLINEPRECURSORIONSELECTION_H
