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

#ifndef OPENMS_ANALYSIS_ID_OFFLINEPRECURSORIONSELECTION_H
#define OPENMS_ANALYSIS_ID_OFFLINEPRECURSORIONSELECTION_H


#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/ANALYSIS/ID/ILPWrapper.h>

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
			 @brief Makes the precursor selection for a given feature map, either feature or scan based.

			 @param features Input feature map
			 @param experiment Input raw data
			 @param ms2 Precursors are added as empty MS2 spectra to this MSExperiment
			 @param charges_set Allowed charge states
			 @param feature_based If true the selection is feature based, if false it is scan based and the highest signals in each spectrum are chosen
		 */
		template <typename InputPeakType>
		void makePrecursorSelectionForKnownLCMSMap(FeatureMap<>& features,MSExperiment< InputPeakType > & experiment,
																							 MSExperiment< InputPeakType > & ms2,std::set<Int>& charges_set,
																							 bool feature_based);

		/**
			 @brief Calculates the mass ranges for each feature and stores them as indices of the raw data.
			 
			 @param features Input feature map
			 @param experiment Input raw data
			 @param indices The boundaries of the features as indices in the raw data
		*/
		template <typename InputPeakType>
		void getMassRanges(FeatureMap<>& features, MSExperiment<InputPeakType>& experiment,
											 std::vector<std::vector<std::pair<Size,Size> > > & indices);

	private:
		/**
			 @brief Calculate the sum of intensities of relevant features for each scan separately.

		 */
		template <typename InputPeakType>
		void calculateXICs_(FeatureMap<> &features,std::vector<std::vector<std::pair<Size,Size> > >& mass_ranges,
												std::vector<std::vector<std::pair<Size,DoubleReal> > >& xics,
												MSExperiment<InputPeakType>& experiment,
												std::set<Int>& charges_set);

		/**
			 @brief Eliminates overlapping peaks.

		 */
		template <typename InputPeakType>		
		void checkMassRanges_(std::vector<std::vector<std::pair<Size,Size> > >& mass_ranges,
													MSExperiment<InputPeakType>&experiment);

		void updateExclusionList_(std::vector<std::pair<Size,Size> >& exclusion_list);
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
		// eliminate nearby peaks
		if(param_.getValue("exclude_overlapping_peaks")=="true") checkMassRanges_(indices,experiment);
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
					weight += experiment[mass_ranges[f][s].first][j].getIntensity();
				}
				// enter xic in the vector for scan 
				xics[mass_ranges[f][s].first].push_back(std::make_pair(f,weight));
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
		DoubleReal rt_dist = 0.;
		if(experiment.size()>1)
		{
				rt_dist = experiment[1].getRT()-experiment[0].getRT();
		}
		
		// feature based selection (e.g. with LC-MALDI)
		if(feature_based)
			{
				// create ILP
				ILPWrapper ilp_wrapper;
			
				std::vector<IndexTriple> variable_indices;
				std::vector<int> solution_indices;
				ilp_wrapper.createAndSolveILPForKnownLCMSMapFeatureBased(features, experiment,variable_indices,
																																 indices,charges_set,
																																 param_.getValue("ms2_spectra_per_rt_bin"),
																																 solution_indices);

				sort(variable_indices.begin(),variable_indices.end(),ILPWrapper::IndexLess());
#ifdef DEBUG_OPS
				std::cout << "best_solution "<<std::endl;
#endif
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
						ms2_spec.setRT(scan->getRT()+rt_dist/2.);
						ms2_spec.setMSLevel(2);
						// link ms2 spectrum with features overlapping its precursor
						// Warning: this depends on the current order of features in the map
						// Attention: make sure to name ALL features that overlap, not only one!
						ms2_spec.setMetaValue("parent_feature_ids", IntList::create(String(feature_index)));
						ms2.push_back(ms2_spec);
						std::cout << " MS2 spectra generated at: " << scan->getRT() << " x " << p.getMZ() << "\n";
					
					}
#ifdef DEBUG_OPS
				std::cout << solution_indices.size() << " out of " << features.size()
									<< " precursors are in best solution.\n";
#endif
			}
		else // scan based selection (take the x highest signals for each spectrum)
			{
#ifdef DEBUG_OPS
				std::cout << "scan based precursor selection"<<std::endl;
#endif
				// if the highest signals for each scan shall be selected we don't need an ILP formulation
				std::vector<std::vector<std::pair<Size,DoubleReal> > > xics;			
				calculateXICs_(features,indices,xics,experiment,charges_set);

				bool dynamic_exclusion = param_.getValue("use_dynamic_exclusion") == "true" ? true : false;
				std::vector<std::pair<Size,Size> > exclusion_list;
				Size exclusion_specs = (Size)(floor((DoubleReal)param_.getValue("exclusion_time") /(DoubleReal) rt_dist));
				// get best x signals for each scan
				for(Size i = 0; i < experiment.size();++i)
					{
						Size max_spec = (Int)param_.getValue("ms2_spectra_per_rt_bin");
#ifdef DEBUG_OPS
						std::cout << "scan "<<experiment[i].getRT() << ":";
#endif
						if(dynamic_exclusion) updateExclusionList_(exclusion_list);
						for(Size j = 0; j < xics[i].size() && j < max_spec; ++j)
							{
								if(dynamic_exclusion)
									{
										// check if feature is currently dynamically excluded
										bool exclude = false;
										for(Size excl_idx = 0; excl_idx < exclusion_list.size();++excl_idx)
											{
												if((xics[i].end()-1-j)->first == exclusion_list[excl_idx].first)
													{
														exclude = true;
														break;
													}
											}
										if(exclude)
											{
												++max_spec;
												continue;
											}
									}
								typename MSExperiment<InputPeakType>::SpectrumType ms2_spec;
								Precursor p;
								std::vector< Precursor > pcs;
								p.setIntensity(features[(xics[i].end()-1-j)->first].getIntensity());
								p.setMZ(features[(xics[i].end()-1-j)->first].getMZ());
								p.setCharge(features[(xics[i].end()-1-j)->first].getCharge());
								pcs.push_back(p);
								ms2_spec.setPrecursors(pcs);
								ms2_spec.setMSLevel(2);
								ms2_spec.setRT(experiment[i].getRT() + (j+1)*rt_dist/(max_spec+1) );
								ms2_spec.setMetaValue("parent_feature_ids", IntList::create(String((xics[i].end()-1-j)->first)));
								ms2.push_back(ms2_spec);
								if(dynamic_exclusion)
									{
										//add feature to exclusion list
										std::pair<Size,Size> pair(std::make_pair((xics[i].end()-1-j)->first,exclusion_specs+1));
										exclusion_list.push_back(pair);
									}
								// link ms2 spectrum with features overlapping its precursor
								// Warning: this depends on the current order of features in the map
								// Attention: make sure to name ALL features that overlap, not only one!
#ifdef DEBUG_OPS
								std::cout << " MS2 spectra generated at: " << experiment[i].getRT() << " x " << p.getMZ()
													<< " int: "<<(xics[i].end()-1-j)->second<< "\n";
#endif
							}
					}
			}

	}

	template <typename InputPeakType>		
	void OfflinePrecursorIonSelection::checkMassRanges_(std::vector<std::vector<std::pair<Size,Size> > >& mass_ranges,
																											MSExperiment<InputPeakType>&experiment)
	{
		std::vector<std::vector<std::pair<Size,Size> > > checked_mass_ranges;
		DoubleReal min_peak_distance = param_.getValue("min_peak_distance");
		checked_mass_ranges.reserve(mass_ranges.size());
		for(Size f = 0; f < mass_ranges.size();++f)
			{
				std::vector<std::pair<Size,Size> > checked_mass_ranges_f;
				for(Size s_idx = 0; s_idx < mass_ranges[f].size();s_idx+=2)
					{
						Size s = mass_ranges[f][s_idx].first;
						bool overlapping_features = false;
						////////////////////////////////////////////////////////////////////////
						// check if other features overlap with this feature in the current scan
						////////////////////////////////////////////////////////////////////////
						InputPeakType & peak_left_border = experiment[s][mass_ranges[f][s_idx].second];
						InputPeakType & peak_right_border = experiment[s][mass_ranges[f][s_idx+1].second];
						for(Size fmr=0; fmr < mass_ranges.size();++fmr)
							{
								if(fmr == f) continue;
								for(Size mr=0; mr < mass_ranges[fmr].size();mr+=2)
									{
										if( mass_ranges[fmr][mr].first ==  s) // same spectrum
											{
												InputPeakType & tmp_peak_left = experiment[s][mass_ranges[fmr][mr].second];
												InputPeakType & tmp_peak_right = experiment[s][mass_ranges[fmr][mr+1].second];
#ifdef DEBUG_OPS
												std::cout << tmp_peak_left.getMZ() << " < "
																	<< peak_left_border.getMZ()-min_peak_distance << " && "
																	<< tmp_peak_right.getMZ() << " < "
																	<< peak_left_border.getMZ()- min_peak_distance<< " ? "
																	<< (tmp_peak_left.getMZ() < peak_left_border.getMZ()-min_peak_distance &&
																			tmp_peak_right.getMZ() < peak_left_border.getMZ()-min_peak_distance)
																	<<" || "
																	<< tmp_peak_left.getMZ() << " > "
																	<< peak_right_border.getMZ()+min_peak_distance << " && "
																	<< tmp_peak_right.getMZ() << " > "
																	<< peak_right_border.getMZ()+ min_peak_distance<< " ? "
																	<< (tmp_peak_left.getMZ() > peak_right_border.getMZ()+min_peak_distance &&
																			tmp_peak_right.getMZ() > peak_right_border.getMZ()+min_peak_distance)
																	<< std::endl;
#endif
												// all other features have to be either completely left or
												// right of the current feature
												if(!((tmp_peak_left.getMZ() < peak_left_border.getMZ()-min_peak_distance &&
															tmp_peak_right.getMZ() < peak_left_border.getMZ()-min_peak_distance) ||
														 (tmp_peak_left.getMZ() > peak_right_border.getMZ()+min_peak_distance &&
															tmp_peak_right.getMZ() > peak_right_border.getMZ()+min_peak_distance)))
													{
#ifdef DEBUG_OPS
														std::cout << "found overlapping peak"<<std::endl;
#endif
														overlapping_features = true;
														break;
													}
											}
									}
							}
						if(!overlapping_features)
							{
#ifdef DEBUG_OPS
								std::cout << "feature in spec ok" << mass_ranges[f][s_idx].second << " in spec "
													<< mass_ranges[f][s_idx].first << std::endl;
#endif
								checked_mass_ranges_f.insert(checked_mass_ranges_f.end(),
																						 mass_ranges[f].begin()+s_idx,
																						 mass_ranges[f].begin()+s_idx+2);
							}
					}
				checked_mass_ranges.push_back(checked_mass_ranges_f);
			}
		mass_ranges.swap(checked_mass_ranges);
	}
	
}

#endif //  OPENMS_ANALYSIS_ID_OFFLINEPRECURSORIONSELECTION_H
