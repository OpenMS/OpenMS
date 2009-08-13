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

	private:
		/**
			 @brief Calculate the sum of intensities of relevant features for each scan separately.

		 */
		template <typename InputPeakType>
		void calculateXICs_(FeatureMap<> &features,std::vector<std::vector<std::pair<Size,Size> > >& mass_ranges,
												std::vector<std::vector<std::pair<Size,DoubleReal> > >& xics,MSExperiment<InputPeakType>& experiment,
												std::set<Int>& charges_set);

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
				std::vector<int> solution_indices;
				ilp_wrapper.createAndSolveILPForKnownLCMSMapFeatureBased(features, experiment,variable_indices,
																																 indices,charges_set,
																																 param_.getValue("ms2_spectra_per_rt_bin"),
																																 param_.getValue("min_peak_distance"),
																																 solution_indices);

				sort(variable_indices.begin(),variable_indices.end(),ILPWrapper::IndexLess());
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
