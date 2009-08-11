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

#ifndef OPENMS_ANALYSIS_ID_ILPWRAPPER_H
#define OPENMS_ANALYSIS_ID_ILPWRAPPER_H
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
//#include <OpenMS/ANALYSIS/ID/OfflinePrecursorIonSelection.h>
#include "coin/CoinModel.hpp"

namespace OpenMS
{


  class OPENMS_DLLAPI ILPWrapper : public DefaultParamHandler
  { 
   
  public:
	
		ILPWrapper();
    virtual ~ILPWrapper();

		/**
			 @brief Struct that holds the indices of the precursors in the feature map and the ilp formulation.
			 
		*/
		struct IndexTriple
		{
			Size feature;
			Size scan;
			Int variable;
			DoubleReal rt_probability;
			DoubleReal signal_weight;
		};

		
    /**
     *	@brief Encode ILP formulation for a given LC-MS map, but unknown protein sample.
     *	
     *	@param features FeatureMap with all possible precursors
		 *  @param 
		 *  @param 
     */
		template <typename InputPeakType>
    void encodeModelForKnownLCMSMapFeatureBased(FeatureMap<>& features, MSExperiment<InputPeakType>& experiment,
																								std::vector<IndexTriple >& variable_indices,
																								std::vector<std::vector<std::pair<Size,Size> > > & mass_ranges,
																								std::set<Int>& charges_set,UInt ms2_spectra_per_rt_bin);

		
    /**
     *	@brief Encode ILP formulation for a given LC-MS map and given ids to determine the optimal set of precursors.
     *	
     *	@param features FeatureMap with all possible precursors
		 *  @param 
		 *  @param protein_precursor_map Vector containing a vector with the precursors for each protein
     */
    void encodeModelForOptimalSolution(FeatureMap<>& features,MSExperiment<>& experiment,
																			 std::vector<std::vector<std::pair<Size,Size> > > & mass_ranges,
																			 std::map<String,std::vector<Size> >& protein_precursor_map,
																			 std::vector<IndexTriple>& variable_indices,
																			 UInt ms2_spectra_per_rt_bin);
    
    /**
     *	@brief Solve the ILP.
     *	
     */
    void solve(std::vector<int>& solution_indices);


		struct IndexLess
			: std::binary_function < IndexTriple , IndexTriple , bool >
		{
			inline bool operator () ( IndexTriple  const & left,
																IndexTriple const & right ) const
			{
				return ( left.variable < right.variable );
			}
		};

		
		struct ScanLess
			: std::binary_function < IndexTriple , IndexTriple , bool >
		{
			inline bool operator () ( IndexTriple  const & left,
																IndexTriple  const & right ) const
			{
				return ( left.scan < right.scan );
			}
		};
		
		struct VariableIndexLess
			: std::binary_function < IndexTriple , IndexTriple , bool >
		{
			inline bool operator () ( IndexTriple  const & left,
																IndexTriple  const & right ) const
			{
				return ( left.variable < right.variable );
			}
		};
		
	private:
		CoinModel model_;
		
		template <typename InputPeakType>
		void getXIC_(std::vector<std::pair<Size,Size> >& end_points,
								 std::vector<DoubleReal>& weights,MSExperiment<InputPeakType>& experiment,bool normalize);

		
  };

	template <typename InputPeakType>
	void ILPWrapper::getXIC_(std::vector<std::pair<Size,Size> >& end_points,
													 std::vector<DoubleReal>& weights,MSExperiment<InputPeakType>& experiment,bool normalize)
	{
		DoubleReal max_weight = 0.;
		weights.clear();
		for(Size i = 0; i < end_points.size();i+=2)
			{
				DoubleReal weight = 0.;
				for(Size j = end_points[i].second;j <= end_points[i+1].second;++j)
					{
						weight += experiment[end_points[i].first][j].getIntensity();
						//					std::cout << " add "<<experiment[end_points[i].first][j].getIntensity()<<std::endl;
					}
				if(weight > max_weight)  max_weight = weight;
				weights.push_back(weight);
			}

		if(normalize)
			{
				// normalize weights
				for(Size i = 0; i < weights.size();++i)
					{
#ifdef DEBUG_OPS
						if(end_points.size()>=i)
							{
								std::cout << "scan "<< end_points[i].first << " "<<weights[i] << " "<<max_weight
													<< " " << weights[i] / max_weight << std::endl;
							}
#endif
						weights[i] /= max_weight;
					}
			}
	}


	template <typename InputPeakType>
	void ILPWrapper::encodeModelForKnownLCMSMapFeatureBased(FeatureMap<>& features, MSExperiment<InputPeakType>& experiment,
																													std::vector<IndexTriple>& variable_indices,
																													std::vector<std::vector<std::pair<Size,Size> > > & mass_ranges,
																													std::set<Int>& charges_set,UInt ms2_spectra_per_rt_bin)
	{

		std::cout << "Feature Based: Build model: first objective"<<std::endl;
		///////////////////////////////////////////////////////////////////////
		// add objective function
		///////////////////////////////////////////////////////////////////////
		model_.setOptimizationDirection(-1); // maximize
		// max \sum_j x_jk * signal_jk
		//                    column_index, feature_index,scan
	
		//
		Int counter = 0;
		for(Size i = 0; i < features.size(); ++i)
			{
				// first check if charge state is allowed
				// charge not in "ChargeFilter" list
#ifdef DEBUG_OPS
				std::cout << "feat: "<<i <<" charge "<<features[i].getCharge() << std::endl;
#endif
				if (charges_set.count(features[i].getCharge())<1) continue;
				if(mass_ranges[i].size()==0) continue;
#ifdef DEBUG_OPS
				if(mass_ranges[i].size() > 0)
					{
						std::cout << "start_scan "<< mass_ranges[i][0].first << " ?= "<<features[i].getQuality(0)
											<< " stop scan "<< (mass_ranges[i].end()-1)->first<< " ?= "	<< features[i].getQuality(1)-1<<std::endl;
					}
#endif

				std::vector<DoubleReal> intensity_weights;
				getXIC_(mass_ranges[i],intensity_weights,experiment,true);
#ifdef DEBUG_OPS
				std::cout << "got xic"<<std::endl;
#endif

				// 			if(intensity_weights.size() != (mass_ranges[i].end()-1)->first - mass_ranges[i][0].first)
				// 				{
				// 					std::cout << "attention: "<<intensity_weights.size() << " != "
				// 										<<(mass_ranges[i].end()-1)->first - mass_ranges[i][0].first<<std::endl;
				// 				}
				Size c = 0;
				// go through all rts of the current feature
				for(Size s = mass_ranges[i][0].first; s <= (mass_ranges[i].end()-1)->first;++s)
					{ 
						model_.setColumnName(counter,(String("x_")+i+","+s).c_str());
#ifdef DEBUG_OPS
						std::cout << "add column "<<counter << std::endl;
#endif
						IndexTriple triple;
						triple.feature = i;
						triple.scan = s;
						triple.variable = counter;
						variable_indices.push_back(triple);
						model_.setColumnUpper(counter,1.);
						model_.setColumnLower(counter,0.);
						model_.setColumnIsInteger(counter,true);
					
#ifdef DEBUG_OPS	
						std::cout << "feat "<<i << " scan "<< s << " intensity_weight "
											<< intensity_weights[c] <<std::endl;
#endif
						model_.setObjective(counter,intensity_weights[c]);
						++counter;
						++c;
					}
			}
	
		///////////////////////////////////////////////////////////////////////
		// add constraints
		///////////////////////////////////////////////////////////////////////
		std::cout << "and now the constraints:"<<std::endl;

		///////////////////////////////////////////////////////////////////////
		// 1: ensure that each precursor is acquired maximally once
		///////////////////////////////////////////////////////////////////////
		std::cout << "first the number of times a precursors is acquired"<<std::endl;
		Size j = 0;
		for(Size i = 0; i < features.size();++i)
			{
				Size start = j;
				while(j < variable_indices.size() && variable_indices[j].feature == i)
					{
#ifdef DEBUG_OPS
						std::cout << j << " "<<variable_indices[j].variable << " "
											<< variable_indices[j].feature << " "
											<< variable_indices[j].scan<<std::endl;
#endif
						++j;
					}

				Size stop = j;
				double* entries = new double[stop-start];
				int* indices = new int[stop-start];
#ifdef DEBUG_OPS
				std::cout << "feature "<<i <<" "<<features[i].getMZ() <<" "<<features[i].getRT()<<" ";
				std::cout << stop-start<<"variables in equation\n";
#endif
				Size c = 0;
				for(Size k = start; k < stop; ++k)
					{
						entries[c] = 1.;
						indices[c] = variable_indices[k].variable;
						//					std::cout << j<<" "<<indices[j]<<std::endl;
						++c;
					}
#ifdef DEBUG_OPS
				std::cout << "\nadd row "<<std::endl;
#endif
				String name = "PREC_ACQU_LIMIT_" + String(i);
			
				model_.addRow((int)(stop-start),indices,entries,-COIN_DBL_MAX,1,name.c_str());
#ifdef DEBUG_OPS
				std::cout << stop-start << " "<<name<<std::endl;
				std::cout << "added row"<<std::endl;
#endif
				delete entries;
				delete indices;
			
			}



		///////////////////////////////////////////////////////////////////////
		// 2: do not exceed rt bin capacity
		///////////////////////////////////////////////////////////////////////
		std::cout << "and now the rt bin capacity"<<std::endl;
		std::cout << ms2_spectra_per_rt_bin << " rt bin capacity"<<std::endl;
		// sort variable_indices according to their scan number
		sort(variable_indices.begin(),variable_indices.end(),ScanLess());
		j = 0;
		for(Size i = 0; i < experiment.size();++i)
			{
				// first determine number of indices:
				Size start = j;
				while(j < variable_indices.size() && variable_indices[j].scan == i)
					{
						++j;
					}
				// no feature occuring in this scan
				if(start == j) continue;

				Size stop = j;
				Size c = 0;			
				double* entries = new double[stop-start];
				int* indices = new int[stop-start];
				for(Size s = start; s < stop; ++s)
					{
						entries[c] = 1.;
						indices[c] = variable_indices[s].variable;
						++c;
					}
#ifdef DEBUG_OPS
				std::cout << "\nadd row "<<std::endl;
#endif
				model_.addRow((int)(stop-start),indices,entries,-COIN_DBL_MAX,ms2_spectra_per_rt_bin,(String("RT_CAP")+i).c_str());
#ifdef DEBUG_OPS
				std::cout << "added row"<<std::endl;
#endif
				delete entries;
				delete indices;
			
			}



#ifdef DEBUG_OPS	
		model_.writeMps("/home/zerck/data/tmp/test_pis_problem.mps",0,0,2,true);
#endif

	}


} // namespace

#endif // OPENMS_ANALYSIS_ID_ILPWRAPPER_H
