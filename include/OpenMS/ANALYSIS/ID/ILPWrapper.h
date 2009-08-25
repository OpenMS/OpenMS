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
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>

class CoinModel;
namespace OpenMS
{

	/**
		 @brief Implements ILP formulation of precursor selection problems
		 
  */
  class OPENMS_DLLAPI ILPWrapper
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
		 *  @param experiment Input raw data
		 *  @param variable_indices Assignment of feature indices and ILP variables
		 *  @param mass_ranges Feature borders as indices in the raw data
		 *  @param charges_set Allowed charge states
		 *  @param ms2_spectra_per_rt_bin Allowed number of precursors per rt bin
		 *  @param solution_indices Indices of ILP variables that are in the optimal solution
     */
		template <typename InputPeakType>
    void createAndSolveILPForKnownLCMSMapFeatureBased(FeatureMap<>& features, MSExperiment<InputPeakType>& experiment,
																											std::vector<IndexTriple >& variable_indices,
																											std::vector<std::vector<std::pair<Size,Size> > > & mass_ranges,
																											std::set<Int>& charges_set,UInt ms2_spectra_per_rt_bin,
																											std::vector<int>& solution_indices);

		
    

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
		
		template <typename InputPeakType>
		void getXIC_(std::vector<std::pair<Size,Size> >& end_points,
								 std::vector<DoubleReal>& weights,MSExperiment<InputPeakType>& experiment,bool normalize);

		/**
     *	@brief Calculates the XICs for all features.
     *	
     */
		template <typename InputPeakType>
		void calculateXICs_(std::vector<std::vector<DoubleReal> >& xics,FeatureMap<>& features,
												MSExperiment<InputPeakType>& experiment,std::vector<std::vector<std::pair<Size,Size> > > & mass_ranges,
												bool normalize);
		
		/**
     *	@brief Creates and solves the ILP.
     *	
     */
    void createAndSolveILP_(FeatureMap<>& features,std::vector<std::vector<DoubleReal> >& intensity_weights,
														std::set<Int>& charges_set,std::vector<std::vector<std::pair<Size,Size> > > & mass_ranges,
														std::vector<IndexTriple>& variable_indices,std::vector<int>& solution_indices,
														UInt ms2_spectra_per_rt_bin,Size number_of_scans);

		
		/**
     *	@brief Solve the ILP.
     *	
     */
    void solveILP_(CoinModel& model,std::vector<int>& solution_indices);


		
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
	void ILPWrapper::calculateXICs_(std::vector<std::vector<DoubleReal> >& xics,FeatureMap<>& features,
																	MSExperiment<InputPeakType>& experiment,
																	std::vector<std::vector<std::pair<Size,Size> > > & mass_ranges,
																	bool normalize)
	{
		xics.clear();
		xics.resize(features.size());
		for(Size i = 0; i < features.size(); ++i)
			{
				getXIC_(mass_ranges[i],xics[i],experiment,normalize);
			}
	}

	template <typename InputPeakType>
	void ILPWrapper::createAndSolveILPForKnownLCMSMapFeatureBased(FeatureMap<>& features, MSExperiment<InputPeakType>& experiment,
																																std::vector<IndexTriple>& variable_indices,
																																std::vector<std::vector<std::pair<Size,Size> > > & mass_ranges,
																																std::set<Int>& charges_set,UInt ms2_spectra_per_rt_bin,
																																std::vector<int>& solution_indices)
	{

		std::vector<std::vector<DoubleReal> > intensity_weights;
		calculateXICs_(intensity_weights,features,experiment,mass_ranges,true);
#ifdef DEBUG_OPS
		std::cout << "got xics"<<std::endl;
#endif
		
		createAndSolveILP_(features,intensity_weights,charges_set,mass_ranges,variable_indices,solution_indices,
											 ms2_spectra_per_rt_bin,experiment.size());
	}
	


} // namespace

#endif // OPENMS_ANALYSIS_ID_ILPWRAPPER_H
