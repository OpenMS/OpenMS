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
		void makePrecursorSelectionForKnownLCMSMap(FeatureMap<>& features,MSExperiment< Peak1D > & experiment,
																							 MSExperiment< Peak1D > & ms2,std::set<Int>& charges_set,
																							 bool feature_based);

		/**
			 @brief Calculates the mass ranges for each feature and stores them as indices of the raw data.
			 
		*/
		void getMassRanges(FeatureMap<>& features, MSExperiment<>& experiment,
											 std::vector<std::vector<std::pair<Size,Size> > > & indices);

		void computeOptimalSolution(std::vector<ProteinIdentification>& prot_ids,
																std::vector<PeptideIdentification>& pep_ids,
																MSExperiment<>& experiment,
																FeatureMap<>& features,
																FeatureMap<>& optimal_set,
																bool filter);
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
		/**
			 @brief Calculate the sum of intensities of relevant features for each scan separately.

		 */
		void calculateXICs_(FeatureMap<> &features,std::vector<std::vector<std::pair<Size,Size> > >& mass_ranges,
												std::vector<std::vector<std::pair<Size,DoubleReal> > >& xics,MSExperiment<>& experiment,
												std::set<Int>& charges_set);

		std::vector<PeptideIdentification> filterPeptideIds_(std::vector<PeptideIdentification>& pep_ids);
		std::map<String,std::vector<Size> > protein_precursor_map_;
		std::vector<std::vector<Int> > rt_bins_;

		
		
  };
}

#endif //  OPENMS_ANALYSIS_ID_OFFLINEPRECURSORIONSELECTION_H
