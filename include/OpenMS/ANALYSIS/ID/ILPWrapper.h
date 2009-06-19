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

#include <OpenMS/ANALYSIS/ID/OfflinePrecursorIonSelection.h>
#include "coin/CoinModel.hpp"

namespace OpenMS
{


  class OPENMS_DLLAPI ILPWrapper : public DefaultParamHandler
  { 
   
  public:
		typedef OfflinePrecursorIonSelection::IndexTriple IndexTriple;

		ILPWrapper();
    virtual ~ILPWrapper();

    /**
     *	@brief Encode ILP formulation for a given LC-MS map, but unknown protein sample.
     *	
     *	@param features FeatureMap with all possible precursors
		 *  @param 
		 *  @param 
     */
    void encodeModelForKnownLCMSMapFeatureBased(FeatureMap<>& features, MSExperiment<>& experiment,
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
																			 UInt ms2_spectra_per_rt_bin);
    
    /**
     *	@brief Solve the ILP.
     *	
     */
    void solve(std::vector<int>& solution_indices);
	private:
		CoinModel model_;

		void getXIC_(std::vector<std::pair<Size,Size> >& end_points,
								 std::vector<DoubleReal>& weights,MSExperiment<>& experiment,bool normalize);

		
  }; 

} // namespace

#endif // OPENMS_ANALYSIS_ID_ILPWRAPPER_H
