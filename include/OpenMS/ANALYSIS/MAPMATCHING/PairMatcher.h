// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_PAIRMATCHER_H
#define OPENMS_ANALYSIS_MAPMATCHING_PAIRMATCHER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/ElementPair.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <cmath>

namespace OpenMS
{
	/**
	 	@brief The PairMatcher allows the matching of labeled features (features with a fixed distance).
		
		Finds feature pairs that have a defined distance in RT and m/z in the same map.
	
		@ref PairMatcher_Parameters are explained on a separate page.

		@ingroup FeatureGrouping
	*/
	class PairMatcher
		: public DefaultParamHandler
	{
		public:
			
			/** @name Type definitions
			*/
			//@{	
			///
      typedef std::vector< ElementPair<Feature> > PairVectorType;
			//@}
			
			/// Default constructor
			PairMatcher();
	
			/// Destructor
			inline virtual ~PairMatcher()
			{
			}
	
			/** 
				@brief Pairing step of the PairMatcher
	
				Return pairs of features that have the same charge and a distance lying within a user-defined range.
				In order to get unique pairs (each feature is contained in at most one pair) use getBestPairs().
			*/
			const PairVectorType& run(const FeatureMap<>& features);
	
			/** 
				@brief Matching step of the PairMatcher
	
				Greedy 2-approximation to extract a set of pairs so that each feature is contained in at most one pair.
			*/
			inline const PairVectorType& getBestPairs()		
			{
				return best_pairs_;
			}
	
			/** @brief Print informations about the pair vector @p pairs to stream @p out
	
				 Print informations (quality, ratio, charge, feature positions, ...)
				 about the pair vector @p pairs to stream @p out
			*/
			static void printInfo(std::ostream& out, const PairVectorType& pairs);

		protected:
			
	    /// Square root of two
	    static const double sqrt2_;
	
			/// all possible pairs (after Pairing)
			PairVectorType pairs_;
	
			/// only the best pairs, no ambiguities (after Matching)
			PairVectorType best_pairs_;
	
			/// Compare to pairs by comparing their qualities
			struct Comparator
			{
				bool operator()(const ElementPair<Feature>& a, const ElementPair<Feature>& b)
				{
					return a.getQuality() > b.getQuality();
				}
			};
	
			/// return the p-value at position x for the bi-Gaussian distribution
			/// with mean @p m and standard deviation @p sig1 (left) and @p sig2 (right)
			inline double PValue_(double x, double m, double sig1, double sig2)
			{
				if (m<x)
				{
					return 1-erf((x-m)/sig2/sqrt2_);
				}
				else
				{
					return 1-erf((m-x)/sig1/sqrt2_);
				}
			}
		
		private:
			/// Copy constructor not implemented => private
			PairMatcher(const PairMatcher& source);

			/// Assignment operator not implemented => private
	    PairMatcher& operator=(const PairMatcher& source);

	}; // end of class PairMatcher

} // end of namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHER_PAIRMATCHER_H
