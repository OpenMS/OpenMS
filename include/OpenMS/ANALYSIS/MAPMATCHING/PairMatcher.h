// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/ANALYSIS/MAPMATCHING/DFeaturePairVector.h>
#include <OpenMS/KERNEL/DimensionDescription.h>
#include <OpenMS/KERNEL/DFeatureMap.h>
#include <OpenMS/KERNEL/KernelTraits.h>
#include <OpenMS/CONCEPT/FactoryProduct.h>
#include <OpenMS/DATASTRUCTURES/QuadTree.h>

#include <cmath>

namespace OpenMS
{
	/**
	 	@brief The PairMatcher allows the matching of labeled features (features with a fixed distance)
		as described in the thesis "Automated LC-MS data analysis for differential quantification
		of MHC ligands using stable isotope labeling".
		
		@note The mz_pair_dist is added while the rt_pair_dist is substracted when searching pairs.
			    This is due to the fact, that the heavier peptide normally elutes earlier!

	  Parameters:
			<table>
			<tr><td>mz_pair_dist</td>
					<td>optimal pair distance in [Th] for features with charge +1</td></tr>
			<tr><td>rt_pair_dist</td>
					<td>optimal pair distance in [min]</td></tr>
			<tr><td>mz_stdev</td>
					<td>standard deviation from optimal m/z distance</td></tr>
			<tr><td>rt_stdev_low</td>
					<td>standard deviation below optimal retention time distance</td></tr>
			<tr><td>rt_stdev_high</td>
					<td>standard deviation above optimal retention time distance</td></tr>
			</table>
	*/
	class PairMatcher
		: public FactoryProduct
	{
		public:
			
			/** @name Type definitions
			*/
			//@{	
			///
			typedef DFeature<2> FeatureType;
			typedef DFeatureMap<2> FeatureMapType;
			typedef DFeaturePair<2> PairType;
			typedef DFeaturePairVector<2> PairVectorType;
			typedef QuadTree< KernelTraits, FeatureType > QuadTreeType;
			//@}
	
		  /// Defines the coordinates of peaks / features.
			enum DimensionId
			{
				RT = DimensionDescription < LCMS_Tag >::RT,
				MZ = DimensionDescription < LCMS_Tag >::MZ
			};
			
			/// Constructor
			PairMatcher(FeatureMapType& features);
	
			/// Copy constructor
			PairMatcher(const PairMatcher& source);
	
			///  Assignment operator
	    virtual PairMatcher& operator = (const PairMatcher& source);
	
			/// Destructor
			virtual ~PairMatcher();
	
			/** @brief Pairing step of the PairMatcher
	
				Return pairs of features that have the same charge and a distance
				lying within a user-defined range.
			*/
			const PairVectorType& run();
	
			/** @brief Matching step of the PairMatcher
	
				Greedy 2-approximation to extract a set of pairs so that each feature
				is contained in at most one pair.
			*/
			const PairVectorType& getBestPairs();
	
			/** @brief Print informations about the pair vector @p pairs to stream @p out
	
				 Print informations (quality, ratio, charge, feature positions, ...)
				 about the pair vector @p pairs to stream @p out
			*/
			static void printInfo(std::ostream& out, const PairVectorType& pairs);
	
			static const String getProductName()
	    {
	      return "PairMatcher";
	    }

		protected:
	    /// Square root of two
	    static const double sqrt2_;
	
			/// features to be paired
			FeatureMapType& features_;
	
			/// all possible pairs (after Pairing)
			PairVectorType pairs_;
	
			/// only the best pairs, no ambiguities (after Matching)
			PairVectorType best_pairs_;
	
			/// Compare to pairs by comparing their qualities
			struct Comparator
			{
				bool operator()(const PairType& a, const PairType& b)
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
				else{
					return 1-erf((m-x)/sig1/sqrt2_);
				}
			}

		private:
			/// constants for accessing feature meta values
			enum Constants
			{
				ID=11,				/**<  used to assocate the feature with its index in the set */
				LOW_QUALITY		/**< used by the greedy approximation */
			};

	}; // end of class PairMatcher

} // end of namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHER_PAIRMATCHER_H
