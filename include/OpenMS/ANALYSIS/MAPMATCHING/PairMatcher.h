// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: PairMatcher.h,v 1.5 2006/06/09 08:19:41 j-joachim Exp $
// $Author: j-joachim $
// $Maintainer: Jens Joachim $
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
	 	@brief The PairMatcher allows the matching of labeled features with a fixed distance.

	  Parameters:
			<table>
			<tr><td>mz_pair_dist</td>
					<td>optimal pair distance in [Th] for features with charge +1</td></tr>
			<tr><td>rt_pair_dist</td>
					<td>optimal pair distance in [min]</td></tr>
			<tr><td>mz_stdev</td>
					<td>standard deviation from optimal m/z distance</td></tr>
			<tr><td>rt_stdev1</td>
					<td>standard deviation below optimal retention time distance</td></tr>
			<tr><td>rt_stdev2</td>
					<td>standard deviation above optimal retention time distance</td></tr>
			</table>
	 **/

	class PairMatcher: public FactoryProduct
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
			RT = DimensionDescription < DimensionDescriptionTagLCMS >::RT,
			MZ = DimensionDescription < DimensionDescriptionTagLCMS >::MZ
		};
		
		/// Constructor
		PairMatcher(FeatureMapType& features);

		/// Copy constructor
		PairMatcher(const PairMatcher& source);

		///  Assignment operator
    virtual PairMatcher& operator = (const PairMatcher& source);

		/// Destructor
		virtual ~PairMatcher();

		/** Pairing step of the PairMatcher

				Return pairs of features that have the same charge and a distance lying within a user-defined range.
		*/
		const PairVectorType& run();

		/** @brief Reduction step of the PairMatcher

				Remove ambiguous pairs by extracting the best pairs so that each feature is associated with only one pair.
		*/
		const PairVectorType& getBestPairs();


		/** @brief Convert pair vector into feature map 

				Convert pair vector into feature map for visualization in TOPPView.
				The pairing is shown as an octagon stored in the convex hull layer of
				the first feature of each pair.
				This feature also contains some of the pairs meta values (quality, ratio, ...).
		*/
		static void fillFeatureMap(FeatureMapType&, const PairVectorType&);

		/** @brief Print informations about the pair vector @p pairs to stream @p out

			 Print informations (quality, ratio, charge, feature positions, ...) 
			 about the pair vector @p pairs to stream @p out
		*/
		static void printInfo(std::ostream& out, const PairVectorType& pairs);

		static const String getName()
    {
      return "PairMatcher";
    }

		protected:
		static const double sqrt2;

		/// feature to be paired
		FeatureMapType& features_;

		/// all possible pairs
		PairVectorType pairs_;

		/// only the best pairs, no ambiguities
		PairVectorType best_pairs_;

		double mz_diff_, rt_min_, rt_max_, mz_pair_dist_, rt_pair_dist_;
		double rt_stdev1_, rt_stdev2_, mz_stdev_;

		/// Compare to pairs by comparing their qualities
		struct Comparator{
			bool operator()(const PairType& a, const PairType& b)
			{
				return a.getQuality() > b.getQuality();
			}
		};

		inline double p_value_(double x, double m, double sig1, double sig2)
		{
			if (m<x)
			{
				return 1-erf((x-m)/sig2/sqrt2);
			}
			else{
				return 1-erf((m-x)/sig1/sqrt2);
			}
		}

		private:
		enum Constants{ ID=11, LOW_QUALITY};

	}; // end of class PairMatcher

} // end of namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHER_PAIRMATCHER_H
