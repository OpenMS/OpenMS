// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Hendrik Weisser $
// $Authors: Steffen Sass, Hendrik Weisser $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_MAPMATCHING_QTCLUSTERFINDER_H
#define OPENMS_ANALYSIS_MAPMATCHING_QTCLUSTERFINDER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/BaseGroupFinder.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/HashGridOld.h>
#include <OpenMS/DATASTRUCTURES/GridFeature.h>
#include <OpenMS/DATASTRUCTURES/QTCluster.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureDistance.h>

namespace OpenMS 
{

/**
   @brief A variant of QT clustering for the detection of feature groups.

	 The algorithm accumulates all features from all input maps, then applies a variant of QT clustering to find groups of corresponding features. In more detail, every feature from every input map is considered as a potential cluster center. For every center, its nearest neighbors from the other input maps are detected and added to the potential cluster. Iteratively, the cluster with the highest quality is extracted and the clustering is updated.

	 <b>Properties affecting the grouping</b>

	 To be included in a particular cluster, a feature has to fulfill the following conditions:
	 @li differences in RT and m/z from the cluster center must be below user-defined thresholds (@p distance_RT:max_difference and @p distance_MZ:max_difference),
	 @li the charge state must match that of the cluster center (unless @p ignore_charge is set),
	 @li if @p use_identifications is set and both the feature and the cluster center are annotated with peptide identifications, the identifications have to match.

	 Every cluster contains at most one feature from each input map - namely the feature closest to the cluster center that meets the criteria and does not belong to a better cluster.

	 The notion of "closeness" for features is defined by the distance function implemented in @ref FeatureDistance, the parameters of which can be set by the user.

	 The quality of a cluster is computed from the number of elements in it and their distances to the cluster center. For more details see QTCluster.

	 <b>Optimization</b>

	 This algorithm includes a number of optimizations to reduce run-time:
	 @li two-dimensional hashing of features,
	 @li a look-up table for feature distances,
	 @li a variant of QT clustering that requires only one round of clustering.

	 @see FeatureGroupingAlgorithmQT

   @htmlinclude OpenMS_QTClusterFinder.parameters

   @ingroup FeatureGrouping
*/

	class OPENMS_DLLAPI QTClusterFinder: public BaseGroupFinder
	{
	private:
		/// Distances between pairs of grid features
		typedef std::map<std::pair<GridFeature*, GridFeature*>, DoubleReal> 
			PairDistances;
	
		/// Number of input maps
		Size num_maps_;

		/// Consider peptide identifications for grouping?
		bool use_IDs_;

		/// Maximum RT difference
		DoubleReal max_diff_rt_;

		/// Maximum m/z difference
		DoubleReal max_diff_mz_;

		/// Feature distance functor
		FeatureDistance feature_distance_;

		/**
			 @brief Distance map.
		 
			 To compute it only once, the distance between two features is accessible by searching for a pair where the first position is the smaller pointer value.
		 */
		PairDistances distances_;

		/**
			 @brief Calculates the distance between two grid features.
		 
			 The distance is looked up in the distance map and only computed (and stored) if it's not already available.
		 */
		DoubleReal getDistance_(GridFeature* left, GridFeature* right);

		/**
			 @brief Checks whether the peptide IDs of a cluster and a neighboring feature are compatible.
			 
			 A neighboring feature without identification is always compatible. Otherwise, the cluster and feature are compatible if the best peptide hits of each of their identifications have the same sequences.
		*/
		bool compatibleIDs_(QTCluster& cluster, const GridFeature* neighbor);
		
		/// Sets algorithm parameters
		void setParameters_(DoubleReal max_intensity);

		/// Generates a consensus feature from the best cluster and updates the clustering
		void makeConsensusFeature_(std::list<QTCluster>& clustering,
															 ConsensusFeature& feature);

		/// Computes an initial QT clustering of the points in the hash grid
		void computeClustering_(HashGridOld& grid, std::list<QTCluster>& clustering);

		/// Runs the algorithm on feature maps or consensus maps
		template <typename MapType> void run_(const std::vector<MapType>& 
																					input_maps, ConsensusMap& result_map);

	protected:
		enum
		{
			RT = Peak2D::RT,
			MZ = Peak2D::MZ
		};

	public:
		/// Constructor
		QTClusterFinder();

		/// Destructor
		virtual ~QTClusterFinder();

		/// Returns the name of the product
		static const String getProductName()
		{
			return "qt";
		}

		/**
			 @brief Runs the algorithm on consensus maps
			 
			 @exception Exception::IllegalArgument is thrown if the input data is not valid.
		*/
		void run(const std::vector<ConsensusMap>& input_maps, 
						 ConsensusMap& result_map);

		/**
			 @brief Runs the algorithm on feature maps
			 
			 @exception Exception::IllegalArgument is thrown if the input data is not valid.
		*/
		void run(const std::vector<FeatureMap<> >& input_maps,
						 ConsensusMap& result_map);

		/// Returns an instance of this class
		static BaseGroupFinder*	create()
		{
			return new QTClusterFinder();
		}
	};

}

#endif /* OPENMS_ANALYSIS_MAPMATCHING_QTCLUSTERFINDER_H */
