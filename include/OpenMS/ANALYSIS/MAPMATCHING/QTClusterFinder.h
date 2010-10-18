// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Steffen Sass $
// $Authors: Steffen Sass, Hendrik Weisser $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_MAPMATCHING_QTCLUSTERFINDER_H
#define OPENMS_ANALYSIS_MAPMATCHING_QTCLUSTERFINDER_H

#include <OpenMS/DATASTRUCTURES/HashGrid.h>
#include <OpenMS/DATASTRUCTURES/GridFeature.h>
#include <OpenMS/DATASTRUCTURES/QTCluster.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/BaseGroupFinder.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

namespace OpenMS 
{

/**
   @brief A variant of QT clustering for the detection of feature groups.

	 The algorithm accumulates all features from all input maps, then applies a variant of QT clustering to find groups of corresponding features. In more detail, every feature from every input map is considered as a potential cluster center. For every center, its nearest neighbors from the other input maps are detected and added to the potential cluster. Iteratively, the cluster with the highest quality is extracted and the clustering is updated.

	 <b>Properties affecting the grouping</b>

	 To be included in a particular cluster, a feature has to fulfill the following conditions:
	 @li distances in RT and m/z from the cluster center must be below user-defined thresholds (@p max_distance:RT and @p max_distance:MZ),
	 @li the charge state must match that of the cluster center,
	 @li if @p use_identifications is set and both the feature and the cluster center are annotated with peptide identifications, the identifications have to match.

	 Every cluster contains at most one feature from each input map - namely the feature closest to the cluster center that meets the criteria and does not belong to a better cluster.

	 The notion of "closeness" for features is defined by the following distance function, the parameters of which can be set by the user.
	 Let \f$\Delta_\textit{RT}\f$ and \f$\Delta_\textit{MZ}\f$ be the absolute values of the RT and m/z differences of two features. Then the distance between the features is:
   \f[
   \big(
   \frac
     {\Delta_\textit{RT}}
     {\textit{max\_distance:RT}}
	 \big)
   ^\textit{diff\_exponent:RT}
   +
	 \big(
   \frac
     {\Delta_\textit{MZ}}
     {\textit{max\_distance:MZ}}
   \big)
   ^\textit{diff\_exponent:MZ}
   \f]

   The parameters @p diff_exponent:RT and @p diff_exponent:MZ control the growth rate of the penalty for differences. By default, the absolute distance in RT (which should not be very susceptible to outliers) and the squared distance in m/z (where small difference occur frequently, but large differences indicate a mismatch) are used.

	 The quality of a cluster is computed from the number of elements in it and their distances to the cluster center. For more details see QTCluster.

	 <b>Optimization</b>

	 This algorithm includes a number of optimizations to reduce run-time:
	 @li two-dimensional hashing of features,
	 @li a look-up table for feature distances,
	 @li a variant of QT clustering that requires only one round of clustering.

	 @see FeatureGroupingAlgorithmQT

   @htmlinclude OpenMS_QTPairFinder.parameters

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

		/// Maximum distance in RT direction
		DoubleReal max_dist_rt_;

		/// Maximum distance in m/z direction
		DoubleReal max_dist_mz_;

		/// Distances in RT are raised to this power
		DoubleReal diff_exp_rt_;

		/// Distances in m/z are raised to this power
		DoubleReal diff_exp_mz_;

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
			 @brief calculates the distance between two features.
		 
			 Returns -1 if the two features have different charge.
		 */
		DoubleReal getDistance_(BaseFeature const& left, BaseFeature const& right) 
			const;

		/**
			 @brief calculates the distance of RT difference and m/z difference.
		 
			 Uses m/z and RT exponent values.
		 */
		DoubleReal getDistance_(DoubleReal pos_diff_rt, DoubleReal pos_diff_mz) 
			const;

		/**
			 @brief Checks whether the peptide IDs of a cluster and a neighboring feature are compatible.
			 
			 A neighboring feature without identification is always compatible. Otherwise, the cluster and feature are compatible if the best peptide hits of each of their identifications have the same sequences.

			 @note A cluster without identification is only compatible to features without identification, because otherwise features with different identifications could end up in the same cluster (we only compare to the cluster center).
		*/
		bool compatibleIDs_(const QTCluster& cluster, const GridFeature* neighbor) 
			const;
		
		/// Sets algorithm parameters
		void setParameters_();

		/// Generates a consensus feature from the best cluster and updates the clustering
		void makeConsensusFeature_(std::list<QTCluster>& clustering,
															 ConsensusFeature& feature);

		/// Computes an initial QT clustering of the points in the hash grid
		void computeClustering_(HashGrid& grid, std::list<QTCluster>& clustering);

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
