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
// $Maintainer: Steffen Sass $
// $Authors: $
// --------------------------------------------------------------------------


#ifndef  OPENMS_COMPARISON_CLUSTERING_HASHCLUSTERING_H
#define  OPENMS_COMPARISON_CLUSTERING_HASHCLUSTERING_H

#include <OpenMS/DATASTRUCTURES/HashGrid.h>
#include <OpenMS/DATASTRUCTURES/DataSubset.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusteringMethod.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{
  class SILACTreeNode;
/**
		@brief Hierarchical clustering based on geometric hashing.

		Only elements in the surrounding are clustered.
		So it is not more necessary to calculate all distances, which leads to a linear runtime instead of quadratic.
		Instead of creating a hierarchical clustering tree of the whole data set, a vector of subtrees of this tree will be created containing only data points, which lie within the range defined by the two thresholds.

		@see ClusterHierarchical
		@ingroup SpectraClustering
	*/

class OPENMS_DLLAPI HashClustering : public ProgressLogger
{
  private:

	/**
	 * @brief current minimal distance
	 */
   DoubleReal min_distance;

	/**
	 * @brief two DataSubsets with current minimal distance
	 */
   std::pair<DataSubset*,DataSubset*> min_distance_subsets;

	/**
	 * @brief set of distances
	 */
   DistanceSet distances;

	/**
	 * @brief method which is used to calculate the distances
	 */
   ClusteringMethod* method;

	/**
	 * @brief the grid for geometric hashing
	 */
   HashGrid grid;

	/**
	 * @brief average silhoutte widths for each subtree
	 */
   std::vector<std::vector<Real> > silhouettes;

	/**
	 * @brief Calculates initial distances
	 */
   void init();

	/**
	 * @brief Merges two DataSubsets
	 */
   void merge();

	/**
	 * @brief Finds the two DataSubsets with minimal distance
	 * @param subset1 first DataSubset
	 * @param subset2 second DataSubset
	 */
   void updateMinElements();

	/**
   * @brief Calculates the distance of two DataSubsets using <i>getDistance</i> of the clustering method
   * @param subset1 first DataSubset
	 * @param subset2 second DataSubset
	 */
   DoubleReal getDistance(DataSubset& subset1, DataSubset& subset2);

	/**
	 * @brief Calculates the distance of two DataPoints using <i>getDistance</i> of the clustering method
	 * @param point1 first DataPoint
	 * @param point2 second DataPoint
	 */
   DoubleReal getDistance(DataPoint& point1, DataPoint& point2);

	/**
	 * @brief Calculates the silhouette values for any possible cluster number
	 * @param tree hierarchical clustering tree
	 */
   std::vector< Real > averageSilhouetteWidth(DataSubset& subset);

	/**
   *@brief Method to calculate a partition resulting from a certain step in clustering given by the number of clusters
   *@param cluster_quantity Size giving the number of clusters
   *@param tree vector of SILACTreeNodes representing the clustering
   *@param clusters vector of vectors holding the clusters
   *@see SILACTreeNode
    after call of this method the argument clusters is filled corresponding to the given @p cluster_quantity with the indices of the elements clustered
   */
   void cut_(const Size cluster_quantity, std::vector< std::vector<DataPoint*> >& clusters, const std::vector<SILACTreeNode>& tree);

public:

	/**
   *@brief Exception thrown if not enough data (<2) is used
		If the set of data to be clustered contains only one data point,
		clustering algorithms would fail for obvious reasons.
	 */
   class OPENMS_DLLAPI InsufficientInput : public Exception::BaseException
   {
     public:
       InsufficientInput(const char* file, int line, const char* function, const char* message= "not enough data points to cluster anything") throw();
       virtual ~InsufficientInput() throw();
   };

	/**
	 * @brief Detailed constructor
	 * @param data this data points will be clustered
	 * @param rt_threshold height of the grid cells
	 * @param mz_threshold width of the grid cells
	 * @param method_ method to use for calculating distances
	 */
   HashClustering(std::vector<DataPoint>& data, DoubleReal rt_threshold, DoubleReal mz_threshold, ClusteringMethod& method_);

	/**
	 * @brief Starts the clustering and returns a vector of subtrees when finished
	 */
   void performClustering();

	/**
	 * @brief Gets the hierarchical clustering subtrees after clustering has been performed. If the data has not been clustered yet, the method returns an empty vector
	 * @param subtrees vector of subtrees, which will be filled after the clustering process
	 */
   void getSubtrees(std::vector<std::vector<SILACTreeNode> >& subtrees);

	/**
   * @brief Extracts the clusters out of the subtrees using average silhoutte widths. If the data has not been clustered yet, the method returns an empty vector
   * @param clusters vector of clusters, which will be filled after the extraction of the clusters
   */
   void createClusters(std::vector<std::vector<DataPoint*> >& clusters);

  /**
   * @brief gets the average silhoutte widths for each subtree
   */
   std::vector<std::vector<Real> > getSilhouetteValues();

};

}

#endif /* HASHCLUSTERING_H_ */
