// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_COMPARISON_CLUSTERING_AVERAGELINKAGE_H
#define OPENMS_COMPARISON_CLUSTERING_AVERAGELINKAGE_H

#include <vector>
#include <set>

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DistanceMatrix.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterFunctor.h>

namespace OpenMS
{
	/**
		@brief AverageLinkage ClusterMethod

		The details of the method can be found in:
		Backhaus, Erichson, Plinke, Weiber Multivariate Analysemethoden, Springer 2000 and
		Ellen M. Voorhees: Implementing agglomerative hierarchic clustering algorithms for use in document retrieval. Inf. Process. Manage. 22(6): 465-476 (1986)
		@see ClusterFunctor

		@ingroup SpectraClustering
	*/
	class OPENMS_DLLAPI AverageLinkage : public ClusterFunctor, public ProgressLogger
	{
		public:

		/// default constructor
		AverageLinkage();

		/// copy constructor
		AverageLinkage(const AverageLinkage& source);

		/// destructor
		virtual ~AverageLinkage();

		/// assignment operator
		AverageLinkage& operator = (const AverageLinkage& source);


		/**
			@brief clusters the indices according to their respective element distances

			@param original_distance DistanceMatrix<Real> containing the distances of the elements to be clustered, will be changed during clustering process, make sure to have a copy or be able to redo
			@param cluster_tree vector< BinaryTreeNode >, represents the clustering, each node contains the next merged clusters (not element indices) and their distance, strict order is kept: left_child < right_child
			@param threshold Real value, the minimal distance from which on cluster merging is considered unrealistic. By default set to 1, i.e. complete clustering until only one cluster remains
			@throw ClusterFunctor::InsufficientInput thrown if input is <2
			The clustering method is average linkage, where the updated distances after merging two clusters are each the average distances between the elements of their clusters. After @p theshold is exceeded, @p cluster_tree is filled with dummy clusteringsteps (children: (0,1), distance: -1) to the root.
			@see ClusterFunctor , BinaryTreeNode
		*/
		void operator () (DistanceMatrix<Real>& original_distance, std::vector<BinaryTreeNode>& cluster_tree, const Real threshold=1) const;

		/// creates a new instance of a AverageLinkage object
		static ClusterFunctor* create()
		{
			return new AverageLinkage();
		}

		/// get the identifier for this object
		static const String getProductName()
		{
			return "AverageLinkage";
		}

	};

}
#endif //OPENMS_COMPARISON_CLUSTERING_AVERAGELINKAGE_H
