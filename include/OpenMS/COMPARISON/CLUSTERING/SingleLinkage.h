// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_COMPARISON_CLUSTERING_SINGLELINKAGE_H
#define OPENMS_COMPARISON_CLUSTERING_SINGLELINKAGE_H

#include <vector>
#include <set>
#include <limits>
#include <algorithm>

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DistanceMatrix.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterFunctor.h>

namespace OpenMS
{
	/**
		@brief SingleLinkage ClusterMethod

		The details of the method can be found in:
		SLINK: An optimally efficient algorithm for the single-link cluster method, The Computer Journal 1973 16(1):30-34; doi:10.1093/comjnl/16.1.30
		@see ClusterFunctor() base class.

		@ingroup SpectraClustering
	*/
	class OPENMS_DLLAPI SingleLinkage : public ClusterFunctor, public ProgressLogger
	{
		public:

			/// default constructor
			SingleLinkage();

			/// copy constructor
			SingleLinkage(const SingleLinkage& source);

			/// destructor
			virtual ~SingleLinkage();

			/// assignment operator
			SingleLinkage& operator = (const SingleLinkage& source);

			/**
				@brief clusters the indices according to their respective element distances

			@param original_distance DistanceMatrix<Real> containing the distances of the elements to be clustered
			@param cluster_tree vector< BinaryTreeNode >, represents the clustering, each node contains the next two clusters merged and their distance, strict order is kept: left_child < right_child
			@param threshold Real value to meet Baseclass interface, will not be used because algorithm used is considerably fast and does not work by growing distances
			@throw ClusterFunctor::InsufficientInput thrown if input is <2
				The clustering method is single linkage, where the updated distances after merging two clusters are each the minimal distance between the elements of their clusters.
			@see ClusterFunctor , BinaryTreeNode
			*/
			void operator () (DistanceMatrix<Real>& original_distance, std::vector<BinaryTreeNode>& cluster_tree, const Real threshold=1) const;

			/// creates a new instance of a SingleLinkage object
			static ClusterFunctor* create()
			{
				return new SingleLinkage();
			}

			/// get the identifier for this object
			static const String getProductName()
			{
				return "SingleLinkage";
			}

	};



}
#endif //OPENMS_COMPARISON_CLUSTERING_SINGLELINKAGE_H
