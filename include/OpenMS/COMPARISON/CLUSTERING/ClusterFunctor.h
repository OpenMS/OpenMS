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
#ifndef OPENMS_COMPARISON_CLUSTERING_CLUSTERFUNCTOR_H
#define OPENMS_COMPARISON_CLUSTERING_CLUSTERFUNCTOR_H

#include <OpenMS/DATASTRUCTURES/DistanceMatrix.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterAnalyzer.h>

#include <vector>

namespace OpenMS
{

	/**
		@brief Base class for cluster functors

		Each cluster functor employs a different method for stepwise merging clusters up to a given threshold, starting from the most elementary partition of data. Elements are represented by indices of a given distance matrix, which also should represent the order of input.

		@ingroup SpectraClustering
	*/
	class OPENMS_DLLAPI ClusterFunctor
	{
		
		public:

			/**
				@brief Exception thrown if not enough data (<2) is used
		
					If the set of data to be clustered contains only one data point,
					clustering algorithms would fail for obvious reasons.
			*/
			class OPENMS_DLLAPI InsufficientInput : public Exception::BaseException
			{
				public:
					InsufficientInput(const char* file, int line, const char* function, const char* message= "not enough data points to cluster anything") throw();
					virtual ~InsufficientInput() throw();
			};
	
	
			/// default constructor
			ClusterFunctor();
	
			/// copy constructor
			ClusterFunctor(const ClusterFunctor& source);
	
			/// destructor
			virtual ~ClusterFunctor();
	
			/// assignment operator
			ClusterFunctor& operator = (const ClusterFunctor& source);
	
			/**
				@brief abstract for clustering the indices according to their respective element distances
	
				@param original_distance DistanceMatrix<Real> containing the distances of the elements to be clustered, will be changed during clustering process, make sure to have a copy or be able to redo
				@param cluster_tree vector< BinaryTreeNode >, represents the clustering, each node contains the next merged clusters (not element indices) and their distance, strict order is kept: left_child < right_child,
				@param threshold Real value, the minimal distance from which on cluster merging is considered unrealistic. By default set to 1, i.e. complete clustering until only one cluster remains
	
				@p original_distance is considered mirrored at the main diagonal, so only entrys up the main diagonal are used.
				The @p threshold can be taken from the maximal distance of two elements considered related and adapted in a way corresponding to the employed clustering method.
				The results are represented by @p cluster_tree, to get the actual clustering (with element indices) from a certain step of the clustering
				@see BinaryTreeNode , ClusterAnalyzer::cut
			*/
			virtual void operator () (DistanceMatrix<Real>& original_distance, std::vector<BinaryTreeNode>& cluster_tree, const Real threshold =1) const= 0 ;
	
			/// registers all derived products
			static void registerChildren();
	
	};

}
#endif // OPENMS_COMPARISON_CLUSTERFUNCTOR_H
