// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Mathias Walzer $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_COMPARISON_CLUSTERING_CLUSTERFUNCTOR_H
#define OPENMS_COMPARISON_CLUSTERING_CLUSTERFUNCTOR_H


#include <OpenMS/DATASTRUCTURES/DistanceMatrix.h>
#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/CONCEPT/FactoryProduct.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <vector>

namespace OpenMS
{

	/**
		@defgroup SpectraClustering classes

		@brief This class contains SpectraClustering classes

		These classes are components for clustering all kinds of data for which a distance relation, normalizable in
		the range of [0,1], is available. Mainly this will be data for which there is a corresponding CompareFunctor
		given (e.g. PeakSpectrum) that is yielding the similarity normalized in the range of [0,1] of such two
		elements, so it can easily converted to the needed distances. @see PeakSpectrumCompareFunctor().
	*/

	/**
		@brief Base class for cluster functors

		Each cluster functor employs a different method for stepwise merging clusters up to a given threshold, starting
		from the most elementary partition of data. Elements are represented by indices of a given distance matrix
		that correspond to the respective position of a element in a vector holding the real elements.

		@ingroup SpectraClustering
	*/
	class ClusterFunctor : public FactoryProduct
	{

  		public:

		/**
		@brief Exception thrown if not enough data (<2) is used

			If the set of data to be clustered contains only one data point,
			clustering algorithms would fail for obvious reasons.
		*/
		class InsufficientInput : public Exception::BaseException
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

			@param original_distance DistanceMatrix<double> containing the distances of the elements to be clustered
			@param actual_distance DistanceMatrix<double> containing the distances of the clusters at current stage
			@param clusters vector< vector<UInt> >, each vector<UInt> represents a cluster, its elements being the indices according to original_dist.
			@param filepath String&, by default empty, when given, a dendrogram will be written in a file created in that path
			@param threshold double value, the minimal distance from which on cluster merging is considered unrealistic. By default set to 1, i.e. complete clustering until only one cluster remains

			original_dist is considered mirrored at the main diagonal, so only entrys up the main diagonal are used.
			The threshold can be taken from the maximal distance of two elements considered related and adapted
			in a way corresponding to the employed clustering method.
			The clusters is emptyed and filled again with the clusters, where a vector of UInt represents
			one cluster containing the indices to the elements.
		*/
		virtual void cluster(const DistanceMatrix<double>& original_distance, DistanceMatrix<double>& actual_distance, std::vector< std::vector<UInt> >& clusters, const String filepath = "", const double threshold =1) const= 0 ;

		/// registers all derived products
		static void registerChildren();


		/// get the identifier for this FactoryProduct
		static const String getProductName()
		{
			return "ClusterFunctor";
		}

	};

}
#endif // OPENMS_COMPARISON_CLUSTERFUNCTOR_H
