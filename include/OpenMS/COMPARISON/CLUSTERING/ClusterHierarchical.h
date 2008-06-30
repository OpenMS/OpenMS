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
// $Maintainer: Mathias Walzer $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_COMPARISON_CLUSTERING_CLUSTERHIERARCHICAL_H
#define OPENMS_COMPARISON_CLUSTERING_CLUSTERHIERARCHICAL_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/DATASTRUCTURES/DistanceMatrix.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterFunctor.h>
#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>
//#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>
//#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrumCompareFunctor.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <vector>


namespace OpenMS
{

	/**
		@brief Hierarchical clustering with generic clustering functions

		ClusterHierarchical clusters objects with the corresponding distancemethod and clusteringmethod.
		
		@param 	threshold_ distance value from which on merging clusters is considered to 'overcluster'
		@ingroup SpectraClustering
	*/
	class ClusterHierarchical 
	{
	private:
	
	  /// the threshold given to the clusterer
	  double threshold_; 
	
	public: 
  
	/// constructor
	ClusterHierarchical(double x=1.0):threshold_(x)
	{
	}

    	/// copy constructor
    	ClusterHierarchical(const ClusterHierarchical& source):threshold_(source.threshold_)
	{
	}
		
    	/// destructor
    	virtual ~ClusterHierarchical()
	{
	}
		
		/**
	   	    @brief clustering function
	    
	   		The generic clustering function up to the given ClusterHierarchical::threshold_. 
	   		The type of the objects to be clustered has to be the first template argument, the 
			similarity functor applicable to this type must be the second template argument, e.g. 
			with a @see PeakSpectrumCompareFunctor for @see PeakSpectrum to be clustered.
	   		The similarity functor must provide the similarity calculation with the ()-operator and 
			yield normalized values in range of [0,1].
	   		
			@param data vector of objects to be clustered	                
			@param comparator similarity functor fitting for types in data	   		
			@param clusterer a clustermethod 
			@param clusters the vector that will hold the index represented clusters @see ClusterFunctor
			

			@ingroup SpectraClustering
		*/
		template <typename Data, typename SimilarityComparator> 
		void clusterForVector(vector<Data>& data, const SimilarityComparator& comparator, const ClusterFunctor& clusterer, vector< vector < UInt > >& clusters)
		{

			//create distancematrix for data with comparator
			DistanceMatrix<double> originalDist_(data.size(),1);
			for(UInt i=0; i < data.size(); i++)
			{
				for(UInt j=i+1; j < data.size(); j++)
				{
					//distance value is 1-similarity value, since similarity is in range of [0,1]
					originalDist_.setValue(i,j,1-comparator(data[i],data[j]));
				}
			}

			//prepare input for clusterer
			DistanceMatrix<double> actualDist_(originalDist_);
						
			//prune clusters vector and fill atomar
			clusters.clear();
			/*
			for (UInt i = 0; i < originalDist_.dimensionsize(); ++i)
			{
				vector<UInt> tmp(1,i);
				clusters.push_back(tmp);
			}
			*/
			
			// create clustering with ClusterMethod, DistanceMatrix and Data
			clusterer.cluster(originalDist_, actualDist_, clusters,"",threshold_);
		}
		
		


		/**
	   	    @brief complete clustering function
	    
	   		The generic function for complete clustering and creating a corresponding dendrogramm. 
	   		The type of the objects to be clustered has to be the first template argument, the
			similarity functor applicable to this type must be the second template argument, e.g. 
			with a @see PeakSpectrumCompareFunctor for @see PeakSpectrum to be clustered.
	   		The similarity functor must provide the similarity calculation with the ()-operator and 
			yield normalized values in range of [0,1].
			In this case any threshold other than 1 will defy the construction of a complete dendrogram.
	   		
			@param data vector of objects to be clustered	                
			@param comparator similarity functor fitting for types in data	   		
			@param clusterer a clustermethod 
			@param clusters a vector to contain the clustering (only indices reffering to data)
			@param filepath the full path for the dendrogramfile to be stored
						
			@ingroup SpectraClustering
	    */
		template <typename Data, typename SimilarityComparator>
		void clusterForDendrogramm( const vector<Data>& data, const SimilarityComparator& comparator, const ClusterFunctor& clusterer, vector< vector < UInt > >& clusters, const String& filepath)
		{
			//create distancematrix for data with comparator
			DistanceMatrix<double> originalDist_(data.size(),1);

			for(UInt i=0; i < data.size(); i++)
			{
				for(UInt j=i+1; j < data.size(); j++)
				{
					//distance value is 1-similarity value, since similarity is in range of [0,1]
					originalDist_.setValue(i,j,1-comparator(data[i],data[j]));
				}
			}

			//prepare input for clusterer			
			DistanceMatrix<double> actualDist_(originalDist_);

			//prune clusters vector and fill atomar
			clusters.clear();
			/*
			for (UInt i = 0; i < originalDist_.dimensionsize(); ++i)
			{
				vector<UInt> tmp(1,i);
				clusters.push_back(tmp);
			}
			*/
			// create Clustering with ClusterMethod, DistanceMatrix
			clusterer.cluster(originalDist_,actualDist_,clusters,filepath,threshold_);
		}
		
		/* *
	   	    @brief clustering function for binned PeakSpectrum
	    
	   		The explicite version of the clustering function up to the given ClusterHierarchical::threshold_ for PeakSpectra 
	   		employing binned similarity methods. From the given PeakSpectrum BinnedSpectrum are generated, 
	   		so the similarity functor @see BinnedSpectrumCompareFunctor can be applied.
	   		
			@param data vector of PeakSpectrum to be clustered	                
			@param comparator a BinnedSpectrumCompareFunctor	   		
			@param clusterer a clustermethod 
			@param sz the desired binsize for the @see BinnedSpectrum s
			@param sp the desired binspread for the @see BinnedSpectrum s
			@param clusters the vector that will hold the index represented clusters @see ClusterFunctor
			@param distances the matrix that will hold the actual distances at each step of clustering (and after exceeding the ClusterHierarchical::threshold_ holding the final cluster distances)
			
			@ingroup SpectraClustering
		*/
		
		/*
		void clusterForVector(const vector<PeakSpectrum>& data, const BinnedSpectrumCompareFunctor& comparator, const ClusterFunctor& clusterer, 
								double sz, UInt sp, vector< vector<UInt> >& clusters)
		{
			
			vector<BinnedSpectrum> binned_data;
			binned_data.reserve(data.size());
			
			// transform each PeakSpectrum to a corresponding BinnedSpectrum with given settings of size and spread
			for(UInt i=0; i < data.size(); i++)
			{
				//double sz(2), UInt sp(1);
				binned_data[i] = BinnedSpectrum(sz,sp,data[i]);
			}

			//create distancematrix for data with comparator
			DistanceMatrix<double> originalDist_(data.size(),1);
			
			for(UInt i=0; i < binned_data.size(); i++)
			{
				for(UInt j=i+1; j < binned_data.size(); j++)
				{
					//distance value is 1-similarity value, since similarity is in range of [0,1]
					originalDist_.setValue(i,j,1-comparator(binned_data[i],binned_data[j]));
				}
			}
			
			// create Clustering with ClusterMethod, DistanceMatrix and Data
			clusterer.cluster(originalDist_, threshold_, clusters);		
		}				
		*/
		
		/* *
			@brief clustering function for binned PeakSpectrum
	    
	   		The explicite version of the complete clustering function and creating a corresponding dendrogramm for PeakSpectra 
	   		employing binned similarity methods. From the given PeakSpectrum BinnedSpectrum are generated, 
	   		so the similarity functor @see BinnedSpectrumCompareFunctor can be applied.
	   		
			@param data vector of PeakSpectrum to be clustered	                
			@param comparator a BinnedSpectrumCompareFunctor	   		
			@param clusterer a clustermethod 
			@param sz the desired binsize for the @see BinnedSpectrum s
			@param sp the desired binspread for the @see BinnedSpectrum s

			@return name String with the tablename holding the dendrogramm and additional output @see ClusterFunctor
			
			@ingroup SpectraClustering
		*/
		
		/*
		void clusterForDendrogram(const vector<PeakSpectrum>& data, const BinnedSpectrumCompareFunctor& comparator, const ClusterFunctor& clusterer, 
									double sz, UInt sp, const String& filepath)
		{
			
			vector<BinnedSpectrum> binned_data;
			binned_data.reserve(data.size());
			
			// transform each PeakSpectrum to a corresponding BinnedSpectrum with given settings of size and spread
			for(UInt i=0; i < data.size(); i++)
			{
				//double sz(2), UInt sp(1);
				binned_data[i] = BinnedSpectrum(sz,sp,data[i]);
			}
			
			//create distancematrix for data with comparator
			DistanceMatrix<double> originalDist_(data.size(),1);
						
			for(UInt i=0; i < binned_data.size(); i++)
			{
				for(UInt j=i; j < binned_data.size(); j++)
				{
					//distance value is 1-similarity value, since similarity is in range of [0,1]
					originalDist_.setValue(i,j,1-comparator(binned_data[i],binned_data[j]));
				}
			}
			
			// create Clustering with ClusterMethod, DistanceMatrix and Data
			clusterer.dendrogramInFile(originalDist_,filepath);	
			
		}
		*/
				
		/// get the Name
		static const String getName()
		{
			return "ClusterHierarchical";
		}
		
		/// get the threshold
		double getThreshold()
		{
			return threshold_;
		}
		
		/// set the threshold
		void setThreshold(double x)
		{
			threshold_= x;
		}

  	};
  	
	/** @brief Exception thrown if clustering is attempted without a normalized compare functor
		
		    due to similarity - distance conversions that are mandatory in some context, compare functors 
		    must return values normalized in the range [0,1] to ensure a clean conversion	
	*/
	class UnnormalizedComparator : public Exception::BaseException
	{	
		public:
		UnnormalizedComparator(const char* file, int line, const char* function, const char* message
          	= "Clustering with unnormalized similarity measurement requested, normalized is mandatory") throw();
		virtual ~UnnormalizedComparator() throw();
	};

	/** @brief class for structs used for keeping the input organized and odered
		...		
	*/
	class InputPieces
	{
		public:
		///
		struct inputPiece
		{
			double precursorMass;
			PeakSpectrum spectrum;
			String file;
		};

		///
		struct comparePiece
		{
			bool operator()(const inputPiece& x, const inputPiece& y) 
			{
				return (x.precursorMass < y.precursorMass);
			}
		};
		
	};
	
}
#endif //OPENMS_COMPARISON_CLUSTERING_CLUSTERHIERARCHICAL_H
