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

#include <OpenMS/COMPARISON/CLUSTERING/SingleLinkage.h>

using namespace std;

namespace OpenMS
{
	SingleLinkage::SingleLinkage()
	  : ClusterFunctor()
	{
			setName(SingleLinkage::getProductName());
	}
	
	SingleLinkage::SingleLinkage(const SingleLinkage& source)
	  : ClusterFunctor(source)
	{
	}
	
	SingleLinkage::~SingleLinkage()
	{
	}
	
	SingleLinkage& SingleLinkage::operator = (const SingleLinkage& source)
	{
		if (this != &source)
		{
			ClusterFunctor::operator = (source);
		}
		return *this;
	}
		
	void SingleLinkage::cluster(const DistanceMatrix<double>& originalDist, DistanceMatrix<double>& actualDist, vector< vector<UInt> >& clusters, const String filepath /*= ""*/, const double threshold /*=1*/) const throw (Exception::UnableToCreateFile,ClusterFunctor::InsufficientInput)
	{
		// attention: clustering process is done by clustering the indices 
		// pointing to elements in inputvector and distances in inputmatrix

		// input MUST have >= 2 elements!
		if(actualDist.dimensionsize()<2) 
		{
			throw ClusterFunctor::InsufficientInput(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Distance matrix to start from only contains one element");
		}
		
		// in case of preclustering actualDist and clusters have to match
		if(clusters.size()>0 && actualDist.dimensionsize()!=clusters.size()) 
		{
			throw ClusterFunctor::InsufficientInput(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Distance matrix to start from does not match the clusters to start from");
		}

		std::ofstream os;
		//make sure writing is successful
		if(filepath.size()!=0)
		{
			os.open(filepath.c_str());
			if (!os)
			{
				throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filepath);
			}
			
			os << "SINGLE LINKAGE" << endl << "dendrogram-file: " << filepath << endl << "timestamp: " << PreciseTime::now() << endl;
			
			//clusters inital state
			if(clusters.size()!=originalDist.dimensionsize())
			{
				os << "preclustering--snip" << endl;
				UInt mergedCounter = 0;
				for(UInt i = 0; i < clusters.size(); ++i)
				{
					if(clusters[i].size() == 2)
					{
						os << "*" << " | " << clusters[i][0]-mergedCounter << " | " << clusters[i][1]-mergedCounter << " | " << "#" << endl;
						++mergedCounter;
					}
				}
				os << "preclustering--snap" << endl;
			}
		}
		
		// when no clustering has been done start atomar
		if(clusters.size()==0) 
		{
			if(actualDist.dimensionsize()!=originalDist.dimensionsize()) 
			{
				throw ClusterFunctor::InsufficientInput(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Distance matrices do not match");
			}
			for (UInt i = 0; i < actualDist.dimensionsize(); ++i)
			{
				vector<UInt> tmp(1,i);
				clusters.push_back(tmp);
			}
		}

		// Initial minimum-distance pair
		actualDist.updateMinElement();
		pair<UInt,UInt> min = actualDist.getMinElementCoordinates();
	
		UInt counter1(0);
		while(actualDist(min.first,min.second) < threshold && actualDist.dimensionsize() > 2)
		{ 
			++counter1;
			
			if(filepath.size()!=0)
			{
				//write in file
				os << counter1 << " | " << min.first << " | " << min.second << " | " << actualDist(min.first,min.second) << endl;
			}

cout << counter1 << " | " << min.first << " | " << min.second << " | " << actualDist(min.first,min.second) << endl;
		//pick minimum-distance pair i,j and merge them
			//pushback elements of second to first (and then erase second)
		 	for (UInt c = 0; c < clusters[min.second].size(); ++c)
		 	{
			 	clusters[min.first].push_back(clusters[min.second][c]);
		 	}

		 	// erase second one
		 	clusters.erase(clusters.begin()+min.second,clusters.begin()+min.second+1);
	
		 	//reduce
			actualDist.reduce(min.second);
	
			//update actualDist matrix (and minimum-distance pair)
			//single linkage: new distcance between clusteres is the minimum distance between elements of each cluster
			for (UInt i = 0; i < min.first; ++i)
		 	{
			 	actualDist.setValueQuick(i,min.first,getMinDist_(i,min.first,clusters,originalDist));
			}
			for (UInt j = min.first+1; j < actualDist.dimensionsize(); ++j)
		 	{
			 	actualDist.setValueQuick(min.first,j,getMinDist_(min.first,j,clusters,originalDist));
			}
			
			//update min
			actualDist.updateMinElement();
				
			//get min-pair from triangular matrix
			min = actualDist.getMinElementCoordinates();
			
		//repeat until only two cluster remains, last step skips matrix operations
		}
		
		if(actualDist(min.first,min.second) < threshold && actualDist.dimensionsize() == 2)
		{
			if(filepath.size()!=0)
			{
				//write in file
				os << counter1 << " | " << min.first << " | " << min.second << " | " << actualDist(min.first,min.second) << endl;
			}				
			//pick minimum-distance pair i,j and merge them
			//pushback elements of second to first (and then erase second)
		 	for (UInt c = 0; c < clusters[min.second].size(); ++c)
		 	{
			 	clusters[min.first].push_back(clusters[min.second][c]);
		 	}
		 	// erase second one
		 	clusters.erase(clusters.begin()+min.second,clusters.begin()+min.second+1);
		 }
		
		if(filepath.size()!=0)
		{
			//done writing
			os.close();
		}
	}

	double SingleLinkage::getMinDist_(UInt& o, UInt x, vector< vector<UInt> >& clusters, const DistanceMatrix<double>& originalDist) const
	{
		double min(DBL_MAX);
		
		#ifdef TO_DEBUG
		cout << "getMinDist_'s clusters.size(): "<< clusters.size() << endl;
		#endif
		
		if(clusters.size() > o && clusters.size() > x )
		{
			
			#ifdef TO_DEBUG
			cout << "looking at "<< o << " and " << x << endl;
			#endif
			
			//look only in cells upper the main diagonal!
			//o, x are indices of clusters
			//wanted is mindist between elements of o to x
			for (UInt i = 0; i < clusters[x].size(); ++i)
			{
							
				for (UInt j = 0; j < clusters[o].size(); ++j)
		 		{
			 				
			 		if(clusters[x][i]<clusters[o][j])
			 		{
				 		#ifdef TO_DEBUG
			 			cout << "comparing element "<< clusters[x][i] << " with element " << clusters[o][j] << endl;
			 			cout << "distance: "<< originalDist.getValue(clusters[x][i],clusters[o][j]) << endl;
			 			#endif
			 			
				 		if(originalDist.getValue(clusters[x][i],clusters[o][j])< min)
				 		{
			 				min = originalDist.getValue(clusters[x][i],clusters[o][j]);
		 				}
				 	}
				 	else
				 	{
					 	#ifdef TO_DEBUG
			 			cout << "comparing element "<< clusters[o][j] << " with element " << clusters[x][i] << endl;
			 			cout << "distance: "<< originalDist.getValue(clusters[o][j],clusters[x][i]) << endl;
			 			#endif
			 			
					 	if(originalDist.getValue(clusters[o][j],clusters[x][i])< min)
			 			{
				 			min = originalDist.getValue(clusters[o][j],clusters[x][i]);
			 			}
				 	}
		 		}
			}
	 		
		}
		return min;
	}
/*
	void SingleLinkage::cluster(const DistanceMatrix<double,1>& originalDist, double& threshold, vector< vector<UInt> >& clusters) const throw (ClusterFunctor::InsufficientInput)
	{
	
		// attention: clustering process is done by clustering the indices 
		// pointing to elements in inputvector and distances in inputmatrix

		// input MUST have >= 2 elements!
		if(originalDist.dimensionsize()<2) 
		{
			throw ClusterFunctor::InsufficientInput(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Distance matrix only contains one element");
		}
				
		// first all indices are in separate(single) Cluster
		clusters.clear();
		for (UInt i = 0; i < originalDist.dimensionsize(); ++i)
		{
			vector<UInt> tmp(1,i);
			clusters.push_back(tmp);
		}
				
		//DistanceMatrix representing the distances between the actual clusters
		DistanceMatrix<double,1> actualDist(originalDist);
		
		// Initial minimum-distance pair
		pair<UInt,UInt> min = actualDist.getMinElementCoordinates();
			
		#ifdef TO_DEBUG 
		cout << "minpair: " << min.first << min.second << endl;
		#endif

		while(actualDist(min.first,min.second) < threshold && actualDist.dimensionsize() > 2)
		{ 		

		//pick minimum-distance pair i,j and merge them
			//pushback elements of second to first (and then erase second)
		 	for (UInt c = 0; c < clusters[min.second].size(); ++c)
		 	{
			 	clusters[min.first].push_back(clusters[min.second][c]);
		 	}
		 	// erase second one
		 	clusters.erase(clusters.begin()+min.second,clusters.begin()+min.second+1);
		 	
		 	#ifdef TO_DEBUG  
		 	for (UInt x = 0; x < clusters.size(); ++x)
  			{
	  			cout << "Cluster " << x << ":" << endl;
	  	  		for (UInt y = 0; y < clusters[x].size(); ++y)
  				{
					cout << clusters[x][y] << endl;
				}
			}
			#endif		 	
			
			//reduce
			actualDist.reduce(min.second);
			
			//update actualDist matrix (and minimum-distance pair)
			//single linkage: new distcance between clusteres is the minimum distance between elements of each cluster
			for (UInt i = 0; i < min.first; ++i)
			{
			 	actualDist.setValueQuick(i,min.first,getMinDist_(i,min.first,clusters,originalDist));						
			}
			for (UInt j = min.first+1; j < actualDist.dimensionsize(); ++j)
			{
			 	actualDist.setValueQuick(min.first,j,getMinDist_(min.first,j,clusters,originalDist));						
			}
			
			//update min
			actualDist.updateMinElement();
							
			#ifdef TO_DEBUG
			cout << actualDist << endl;
			#endif
		
			//get min-pair from actual distance matrix
			min = actualDist.getMinElementCoordinates();
			
			#ifdef TO_DEBUG 
			cout << "minpair: " << min.first << min.second << endl;
			cout << "dist.dimensionsize: " << actualDist.dimensionsize() << endl;
			#endif
			
		//repeat until threshhold distance between all clusters is exceeded
		//or actual distance matrix dimensionsize will be 2, i.e. next step there 
		//will be only one cluster, last step skips matrix operations
		}

		// if threshold still not exceeded and only two clusters left to merge
		if(actualDist(min.first,min.second) < threshold && actualDist.dimensionsize() == 2)
		{
			for (UInt c = 0; c < clusters[min.second].size(); ++c)
		 	{
			 	clusters[min.first].push_back(clusters[min.second][c]);
		 	}
		 	// erase second one
		 	clusters.erase(clusters.begin()+min.second,clusters.begin()+min.second+1);
		}		
	}
*/		

}
