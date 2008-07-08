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

#include <OpenMS/COMPARISON/CLUSTERING/CompleteLinkage.h>

using namespace std;

namespace OpenMS
{
	CompleteLinkage::CompleteLinkage()
	  : ClusterFunctor()
	{
			setName(CompleteLinkage::getProductName());
	}
	
	CompleteLinkage::CompleteLinkage(const CompleteLinkage& source)
	  : ClusterFunctor(source)
	{
	}
	
	CompleteLinkage::~CompleteLinkage()
	{
	}
	
	CompleteLinkage& CompleteLinkage::operator = (const CompleteLinkage& source)
	{
		if (this != &source)
		{
			ClusterFunctor::operator = (source);
		}
		return *this;
	}
	
	void CompleteLinkage::cluster(const DistanceMatrix<double>& original_distance, DistanceMatrix<double>& actual_distance, vector< vector<UInt> >& clusters, const String filepath /*= ""*/, const double threshold /*=1*/) const throw (Exception::UnableToCreateFile,ClusterFunctor::InsufficientInput)
	{
		// attention: clustering process is done by clustering the indices 
		// pointing to elements in inputvector and distances in inputmatrix
		
		// input MUST have >= 2 elements!
		if(actual_distance.dimensionsize()<2) 
		{
			throw ClusterFunctor::InsufficientInput(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Distance matrix to start from only contains one element");
		}
		
		// in case of preclustering actual_distance and clusters have to match
		if(clusters.size()>0 && actual_distance.dimensionsize()!=clusters.size()) 
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
			
			os << "COMPLETE LINKAGE" << endl << "dendrogram-file: " << filepath << endl << "timestamp: " << PreciseTime::now() << endl;
			
			//clusters inital state
			if(clusters.size()!=original_distance.dimensionsize())
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
			if(actual_distance.dimensionsize()!=original_distance.dimensionsize()) 
			{
				throw ClusterFunctor::InsufficientInput(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Distance matrices do not match");
			}
			for (UInt i = 0; i < actual_distance.dimensionsize(); ++i)
			{
				vector<UInt> tmp(1,i);
				clusters.push_back(tmp);
			}
		}
		
		// Initial minimum-distance pair
		actual_distance.updateMinElement();
		pair<UInt,UInt> min = actual_distance.getMinElementCoordinates();
				
		UInt clustersteps_counter(0);
		while(actual_distance(min.first,min.second) < threshold && actual_distance.dimensionsize() > 2)
		{ 
			++clustersteps_counter;
			
			if(filepath.size()!=0)
			{
				//write in file
				os << clustersteps_counter << " | " << min.first << " | " << min.second << " | " << actual_distance(min.first,min.second) << endl;
			}
				
		//pick minimum-distance pair i,j and merge them
			//pushback elements of second to first (and then erase second)
			for (UInt c = 0; c < clusters[min.second].size(); ++c)
			{
				clusters[min.first].push_back(clusters[min.second][c]);
			}
			// erase second one
			clusters.erase(clusters.begin()+min.second,clusters.begin()+min.second+1);
					
			//reduce
			actual_distance.reduce(min.second);
					 	 	
			//update actual_distance matrix (and minimum-distance pair)
			//complete linkage: new distcance between clusteres is the minimum distance between elements of each cluster
			for (UInt i = 0; i < min.first; ++i)
			{
			 	actual_distance.setValueQuick(i,min.first,getMaxDist_(i,min.first,clusters,original_distance));
			}
			for (UInt j = min.first+1; j < actual_distance.dimensionsize(); ++j)
			{
			 	actual_distance.setValueQuick(min.first,j,getMaxDist_(min.first,j,clusters,original_distance));
			}
			
			//update min
			actual_distance.updateMinElement();
				
			//get min-pair from triangular matrix
			min = actual_distance.getMinElementCoordinates();
			
		//repeat until only two cluster remains, last step skips matrix operations
		}
		
		if(actual_distance(min.first,min.second) < threshold && actual_distance.dimensionsize() == 2)
		{
			if(filepath.size()!=0)
			{
				//write in file
				os << clustersteps_counter << " | " << min.first << " | " << min.second << " | " << actual_distance(min.first,min.second) << endl;
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

	
	double CompleteLinkage::getMaxDist_(UInt& o, UInt x, vector< vector<UInt> >& clusters, const DistanceMatrix<double>& original_dist) const
	{
		double max(DBL_MIN);
		
		if(clusters.size() > o && clusters.size() > x )
		{
			
			//look only in cells upper the main diagonal!
			//o, x are indices of clusters
			//wanted is mindist between elements of o to x
			for (UInt i = 0; i < clusters[x].size(); ++i)
			{
							
				for (UInt j = 0; j < clusters[o].size(); ++j)
				{
							
					if(clusters[x][i]<clusters[o][j])
					{
						#ifdef OPENMS_DEBUG
						cout << "comparing element "<< clusters[x][i] << " with element " << clusters[o][j] << endl;
						cout << "distance: "<< original_dist.getValue(clusters[x][i],clusters[o][j]) << endl;
						#endif
						
						if(original_dist.getValue(clusters[x][i],clusters[o][j])> max)
							max = original_dist.getValue(clusters[x][i],clusters[o][j]);
					}
					else
					{
						#ifdef OPENMS_DEBUG
						cout << "comparing element "<< clusters[o][j] << " with element " << clusters[x][i] << endl;
						cout << "distance: "<< original_dist.getValue(clusters[o][j],clusters[x][i]) << endl;
						#endif
						
						if(original_dist.getValue(clusters[o][j],clusters[x][i])> max)
							max = original_dist.getValue(clusters[o][j],clusters[x][i]);
					}
				}
			}
			
		}
		return max;
	}

}
