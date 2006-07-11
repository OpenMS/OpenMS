// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer:  $
// --------------------------------------------------------------------------
//
#include <OpenMS/COMPARISON/CLUSTERING/ClusterCompareFunctor.h>

#include <OpenMS/COMPARISON/CLUSTERING/ClusterExperiment.h>

#include <cmath>

using namespace std;

namespace OpenMS
{
  ClusterCompareFunctor::ClusterCompareFunctor()
    : AnalysisFunctor()
  {
		name_ = ClusterCompareFunctor::getName();
    needsAdapter_ = false;
    
    /// the reference Clustering is taken from the other ClusterRun <br>  
    needsClusterrun_ = true;
  }

  ClusterCompareFunctor::ClusterCompareFunctor(const ClusterCompareFunctor& source)
    :AnalysisFunctor(source)
  {
  }

  ClusterCompareFunctor& ClusterCompareFunctor::operator=(const ClusterCompareFunctor& source)
  {
    AnalysisFunctor::operator=(source);
    return *this;
  }

  ClusterCompareFunctor::~ClusterCompareFunctor()
  {
  }

  map<String,double> ClusterCompareFunctor::operator()(const map<int,ClusterNode*>& clustering)
  {
    
    canRun();
    
    map<String,double> result;
    
    //which id is in what cluster in the referenceclustering 
    vector<pair<int,int> > refdist;
    //same for the argument ("clustering")
    map<int,int> argdist;

    //calculate the refdist
    for (map<int,ClusterNode*>::const_iterator cmit = clusterrunp_->getClustering().begin(); cmit != clusterrunp_->getClustering().end(); ++cmit)
    {
      for (list<SignedInt>::const_iterator cit = cmit->second->children().begin(); cit != cmit->second->children().end(); ++cit )
      {
        refdist.push_back(make_pair(*cit,cmit->first));
      }
    }
    for (map<int,ClusterNode*>::const_iterator cmit = clustering.begin(); cmit != clustering.end(); ++cmit)
    {
      for (list<SignedInt>::const_iterator cit = cmit->second->children().begin(); cit != cmit->second->children().end(); ++cit )
      {
        argdist.insert(make_pair(*cit,cmit->first));
      }
    }

    //to count co-clusterdness 
    long a = 0;
    long b = 0;
    long c = 0;
    long d = 0;
    
    //compare the two distributions
    for (uint i = 0; i < refdist.size(); ++i)
    {
      for (uint j = i+1; j < refdist.size(); ++j) 
      {
        int id1 = refdist[i].first;
        int id2 = refdist[j].first;
        int rcluster1 = refdist[i].second;
        int rcluster2 = refdist[j].second;
        int acluster1 = argdist[id1];
        int acluster2 = argdist[id2];
        if ( rcluster1 == rcluster2 && acluster1 == acluster2 )
        {
          a++;
        }
        else if (rcluster1 == rcluster2 && acluster1 != acluster2)
        {
          b++;
        }
        else if (rcluster1 != rcluster2 && acluster1 == acluster2 )
        {
          c++;
        }
        else if (rcluster1 != rcluster2 && acluster1 != acluster2 )
        {
          d++;
        }
      }
    }
    result["a"] = a;
    result["b"] = b;
    result["c"] = c;
    result["d"] = d;
    result["jackard"] = (double)a/(a+b+c);
    result["folkes_mallows"] = sqrt( ((double)a/(a+b)) * ((double)a/(a+c)) );
    result["rand_statistics"] = (double)(a+d)/(a+b+c+d);
    result["hubert"] = (double)(a*d-b*c)/(sqrt((double)((a+b)*(c+d)*(a+c)*(b+d))));
    return result;
  }
}
