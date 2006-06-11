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
// $Id: LinkageCluster.C,v 1.5 2006/03/29 12:30:29 andreas_bertsch Exp $
// $Author: andreas_bertsch $
// $Maintainer:  $
// --------------------------------------------------------------------------
//
#include <OpenMS/COMPARISON/CLUSTERING/LinkageCluster.h>
#include <OpenMS/COMPARISON/SPECTRA/CompareFunctor.h>

#include <OpenMS/CONCEPT/Exception.h>

#include <cmath>

using namespace std;

namespace OpenMS
{
  
  const String LinkageCluster::info_ = "complete,average and single Linkage";
  
  LinkageCluster::LinkageCluster()
    :ClusterFunctor()
  {
		name_ = LinkageCluster::getName();
    defaults_.setValue("cutoff", 0.5);
    defaults_.setValue("linkage", 1); // 0 = complete, 1 = average, 2 = single
    defaults_.setValue("massdiff", 2.3);
    defaults_.setValue("split", 0);
    defaults_.setValue("stepsize", 0);
		param_ = defaults_;
  }

  LinkageCluster::LinkageCluster(const LinkageCluster& source)
    : ClusterFunctor(source)
  {
  }
  
  LinkageCluster& LinkageCluster::operator=(const LinkageCluster& source)
  {
    ClusterFunctor::operator=(source);
    return *this;
  }
 
  map<int,OpenMS::ClusterNode*> LinkageCluster::operator()(const ClusterExperiment::ClusterRun* clusterrunp)
  {
    double massdiff = (double)param_.getValue("massdiff");
    double cutoff = (double)param_.getValue("cutoff");
    int linkage = (int)param_.getValue("linkage");
    double low_cutoff  = cutoff/4;  
    uint stepsize = (unsigned int)param_.getValue("stepsize");
    
    uint splitstart = 0;
    uint splitend = 0;
    if ( fabs((double)param_.getValue("split")) > 1e-8 )
    {
      splitstart = 1;
      splitend = 3;
    }
    
    if ( linkage == 2 )
    {
      low_cutoff = cutoff;
    }
    map<int,ClusterNode*> result;
   
    bool overlap = 0;
   
    for ( uint charge = splitstart; charge <= splitend ; ++charge )
    {

      bool finished = 0;
      double startmz = 0;
      while ( !finished || overlap )
      {

        simmatrix* matrixp = 0;
        vector<ClusterSpectrum*>* cspectrap = 0;
        if ( overlap )
        {
          cspectrap = clusterrunp->getOverlap(result);
        }
        else
        {
          cspectrap = clusterrunp->getSpectra(charge,&finished,stepsize,&startmz);
        }
        vector<pair<int,int> > comparerow(cspectrap->size());
        matrixp = new simmatrix(cspectrap->size());
        for ( uint i = 0; i < matrixp->size(); ++i )
        {
          (*matrixp)[i] = new SparseVector();
        }
        vector<double> autocorr(cspectrap->size());
        uint startcompare = 0;
        for (uint i = 0; i < cspectrap->size(); ++i)
        {
          for (uint j = 0; j < cspectrap->size(); ++j)
          {
            if ( fabs(autocorr[j]) < 1e-8 )
            {
              clusterrunp->preprocess((*cspectrap)[j]->spec());
              autocorr[j] = (*clusterrunp->getSimFunc())(*(*cspectrap)[j],*(*cspectrap)[j]);
              // by convention all autocorrs of useful spectra are > 0 
              if ( autocorr[j] < 1e-8 ) autocorr[j] = -1;
            }
            if ( (*cspectrap)[j]->getParentMass() < (*cspectrap)[i]->getParentMass() - massdiff ) 
            {
              //startcompare is the first spectrum which has a compatible size
              //the spectra before that surely have similarity 0
              if ( j < startcompare ) 
              {
                j = startcompare;
              }
              else 
              {
                (*cspectrap)[startcompare++]->strip();
              }
              continue;
            }
            else if ( (*cspectrap)[j]->getParentMass() > (*cspectrap)[i]->getParentMass() + massdiff )
            {
              //no need to compare with the rest, they are all too big
              comparerow[i] = make_pair(startcompare,j);
              j = cspectrap->size() -1;
              continue;
            }
            else if ( j == cspectrap->size() -1)
            {
              comparerow[i] = make_pair(startcompare,j);
            }
            
            {
              // has already been calculated, but clustering depends on 
              // full matrix
              if ( j < i ) 
              {
                (*(*matrixp)[i])[j] = (*(*matrixp)[j])[i];
              }
              else 
              {
                (*(*matrixp)[i])[j] = clusterrunp->similarity(*(*cspectrap)[i],*(*cspectrap)[j],autocorr[i],autocorr[j]);
              }
            }
          }
        }


        multimap<double,pair<int,int> > correlation;
    
        for ( uint i = 0; i < matrixp->size(); ++i )
        {
          for ( uint j = 0; j < matrixp->size(); ++j )
          {
            if ( i != j && (*(*matrixp)[i])[j] > cutoff )
            {
              correlation.insert(std::make_pair((*(*matrixp)[i])[j] ,std::make_pair(i,j))); 
            }
          }
        }
       
        //the top level clusters, i.e. the roots
        map<int,ClusterNode*> clusterfront;
        
        int countSpectra = matrixp->size();
        
        //these provide a mapping from clusternr to rownumber
        vector<int> clnr2row(countSpectra);
        vector<int> row2clnr(countSpectra);
        
        //create a new cluster for every ClusterSpectrum
        for ( int i = 0; i < countSpectra; ++i)
        {
          clusterfront[i] = new ClusterNode((*cspectrap)[i]->id(),(*cspectrap)[i]->getParentMass());
          clnr2row[i] = i;
          row2clnr[i] = i;
        }
        
        
        //saves which clusters are inside other clusters
        map<int,int> incluster;
        
        int k = cspectrap->size();
        for ( uint di = 0; di < cspectrap->size(); ++di )
        {
          delete (*cspectrap)[di];
        }
        delete cspectrap;
       
        //while we have clusterpairs that are more similar than cutoff
        while ( correlation.size() > 0  && correlation.rbegin()->first > cutoff ) 
        //while ( correlation_list.size() > 0  && correlation_list.rbegin()->first > cutoff ) 
        {
         
          //i,j = clusters with similarity = correlation.rbegin()->first
          int i = correlation.rbegin()->second.second;
          int j = correlation.rbegin()->second.first;

          //rownrs of i and j
          int ri = clnr2row.at(correlation.rbegin()->second.second); 
          int rj = clnr2row.at(correlation.rbegin()->second.first);
          
          //there is no use in clustering something with itself
          if ( i != j ) 
          {
            
            //if the clusters i,j are already inside other clusters, ignore them
            if (incluster.find(i) != incluster.end()  || incluster.find(j) != incluster.end()) 
            {
            }
            else 
            {
              
              //reuse entries from cluster i for new cluster k 
              //create new cluster (= ClusterNode)
              uint iclsize = clusterfront[i]->size();
              uint jclsize = clusterfront[j]->size();
              clusterfront[k] = new ClusterNode(clusterfront[i], clusterfront[j]);
              
              uint mincompare = min(comparerow[ri].first,comparerow[rj].first);
              uint maxcompare = max(comparerow[ri].second,comparerow[rj].second);
              comparerow[ri] = make_pair(mincompare,maxcompare);
              
              //calculate distances to new Cluster
              for ( int l = comparerow[ri].first; l <= comparerow[ri].second; l++)
              {
                if ( incluster.find(row2clnr[l]) != incluster.end() )
                {
                  continue;
                }
                switch (linkage)
                {
                  case 2:
                    {
                      (*(*matrixp)[l])[ri] = min((*(*matrixp)[l])[ri],(*(*matrixp)[l])[rj]);
                    }
                    break;
                  case 0:
                    {
                      (*(*matrixp)[l])[ri] = max((*(*matrixp)[l])[ri],(*(*matrixp)[l])[rj]);
                    }
                    break;
                  case 1:
                    {
                      (*(*matrixp)[l])[ri] = ((*(*matrixp)[l])[ri]*iclsize+(*(*matrixp)[l])[rj]*jclsize)/(iclsize+jclsize);
                    }
                    break;
                 }
                (*(*matrixp)[ri])[l] = (*(*matrixp)[l])[ri];
                
                //dont save zeros
                if ( (*(*matrixp)[ri])[l] > low_cutoff ) 
                {
                  correlation.insert(std::make_pair((*(*matrixp)[ri])[l],std::make_pair(row2clnr[l],k)));
                }
              }
              delete (*matrixp)[rj];
              (*matrixp)[rj] = 0;
              clusterfront.erase(i);
              clusterfront.erase(j);
              incluster[i] = j;
              incluster[j] = i;
              clnr2row.push_back(ri);
              row2clnr.at(ri) = k;
              k++;
            }
          }
          correlation.erase(--correlation.end());
        }
        for ( uint i = 0; i < matrixp->size(); ++i )
        {
          delete (*matrixp)[i];
        }
        uint size = 0;
        if ( result.size() ) size = result.rbegin()->first;
        for ( map<int,ClusterNode*>::const_iterator cmit = clusterfront.begin(); cmit != clusterfront.end(); ++cmit )
        {
          result.insert(make_pair(++size,cmit->second));
        }
        if ( overlap )
        {
          overlap = false;
        }
        else
        {
          overlap = true;
        }
      }
    }
    return result;
  }

  String LinkageCluster::info() const
  {
    return info_;
  }

}
