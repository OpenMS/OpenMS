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
// $Id: DistanceAnalyzer.C,v 1.6 2006/03/29 12:30:29 andreas_bertsch Exp $
// $Author: andreas_bertsch $
// $Maintainer:  $
// --------------------------------------------------------------------------
//
#include <OpenMS/COMPARISON/CLUSTERING/DistanceAnalyzer.h>
#include <OpenMS/MATH/STATISTICS/ROCCurve.h>
#include <OpenMS/COMPARISON/SPECTRA/CompareFunctor.h>
#include <OpenMS/FORMAT/DBAdapter.h>

#include <cmath>
#include <cassert>
#include <sstream>

using namespace std;

namespace OpenMS
{
  const String DistanceAnalyzer::info_ = "returns the distances in and between clusters";

  DistanceAnalyzer::DistanceAnalyzer()
    : AnalysisFunctor()
  {
		name_ = DistanceAnalyzer::getName();
    defaults_.setValue("clustersize_cutoff", 2);
    defaults_.setValue("mass_cutoff", 2.2);
    defaults_.setValue("resolution", 10);
    // for comparefunctor
    needsClusterrun_ = true;
    needsAdapter_ = true;
		param_ = defaults_;
  }

  DistanceAnalyzer::DistanceAnalyzer(const DistanceAnalyzer& source)
    :AnalysisFunctor(source)
  {
  }

  DistanceAnalyzer& DistanceAnalyzer::operator=(const DistanceAnalyzer& source)
  {
    AnalysisFunctor::operator=(source);
		needsClusterrun_ = source.needsClusterrun_;
		needsAdapter_ = source.needsAdapter_;
    return *this;
  }

  DistanceAnalyzer::~DistanceAnalyzer()
  {
  }
    
  map<String,double> DistanceAnalyzer::operator()(const map<int,ClusterNode*>& clustering)
  {
    canRun();
    map<String,double> result;
    ostringstream ss;
    double mass_cutoff = (double)param_.getValue("mass_cutoff");
    uint resolution = (unsigned int)param_.getValue("resolution");
    ROCCurve roc; 
    map<int,vector<double> > autocorr;
    // todo tmp
    map<int,vector<ClusterSpectrum>* > clusters;
    vector<uint> intracluster(resolution);
    vector<uint> intercluster(resolution);
    for(map<int,ClusterNode*>::const_iterator cmit = clustering.begin(); cmit != clustering.end(); ++cmit)
    {
      if ( cmit->second->size() < (int)param_.getValue("clustersize_cutoff")) continue; 
      clusters[cmit->first] = new vector<ClusterSpectrum>();
      clusters[cmit->first]->reserve(cmit->second->size());
      for(list<int>::const_iterator cit = cmit->second->children().begin(); cit != cmit->second->children().end(); ++cit )
      {
        clusters[cmit->first]->push_back(ClusterSpectrum(*cit,adapterp_,clusterrunp_->getBinSize(),clusterrunp_->getBinSpread()));
        if ( clusterrunp_ ) clusterrunp_->preprocess(clusters[cmit->first]->rbegin()->spec());
        autocorr[cmit->first].push_back((*clusterrunp_->getSimFunc())(*clusters[cmit->first]->rbegin()));
      }
    }
    for(map<int,ClusterNode*>::const_iterator cmit = clustering.begin(); cmit != clustering.end(); ++cmit)
    {
      if ( cmit->second->size() < (int)param_.getValue("clustersize_cutoff") ) continue; 
      //calculate intracluster similarities
      vector<ClusterSpectrum>& cluster1 = *clusters[cmit->first];
      for(int i = 0; i < cmit->second->size(); ++i)
      {
        for(int j = 0; j < cmit->second->size(); ++j)
        {
          if ( i < j ) 
          { 
            double score = clusterrunp_->similarity(cluster1[i],cluster1[j],autocorr[cmit->first][i],autocorr[cmit->first][j]);
            roc.insertPair(score,1);
            uint pos = (int)(score*resolution);
            if ( pos == resolution ) --pos;
            intracluster[pos]++;
          }
        }
      }
      for(map<int,ClusterNode*>::const_iterator cmit2 = cmit; cmit2 != clustering.end(); ++cmit2)
      {
        if ( cmit2->second->size() < (int)param_.getValue("clustersize_cutoff")) continue; 
        if (//cl2 is included in cl 
            cmit->second->getMaxParentMass() > cmit2->second->getMaxParentMass() &&
              cmit->second->getMinParentMass() < cmit2->second->getMinParentMass() ||
            //cl in included in lc2
            cmit2->second->getMaxParentMass() > cmit->second->getMaxParentMass() &&
              cmit2->second->getMinParentMass() < cmit->second->getMinParentMass() ||
            //cl bigger, but overlaps
            cmit->second->getMinParentMass() < cmit2->second->getMaxParentMass() + mass_cutoff &&
              cmit->second->getMaxParentMass() > cmit2->second->getMaxParentMass() ||
            //cl2 bigger, but overlaps
            cmit2->second->getMinParentMass() < cmit->second->getMaxParentMass() + mass_cutoff &&
              cmit2->second->getMaxParentMass() > cmit->second->getMaxParentMass() 
              )
        {
          vector<ClusterSpectrum>& cluster2 = *clusters[cmit2->first];
          //calculate intercluster similarities  
          for(int i = 0; i < cmit->second->size(); ++i)
          {
            //todo use symetry
            for(int j = 0; j < cmit2->second->size(); ++j)
            {
              if ( i < j ) 
              { 
                // ignore scores from peptides that have mass differences
                // above mass_cutoff, because scores arising from the mass 
                // filter in the similarity functions dont interest us here
                if ( fabs( cluster1[i].getParentMass() - cluster2[j].getParentMass() ) > mass_cutoff ) continue;
                //todo use norm
                double score = clusterrunp_->similarity(cluster1[i],cluster2[j],autocorr[cmit->first][i],autocorr[cmit2->first][j]);
                roc.insertPair(score,0);
                uint pos = (int)(score*resolution);
                if ( pos == resolution ) --pos;
                intercluster[pos]++;
              }
            }
          }
        }
      }
    }
    for ( map<int, vector<ClusterSpectrum>*>::iterator mit = clusters.begin(); mit != clusters.end(); ++mit )
    {
      delete mit->second;
    }
    result.insert(make_pair("ROC_AUC",roc.AUC()));
    vector<pair<double,double> > curvepoints = roc.curve(resolution);
    result.insert(make_pair("cutoff for 0.95 pos",roc.cutoffPos(0.95)));
    result.insert(make_pair("cutoff for 0.95 neg",roc.cutoffNeg(0.95)));
    ss.str("");
    ss << "ROCPoints:\n";
    for ( vector<pair<double,double> >::const_iterator cvit = curvepoints.begin(); cvit != curvepoints.end(); ++cvit )
    {
      ss << cvit->first << "\t" << cvit->second << "\n";
    }
    result.insert(make_pair(ss.str(),0));
    
    // intracluster and intercluster contain all similarities that are < 1/resolution*pos+1 and > 1/resolution
    for ( uint i = 0; i < resolution ; ++i )
    {
      ss.str("");
      ss << "intracluster " << (double)(i+1)/resolution;
      result.insert(make_pair(ss.str(),intracluster[i]));
      ss.str("");
      ss << "intercluster " << (double)(i+1)/resolution ;
      result.insert(make_pair(ss.str(),intercluster[i]));
    }
    return result;
  }

  String DistanceAnalyzer::info() const
  {
    return info_;
  }
}
