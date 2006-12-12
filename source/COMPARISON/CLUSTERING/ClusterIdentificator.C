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

#include <OpenMS/COMPARISON/CLUSTERING/ClusterIdentificator.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterSpectrum.h>
#include <OpenMS/COMPARISON/CLUSTERING/SpectrumGenerator.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/COMPARISON/CLUSTERING/Cluster.h>

#include <algorithm>
#include <sstream>

using namespace std;

namespace OpenMS
{
  ClusterIdentificator::ClusterIdentificator(DBAdapter* adapterp)
    : adapterp_(adapterp), clusterrunp_(0), masstolerance_(2.3),clustermass2Clusterid_(),clusterid2Clusterp_(),usemedian_(0)
  {
  }

  ClusterIdentificator::ClusterIdentificator(const ClusterIdentificator& source )
    : adapterp_(source.adapterp_), clusterrunp_(source.clusterrunp_), masstolerance_(source.masstolerance_),clustermass2Clusterid_(source.clustermass2Clusterid_),usemedian_(source.usemedian_)
  {
  }

  ClusterIdentificator::~ClusterIdentificator()
  {
    for (map<int,Cluster*>::const_iterator cmit = clusterid2Clusterp_.begin();
        cmit != clusterid2Clusterp_.end(); ++cmit )
    {
      delete cmit->second;
    }
  }
  
  ClusterIdentificator& ClusterIdentificator::operator=(const ClusterIdentificator& source)
  {
    adapterp_ = source.adapterp_;
    clusterrunp_ = source.clusterrunp_;
    masstolerance_ = source.masstolerance_;
    clustermass2Clusterid_ = source.clustermass2Clusterid_;
    for (map<int,Cluster*>::const_iterator cmit = clusterid2Clusterp_.begin();
        cmit != clusterid2Clusterp_.end(); ++cmit )
    {
      delete cmit->second;
    }
    for (map<int,Cluster*>::const_iterator cmit = source.clusterid2Clusterp_.begin();
        cmit != source.clusterid2Clusterp_.end(); ++cmit )
    {
      clusterid2Clusterp_.insert(make_pair(cmit->first,new Cluster(*cmit->second)));
    }
    usemedian_ = source.usemedian_;
    return *this;
  }

  uint ClusterIdentificator::existsAsReference(String sequence,uint charge)
  {
    //SequestCompareFunctor scf;
    double mass = SpectrumGenerator::instance()->getPeptidemass(sequence)+charge;
    vector<const Cluster*> candidates = getCandidateClusters_(mass,charge);
    uint result = 0;
    for (vector<const Cluster*>::const_iterator cvit = candidates.begin(); cvit != candidates.end(); ++cvit )
    {
      //const Cluster* clusterp = *cvit;
      //if ( scf.matchIsobaric(clusterp->sequence(),sequence) > 0.5 && clusterp->median().getPrecursorPeak().getCharge() ==  (int) charge ) ++result;
    }
    return result;
  }
  
  void ClusterIdentificator::setClusterRun( const ClusterExperiment::ClusterRun* crp )
  {
    for (map<int,Cluster*>::const_iterator cmit = clusterid2Clusterp_.begin();
        cmit != clusterid2Clusterp_.end(); ++cmit )
    {
      delete cmit->second;
    }
    clusterid2Clusterp_.clear();
    clustermass2Clusterid_.clear();
    clusterrunp_ = crp;
    Cluster* clusterp = 0;
    for ( map<int,ClusterNode*>::const_iterator cmit = crp->getClustering().begin();
        cmit != crp->getClustering().end(); ++cmit )
    {
      clusterp = new Cluster(*cmit->second,clusterrunp_,adapterp_);
      uint clusterid = clusterp->getPersistenceId();
      // if cluster has no id (i.e. is not handled by a PersistenceManager) use other number
      if ( !clusterid ) clusterid = cmit->first;
      clustermass2Clusterid_.insert(make_pair(clusterp->getMinParentMass(),clusterid));
      clustermass2Clusterid_.insert(make_pair(clusterp->getMaxParentMass(), clusterid));
      clusterid2Clusterp_[clusterid] = clusterp;
    }
  }

  void ClusterIdentificator::setMassTolerance(double tol)
  {
    masstolerance_ = tol;
  }
 
  vector<const Cluster*> ClusterIdentificator::getCandidateClusters_(double mass, uint charge)
  {
    vector<const Cluster*>result;
    multimap<double,int>::const_iterator cmmit = clustermass2Clusterid_.lower_bound(mass - masstolerance_);
    for (; cmmit != clustermass2Clusterid_.end() && cmmit->first < mass + masstolerance_ ; ++cmmit)
    {
      Cluster* clusterp = clusterid2Clusterp_[cmmit->second];
      if ( clusterp )
      {
        uint clcharge = clusterp->centroid().getPrecursorPeak().getCharge();
        if ( clcharge == charge )
        {
          result.push_back(clusterp);
        }
      }
      else
      {
        // should not happen
        cerr << "cluster " << cmmit->second << " is not cached\n";
      }
    }
    sort(result.begin(),result.end());
    vector<const Cluster*>::iterator enduniquerange = unique(result.begin(),result.end());
    result.erase(enduniquerange,result.end());
    return result;
  }
  
  Identification ClusterIdentificator::id(const ClusterSpectrum& cspec)
  {
    if (!clusterrunp_)
    {
      throw Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__,"no ClusterRun in ClusterIdentificator","set with setClusterRun or use quickid");
    }
    Identification result;
    stringstream ss;
    vector<const Cluster*> candidates = getCandidateClusters_(cspec.getParentMass(),cspec.getParentionCharge());
    for ( vector<const Cluster*>::const_iterator cvit = candidates.begin(); 
        cvit != candidates.end(); ++cvit )
    {
      ClusterSpectrum cls;
      if ( usemedian_ ) cls = ClusterSpectrum((**cvit).median(),adapterp_,clusterrunp_->getBinSize(),clusterrunp_->getBinSpread());
      else cls = ClusterSpectrum((**cvit).centroid(),adapterp_,clusterrunp_->getBinSize(),clusterrunp_->getBinSpread());
      clusterrunp_->preprocess(cls.spec());
      double score = clusterrunp_->similarity(cspec,cls);
      String scoretype;
      String sequence;
      ss.str("");
      ss << (*cvit)->size() << " " << cls.id();
      String information = String(ss.str());
      sequence = (*cvit)->sequence();
      ss.str("");
      //ss << clusterrunp_->getSimFunc()->getName() << "_similarity"; 
      scoretype = ss.str();
      result.insertPeptideHit(PeptideHit(score,scoretype,0,sequence));
    }
    result.sort();
    result.assignRanks();
    return result;
  }

  Identification ClusterIdentificator::quickid(const ClusterExperiment::ClusterRun& /*cr*/ ,const ClusterSpectrum& /*cspec*/, uint /*clusterrunid*/)
  {
    Identification result;
//TODO Persistence
//    double mass = cspec.getParentMass();
//    uint charge = cspec.getParentionCharge();
//    Stringstream ss;
//    if ( usemedian_ )
//    {
//      ss << "SELECT Cluster.id, MSMSFraction.precursor_peak_charge FROM Cluster JOIN MSMSFraction on Cluster.median = MSMSFraction.peak_list WHERE ";
//    }
//    else
//    {
//      ss << "SELECT Cluster.id, MSMSFraction.precursor_peak_charge FROM Cluster JOIN MSMSFraction on Cluster.centroid = MSMSFraction.peak_list WHERE ";
//    }
//     ss << " Cluster.minmass - " << masstolerance_ << " < " << mass << " AND " 
//      << " Cluster.maxmass + " << masstolerance_ << " > " << mass << " AND "
//      << " Cluster.clusterrun_id = " << clusterrunid;
//    vector<const Cluster*> candidates;
//    adapterp_->executeQuery(ss.str(),false);
//    QSqlQuery sqlres = adapterp_->lastResult();
//    vector<int> clusters;
//    while ( sqlres.next() )
//    {
//      if ( sqlres.value(1).toInt() == (int) charge )
//      {
//        clusters.push_back(sqlres.value(0).toInt());
//      }
//    }
//    for ( vector<int>::const_iterator cvit = clusters.begin();
//        cvit != clusters.end(); ++cvit )
//    {
//      //Cluster* clp = dynamic_cast<Cluster*>(adapterp_->createObject(*cvit)); //TODO Persistence
//      Cluster* clp;
//      if (clp)
//      {
//        candidates.push_back(clp);
//        if ( usemedian_ )
//        {
//          clusterrunp_->preprocess(clp->median());
//        }
//        else
//        {
//          clusterrunp_->preprocess(clp->centroid());
//        }
//      }
//    }
//    for ( vector<const Cluster*>::const_iterator cvit = candidates.begin(); 
//        cvit != candidates.end(); ++cvit )
//    {
//      ClusterSpectrum cls;
//      if ( usemedian_ ) cls = ClusterSpectrum((**cvit).median(),adapterp_,clusterrunp_->getBinSize(),clusterrunp_->getBinSpread());
//      else cls = ClusterSpectrum((**cvit).centroid(),adapterp_,clusterrunp_->getBinSize(),clusterrunp_->getBinSpread());
//      double score = cr.similarity(cspec,cls);
//      String scoretype;
//      String sequence;
//      ss.str("");
//      ss << cls.id();
//      String information = String(ss.str());
//      sequence = (*cvit)->sequence();
//      ss.str("");
//      ss << clusterrunp_->getSimFunc()->name() << "_similarity";
//      scoretype = ss.str();
//      result.insertPeptideHit(PeptideHit(score,scoretype,0, sequence));
//    }
//    result.sort();
//    result.assignRanks();
//    for ( vector<const Cluster*>::iterator vit = candidates.begin();
//        vit != candidates.end(); ++vit )
//    {
//      delete *vit;
//    }
    return result;
  }

}
