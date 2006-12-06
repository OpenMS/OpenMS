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
#include <OpenMS/COMPARISON/CLUSTERING/Cluster.h>
#include <OpenMS/FORMAT/DataSetInfo.h>
#include <OpenMS/FORMAT/DBAdapter.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumCheapDPCorr.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterSpectrum.h>

#include <iostream>
#include <cmath>
#include <sstream>

using namespace std;

namespace OpenMS
{

  Cluster::Cluster()
    : PersistentObject(),
    datasetp_(0),medianp_(0),centroidp_(0),minmass_(0),maxmass_(0),size_(0),sequencecount_(0),sequence_("")
  {
  }
  
  Cluster::Cluster( ClusterNode* clusterp, const ClusterExperiment::ClusterRun* crp,DBAdapter* adapterp)
    :PersistentObject(),
    minmass_(-1),maxmass_(-1),sequencecount_(0),sequence_("")
  {
    datasetp_ = clusterp;
    centroidp_ = findcentroid_( crp,adapterp );
    medianp_ = findmedian_( crp,adapterp );
    updatesizes_(adapterp);
  }

  Cluster::Cluster(const ClusterNode& clusterr, const ClusterExperiment::ClusterRun* crp, DBAdapter* adapterp)
    :PersistentObject(),
    minmass_(-1),maxmass_(-1),sequencecount_(0),sequence_("")
  {
    datasetp_ = new ClusterNode(clusterr);
    centroidp_ = findcentroid_ (crp,adapterp);
    medianp_ = findmedian_( crp,adapterp );
    updatesizes_(adapterp);
  }

  Cluster::Cluster( const Cluster& source )
    :PersistentObject(source),
    	datasetp_(0),medianp_(0),centroidp_(0),sequencecount_(source.sequencecount_),sequence_(source.sequence_)
  {
    if ( source.datasetp_) datasetp_ = new ClusterNode( *source.datasetp_);
    if ( source.medianp_ ) 
    {
      medianp_ = new MSSpectrum< DPeak<1> >(*source.medianp_);
    }
    if ( source.centroidp_ ) 
    {
      centroidp_ = new MSSpectrum< DPeak<1> >(*source.centroidp_);
    }
  }

  Cluster::~Cluster()
  {
    delete datasetp_;
    delete centroidp_;
    delete medianp_;
  }

  MSSpectrum< DPeak<1> >& Cluster::median() 
  {
    return *medianp_;
  }
  
  MSSpectrum< DPeak<1> >& Cluster::centroid()
  {
    return *centroidp_;
  }
 
  const MSSpectrum< DPeak<1> >& Cluster::median() const
  {
    return *medianp_;
  }
  
  const MSSpectrum< DPeak<1> >& Cluster::centroid() const
  {
    return *centroidp_;
  }
 
  void Cluster::loadClusterNode(DBAdapter* /*adapterp*/)
  {
    if (datasetp_) return;
    else
    {
//TODO Persistence
//      stringstream ss;
//      ss << "SELECT DataSetInfo.id FROM DataSetInfo JOIN Cluster ON DataSetInfo.dataset_id = Cluster.dataset WHERE  Cluster.id = " << persistence_id_ ;
//      adapterp->executeQuery(ss.str());
//      QSqlQuery sqlres = adapterp->lastResult();
//      //int dsid = sqlres.value(0).toInt(); //TODO Persistence
//      //DataSetInfo* dsi = dynamic_cast<DataSetInfo*>(adapterp->createObject(dsid)); //TODO Persistence
//      DataSetInfo* dsi ;
//      datasetp_ = new ClusterNode();
//      for ( vector<int>::const_iterator cvit = dsi->contents("PeakList").begin(); cvit != dsi->contents("PeakList").end(); ++cvit )
//      {
//        datasetp_->insert(new ClusterNode(*cvit));
//      }
    }
  }
  
//  void Cluster::serialize(PersistenceManager& f)
//  {
//    // warning only works with DBAdapter!
//    DBAdapter& adapter = dynamic_cast<DBAdapter&>(f);
//    //updatesizes_(&adapter);
//    f.writeAttributeString("centroid",(int)centroidp_->id());
//
//    // begin ugly hack
//    if ( size() > 1 )
//    {
//      adapter.preWritePointer("PEAKLIST1D",medianp_);
//      adapter.writeConnectedPointer("info",&medianp_->,(int)medianp_);
//      for ( MSSpectrum< DPeak<1> >::iterator dit = medianp_->getContainer().begin();
//          dit != medianp_->getContainer().end(); ++dit )
//      {
//        adapter.writeConnectedPointer("peak",&*dit,(int)medianp_);
//      }
//    }
//    // end
//    f.writeAttributeString("median",(int)medianp_->id());
//    f.writeAttributeString("minmass",getMinParentMass());
//    f.writeAttributeString("maxmass",getMaxParentMass());
//    f.writeAttributeString("sequence",sequences(&adapter).rbegin()->second);
//    f.writeAttributeString("size",(int)size());
//    f.writeAttributeString("sequencecount",(int)sequencecount());
//  }
  
  Cluster& Cluster::operator=( const Cluster& source)
  {
    PersistentObject::operator=(source);
    medianp_ = source.medianp_;
    centroidp_ = source.centroidp_;
    datasetp_ = new ClusterNode(*source.datasetp_);
    minmass_ = source.minmass_;
    maxmass_ = source.maxmass_;
    return *this;
  }
  
  void Cluster::updatesizes_( DBAdapter* adapterp )
  {
    if ( !datasetp_->size() ) return;
    size_ = 0;
    for(list<int>::const_iterator cit = datasetp_->children().begin(); cit != datasetp_->children().end(); ++cit )
    {
      ++size_;
      //Spectrum< DPeak<1> >* spec = dynamic_cast<MSSpectrum< DPeak<1> >*>(adapterp->createObject(*cit));//TODO Persistence
      MSSpectrum< DPeak<1> >* spec = 0;
      
      if ( !spec ) 
      {
        stringstream ss;
   //TODO Persistence
   //     ss << "DataSet " << datasetp_->id() << " contains id " << *cit << " which is no PeakList in DB " << adapterp->dbname();
        throw Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__,"DataSet Content Problems",ss.str().c_str() );
      }
      
      if ( minmass_ < 0 ) 
      {
        minmass_ = spec->getPrecursorPeak().getPosition()[0];
      }
      
      else 
      {
        minmass_ = min(minmass_,spec->getPrecursorPeak().getPosition()[0]);
      }
      if ( maxmass_ < 0 ) 
      {
        maxmass_ = spec->getPrecursorPeak().getPosition()[0];
      }
      else
      {
        maxmass_ = max(maxmass_,spec->getPrecursorPeak().getPosition()[0]);
      }
      delete spec;
    }
    map<int, String> seqs = sequences(adapterp);
    if ( seqs.size() ) 
    {
      sequence_ = seqs.rbegin()->second;
      sequencecount_ = seqs.rbegin()->first;
    }
  }
  
  MSSpectrum< DPeak<1> >* Cluster::findcentroid_( const ClusterExperiment::ClusterRun* crp, DBAdapter* adapterp )
  {
    vector<ClusterSpectrum> spectra;
    loadClusterNode(adapterp);
    for(list<int>::const_iterator cit = datasetp_->children().begin(); cit != datasetp_->children().end(); ++cit )
    {
      ClusterSpectrum spectrum(*cit,adapterp,crp->getBinSize(),crp->getBinSpread());
      crp->preprocess(spectrum.spec());
      spectra.push_back(spectrum);
    }
    vector<vector<double> > simmatrix = vector<vector<double> >(spectra.size(),vector<double>(spectra.size()));
    vector<double> autocorr;
    for ( uint i = 0; i < spectra.size(); ++i )
    {
      autocorr.push_back((*crp->getSimFunc())(spectra[i]));
    }
    
    for ( uint i = 0; i < spectra.size(); ++i )
    {
      for ( uint j = i ; j < spectra.size() ; ++j )
      {
        simmatrix[i][j] = crp->similarity(spectra[i],spectra[j],autocorr[i],autocorr[j]);
        simmatrix[j][i] = simmatrix[i][j];
      }
    }

    vector<double> sumsim;
    
    double max = 0;
    uint maxidx = 0;
    for ( uint i = 0; i < spectra.size(); ++i )
    {
      double sum = 0;
      for ( uint j = 0; j < spectra.size(); ++j )
      {
        sum += simmatrix[i][j];
      }
      if ( sum > max ) 
      {
        max = sum;
        maxidx = i;
      }
    }
    
    return new MSSpectrum< DPeak<1> >(spectra[maxidx].getSpec());

  }
  
  MSSpectrum< DPeak<1> >* Cluster::findmedian_( const ClusterExperiment::ClusterRun* crp, DBAdapter* adapterp )
  {
    vector<ClusterSpectrum> spectra;
    loadClusterNode(adapterp);
    for(list<int>::const_iterator cit = datasetp_->children().begin(); cit != datasetp_->children().end(); ++cit )
    {
      ClusterSpectrum spectrum(*cit,adapterp,crp->getBinSize(),crp->getBinSpread());
      crp->preprocess(spectrum.spec());
      spectra.push_back(spectrum);
    }
    SpectrumCheapDPCorr merger;
    MSSpectrum< DPeak<1> >* curspecp = new MSSpectrum< DPeak<1> >(spectra[0].getSpec());
    for ( uint i = 1; i < spectra.size(); ++i )
    {
      merger.setFactor(1.0/(i+1));
      merger(*curspecp,spectra[i].getSpec());
      delete curspecp;
      curspecp = new MSSpectrum< DPeak<1> >(merger.lastconsensus());
    }
    return curspecp;
  }
    
  map<int, String> Cluster::sequences(DBAdapter* adapterp)
  {
    map<String,int> sequencecount;
    loadClusterNode(adapterp); 
    for(list<int>::const_iterator cit = datasetp_->children().begin(); cit != datasetp_->children().end(); ++cit )
    {
      ClusterSpectrum tmp(*cit,adapterp);
      sequencecount[tmp.getTophit().getSequence()]++;
    }

    map<int, String> result;
    
    for ( map<String,int>::const_iterator cmit = sequencecount.begin(); cmit != sequencecount.end(); ++cmit )
    {
      result.insert(make_pair(cmit->second,cmit->first));
    }

    return result;
  }

  uint Cluster::size() const
  {
    return size_;
  }
  
  uint Cluster::sequencecount() const
  {
    return sequencecount_;
  }
 
  const String& Cluster::sequence() const
  {
    return sequence_;
  }

//	void Cluster::persistentWrite(PersistenceManager& pm, const char* name) const throw (Exception::Base)
//	{
//		pm.writeObjectHeader(this,name);
//		//TODO Persistence
//		pm.writeObjectTrailer(name);
//	}
//	
//	void Cluster::persistentRead(PersistenceManager& pm) throw (Exception::Base)
//	{
//		//TODO Persistence
//		int dummy;
//		pm.readPrimitive(dummy,"dummy_");
//	}

}

