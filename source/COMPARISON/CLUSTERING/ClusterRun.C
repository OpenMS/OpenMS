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
#include <OpenMS/COMPARISON/CLUSTERING/ClusterExperiment.h>
#include <fstream>
#include <OpenMS/SYSTEM/StopWatch.h>

#include <OpenMS/COMPARISON/CLUSTERING/ClusterFunctor.h>
#include <OpenMS/COMPARISON/SPECTRA/CompareFunctor.h>
#include <OpenMS/FILTERING/TRANSFORMERS/PreprocessingFunctor.h>
#include <OpenMS/FORMAT/DataSetInfo.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterFactory.h>
#include <OpenMS/COMPARISON/CLUSTERING/AnalysisFunctor.h>
#include <OpenMS/COMPARISON/CLUSTERING/Cluster.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <OpenMS/COMPARISON/CLUSTERING/helper.h>

using namespace std;

namespace OpenMS
{

  ClusterExperiment::ClusterRun::ClusterRun()
    :PersistentObject(),parentp_(0),  clusters_(),binsize_(1),binspread_(1),didrun_(0),sim_funcp_(0), preprocess_queue_(), cluster_funcp_(0), analysis_queue_(),norm_(geometric)
  {
  }
  
  ClusterExperiment::ClusterRun::ClusterRun( const ClusterExperiment& parent)
    :PersistentObject(),parentp_(&parent),  clusters_(),binsize_(1),binspread_(1),didrun_(0),sim_funcp_(0), preprocess_queue_(), cluster_funcp_(0), analysis_queue_(),norm_(geometric)
  {
  }

  ClusterExperiment::ClusterRun::ClusterRun(const ClusterExperiment::ClusterRun& source)
    :PersistentObject(source),parentp_(source.parentp_),clusters_(),binsize_(source.binsize_),binspread_(source.binspread_),didrun_(source.didrun_),sim_funcp_(0),preprocess_queue_(), cluster_funcp_(0),analysis_queue_(source.analysis_queue_),norm_(source.norm_)
  {
    PersistentObject::operator=(source);
    ClusterFactory* fp = ClusterFactory::instance();
    for ( map<int,ClusterNode*>::const_iterator cmit = source.clusters_.begin();
        cmit != source.clusters_.end(); ++cmit )
    {
      clusters_.insert(make_pair(cmit->first,new ClusterNode(*cmit->second)));
    }
    if ( source.sim_funcp_ ) 
    {
      sim_funcp_ = dynamic_cast<CompareFunctor*>(fp->duplicate(source.sim_funcp_));
    }
    if ( source.cluster_funcp_ ) 
    {
      cluster_funcp_ = dynamic_cast<ClusterFunctor*>(fp->duplicate(source.cluster_funcp_));
    }
    for ( vector<PreprocessingFunctor*>::const_iterator cvit = source.preprocess_queue_.begin();
        cvit != source.preprocess_queue_.end(); ++cvit )
    {
      PreprocessingFunctor* mfp = dynamic_cast<PreprocessingFunctor*>(fp->duplicate(*cvit));
      preprocess_queue_.push_back(mfp);
    }
  }

  ClusterExperiment::ClusterRun& ClusterExperiment::ClusterRun::operator = ( const ClusterExperiment::ClusterRun& source)
  {
    PersistentObject::operator=(source);
    parentp_=source.parentp_;
    didrun_ = source.didrun_;
    binsize_=source.binsize_; 
    binspread_=source.binspread_; 
    analysis_queue_=source.analysis_queue_;
    norm_ = source.norm_;
    ClusterFactory* fp = ClusterFactory::instance();
    for ( map<int,ClusterNode*>::const_iterator cmit = source.clusters_.begin();
        cmit != source.clusters_.end(); ++cmit )
    {
      clusters_.insert(make_pair(cmit->first,new ClusterNode(*cmit->second)));
    }
    if ( source.sim_funcp_ ) 
    {
      sim_funcp_ = dynamic_cast<CompareFunctor*>(fp->duplicate(source.sim_funcp_));
    }
    if ( source.cluster_funcp_ ) 
    {
      cluster_funcp_ = dynamic_cast<ClusterFunctor*>(fp->duplicate(source.cluster_funcp_));
    }
    for ( vector<PreprocessingFunctor*>::const_iterator cvit = source.preprocess_queue_.begin();
        cvit != source.preprocess_queue_.end(); ++cvit )
    {
      PreprocessingFunctor* mfp = dynamic_cast<PreprocessingFunctor*>(fp->duplicate(*cvit));
      preprocess_queue_.push_back(mfp);
    }
    return *this;
  }
 
  void ClusterExperiment::ClusterRun::setNorm(Norm norm)
  {
    norm_ = norm;
  }
  
  void ClusterExperiment::ClusterRun::deleteContents()
  {
    didrun_ = 0;
    for (map<int,ClusterNode*>::iterator mit = clusters_.begin(); mit != clusters_.end(); ++mit)
    {
      delete mit->second;
    }
    delete sim_funcp_;
    for (vector<PreprocessingFunctor*>::iterator lit = preprocess_queue_.begin(); lit != preprocess_queue_.end() ;++lit)
    {
      delete *lit;
    }
    delete cluster_funcp_;
    deleteCachedClusters_();
  }

  ClusterExperiment::ClusterRun::~ClusterRun()
  {
    deleteContents();
  }

  void ClusterExperiment::ClusterRun::setBinSize(double size)
  {
    didrun_ = 0;
    binsize_ = size;
  }

  void ClusterExperiment::ClusterRun::setBinSpread(uint spread)
  {
    didrun_ = 0;
    binspread_ = spread;
  }

  int ClusterExperiment::ClusterRun::addMower(PreprocessingFunctor* mower)
  {
    didrun_ = 0;
    preprocess_queue_.push_back(mower);
    return preprocess_queue_.size() -1;
  }

  void ClusterExperiment::ClusterRun::setSimFunc(CompareFunctor* sim_func)
  {
    didrun_ = 0;
    sim_funcp_ = sim_func;
  }

  void ClusterExperiment::ClusterRun::setClusterFunc(ClusterFunctor* clfunc)
  {
    didrun_ = 0;
    cluster_funcp_ = clfunc;
  }

  void ClusterExperiment::ClusterRun::setClustering(const map<int,ClusterNode*>& clusters)
  {
    clusters_ = clusters;
    didrun_ = 1;
  }

  int ClusterExperiment::ClusterRun::addAnalysisFunctor(AnalysisFunctor* anfunc)
  {
    didrun_ = 0;
    analysis_queue_.push_back(Analysis(anfunc));
    return analysis_queue_.size() -1;
  }

  void ClusterExperiment::ClusterRun::preprocess(MSSpectrum< DPeak<1> >& spec) const
  {
    for (vector<PreprocessingFunctor*>::const_iterator it = preprocess_queue_.begin(); it != preprocess_queue_.end(); ++it)
    {
      (**it).apply(spec);
    }
  }

  bool ClusterExperiment::ClusterRun::canRun() const
  {
    if (isComplete() == 1)
    {
      return 1;
    }
    else return 0;
  }

  // todo didnt complain without clusterfuncp
  //more verbose alternative to canRun
  // 1 : all is well
  //-1 : DBAdapter not working
  //-2 : DBAdapter not there
  //-3 : no DataSet
  //-4 : no Similarity Function
  //-5 : no Clustering Function
  int ClusterExperiment::ClusterRun::isComplete() const
  {
    //debug
    return 1;
//TODO Persistence
//    if (parentp_->adapterp_)
//    {
//      try 
//      {
//        parentp_->adapterp_->executeQuery("select 1",false);
//      }
//      //todo do i catch all?
//      catch (DBAdapter::NotConnected)
//      {
//        return -1;
//      }
//      catch (Exception::Base& e)
//      {
//        cerr << "unexpected exception " << e.what();
//        throw (e);
//      }
//    }
//    else return -2;
//    if ( parentp_->datasetname_.size() == 0 ) return -3;
//    if ( !sim_funcp_ ) return -4;
//    if ( !cluster_funcp_ ) return 5;
//    else return 1;
  }

  void ClusterExperiment::ClusterRun::cluster()
  {
    clusters_ = (*cluster_funcp_)(this);
  }

  void ClusterExperiment::ClusterRun::analyze()
  {
    for (vector<Analysis>::iterator it = analysis_queue_.begin(); it != analysis_queue_.end(); ++it)
    {
      it->setAdapter(parentp_->adapterp_);
      it->run(clusters_);
    }
  }
  
  void ClusterExperiment::ClusterRun::run() throw(CanNotRun)
  {
    StopWatch stopwatch;
    if (didrun_  == 1) 
    {
      return;
    }
    
    if ( isComplete() != 1 ) 
    {
      switch (isComplete())
      {
        case -1:
          throw(CanNotRun("DBAdapter not working",__FILE__, __LINE__, __PRETTY_FUNCTION__));
        case -2:
          throw(CanNotRun("DBAdapter missing",__FILE__, __LINE__, __PRETTY_FUNCTION__));
        case -3:
          throw(CanNotRun("DataSet missing",__FILE__, __LINE__, __PRETTY_FUNCTION__));
        case -4:
          throw(CanNotRun("CompareFunctor missing",__FILE__, __LINE__, __PRETTY_FUNCTION__));
        case -5:
          throw(CanNotRun("ClusterFunctor missing",__FILE__, __LINE__, __PRETTY_FUNCTION__));
      }
    }

    stopwatch.start();
    
    cluster();
    
    stopwatch.reset();
    
    //analyze
    analyze();
    stopwatch.reset();
    didrun_ = 1;
  }

//  void ClusterExperiment::ClusterRun::serialize(PersistenceManager& f)
//  {
//    DBAdapter* adapterp = dynamic_cast<DBAdapter*>(&f);
//    if ( !adapterp )
//    {
//      throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
//    }
//    f.writeAttributeString("binsize",(double)getBinSize());
//    f.writeAttributeString("binspread",(int)getBinSpread());
//    if ( getNorm() == arithmetic ) f.writeAttributeString("norm","arithmetic");
//    else if ( getNorm() == geometric ) f.writeAttributeString("norm","geometric");
//    else if ( getNorm() == none ) f.writeAttributeString("norm","none");
//    for (vector<PreprocessingFunctor*>::const_iterator lit = preprocess_queue_.begin(); lit != preprocess_queue_.end(); ++lit)
//    {
//      f.writePointer("FactoryProduct",*lit);
//    }
//    if ( sim_funcp_ != 0 )
//    {
//      f.writePointer("FactoryProduct",sim_funcp_);
//    }
//    if ( cluster_funcp_ != 0 )
//    {
//      f.writePointer("FactoryProduct",cluster_funcp_);
//    }
//    if ( clusters_.size() > 0 )
//    {
//      deleteCachedClusters_();
//      for (map<int,ClusterNode*>::const_iterator mit = clusters_.begin(); mit != clusters_.end(); ++mit)
//      {
//        // create Clusters
//        Cluster* newcluster = new Cluster(*mit->second,this,adapterp);
//        clustersToDelete_.push_back(newcluster);
//        f.writePointer("Cluster",newcluster);
//      }
//    } 
//    if ( analysis_queue_.size() > 0 ) 
//    {
//      for (vector<Analysis>::iterator vit = analysis_queue_.begin(); vit != analysis_queue_.end(); ++vit)
//      {
//        if ( parentp_ ) vit->setDataSet(parentp_->getDataSetInfo()->id());
//        f.writePointer("Analysis",&*vit);
//      }
//    }
//
//  }
  
  void ClusterExperiment::ClusterRun::deleteCachedClusters_()
  {
//TODO Persistence
//    for (vector<Cluster*>::const_iterator cvit = clustersToDelete_.begin(); 
//        cvit != clustersToDelete_.end(); ++cvit )
//    {
//      delete *cvit;
//    }
//    clustersToDelete_.clear();
	}

  void ClusterExperiment::ClusterRun::save(ostream&  document, int& ind)
  {
    document << "<ClusterRun finished=\"" << didrun_ << "\" >\n";
    document << indent(++ind);
   
    if ( preprocess_queue_.size() > 0 )
    {
      document << "<Preprocessing>\n";
      document << indent(++ind);

      for (vector<PreprocessingFunctor*>::const_iterator lit = preprocess_queue_.begin(); lit != preprocess_queue_.end(); ++lit)
      {
        document << "<FilterFunc ";
        // TODO (*lit)->save(document ,ind);
        document << indent(--ind);
        document << "</FilterFunc>\n";
      }
      document << indent(--ind);
      document << "</Preprocessing>\n";
    }
    if ( sim_funcp_ != 0 )
    {
      if ( norm_ == arithmetic )
      {
        document << "<Norm mean = \"arithmetic\"/>\n";
      }
      else if ( norm_ == geometric )
      {
        document << "<Norm mean = \"geometric\"/>\n";
      }
      else if ( norm_ == none )
      {
        document << "<Norm mean = \"none\"/>\n";
      }
      if ( sim_funcp_->usebins() )
      {
        document << "<Bins size = \"" << binsize_  << "\" spread = \"" << binspread_<< "\"/>\n"; 
      }
      document << "<SimFunc " ;
      // TODO sim_funcp_->save(document,ind);
      document << indent(--ind);
      document << "</SimFunc>\n";
    }
    if ( cluster_funcp_ != 0 )
    {
      document << "<ClustFunc ";
      // TODO cluster_funcp_->save(document,ind);
      document << indent(--ind);
      document << "</ClustFunc>\n";
    }
    if ( clusters_.size() > 0 )
    {
      document << "<Clustering>\n";
      document << indent(++ind);
      for (map<int,ClusterNode*>::const_iterator mit = clusters_.begin(); mit != clusters_.end(); ++mit)
      {
        document << "<Cluster id =\"" << mit->first << "\" size =\"" << mit->second->size() << "\" mininum_parent_mass =\"" << mit->second->getMinParentMass() << "\" maximum_parent_mass =\"" << mit->second->getMaxParentMass() << "\">\n";
        document << indent(++ind);
        for (list<SignedInt>::const_iterator cit = mit->second->children().begin(); cit != mit->second->children().end(); ++cit )
        {
          document << "<member_id>" << *cit << "</member_id>\n";
        }
        //todo? median id 
        document << indent(--ind);
        document << "</Cluster>\n";
      }
      document << indent(--ind);
      document << "</Clustering>\n";
    } 
    if ( analysis_queue_.size() > 0 ) 
    {
      document << "<Evaluation>\n";
      document << indent(++ind);
      for (vector<Analysis>::iterator vit = analysis_queue_.begin(); vit != analysis_queue_.end(); ++vit)
      {
        document << "<Analysis>\n";
        document << indent(++ind);
        vit->save(document,ind);
        document << indent(--ind);
        document << "</Analysis>\n";
      }
      document << indent(--ind);
      document << "</Evaluation>\n";
    }
    document << indent(--ind);
    document << "</ClusterRun>\n";
  }

  const ClusterExperiment::Analysis& ClusterExperiment::ClusterRun::operator[] (uint pos)const throw(Exception::IndexOverflow)
  {
    if (pos >= analysis_queue_.size())
    {
      throw(Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__,pos,analysis_queue_.size()));
    }
    else
    {
      return analysis_queue_[pos];
    }
  }
  
  /**
  divides dataset into smaller subsets for saving space during clustering<br>
  */
  vector<ClusterSpectrum*>* ClusterExperiment::ClusterRun::getSpectra(uint /*charge*/, bool* /*finished*/, uint /*stepsize*/, double* /*startmz*/ ) const
  {
    vector<ClusterSpectrum*>* cspectra = new vector<ClusterSpectrum*>();
//TODO Persistence
//    stringstream ss;
//
//    // for m/z sorting of spectra
//    const DataSetInfo* dsip = parentp_->getDataSetInfo();
//    if ( !dsip ) cerr << "dsip problems\n";
//    ss.str(""); 
//    ss << "CREATE TEMPORARY TABLE tmpDataSet ( id BIGINT, member BIGINT ) ";
//    parentp_->adapterp_->executeQuery(ss.str(),false);
//    for ( vector<int>::const_iterator cvit = dsip->contents("PeakList").begin(); cvit != dsip->contents("PeakList").end(); ++cvit )
//    {
//      ss.str("");
//      ss << "INSERT INTO tmpDataSet VALUES(" << dsip->dataset_id() << ","<< *cvit << ")";
//      parentp_->adapterp_->executeQuery(ss.str(),false);
//    }
//    
//    // get spectra
//    ss.str("");
//    ss << "SELECT tmpDataSet.member FROM tmpDataSet JOIN MSMSFraction ON tmpDataSet.member = MSMSFraction.peak_list WHERE";
//    if ( charge > 0 ) ss << " MSMSFraction.precursor_peak_charge = " << charge << " AND ";
//    if ( startmz ) ss << " MSMSFraction.target_m_to_z >= " << *startmz << " AND ";
//    ss << " tmpDataSet.id = " << dsip->dataset_id()  << " ORDER BY MSMSFraction.target_m_to_z";
//    if ( stepsize > 0 ) ss << " LIMIT " << stepsize;
//    
//    parentp_->adapterp_->executeQuery(ss.str(),false);
//    QSqlQuery sqlres = parentp_->adapterp_->lastResult();
//    cspectra->reserve(sqlres.size());
//    
//    while ( sqlres.next() )
//    {
//      cspectra->push_back(new ClusterSpectrum(sqlres.value(0).toInt(),parentp_->adapterp_,binsize_,binspread_));
//    }
//
//    if ( finished )
//    {
//      ss.str("");
//      ss << "SELECT count(tmpDataSet.member) FROM tmpDataSet JOIN MSMSFraction ON tmpDataSet.member = MSMSFraction.peak_list WHERE";
//      if ( charge > 0 ) ss << " MSMSFraction.precursor_peak_charge = " << charge << " AND ";
//      if ( startmz > 0 ) ss << " MSMSFraction.target_m_to_z >= " << *startmz << " AND ";
//      ss << " tmpDataSet.id = " << dsip->dataset_id() << " ORDER BY MSMSFraction.target_m_to_z";
//      parentp_->adapterp_->executeQuery(ss.str(),false);
//      QSqlQuery sqlres2 = parentp_->adapterp_->lastResult();
//      sqlres2.first();
//      int size = sqlres2.value(0).toInt();
//      if ( size == sqlres.size() ) 
//      {
//        *finished = 1;
//      }
//      else
//      {
//        *finished = 0;
//      }
//    }
//    
//    double overlapsize = 10;
//    ss.str("");
//    if ( cspectra->size() ) 
//    {
//      ss << "SELECT target_m_to_z FROM MSMSFraction WHERE peak_list = " << (*cspectra->rbegin())->id();
//      parentp_->adapterp_->executeQuery(ss.str(),false);
//      sqlres = parentp_->adapterp_->lastResult();
//      sqlres.first();
//   
//      if ( !finished && sqlres.value(0).toDouble() - overlapsize <= *startmz )
//      {
//        throw Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__,"too small steps","with the current ovlerlapsize and stepsize the clustering will never be complete");
//      }
//      else
//      {
//        *startmz = sqlres.value(0).toDouble() - overlapsize;
//      }
//    }
//    ss.str("");
//    ss << "DROP TABLE tmpDataSet";
//    parentp_->adapterp_->executeQuery(ss.str(),false);

    return cspectra;
  }
  
  /**
  remove clusters that overlap and return the spectra for reevaluation
  */
  std::vector<ClusterSpectrum*>* ClusterExperiment::ClusterRun::getOverlap(map<int,ClusterNode*>& clusters) const
  {

    StopWatch watch;
    watch.start();

    vector<int> contents;
    for ( map<int,ClusterNode*>::const_iterator cmit = clusters.begin(); cmit != clusters.end(); ++cmit )
    {
      for ( std::list<int>::const_iterator cit = cmit->second->children().begin(); cit != cmit->second->children().end(); ++cit )
      {
        contents.push_back(*cit);
      }
    }

    sort(contents.begin(), contents.end());
    
    vector<int> duplicate;
   
    int oldid = 0;
    for ( vector<int>::const_iterator cvit = contents.begin(); cvit != contents.end(); ++cvit )
    {
      if ( *cvit == oldid )
      {
        duplicate.push_back(*cvit);
      }
      oldid = *cvit;
    }
    map<int,int> overlap;

    for ( map<int,ClusterNode*>::iterator mit = clusters.begin(); mit != clusters.end(); ++mit )
    {
      bool found = false;
      for ( std::list<int>::const_iterator cit = mit->second->children().begin(); cit != mit->second->children().end(); ++cit )
      {
        if ( find(duplicate.begin(),duplicate.end(),*cit) != duplicate.end()) found = true;
      }
      if ( found )
      {
        for ( std::list<int>::const_iterator cit = mit->second->children().begin(); cit != mit->second->children().end(); ++cit )
        {
          overlap.insert(make_pair(*cit,1));
        }
        delete mit->second;
        clusters.erase(mit);
      }
    }
    
//TODO Persistence
//    stringstream ss;
//    ss << "CREATE TEMPORARY TABLE tmpOverlap ( id bigint ) ";
//    parentp_->adapterp_->executeQuery(ss.str(),false);
//    for ( map<int,int>::const_iterator cmit = overlap.begin(); cmit != overlap.end(); ++cmit )
//    {
//      ss.str("");
//      ss << "INSERT INTO tmpOverlap VALUES ( " << cmit->first << ")";
//      parentp_->adapterp_->executeQuery(ss.str(),false);
//    }
//
    vector<ClusterSpectrum*>* result = new vector<ClusterSpectrum*>();
//TODO Persistence
//
//    ss.str("");
//    ss << "SELECT MSMSFraction.peak_list FROM MSMSFraction JOIN tmpOverlap on tmpOverlap.id = MSMSFraction.peak_list ORDER BY MSMSFraction.target_m_to_z";
//    parentp_->adapterp_->executeQuery(ss.str(),false);
//    QSqlQuery sqlres = parentp_->adapterp_->lastResult();
//    
//    
//    while ( sqlres.next() )
//    {
//      result->push_back(new ClusterSpectrum(sqlres.value(0).toInt(),parentp_->adapterp_,binsize_,binspread_));
//    }
//
//    ss.str("");
//    ss << "DROP TABLE tmpOverlap";
//    parentp_->adapterp_->executeQuery(ss.str(),false);
//    
    return result;
  }
  
  double ClusterExperiment::ClusterRun::similarity(const ClusterSpectrum& a, const ClusterSpectrum& b) const
  {
    double autocorr_a = (*sim_funcp_)(a,a);
    double autocorr_b = (*sim_funcp_)(b,b);
    return similarity(a,b,autocorr_a,autocorr_b);
  }
  
  double ClusterExperiment::ClusterRun::similarity(const ClusterSpectrum& a, const ClusterSpectrum& b,double autocorr_a, double autocorr_b) const
  {
    if ( norm_ == none ) return (*sim_funcp_)(a,b);
    if ( autocorr_a < 1e-8 )
    {
      cerr << a.id() << " has self similarity == 0 ! Thats unlikely\n";
      return 0;
    }
    if ( autocorr_b < 1e-8 )
    {
      cerr << b.id() << " has self similarity == 0 ! Thats unlikely\n";
      return 0;
    }
    double result;
    if ( norm_ == arithmetic )
    {
      result = (*sim_funcp_)(a,b)/(autocorr_a/2+autocorr_b/2);
    }
    else if ( norm_ == geometric )
    {
      result = (*sim_funcp_)(a,b)/sqrt(autocorr_a*autocorr_b);
    }
    else
    {
      throw Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__,"unknown mean","dont know what mean to use for similarity");
    }
    if ( isnan(result) ) 
    {
      cerr << "result = " << result << " autocorr_a = " << autocorr_a << " autocorr_b " << autocorr_b << endl;
      return 0;
    }

    return result;
  }

	void ClusterExperiment::ClusterRun::persistentWrite(PersistenceManager& pm, const char* name) const throw (Exception::Base)
	{
		pm.writeObjectHeader(this,name);
		//TODO Persistence
		pm.writeObjectTrailer(name);
	}
	
	void ClusterExperiment::ClusterRun::persistentRead(PersistenceManager& pm) throw (Exception::Base)
	{
		//TODO Persistence
		int dummy;
		pm.readPrimitive(dummy,"dummy_");
	}

}
