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
// $Id: ClusterExperiment.h,v 1.17 2006/05/30 15:46:38 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:  $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_COMPARISON_CLUSTERING_CLUSTEREXPERIMENT_H
#define OPENMS_COMPARISON_CLUSTERING_CLUSTEREXPERIMENT_H

#include <vector>
#include <map>
#include <list>
#include <fstream>
#include <cmath>
#include <iostream>

//these are Clustered
#include <OpenMS/COMPARISON/CLUSTERING/ClusterSpectrum.h>

//the Clusters are represented as ClusterNodes
#include <OpenMS/COMPARISON/CLUSTERING/ClusterNode.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/DataSetInfo.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

namespace OpenMS
{

  // forward declarations
  class FactoryProduct;
  class ClusterFunctor;
  class AnalysisFunctor;
  class CompareFunctor;
  class PreprocessingFunctor;
  class DBAdapter;
  class Cluster;
  
  enum Norm { arithmetic, geometric, none };
  /**
  a ClusterExperiment consists of a DataSet and several (at least one)
  ClusterRun (a sublcass of ClusterExperiment)<br>
  The rationale behind this is to allow to test many different parameter settings
  and/or different Functors (Similarity, Clustering, Preprocessing...) in a easy to manage way<br>
  ClusterExperiments can be saved to xml
  this allows the user to create and run the Object on a fast Computer and analyze the Results
  later, possibly graphically<br>
  */
  class ClusterExperiment
  {
  public:
    //Exceptions

    /**
    CanNotRun is thrown if the Cluster is not ready to run
    this happens if some vital Parts miss, like DBAdapter ,
    CompareFunctor, ClusterFunctor or Data
    */
    class CanNotRun
      :public OpenMS::Exception::Base
    {
    public:
      /** @name constructors and destructors */
      //@{
      CanNotRun(const char* file, int line, const char* function) throw();
      CanNotRun(const char* message, const char* file, int line, const char* function) throw();
      ~CanNotRun() throw();
      //@}
    };

    /**
    NoClusterRun is thrown if an operation was performed on
    the ClusterExperiment that needs a ClusterRun but the
    ClusterExperiment contains no clusterrun<br>
    use ClusterExperiment::createrun() to create runs
    */
    class NoClusterRun
      :public OpenMS::Exception::Base
    {
      public:
        /** @brief constructor <br> */
        NoClusterRun(const char* file, int line, const char* function) throw();
	/** @brief destructor <br> */
        ~NoClusterRun() throw();
    };

    /**
    //represents a whole Analysis, consisting of a AnalysisFunctor and the corresponding result
    */
    class Analysis : public PersistentObject
    {
      friend class ClusterExperiment;
      friend class ClusterExperimentXMLHandler;
      friend class DBAdapter;
    public:
        /** @name constructors and destructors */
	//@{
        Analysis(AnalysisFunctor* func);
        Analysis(const Analysis& source);
        ~Analysis();
	//@}

	/** @brief assignment operator <br> */
        Analysis& operator = (const Analysis& source);

	/** @brief needs the analysis to be run? <br> */
        bool done() const;

	/** @brief perform analysis <br> */
        void run(std::map<int,ClusterNode*> clusters);

	/** @brief name of AnalysisFunctor <br> */
        String name() const;

	/** @brief read access to Analysis Functor <br> */
        const FactoryProduct* anafuncp() const;

	/** @brief results of AnalysisFunctor <br> */
        const std::map<String,double>& results() const {return result_;}

	/** @brief set DBAdapter for AnalysisFunctor <br> */
        void setAdapter(DBAdapter* adapterp);

	/** @brief save to xml <br> */
        void save(std::ostream& document, int& ind);

        /** @brief for internal use <br> */
        void setDataSet(uint);

			// Docu in base class
			virtual void persistentWrite(PersistenceManager& pm, const char* name=0) const throw (Exception::Base);
			
			// Docu in base class
			virtual void persistentRead(PersistenceManager& pm) throw (Exception::Base);

			Analysis() {} //TODO Persistence : make private again
				
		protected:
			// Docu in base class
	    virtual void clearChildIds_()
	    {
	    	//TODO Persistence	
	    };
	    

    private:

        
        AnalysisFunctor* anafuncp_;
        std::map<String,double> result_;
        bool analyzed_;
        uint dataset_;
    };
    /**
    represents a single ClusterRun <br>
    The ClusterRuns in one ClusterExperiment all use the
    same DataSet and the same DataBase
    everything else can be varied, like the clustering algorithm,
    the representation of the spectra ( BinnedRep or not, the size of
    bins etc), Preprocessing and the Analysis used on the resulting
    Clustering
    */
    class ClusterRun: public PersistentObject
    {
      friend class DBAdapter;
      friend class ClusterExperimentXMLHandler;
    public:
      /** @name constructors, destructor, assigment operator <br> */
      //@{
      ClusterRun();
      ClusterRun(const ClusterExperiment& parent);
      ClusterRun(const ClusterRun& source);
      virtual ~ClusterRun();
      ClusterRun& operator = (const ClusterRun& source);
      //@}
      /** @brief clear ClusterRun <br> */
      void deleteContents();

      /** @name write accessors */
      //@{
      void setBinSize(double);
      void setBinSpread(uint);
      void setNorm(Norm);
      int addMower(PreprocessingFunctor*);
      void setSimFunc(CompareFunctor* );
      void setClusterFunc(ClusterFunctor*);
      int addAnalysisFunctor(AnalysisFunctor*);
      //@}

      /** @name read accessors <br> */
      //@{
      const double& getBinSize() const {return binsize_;}
      const uint& getBinSpread() const { return binspread_;}
      const Norm& getNorm() const { return norm_;}
      const std::map<int,ClusterNode*>& getClustering() const { return clusters_;}
      inline const CompareFunctor* getSimFunc() const {return sim_funcp_;}
      const std::vector<PreprocessingFunctor*>& getPreprocessqueue() const {return preprocess_queue_;}
      const ClusterFunctor* getClusterFunc() const { return cluster_funcp_;}
      /** @brief read access to the underlying ClusterExperiment::Analysis <br> */
      const Analysis& operator[](uint pos) const throw(Exception::IndexOverflow);

      //@}
     
      /** @brief set clusters manually */
      void setClustering(const std::map<int,ClusterNode*>& clusters);
      
      /** @brief preprocess stick spectrum <br> */
      void preprocess(MSSpectrum< DPeak<1> >& ) const;

      /** @brief get spectra to cluster <br> */
      std::vector<ClusterSpectrum*>* getSpectra(uint charge = 0, bool* finished = 0, uint stepsize = 0 ,double* startmz = 0 ) const;

      /** @brief get ambigously clustered spectra <br> */
      std::vector<ClusterSpectrum*>* getOverlap(std::map<int,ClusterNode*>&) const;

      /** @brief return similarity between 0 and 1<br> */
      double similarity(const ClusterSpectrum& a, const ClusterSpectrum& b) const;

      /** @brief return similarity between 0 and 1<br> */
      double similarity(const ClusterSpectrum& a, const ClusterSpectrum& b,double autocorr_a, double autocorr_b) const;

      /** @brief number of ClusterExperiment::Analysis <br> */
      uint size() const { return analysis_queue_.size();}

      /** @brief write to xml <br> */
      void save(std::ostream& document, int& ind );

      /** @brief check if ClusterRun has enough information to run <br> */
      bool canRun() const;

      /** @brief more verbose version of canRun() <br> */
      int isComplete() const;

      /** @brief cluster spectra <br> */
      void cluster();

      /** @brief run Analysis <br> */
      void analyze();

      /** @brief cluster() then analyze() <br> */
      void run() throw(CanNotRun);

      //temp todo
      void resetnr()const {nr_ = 0;}
			
			// Docu in base class
			virtual void persistentWrite(PersistenceManager& pm, const char* name=0) const throw (Exception::Base);
			
			// Docu in base class
			virtual void persistentRead(PersistenceManager& pm) throw (Exception::Base);
	
		protected:
			// Docu in base class
	    virtual void clearChildIds_()
	    {
	    	//TODO Persistence	
	    };
    
    private:

      void deleteCachedClusters_();

      const ClusterExperiment* parentp_;

      std::map<int,ClusterNode*> clusters_;

      /** @name bin dimensions <br> */
      //@{
      double binsize_;
      uint binspread_;
      //@}

      bool didrun_;
      CompareFunctor* sim_funcp_;
      std::vector<PreprocessingFunctor*> preprocess_queue_;
      ClusterFunctor* cluster_funcp_;
      std::vector<Analysis> analysis_queue_;
      Norm norm_;

      //temp todo
      mutable uint nr_;
    };

    /**
    Predicate Object to compare ClusterRuns based on some Analysis
    */
    class ClusterRunAnalysisLess
    {
    public:
      /** @brief constructor <br> */
      ClusterRunAnalysisLess(String cfig, String param);

      /** @brief destructor <br> */
      ~ClusterRunAnalysisLess();

      /** @brief comparison <br> */
      bool operator()(const ClusterRun* ap, const ClusterRun* bp);

      /** @brief */
      void setRequirement(String param, double value);
    private:
      String configurablename_;
      String paramname_;
      std::map<String, double> requirements_;
    };

    friend class ClusterExperimentXMLHandler;
    friend class ClusterRun;
  public:
    /** @brief standard constructor <br> */
    ClusterExperiment();

    /** @brief copy constructor <br> */
    ClusterExperiment(const ClusterExperiment& source);

    /** @brief destructor <br> */
    ~ClusterExperiment();

    /** @brief assignment operator <br> */
    ClusterExperiment& operator = (const ClusterExperiment& source);

    /** @brief clear <br> */
    void deleteContents();

    /** @name write accessors for the current run <br> */
    //@{
    void setBinSize(double size, int pos = -1);
    void setBinSpread(uint spread, int pos = -1);
    void setNorm(Norm norm, int pos = -1);
    int addMower(PreprocessingFunctor* funcp, int pos = -1);
    void setClusterFunc(ClusterFunctor* funcp, int pos = -1);
    void setSimFunc(CompareFunctor* funcp, int pos = -1);
    int addAnalysisFunctor(AnalysisFunctor* func, int pos = -1);
    //@}

    /** @brief create new run <br> */
    int createrun();

    /** @brief return number of ClusterRuns <br> */
    uint size() const;

    /** @brief read access to ClusterRuns <br> */
    const ClusterRun& operator[](uint pos) const throw(Exception::IndexOverflow);

    ClusterRun* getClusterRun(uint pos) throw(Exception::IndexOverflow);

    /** @brief save to xml <br> */
    void save(String);

    /** @brief load from xml <br> */
    void load(String);

    /** @name write accessors for dataset <br> */
    //@{
    void setDBAdapter(DBAdapter*);
    void setDataSetInfo(DataSetInfo* dsip);
    void setInfoUser(String name);
    void setInfoComment(String comment);
    //@}

    /** @name read accessors for dataset <br> */
    //@{
    const String& infouser() const ;
    const String& infocomment() const ;
    const String& infodate() const;
    String datasetname() const;
    long datasetid() const;
    uint datasetsize() ; //const?
    //@}

    /** @brief relayed to ClusterRun::canRun at position <i>nr</i> <br> */
    bool canRun(int nr = -1) const;

    /** @brief relayed to ClusterRun::run at position <i>nr</i> <br> */
    void run(int nr = -1);

    /** @brief relayed to ClusterRun::cluster at position <i>nr</i> <br> */
    void cluster( int nr = -1 );

    /** @brief relayed to ClusterRun::analyze at position <i>nr</i> <br> */
    void analyze ( int nr = -1 );

    /** @brief sort by Analysis result <br> */
    void sortbyResult(String name, String param, bool desc = true);

  private:

    const DataSetInfo* getDataSetInfo() const;

    void clear_();

    mutable DBAdapter* adapterp_;

    //data
    String datasetname_;
    mutable DataSetInfo* dsip_;

    //info
    String info_user_;
    String info_date_;
    String info_comment_;

    //runs
    std::vector<ClusterRun*> runps_;
    mutable int currentrun_;
  };

}
#endif //OPENMS_COMPARISON_CLUSTERING_CLUSTEREXPERIMENT_H
