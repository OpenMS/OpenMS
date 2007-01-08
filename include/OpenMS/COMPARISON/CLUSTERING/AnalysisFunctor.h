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

#ifndef OPENMS_COMPARISON_CLUSTERING_ANALYSISFUNCTOR_H
#define OPENMS_COMPARISON_CLUSTERING_ANALYSISFUNCTOR_H

#include <vector>
#include <map>

#include <OpenMS/CONCEPT/FactoryProduct.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterExperiment.h>

namespace OpenMS
{
	class DBAdapter;
	class ClusterNode;
	
  /**
	  @brief AnalysisFunctor classes analyze the output of a ClusterRun
	  
	  @todo write test
  */
  class AnalysisFunctor : public FactoryProduct
  {
  public:
    /** @brief standard constructor <br>*/
    AnalysisFunctor();
    /** @brief copy constructor <br>*/
    AnalysisFunctor(const AnalysisFunctor& source);
    /** @brief assignment operator <br>*/
    AnalysisFunctor& operator=(const AnalysisFunctor& source);
    /** @brief destructor <br>*/
    virtual ~AnalysisFunctor();

    /** @brief function call operator <br>*/
    /**
    needs to be implemented in subclasses<br>
    the actual analysis is performed here<br>
    */
    virtual std::map<String,double> operator()(const std::map<int,ClusterNode*>&)  = 0;
    /** @brief write accessor for the DBAdapter* <br>*/
    void setDBAdapter(DBAdapter* adapter);

    /** @brief set an additional ClusterRun, if needed <br> */
    void setClusterRun(const ClusterExperiment::ClusterRun* crp);

    /**
    check if all prerequisited ( adapter, clusterun ) are met
    */
    virtual void canRun() const;

		static void registerChildren();
		
    /** @brief return if the class needs a connection to the DB <br> */
    bool needsDBAdapter() const {return needsAdapter_;}

    /** @brief return if the class needs an additional ClusterRun<br> */
    bool needsClusterRun() const {return needsClusterrun_;}

    /** @brief return the additional ClusterRun <br>*/
    const ClusterExperiment::ClusterRun* reference() const {return clusterrunp_;}
    
  public:
    /** TRUE if <b>this</b> needs an DBAdapter* <br>*/
    bool needsAdapter_;
    /** access to the DataBase <br>*/
    DBAdapter* adapterp_;

    bool needsClusterrun_;

    const ClusterExperiment::ClusterRun* clusterrunp_;
  };

}
#endif //OPENMS_COMPARISON_CLUSTERING_ANALYSISFUNCTOR_H
