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
#include <OpenMS/COMPARISON/CLUSTERING/AnalysisFunctor.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterCompareFunctor.h>
#include <OpenMS/COMPARISON/CLUSTERING/DistanceAnalyzer.h>
#include <OpenMS/CONCEPT/Factory.h>

#include <sstream>

using namespace std;

namespace OpenMS
{

  AnalysisFunctor::AnalysisFunctor()
    : FactoryProduct(),
			needsAdapter_(0),
			adapterp_(0),
			needsClusterrun_(0),
			clusterrunp_(0)
  {
		name_ = AnalysisFunctor::getName();
  }

  AnalysisFunctor::AnalysisFunctor(const AnalysisFunctor& source)
  	: FactoryProduct(source),
			needsAdapter_(source.needsAdapter_),
			adapterp_(source.adapterp_),
			needsClusterrun_(source.needsClusterrun_),
			clusterrunp_(source.clusterrunp_)
  {
  }
  
  AnalysisFunctor& AnalysisFunctor::operator=(const AnalysisFunctor& source)
  {
    FactoryProduct::operator=(source);
    needsAdapter_ = source.needsAdapter_;
    adapterp_ = source.adapterp_;
    needsClusterrun_ = source.needsClusterrun_;
    clusterrunp_ = source.clusterrunp_;
    return *this;
  }
  
  AnalysisFunctor::~AnalysisFunctor()
  {
  }
 
  void AnalysisFunctor::setDBAdapter(DBAdapter* adapterp)
  {
    adapterp_ = adapterp;
  }

  /**
   some example usages of uses for an additional ClusterRun: <br>
   - compare two Clusterings (see ClusterCompareFunctor) <br>
   - analyze one clustering with anothers preprocessing,compare settings 
     (see DistanceAnalyzer) <br>
   */
  void AnalysisFunctor::setClusterRun(const ClusterExperiment::ClusterRun* crp)
  {
    clusterrunp_ = crp;
  }
 
	void AnalysisFunctor::registerChildren()
	{
		Factory<AnalysisFunctor>::registerProduct(ClusterCompareFunctor::getName(),&ClusterCompareFunctor::create);
		Factory<AnalysisFunctor>::registerProduct(DistanceAnalyzer::getName(),&DistanceAnalyzer::create);

	}
	
  void AnalysisFunctor::canRun() const
  {
    ostringstream oss;
    if ( needsAdapter_ )
    {
      if ( !adapterp_ )
      {
        oss << this->getName() << " needs a DBAdapter";
        throw Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__,"no DBAdapter",oss.str().c_str());
      }
    }
    if ( needsClusterrun_ )
    {
      if ( !clusterrunp_ )
      {
        oss << this->getName() << " needs a ClusterRun";
        throw Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__,"no ClusterRun",oss.str().c_str());
      }
    }
  }
 
}
