//-*- Mode: C++; tab-width: 2; -*-
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
#include <OpenMS/COMPARISON/CLUSTERING/ClusterExperimentXMLHandler.h>

#include <iostream>
#include <cassert>

#include <OpenMS/COMPARISON/CLUSTERING/ClusterFactory.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterFunctor.h>
#include <OpenMS/COMPARISON/CLUSTERING/AnalysisFunctor.h>
#include <OpenMS/FILTERING/TRANSFORMERS/PreprocessingFunctor.h>
#include <OpenMS/COMPARISON/SPECTRA/CompareFunctor.h>

using namespace std;

namespace OpenMS
{
  ClusterExperimentXMLHandler::ClusterExperimentXMLHandler(ClusterExperiment& crun)
    : XMLHandler(),crun_(crun), forwardconfigurablep_(0),isclusterexperiment_(0),tmpclusterp_(0),ismemberid_(0),iscomment_(0),clid_(-1),cl_min_mass_(0),cl_max_mass_(0)
  {
		file_ = __FILE__;
  }

  ClusterExperimentXMLHandler::~ClusterExperimentXMLHandler()
  {
  }

  bool ClusterExperimentXMLHandler::startElement(const QString & /*uri*/, const QString & /*local_name*/, 
																const QString & qname, const QXmlAttributes & attributes )
  {
    // if a FactoryProduct is found its parameters are set here
    if ( forwardconfigurablep_)
    {
      if (qname == "param")
      {
        const char* name = attributes.value("name").ascii();
        double value = asDouble_(attributes.value("value"));
        //should not happen since attributes in param are required
        //but in case validity is not checked (it is not)
        if (name != 0)
        {
          forwardconfigurablep_->getParam().setValue(name,value);
        }
      }
    }
    //simple way of checking if it is a ClusterExperiment xml file
    //real error checking should (hopefully) not be needed since the files are only written by
    //ClusterExperiment.save()
    else if (qname == "ClusterExperiment" )
    {
      isclusterexperiment_ = 1;
    }
    else if (!isclusterexperiment_) return false;
    //top elements
    else if (qname == "Data" )
    {
      const char* name = attributes.value("Name").ascii();
      if ( name != 0 )
      {
        crun_.datasetname_ = String(name);
      }
    }
    else if (qname == "Info")
    {
      const char* user = attributes.value("User").ascii();
      const char* date = attributes.value("Date").ascii();
      if ( user != 0 )
      {
        crun_.info_user_ = String(user);
      }
      if ( date != 0 )
      {
        crun_.info_date_ = String(date);
      }
    }
    else if (qname == "Norm")
    {
      if ( attributes.value("mean") == "arithmetic" )
      {
        crun_.setNorm(arithmetic);
      }
      else if ( attributes.value("mean") == "none" )
      {
        crun_.setNorm(none);
      }
      else if ( attributes.value("mean") == "geometric" )
      {
        crun_.setNorm(geometric);
      }
    }
    else if (qname == "Bins")
    {
      double size = 0;
      int spread = -1;
      size = asDouble_(attributes.value("size")); 
      spread = asUnsignedInt_(attributes.value("spread"));
      if ( size > 1e-8 )
      {
        crun_.setBinSize(size);
        crun_.setBinSpread(spread);
      }
    }
    else if (qname == "Comment")
    {
      iscomment_ = 1;
    }
    else if (qname == "ClusterRun")
    {
      crun_.createrun();
      didrun_ = asSignedInt_(attributes.value("finished"));
    }
    else if (qname == "Preprocessing")
    {
      //noop
      //FilterFunc is a child of Preprocessing for cosmetic reasons
    }
    else if (qname == "FilterFunc" )
    {
      const char* name = attributes.value("name").ascii();
      if (name != 0 )
      {
        forwardconfigurablep_ = dynamic_cast<PreprocessingFunctor*>(ClusterFactory::instance()->create(name));
        if (!forwardconfigurablep_)
        {
          throw Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__,"unknown FactoryProduct",name);
        }
      }
    }
    else if (qname == "SimFunc")
    {
      const char* name = attributes.value("name").ascii();
      if (name != 0)
      {
        forwardconfigurablep_ = dynamic_cast<CompareFunctor*>(ClusterFactory::instance()->create(name));
        if (!forwardconfigurablep_)
        {
          throw Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__,"unknown FactoryProduct",name);
				}
      }
    }
    else if (qname == "ClustFunc" )
    {
      const char* name = attributes.value("name").ascii();
      if (name != 0 )
      {
        forwardconfigurablep_ = dynamic_cast<ClusterFunctor*>(ClusterFactory::instance()->create(name));
        if (!forwardconfigurablep_)
        {
          throw Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__,"unknown FactoryProduct",name);
        }
      }
    }
    else if (qname == "Clustering")
    {
      //noop, just for cosmetical purposes
    }
    else if (qname == "Cluster")
    {
      for ( int i = 0 ; i < attributes.length(); ++i )
      {
        QString attributesValue =  attributes.value(i) ;
        QString attributesName = attributes.qName(i);
        if( attributesName == "id" )
        {
          clid_ = asSignedInt_(attributesValue);
        }
        else if( attributesName == "mininum_parent_mass" )
        {
          cl_min_mass_ = asDouble_(attributesValue);
        }
        else if( attributesName == "maximum_parent_mass")
        {
          cl_max_mass_ = asDouble_(attributesValue);
        }    
      }
    }
    else if (qname == "member_id" )
    {
      ismemberid_ = 1;  
    }
    else if (qname == "Evaluation")
    {
      //noop, just for cosmetical purposes
    }
    else if (qname == "Analysis")
    {
    }
    else if (qname == "AnaFunc")
    {
      const char* name = attributes.value("name").ascii();
      if (name != 0 )
      {
        forwardconfigurablep_ = dynamic_cast<AnalysisFunctor*>(ClusterFactory::instance()->create(name));
        if (!forwardconfigurablep_)
        {
          throw Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__,"unknown FactoryProduct",name);
        }
      }
    }
    else if (qname == "Results")
    {
    }
    else if ( qname == "Result")
    {
      const char* rString = attributes.value("string").ascii();
      double rdouble = asDouble_(attributes.value("double"));
			if ( rString != 0 )
      {
        crun_.runps_[crun_.currentrun_]->analysis_queue_.back().result_.insert(make_pair(rString,rdouble));
      }
    }
    else if ( qname == "Reference" )
    {
      uint clusterrun = asUnsignedInt_(attributes.value("ClusterRunNr"));
      uint ananr = asUnsignedInt_(attributes.value("AnalysisFunctorNr")); 
      uint reference = asUnsignedInt_(attributes.value("references"));
      if ( clusterrun >= crun_.runps_.size() || crun_.runps_[clusterrun]->size() <= ananr || reference >= crun_.runps_.size() )
      {
        cerr << "AnalysisFunctorReferences inconsistent, ignored \n";
        return false;
      }
      (*crun_.runps_[clusterrun])[ananr].anafuncp_->setClusterRun(crun_.runps_[reference]);
    }
		return no_error_;
  }

  bool ClusterExperimentXMLHandler::endElement
	(const QString & /*uri*/, const QString & /*local_name*/, const QString & qname )
	{
    if (qname == "ClusterExperiment" )
    {
      isclusterexperiment_ = 0;
    }
    else if (!isclusterexperiment_) return false;
    //top elements
    else if (qname == "FilterFunc" )
    {
      if ( forwardconfigurablep_)
      {
        //if configurablep_ is really a PreprocessingFunctor* has already been checked in startElement
        PreprocessingFunctor* tmp = dynamic_cast<PreprocessingFunctor*>(forwardconfigurablep_);
        crun_.addMower(tmp);
        forwardconfigurablep_ = 0;
      }
    }
    // todo: why is it just called once? ( xerces-2.6 )
    else if (qname == "ClusterRun")
    {
      (*crun_.runps_.rbegin())->didrun_ = didrun_;
      didrun_ = 0;
    }
    else if (qname == "SimFunc")
    {
      if (forwardconfigurablep_)
      {
        CompareFunctor* tmp = dynamic_cast<CompareFunctor*>(forwardconfigurablep_);
        if (tmp) crun_.setSimFunc(tmp);
        forwardconfigurablep_ = 0;
      }
    }
    else if (qname == "ClustFunc" )
    {
      if ( forwardconfigurablep_)
      {
        //if configurablep_ is really a PreprocessingFunctor* has already been checked in startElement
        ClusterFunctor* tmp = dynamic_cast<ClusterFunctor*>(forwardconfigurablep_);
        crun_.setClusterFunc(tmp);
        forwardconfigurablep_ = 0;
      }
    }
    else if (qname == "Clustering")
    {
      //noop, just for cosmetical purposes
    }
    else if (qname == "Cluster")
    {
      //todo I need a generalized ClusterNode
      //or maybe i need to save the structure in the ClusterNodes
      //todo median id
      assert(clid_ >= 0);
      crun_.runps_[crun_.currentrun_]->clusters_.insert(make_pair(clid_ ,tmpclusterp_));
      tmpclusterp_->min_parent_mass_ = cl_min_mass_;
      tmpclusterp_->max_parent_mass_ = cl_max_mass_;
      tmpclusterp_ =  0;
      clid_ = -1;
      cl_min_mass_ = 0;
      cl_max_mass_ = 0;
    }
    else if (qname == "member_id" )
    {
      ismemberid_ = 0;
    }
    else if (qname == "AnaFunc")
    {
      if ( forwardconfigurablep_)
      {
        //if configurablep_ is really a PreprocessingFunctor* has already been checked in startElement
        AnalysisFunctor* tmp = dynamic_cast<AnalysisFunctor*>(forwardconfigurablep_);
        crun_.addAnalysisFunctor(tmp);
        forwardconfigurablep_ = 0;
      }
    }
		return no_error_;
  }

  bool ClusterExperimentXMLHandler::characters(const QString & chars )
  {
    if (ismemberid_ )
    {
      if (tmpclusterp_)
      {
        tmpclusterp_ = new ClusterNode(tmpclusterp_, new ClusterNode(asSignedInt_(chars)));
      }
      else 
      {
        tmpclusterp_ = new ClusterNode(asSignedInt_(chars));
      }
    }
    else if (iscomment_)
    {
      crun_.info_comment_ = chars.ascii();
      iscomment_ = 0;
    }
		return no_error_;
  }

}
