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

#include <OpenMS/COMPARISON/CLUSTERING/ClusterExperimentXMLHandler.h>

#include <iostream>
#include <cassert>

#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/COMPARISON/CLUSTERING/AnalysisFunctor.h>
#include <OpenMS/FILTERING/TRANSFORMERS/PreprocessingFunctor.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterFunctor.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterNode.h>

#include <xercesc/sax2/Attributes.hpp>

using namespace xercesc;
using namespace std;

namespace OpenMS
{
  ClusterExperimentXMLHandler::ClusterExperimentXMLHandler(ClusterExperiment& crun, const String& filename)
    : XMLHandler(filename),
    	crun_(crun),
    	forwardconfigurablep_(0),
    	isclusterexperiment_(0),
    	tmpclusterp_(0),
    	ismemberid_(0),
    	iscomment_(0),
    	clid_(-1),
    	cl_min_mass_(0),
    	cl_max_mass_(0)
  {
		file_ = __FILE__;
  }

  ClusterExperimentXMLHandler::~ClusterExperimentXMLHandler()
  {
  }

  void ClusterExperimentXMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const Attributes& attributes)
  {
  	String qname_string = XMLString::transcode(qname);
  	
    // if a FactoryProduct is found its parameters are set here
    if (forwardconfigurablep_)
    {
      if (qname_string =="param")
      {
        //should not happen since attributes in param are required
        //but in case validity is not checked (it is not)
        if (attributes.getIndex(XMLString::transcode("name"))!=-1)
        {
          forwardconfigurablep_->getParam().setValue(XMLString::transcode(attributes.getValue(XMLString::transcode("name"))),asDouble_(XMLString::transcode(attributes.getValue(XMLString::transcode("value")))));
        }
      }
    }
    //simple way of checking if it is a ClusterExperiment xml file
    //real error checking should (hopefully) not be needed since the files are only written by
    //ClusterExperiment.save()
    else if (qname_string =="ClusterExperiment" )
    {
      isclusterexperiment_ = 1;
    }
    else if (!isclusterexperiment_)
    {
			const Locator* loc = 0;
			setDocumentLocator(loc);
			String message = String("Invalid file or not a ClusterExperiment file!");
			error(SAXParseException(XMLString::transcode(message.c_str()), *loc ));
    }
    //top elements
    else if (qname_string =="Data" )
    {
      const char* name = XMLString::transcode(attributes.getValue(XMLString::transcode("Name")));
      if ( name != 0 )
      {
        crun_.datasetname_ = String(name);
      }
    }
    else if (qname_string =="Info")
    {
      const char* user = XMLString::transcode(attributes.getValue(XMLString::transcode("User")));
      const char* date = XMLString::transcode(attributes.getValue(XMLString::transcode("Date")));
      if ( user != 0 )
      {
        crun_.info_user_ = String(user);
      }
      if ( date != 0 )
      {
        crun_.info_date_ = String(date);
      }
    }
    else if (qname_string =="Norm")
    {
      if ( String(XMLString::transcode(attributes.getValue(XMLString::transcode("mean")))) == "arithmetic" )
      {
        crun_.setNorm(arithmetic);
      }
      else if ( String(XMLString::transcode(attributes.getValue(XMLString::transcode("mean")))) == "none" )
      {
        crun_.setNorm(none);
      }
      else if ( String(XMLString::transcode(attributes.getValue(XMLString::transcode("mean")))) == "geometric" )
      {
        crun_.setNorm(geometric);
      }
    }
    else if (qname_string =="Bins")
    {
      double size = 0;
      int spread = -1;
      size = asDouble_(XMLString::transcode(attributes.getValue(XMLString::transcode("size")))); 
      spread = asUnsignedInt_(XMLString::transcode(attributes.getValue(XMLString::transcode("spread"))));
      if ( size > 1e-8 )
      {
        crun_.setBinSize(size);
        crun_.setBinSpread(spread);
      }
    }
    else if (qname_string =="Comment")
    {
      iscomment_ = 1;
    }
    else if (qname_string =="ClusterRun")
    {
      crun_.createrun();
      didrun_ = asSignedInt_(XMLString::transcode(attributes.getValue(XMLString::transcode("finished"))));
    }
    else if (qname_string =="Preprocessing")
    {
      //noop
      //FilterFunc is a child of Preprocessing for cosmetic reasons
    }
    else if (qname_string =="FilterFunc" )
    {
      const char* name = XMLString::transcode(attributes.getValue(XMLString::transcode("name")));
      if (name != 0 )
      {
        forwardconfigurablep_ = Factory<PreprocessingFunctor>::instance()->create(name);
        if (!forwardconfigurablep_)
        {
					const Locator* loc = 0;
					setDocumentLocator(loc);
					String message = String("unknown FactoryProduct \"") + name + "\"";
					error(SAXParseException(XMLString::transcode(message.c_str()), *loc ));
        }
      }
    }
    else if (qname_string =="SimFunc")
    {
      const char* name = XMLString::transcode(attributes.getValue(XMLString::transcode("name")));
      if (name != 0)
      {
        //forwardconfigurablep_ = Factory<CompareFunctor>::instance()->create(name);
        if (!forwardconfigurablep_)
        {
					const Locator* loc = 0;
					setDocumentLocator(loc);
					String message = String("unknown FactoryProduct \"") + name + "\"";
					error(SAXParseException(XMLString::transcode(message.c_str()), *loc ));
				}
      }
    }
    else if (qname_string =="ClustFunc" )
    {
      const char* name = XMLString::transcode(attributes.getValue(XMLString::transcode("name")));
      if (name != 0 )
      {
        //forwardconfigurablep_ = Factory<CompareFunctor>::instance()->create(name);
        if (!forwardconfigurablep_)
        {
					const Locator* loc = 0;
					setDocumentLocator(loc);
					String message = String("unknown FactoryProduct \"") + name + "\"";
					error(SAXParseException(XMLString::transcode(message.c_str()), *loc ));
        }
      }
    }
    else if (qname_string =="Clustering")
    {
      //noop, just for cosmetical purposes
    }
    else if (qname_string =="Cluster")
    {
      for ( UnsignedInt i = 0 ; i < attributes.getLength(); ++i )
      {
        String attributesValue =  XMLString::transcode(attributes.getValue(i)) ;
        String attributesName = XMLString::transcode(attributes.getQName(i));
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
    else if (qname_string =="member_id" )
    {
      ismemberid_ = 1;  
    }
    else if (qname_string =="Evaluation")
    {
      //noop, just for cosmetical purposes
    }
    else if (qname_string =="Analysis")
    {
    }
    else if (qname_string =="AnaFunc")
    {
      const char* name = XMLString::transcode(attributes.getValue(XMLString::transcode("name")));
      if (name != 0 )
      {
        forwardconfigurablep_ = Factory<AnalysisFunctor>::instance()->create(name);
        if (!forwardconfigurablep_)
        {
					const Locator* loc = 0;
					setDocumentLocator(loc);
					String message = String("unknown FactoryProduct \"") + name + "\"";
					error(SAXParseException(XMLString::transcode(message.c_str()), *loc ));
        }
      }
    }
    else if (qname_string =="Results")
    {
    }
    else if ( qname_string =="Result")
    {
      const char* rString = XMLString::transcode(attributes.getValue(XMLString::transcode("string")));
      double rdouble = asDouble_(XMLString::transcode(attributes.getValue(XMLString::transcode("double"))));
			if ( rString != 0 )
      {
        crun_.runps_[crun_.currentrun_]->analysis_queue_.back().result_.insert(make_pair(rString,rdouble));
      }
    }
    else if ( qname_string =="Reference" )
    {
      uint clusterrun = asUnsignedInt_(XMLString::transcode(attributes.getValue(XMLString::transcode("ClusterRunNr"))));
      uint ananr = asUnsignedInt_(XMLString::transcode(attributes.getValue(XMLString::transcode("AnalysisFunctorNr")))); 
      uint reference = asUnsignedInt_(XMLString::transcode(attributes.getValue(XMLString::transcode("references"))));
      if ( clusterrun >= crun_.runps_.size() || crun_.runps_[clusterrun]->size() <= ananr || reference >= crun_.runps_.size() )
      {
				const Locator* loc = 0;
				setDocumentLocator(loc);
				String message = String("AnalysisFunctorReferences inconsistent!");
				error(SAXParseException(XMLString::transcode(message.c_str()), *loc ));
      }
      (*crun_.runps_[clusterrun])[ananr].anafuncp_->setClusterRun(crun_.runps_[reference]);
    }
  }

  void ClusterExperimentXMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
	{
		String qname_string = XMLString::transcode(qname);
		
    if (qname_string =="ClusterExperiment" )
    {
      isclusterexperiment_ = 0;
    }
    else if (!isclusterexperiment_)
    {
			const Locator* loc = 0;
			setDocumentLocator(loc);
			String message = String("Invalid file or not a ClusterExperiment file!");
			error(SAXParseException(XMLString::transcode(message.c_str()), *loc ));
    }
    //top elements
    else if (qname_string =="FilterFunc" )
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
    else if (qname_string =="ClusterRun")
    {
      (*crun_.runps_.rbegin())->didrun_ = didrun_;
      didrun_ = 0;
    }
    else if (qname_string =="SimFunc")
    {
      if (forwardconfigurablep_)
      {
        PeakSpectrumCompareFunctor* tmp = 0/*dynamic_cast<CompareFunctor*>(forwardconfigurablep_)*/;
        if (tmp) crun_.setSimFunc(tmp);
        forwardconfigurablep_ = 0;
      }
    }
    else if (qname_string =="ClustFunc" )
    {
      if ( forwardconfigurablep_)
      {
        //if configurablep_ is really a PreprocessingFunctor* has already been checked in startElement
        ClusterFunctor* tmp = dynamic_cast<ClusterFunctor*>(forwardconfigurablep_);
        crun_.setClusterFunc(tmp);
        forwardconfigurablep_ = 0;
      }
    }
    else if (qname_string =="Clustering")
    {
      //noop, just for cosmetical purposes
    }
    else if (qname_string =="Cluster")
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
    else if (qname_string =="member_id" )
    {
      ismemberid_ = 0;
    }
    else if (qname_string =="AnaFunc")
    {
      if ( forwardconfigurablep_)
      {
        //if configurablep_ is really a PreprocessingFunctor* has already been checked in startElement
        AnalysisFunctor* tmp = dynamic_cast<AnalysisFunctor*>(forwardconfigurablep_);
        crun_.addAnalysisFunctor(tmp);
        forwardconfigurablep_ = 0;
      }
    }
  }

  void ClusterExperimentXMLHandler::characters(const XMLCh* const chars, const unsigned int /*length*/)
  {
    if (ismemberid_ )
    {
      if (tmpclusterp_)
      {
        tmpclusterp_ = new ClusterNode(tmpclusterp_, new ClusterNode(asSignedInt_(XMLString::transcode(chars))));
      }
      else 
      {
        tmpclusterp_ = new ClusterNode(asSignedInt_(XMLString::transcode(chars)));
      }
    }
    else if (iscomment_)
    {
      crun_.info_comment_ = XMLString::transcode(chars);
      iscomment_ = 0;
    }
  }

}
