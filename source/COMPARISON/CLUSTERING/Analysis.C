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
// $Id: Analysis.C,v 1.6 2006/03/29 12:30:29 andreas_bertsch Exp $
// $Author: andreas_bertsch $
// $Maintainer:  $
// --------------------------------------------------------------------------
//
#include <OpenMS/COMPARISON/CLUSTERING/ClusterExperiment.h>

#include <OpenMS/COMPARISON/CLUSTERING/AnalysisFunctor.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterFactory.h>

#include <OpenMS/COMPARISON/CLUSTERING/helper.h>

using namespace std;

namespace OpenMS
{
  ClusterExperiment::Analysis::Analysis(AnalysisFunctor* func)
    :PersistentObject(), anafuncp_(func),result_(), analyzed_(0),dataset_(0)
  {
  }
  
  ClusterExperiment::Analysis::Analysis(const Analysis& source)
    :PersistentObject(source),anafuncp_(0),result_(source.result_),analyzed_(source.analyzed_),dataset_(source.dataset_)
  {
    if ( source.anafuncp_) 
    {
      anafuncp_ = dynamic_cast<AnalysisFunctor*>(ClusterFactory::instance()->duplicate(source.anafuncp_));
    }
    else anafuncp_ = 0;
  }

  ClusterExperiment::Analysis::~Analysis()
  {
    delete anafuncp_;
  }

  /**
  dataset_ is only used for serialize <br> 
  */
  void ClusterExperiment::Analysis::setDataSet(uint dsid)
  {
    dataset_ = dsid;
  }
  
  ClusterExperiment::Analysis& ClusterExperiment::Analysis::operator = (const ClusterExperiment::Analysis& source)
  {
    PersistentObject::operator=(source);
    if ( source.anafuncp_) 
    {
      anafuncp_ = dynamic_cast<AnalysisFunctor*>(ClusterFactory::instance()->duplicate(source.anafuncp_));
    }
    else anafuncp_ = 0;
    result_ = source.result_;
    analyzed_ = source.analyzed_;
    dataset_ = source.dataset_;
    return *this;
  }
 
  const FactoryProduct* ClusterExperiment::Analysis::anafuncp() const 
  {
    return anafuncp_;
  }
  
  void ClusterExperiment::Analysis::run(map<int,ClusterNode*> clusters)
  {
    if (!analyzed_)
    {
      result_ = (*anafuncp_)(clusters);
    }
  }

  bool ClusterExperiment::Analysis::done() const
  {
    return analyzed_;
  }

  String ClusterExperiment::Analysis::name() const
  {
    return anafuncp_->getName();
  }
 
  void ClusterExperiment::Analysis::setAdapter(DBAdapter* adapterp)
  {
    anafuncp_->setDBAdapter(adapterp);
  }

//  void ClusterExperiment::Analysis::serialize(PersistenceManager& f)
//  {
//    f.writePointer("FactoryProduct",anafuncp_);
//    f.writeAttributeString("dataset",(int)dataset_);
//  }
  
  //writes formatted xml with the use of indent (see helper.h/C)
  void ClusterExperiment::Analysis::save(ostream& document, int& ind)
  {
    document << "<AnaFunc ";
		// TODO Persistence
    //anafuncp_->save(document,ind);
    document << indent(--ind);
    document << "</AnaFunc>\n";
    document << "<Results>\n";
    document << indent(++ind);
    for (map<String,double>::iterator mit = result_.begin(); mit != result_.end(); ++mit)
    {
      document << "<Result String = \"" << mit->first << "\" double = \"" << mit->second << "\"/>\n";
    }
    document << indent(--ind);
    document << "</Results>\n";
  }

	void ClusterExperiment::Analysis::persistentWrite(PersistenceManager& pm, const char* name) const throw (Exception::Base)
	{
		pm.writeObjectHeader(this,name);
		//TODO Persistence
		pm.writeObjectTrailer(name);
	}
	
	void ClusterExperiment::Analysis::persistentRead(PersistenceManager& pm) throw (Exception::Base)
	{
		//TODO Persistence
		int dummy;
		pm.readPrimitive(dummy,"dummy_");
	}
 
}
