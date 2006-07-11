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
#ifndef OPENMS_COMPARISON_CLUSTERING_CLUSTEREXPERIMENTXMLHANDLER_H
#define OPENMS_COMPARISON_CLUSTERING_CLUSTEREXPERIMENTXMLHANDLER_H

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/CONCEPT/FactoryProduct.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterExperiment.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterNode.h>

namespace OpenMS
{

  class ClusterExperimentXMLHandler : public Internal::XMLHandler
  {
  public:
    ClusterExperimentXMLHandler(ClusterExperiment&);

    virtual ~ClusterExperimentXMLHandler();

		bool characters( const QString & chars );
		
		bool startElement(const QString & uri, const QString & local_name, 
											const QString & qname, const QXmlAttributes & attributes );
		
		bool endElement( const QString & uri, const QString & local_name,
										 const QString & qname ); 
 
  private:
    //not supported by xerces (and not useful)
    ClusterExperimentXMLHandler( const ClusterExperimentXMLHandler& source);

    ClusterExperimentXMLHandler& operator=(const ClusterExperimentXMLHandler& source);
    //not possible because of the reference
    //and not necessary
    ClusterExperimentXMLHandler();
    ClusterExperiment& crun_;
    FactoryProduct* forwardconfigurablep_;
    bool isclusterexperiment_;
    ClusterNode* tmpclusterp_;
    bool ismemberid_;
    bool iscomment_;
    bool isreference_;
    bool didrun_;
    int clid_;
    double cl_min_mass_;
    double cl_max_mass_;
  };

}
#endif //OPENMS_COMPARISON_CLUSTERING_CLUSTEREXPERIMENTXMLHANDLER_H
