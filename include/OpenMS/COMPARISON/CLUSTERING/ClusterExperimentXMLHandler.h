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

#ifndef OPENMS_COMPARISON_CLUSTERING_CLUSTEREXPERIMENTXMLHANDLER_H
#define OPENMS_COMPARISON_CLUSTERING_CLUSTEREXPERIMENTXMLHANDLER_H

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>

namespace OpenMS
{
	class FactoryProduct;
	class ClusterExperiment;
	class ClusterNode;

  class ClusterExperimentXMLHandler : public Internal::XMLHandler
  {
	  public:
	  	//Constructor
	    ClusterExperimentXMLHandler(ClusterExperiment&, const String& filename);
			
			//Destructor
	    virtual ~ClusterExperimentXMLHandler();

			// Docu in base class
      virtual void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);
			
			// Docu in base class
      virtual void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);
			
			// Docu in base class
      virtual void characters(const XMLCh* const chars, const unsigned int length);
 
  private:
    //not supported by xerces (and not useful)
    ClusterExperimentXMLHandler( const ClusterExperimentXMLHandler& source);

    //not possible because of the reference
    //and not necessary
    ClusterExperimentXMLHandler& operator=(const ClusterExperimentXMLHandler& source);

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
