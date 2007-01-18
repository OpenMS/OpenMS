// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_PARAMXMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_PARAMXMLHANDLER_H

#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>

#include <vector>
#include <map>

namespace OpenMS
{
	namespace Internal
	{
	/**
		@brief XML Handler for Param files.

	*/
  class ParamXMLHandler
  	: public XMLHandler
  {
    public:
    	/// Default constructor
      ParamXMLHandler(std::map<std::string, DataValue>& values, const String& filename);
			/// Destructor
      virtual ~ParamXMLHandler();

			// Docu in base class
      virtual void endElement( const XMLCh* const uri, const XMLCh* const local_name, const XMLCh* const qname);
			
			// Docu in base class
      virtual void startElement(const XMLCh* const uri, const XMLCh* const local_name, const XMLCh* const qname, const xercesc::Attributes& attributes);

    protected:
    	/// Not implemented => protected
    	ParamXMLHandler();
    	
      std::vector<std::string> nodes_;
      std::string path_;
      std::map<std::string, DataValue>& values_;
  };

	} // namespace Internal
} // namespace OpenMS

#endif // OPENMS_FORMAT_HANDLERS_PARAMXMLHANDLER_H
