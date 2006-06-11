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
// $Id: ParamXMLHandler.h,v 1.7 2006/05/30 16:41:13 marc_sturm Exp $
// $Author: marc_sturm $
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
      ParamXMLHandler(std::map<std::string, DataValue>& values);

      virtual ~ParamXMLHandler();

      virtual bool startElement(const QString & namespaceURI, const QString & localName, 
																const QString & qName, const QXmlAttributes & atts );

      virtual bool endElement( const QString & namespaceURI, const QString & localName,
															 const QString & qName );

    protected :
      std::vector<std::string> nodes_;
      std::string path_;
      std::map<std::string, DataValue>& values_;
  };

	} // namespace Internal
} // namespace OpenMS

#endif // OPENMS_FORMAT_HANDLERS_PARAMXMLHANDLER_H
