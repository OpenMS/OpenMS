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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/ParamXMLHandler.h>

#include <iostream>

using namespace std;

namespace OpenMS
{
	namespace Internal
	{

	ParamXMLHandler::ParamXMLHandler(map<string,DataValue>& values): 	values_(values)
	{
		file_ = __FILE__;
	}
	
	ParamXMLHandler::~ParamXMLHandler()
	{
	}
	
	bool ParamXMLHandler::startElement(const QString & /*uri*/, const QString & /*local_name*/, 
																const QString & qname, const QXmlAttributes & attributes )
	{
		if ("ITEM" == qname)
		{
			QString name, type, value;
			for(int n = 0 ; n < attributes.length() ; n++ )
			{
			  QString attributesValue = attributes.value(n);
			  QString attributesName  = attributes.qName(n);
		
			  if(attributesName == "name")
			  {
			    name = attributesValue;
				}
			  if(attributesName == "value")
			  {
			    value = attributesValue;
				}
			  if(attributesName == "type")
			  {
			    type = attributesValue;
				}
			}
			if (type == "int")
			{
				values_[path_+name.ascii()]=DataValue(asSignedInt_(value.ascii()));
			}
			if (type == "string")
			{
				values_[path_+name.ascii()]=DataValue(value.ascii());
			}
			if (type == "float")
			{
				values_[path_+name.ascii()]=DataValue(asFloat_(value));
			}
		}
		
		if ("NODE" == qname)
		{
			QString tmp;
			for(int n = 0 ; n < attributes.length() ; n++ )
			{
			  QString attributesValue = attributes.value(n);
			  QString attributesName = attributes.qName(n);
		
			  if( attributesName == "name" )
			  {
			    nodes_.push_back(attributesValue.ascii());
					tmp += attributesValue + ":";
				}
			}
			path_ += tmp.ascii();
		}

		return true;
	}

	bool ParamXMLHandler::endElement( const QString & /*uri*/, const QString & /*local_name*/, const QString & qname )
	{
		if ("NODE" == qname)
		{
			nodes_.pop_back();
			//renew path
			path_ = "";
			for (vector<string>::iterator it = nodes_.begin(); it != nodes_.end();++it)
			{
				path_ += *it+":";
			}
			
		}		
		return true;
	}


	} // namespace Internal
} // namespace OpenMS
