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
// $Id: XMLHandler.C,v 1.8 2006/03/28 18:51:06 j-joachim Exp $
// $Author: j-joachim $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <qxml.h>

#include <iostream>
#include <vector>
#include <string>
#include <netinet/in.h> //network format

using namespace std;

namespace OpenMS
{
	namespace Internal
	{

	XMLHandler::XMLHandler()
	: error_message_("error triggered by consumer"),
		no_error_(true),
		use_warnings_(true)
	{
	}
	
	XMLHandler::~XMLHandler()
	{
	}

	bool XMLHandler::error(const QXmlParseException& exception)
	{
		cerr << "Error: " << exception.message();
		if (exception.lineNumber()!=-1)	cerr << " at line " << exception.lineNumber();
		cerr << endl;
		return false;
	}
	
	bool XMLHandler::fatalError(const QXmlParseException& exception)
	{
		QString error = exception.message();
		if (exception.lineNumber()!=-1)	error.append( QString(" at line %1").arg(exception.lineNumber()));
		throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__,file_.ascii(),error.ascii());
		return false;
	}
	
	bool XMLHandler::warning(const QXmlParseException& exception)
	{
		cout << "Warning: "<< exception.message();
		if (exception.lineNumber()!=-1)	cerr << " at line " << exception.lineNumber();
		cout << endl;
		return false;
	}
	
	bool XMLHandler::characters( const QString & /*chars*/ )
	{
		return true;
	}
	
	bool XMLHandler::startElement(const QString & /*uri*/, const QString & /*local_name*/, 
																const QString & /*qname*/, const QXmlAttributes & /*attributes*/ ) 
	{
		return true;
	}

	bool XMLHandler::endElement( const QString & /*uri*/, const QString & /*local_name*/, const QString & /*qname*/ ) 
	{
		return true;
	}

	QString XMLHandler::errorString()
	{
		return error_message_;
	}

	void XMLHandler::setUseWarnings(bool doUse)
	{
		use_warnings_ = doUse;
	}

	bool XMLHandler::useWarnings()
	{
		return use_warnings_;
	}

	} // namespace Internal
} // namespace OpenMS
