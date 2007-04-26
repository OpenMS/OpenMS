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

#include <OpenMS/FORMAT/Schema2CV.h>
#include <OpenMS/SYSTEM/File.h>

#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>

#include <fstream>
#include <iostream>
#include <utility>
#include <algorithm>

using namespace std;

namespace OpenMS 
{
	
	Schema2CV::Schema2CV()
	{
		
	}

	Schema2CV::~Schema2CV()
	{
		
	}
	
	void Schema2CV::loadFromFile(const String& filename) throw (Exception::FileNotFound, Exception::ParseError)
	{
	    //try to open file
	    if (!File::exists(filename))
	    {
	      throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
	    }
	
	    // initialize parser
	    try
	    {
	      xercesc::XMLPlatformUtils::Initialize();
	    }
	    catch (const xercesc::XMLException& toCatch)
	    {
	      throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("Error during initialization: ") + xercesc::XMLString::transcode(toCatch.getMessage()) );
	    }
	
	    xercesc::SAX2XMLReader* parser = xercesc::XMLReaderFactory::createXMLReader();
	    parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpaces,false);
	    parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpacePrefixes,false);
	    Internal::Schema2CVHandler handler(filename, *this);
	    parser->setContentHandler(&handler);
	    parser->setErrorHandler(&handler);
	
	    xercesc::LocalFileInputSource source( xercesc::XMLString::transcode(filename.c_str()) );
	    try
	    {
	      parser->parse(source);
	      delete(parser);
	    }
	    catch (const xercesc::XMLException& toCatch)
	    {
	      throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("XMLException: ") + xercesc::XMLString::transcode(toCatch.getMessage()) );
	    }
	    catch (const xercesc::SAXException& toCatch)
	    {
	      throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("SAXException: ") + xercesc::XMLString::transcode(toCatch.getMessage()) );
	    }
	}
		
	const std::vector<Schema2CV::CVDesc>& Schema2CV::getCVs() const
	{
		return cvs_;	
	}
	
	const std::vector<Schema2CV::LocDesc>& Schema2CV::getLocations() const
	{
		return locs_;	
	}

	std::ostream& operator << (std::ostream& os, const Schema2CV& mapping)
	{
		for( std::vector<Schema2CV::CVDesc>::const_iterator it =mapping.cvs_.begin(); it != mapping.cvs_.end(); ++it)
		{
			
		}
		
		return os;
	}
	
	namespace Internal
	{
		Schema2CVHandler::Schema2CVHandler(const String& filename, Schema2CV& mapping)
			:	XMLHandler(filename),
				mapping_(mapping)
		{
			
		}
		
		void Schema2CVHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*localname*/, const XMLCh* const qname, const xercesc::Attributes& attrs)
		{
			char* tag = xercesc::XMLString::transcode(qname);
			
			//<CVSource name="GelCV" version="0.1" uri="psidev.sf.net/gel/cv" cvIdentifier="gel" cvFormat="OBO"/>
			if (xercesc::XMLString::equals(tag,"CVSource"))
			{
				Schema2CV::CVDesc tmp;
				tmp.name =  xercesc::XMLString::transcode(attrs.getValue(xercesc::XMLString::transcode("name")));
				tmp.version =  xercesc::XMLString::transcode(attrs.getValue(xercesc::XMLString::transcode("version")));
				tmp.uri =  xercesc::XMLString::transcode(attrs.getValue(xercesc::XMLString::transcode("uri")));
				tmp.id =  xercesc::XMLString::transcode(attrs.getValue(xercesc::XMLString::transcode("cvIdentifier")));
				tmp.format =  xercesc::XMLString::transcode(attrs.getValue(xercesc::XMLString::transcode("cvFormat")));
				mapping_.cvs_.push_back(tmp);
			}
			//<ModelElementMap elementPath="//GelRange/RangeType/OntologyTerm" requirementLevel="MAY">
			else if (xercesc::XMLString::equals(tag,"ModelElementMap"))
			{
				Schema2CV::LocDesc tmp;
				tmp.location =  xercesc::XMLString::transcode(attrs.getValue(xercesc::XMLString::transcode("elementPath")));
				if (xercesc::XMLString::equals(xercesc::XMLString::transcode(attrs.getValue(xercesc::XMLString::transcode("requirementLevel"))),"MUST"))
				{
					tmp.strict = true;
				}
				else
				{
					tmp.strict = false;
				}
				mapping_.locs_.push_back(tmp);
			}
			//<CVTerm termName="GradientRange" termAccession="SEP:to_come" useTerm="false" allowChildren="true" cvRef="gel"/>
			else if (xercesc::XMLString::equals(tag,"CVTerm"))
			{
				Schema2CV::TermDesc tmp;
				tmp.accession = xercesc::XMLString::transcode(attrs.getValue(xercesc::XMLString::transcode("termAccession")));
				tmp.cv = xercesc::XMLString::transcode(attrs.getValue(xercesc::XMLString::transcode("cvRef")));	
				if (xercesc::XMLString::equals(xercesc::XMLString::transcode(attrs.getValue(xercesc::XMLString::transcode("useTerm"))),"true"))
				{
					tmp.allowSelf = true;
				}
				else
				{
					tmp.allowSelf = false;
				}
				if (xercesc::XMLString::equals(xercesc::XMLString::transcode(attrs.getValue(xercesc::XMLString::transcode("allowChildren"))),"true"))
				{
					tmp.allowChildren = true;
				}
				else
				{
					tmp.allowChildren = false;
				}
				tmp.repeatable = true;
				mapping_.locs_.back().terms.push_back(tmp);
			}
		}
	}
	
} // namespace OpenMS

