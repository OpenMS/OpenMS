// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/UnimodXMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/UnimodXMLHandler.h>
#include <OpenMS/SYSTEM/File.h>

#include <OpenMS/CHEMISTRY/ResidueModification.h>

#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>

using namespace xercesc;
using namespace std;

namespace OpenMS 
{

	UnimodXMLFile::UnimodXMLFile()
	{
	  	
	}
	
	UnimodXMLFile::~UnimodXMLFile()
	{
	}
	
  void UnimodXMLFile::load(const String& filename, vector<ResidueModification*>& modifications) const
  {
		String file = File::find(filename);
						
  	//try to open file
		if (!File::exists(file))
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }
		
		// initialize parser
		try 
		{
			XMLPlatformUtils::Initialize();
		}
		catch (const XMLException& toCatch) 
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("Error during initialization: ") + Internal::StringManager().convert(toCatch.getMessage()) );
	  }

		SAX2XMLReader* parser = XMLReaderFactory::createXMLReader();
		parser->setFeature(XMLUni::fgSAX2CoreNameSpaces,false);
		parser->setFeature(XMLUni::fgSAX2CoreNameSpacePrefixes,false);

		Internal::UnimodXMLHandler handler(modifications, file);
		
		parser->setContentHandler(&handler);
		parser->setErrorHandler(&handler);
		
		LocalFileInputSource source(Internal::StringManager().convert(file.c_str()));
		try 
    {
    	parser->parse(source);
    	delete(parser);
    }
    catch (const XMLException& toCatch) 
    {
      throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("XMLException: ") + Internal::StringManager().convert(toCatch.getMessage()) );
    }
    catch (const SAXException& toCatch) 
    {
      throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("SAXException: ") + Internal::StringManager().convert(toCatch.getMessage()) );
    }
  }  					 
  					 
} // namespace OpenMS
