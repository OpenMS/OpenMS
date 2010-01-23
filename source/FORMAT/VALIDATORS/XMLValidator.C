// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/VALIDATORS/XMLValidator.h>
#include <OpenMS/SYSTEM/File.h>

#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
#include <xercesc/validators/common/Grammar.hpp>

using namespace xercesc;
using namespace std;

namespace OpenMS
{
	XMLValidator::XMLValidator()
		: valid_(true),
			os_(0)
	{
	}

	bool XMLValidator::isValid(const String& filename, const String& schema,  std::ostream& os)
	{
		filename_ = filename;
		os_ = &os;
		
		//try to open file
		if (!File::exists(filename))
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
		parser->setFeature(XMLUni::fgSAX2CoreNameSpaces, true);
		parser->setFeature(XMLUni::fgSAX2CoreValidation, true);
		parser->setFeature(XMLUni::fgXercesDynamic, false);
		parser->setFeature(XMLUni::fgXercesSchema, true);
		parser->setFeature(XMLUni::fgXercesSchemaFullChecking, true);
		
		//set this class as error handler
		parser->setErrorHandler(this);
		parser->setContentHandler(NULL);
		parser->setEntityResolver(NULL);
		
		//load schema
		LocalFileInputSource schema_file(Internal::StringManager().convert(schema));
		parser->loadGrammar(schema_file, Grammar::SchemaGrammarType, true);
		parser->setFeature(XMLUni::fgXercesUseCachedGrammarInParse, true);
		
		// try to parse file
		LocalFileInputSource source(Internal::StringManager().convert(filename.c_str()));
			
		try 
		{
			parser->parse(source);
			delete(parser);
		}
		catch (...) 
		{
			/// nothing to do here
		}
		
		return valid_;
	}
  
	void XMLValidator::warning(const SAXParseException& exception)
	{
		char* message = XMLString::transcode(exception.getMessage());
		String error_message = String("Validation warning in file '") + filename_ + "' line " + (UInt) exception.getLineNumber() + " column " + (UInt) exception.getColumnNumber() + ": " + message;
		(*os_) << error_message << endl;
		valid_ = false;
		XMLString::release(&message);
	}

	void XMLValidator::error(const SAXParseException& exception)
	{
		char* message = XMLString::transcode(exception.getMessage());
		String error_message = String("Validation error in file '") + filename_ + "' line " + (UInt) exception.getLineNumber() + " column " + (UInt) exception.getColumnNumber() + ": " + message;
		(*os_) << error_message << endl;
		valid_ = false;
		XMLString::release(&message);
	}

	void XMLValidator::fatalError(const SAXParseException& exception)
	{
		char* message = XMLString::transcode(exception.getMessage());
		String error_message = String("Validation error in file '") + filename_ + "' line " + (UInt) exception.getLineNumber() + " column " + (UInt) exception.getColumnNumber() + ": " + message;
		(*os_) << error_message << endl;
		valid_ = false;
		XMLString::release(&message);
	}

	void XMLValidator::resetErrors()
	{
		valid_ = true;
	}

} // namespace OpenMS

