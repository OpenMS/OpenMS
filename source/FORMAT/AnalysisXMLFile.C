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
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/AnalysisXMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/AnalysisXMLHandler.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>

#include <iostream>

using namespace std;

namespace OpenMS 
{

	AnalysisXMLFile::AnalysisXMLFile()
	{
	  	
	}
	
	AnalysisXMLFile::~AnalysisXMLFile()
	{
	  	
	}

  void AnalysisXMLFile::load(const String& filename, 
  					 								 vector<ProteinIdentification>* protein_identifications,
  													 vector<Identification>* identifications, 
  													 vector<float>* precursor_retention_times,
  													 vector<float>* precursor_mz_values,
  													 ContactPerson* contact_person)
  	const throw (Exception::FileNotFound, Exception::ParseError)
  {
  	//try to open file
		if (!TextFile::exists(filename))
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

		protein_identifications->clear();
		identifications->clear();
		precursor_retention_times->clear();
		precursor_mz_values->clear();
		*contact_person = ContactPerson();

		Internal::AnalysisXMLHandler handler(protein_identifications,
															 identifications,
															 precursor_retention_times, 
															 precursor_mz_values,
															 contact_person);

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
    catch (const xercesc::SAXParseException& toCatch) 
    {
      throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("SAXParseException: ") + xercesc::XMLString::transcode(toCatch.getMessage()) );
    }
    catch (...) 
    {
      throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("Unexpexted parse exception!"));
    }
  }
  					 
  void AnalysisXMLFile::load(const String& filename, 
  					 								 vector<ProteinIdentification>* protein_identifications,
  													 vector<Identification>* identifications, 
  													 vector<float>* precursor_retention_times,
  													 vector<float>* precursor_mz_values,
  													 ContactPerson* contact_person,
      											 std::map<String, double>* predicted_retention_times,
      											 DoubleReal* predicted_sigma)
  	const throw (Exception::FileNotFound, Exception::ParseError)
  {
  	//try to open file
		if (!TextFile::exists(filename))
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

		/// clear information
		protein_identifications->clear();
		identifications->clear();
		precursor_retention_times->clear();
		precursor_mz_values->clear();
		*contact_person = ContactPerson();
		predicted_retention_times->clear();
		*predicted_sigma = 0;

		Internal::AnalysisXMLHandler handler(protein_identifications,
															 identifications,
															 precursor_retention_times, 
															 precursor_mz_values,
															 contact_person,
															 predicted_retention_times,
															 predicted_sigma);

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
    catch (const xercesc::SAXParseException& toCatch) 
    {
      throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("SAXParseException: ") + xercesc::XMLString::transcode(toCatch.getMessage()) );
    }
    catch (...) 
    {
      throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("Unexpexted parse exception!"));
    }
  }
  					 
  void AnalysisXMLFile::store(String filename, 
  					 									const vector<ProteinIdentification>& protein_identifications,
  					 									const vector<Identification>& identifications, 
  					 									const vector<float>& precursor_retention_times,
  					 									const vector<float>& precursor_mz_values) const throw (Exception::UnableToCreateFile)
  {
		std::ofstream os(filename.c_str());
		if (!os)
		{
			throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
		}

		//read data and close stream
		Internal::AnalysisXMLHandler handler(protein_identifications,
															 					 identifications, 
															 					 precursor_retention_times, 
															 					 precursor_mz_values);
		handler.writeTo(os);
		os.close();

  } 

  void AnalysisXMLFile::store(String filename, 
  					 									const vector<ProteinIdentification>& protein_identifications,
  					 									const vector<Identification>& identifications, 
  					 									const vector<float>& precursor_retention_times,
  					 									const vector<float>& precursor_mz_values,
  					 									const ContactPerson& contact_person) const throw (Exception::UnableToCreateFile)
  {
		std::ofstream os(filename.c_str());
		if (!os)
		{
			throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
		}

		//read data and close stream
		Internal::AnalysisXMLHandler handler(protein_identifications,
															 identifications,
															 precursor_retention_times, 
															 precursor_mz_values, 
															 contact_person);
		handler.writeTo(os);
		os.close();

  } 

  void AnalysisXMLFile::store(String filename, 
  					 									const vector<ProteinIdentification>& protein_identifications,
  					 									const vector<Identification>& identifications, 
  					 									const vector<float>& precursor_retention_times,
  					 									const vector<float>& precursor_mz_values,
  					 									const ContactPerson& contact_person,
  					 									const map<String, double>& predicted_retention_times,
  					 									DoubleReal predicted_sigma) 
  	const throw (Exception::UnableToCreateFile)
  {
		std::ofstream os(filename.c_str());
		if (!os)
		{
			throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
		}

		//read data and close stream
		Internal::AnalysisXMLHandler handler(protein_identifications,
															 identifications, 
															 precursor_retention_times, 
															 precursor_mz_values, 
															 contact_person,
															 predicted_retention_times,
															 predicted_sigma);
		handler.writeTo(os);
		os.close();

  } 


} // namespace OpenMS
