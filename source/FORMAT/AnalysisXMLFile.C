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
// $Id: AnalysisXMLFile.C,v 1.10 2006/06/09 23:47:35 nicopfeifer Exp $
// $Author: nicopfeifer $
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/AnalysisXMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/AnalysisXMLHandler.h>

#include <qfileinfo.h>
#include <qfile.h>

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
  	const throw (Exception::FileNotFound, 
  							 Exception::FileNotReadable, 
  							 Exception::FileEmpty,
  							 Exception::ParseError)
  {
  	//try to open file
		QFileInfo file_info;
		file_info.setFile(filename.c_str());
    if (!file_info.exists())
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }

    if (!file_info.isReadable())
    {
      throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }
    if (file_info.size() == 0)
    {
      throw Exception::FileEmpty(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }
		QFile file(filename.c_str());
		QXmlSimpleReader parser;
		srand(static_cast<unsigned>(time(0)));
		parser.setFeature("http://xml.org/sax/features/namespaces",false);
		parser.setFeature("http://xml.org/sax/features/namespace-prefixes", false);


		*protein_identifications = vector<ProteinIdentification>(); // clear information
		*identifications = vector<Identification>();  					 		// clear information
		*precursor_retention_times = vector<float>();  							// clear information
		*precursor_mz_values = vector<float>();  			 							// clear information
		*contact_person = ContactPerson();  					 							// clear information

		Internal::AnalysisXMLHandler handler(protein_identifications,
															 identifications,
															 precursor_retention_times, 
															 precursor_mz_values,
															 contact_person);
		parser.setContentHandler(&handler);
		parser.setErrorHandler(&handler);

		QXmlInputSource source( file );
		try
		{
			parser.parse(source);
  	}
		catch(exception e)
		{
			cout << e.what() << endl;
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename, "Parse error");
		}  
		catch(...)
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename, "Parse error");
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
  	const throw (Exception::FileNotFound, 
  							 Exception::FileNotReadable, 
  							 Exception::FileEmpty,
  							 Exception::ParseError)
  {
  	//try to open file
		QFileInfo file_info;
		file_info.setFile(filename.c_str());
    if (!file_info.exists())
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }

    if (!file_info.isReadable())
    {
      throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }
    if (file_info.size() == 0)
    {
      throw Exception::FileEmpty(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }
		QFile file(filename.c_str());
		QXmlSimpleReader parser;
		srand(static_cast<unsigned>(time(0)));
		parser.setFeature("http://xml.org/sax/features/namespaces",false);
		parser.setFeature("http://xml.org/sax/features/namespace-prefixes", false);

		/// clear information
		*protein_identifications = vector<ProteinIdentification>();
		*identifications = vector<Identification>();  					 
		*precursor_retention_times = vector<float>();  
		*precursor_mz_values = vector<float>();  			 
		*contact_person = ContactPerson();  					
		*predicted_retention_times = std::map<String, double>();
		*predicted_sigma = 0;

		Internal::AnalysisXMLHandler handler(protein_identifications,
															 identifications,
															 precursor_retention_times, 
															 precursor_mz_values,
															 contact_person,
															 predicted_retention_times,
															 predicted_sigma);
		parser.setContentHandler(&handler);
		parser.setErrorHandler(&handler);

		QXmlInputSource source( file );
		try
		{
			parser.parse(source);
  	}
		catch(exception e)
		{
			cout << e.what() << endl;
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename, "Parse error");
		}  
		catch(...)
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename, "Parse error");
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
