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

#include <OpenMS/FORMAT/MascotXMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/MascotXMLHandler.h>

#include <qfileinfo.h>
#include <qfile.h>

#include <iostream>

using namespace std;

namespace OpenMS 
{

	MascotXMLFile::MascotXMLFile()
	{
	  	
	}
	
	MascotXMLFile::~MascotXMLFile()
	{
	  	
	}

  void MascotXMLFile::load(const String& filename,
  												 ProteinIdentification* protein_identification, 
  												 vector<Identification>* identifications, 
  												 vector<float>* precursor_retention_times,
  												 vector<float>* precursor_mz_values)
  	const throw (Exception::FileNotFound, 
  							 Exception::FileNotReadable, 
  							 Exception::FileEmpty,
  							 Exception::ParseError)
  {
  	vector<PeptideHit> 									peptide_hits;
  	vector<float>::iterator 						rt_it;
  	vector<float>::iterator							mz_it;
  	vector<Identification>::iterator 		id_it;
  	
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


		*protein_identification = ProteinIdentification();
		*identifications = vector<Identification>();  					 // clear information
		*precursor_retention_times = vector<float>();  // clear information
		*precursor_mz_values = vector<float>();  			 // clear information

		Internal::MascotXMLHandler handler(protein_identification, 
																			 identifications,
															 				 precursor_retention_times, 
															 				 precursor_mz_values);
		parser.setContentHandler(&handler);
		parser.setErrorHandler(&handler);

		QXmlInputSource source( file );
		try
		{
			parser.parse(source);
  	}
		catch(exception& e)
		{
			cout << e.what();
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename, "Parse error");
		}
		
		id_it = identifications->begin();
		rt_it = precursor_retention_times->begin();
		mz_it = precursor_mz_values->begin();
		
		/// Since the mascot xml can contain "peptides" without sequences the identifications 
		/// without any real peptide hit are removed as well as the corresponding mz and rt values.
		while(id_it != identifications->end() 
					&& rt_it != precursor_retention_times->end()
					&& mz_it != precursor_mz_values->end())
		{
			peptide_hits = id_it->getPeptideHits();
			if (peptide_hits.size() == 0
					|| (peptide_hits.size() == 1
						&& peptide_hits[0].getSequence() == ""))
			{
				id_it = identifications->erase(id_it);
				rt_it = precursor_retention_times->erase(rt_it);
				mz_it = precursor_mz_values->erase(mz_it);
			}
			else
			{
				id_it++;
				rt_it++;
				mz_it++;
			}
		}  
  }  					 
  					 
} // namespace OpenMS
