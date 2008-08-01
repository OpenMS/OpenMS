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
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MascotXMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/MascotXMLHandler.h>
#include <OpenMS/SYSTEM/File.h>

#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>

using namespace xercesc;
using namespace std;

namespace OpenMS 
{

	MascotXMLFile::MascotXMLFile()
	{
	  	
	}

  void MascotXMLFile::load(const String& filename, 
						      					ProteinIdentification& protein_identification, 
						      					vector<PeptideIdentification>& id_data) const
  {
  	map<String, vector<AASequence> > peptides;
  	
  	load(filename, protein_identification, id_data, peptides);      
  }  					 
  					 
  void MascotXMLFile::load(const String& filename, 
						      					ProteinIdentification& protein_identification, 
						      					vector<PeptideIdentification>& id_data,
						      					map<String, vector<AASequence> >& peptides) const
  {
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

		protein_identification = ProteinIdentification();
		id_data.clear();

		SAX2XMLReader* parser = XMLReaderFactory::createXMLReader();
		parser->setFeature(XMLUni::fgSAX2CoreNameSpaces,false);
		parser->setFeature(XMLUni::fgSAX2CoreNameSpacePrefixes,false);

		Internal::MascotXMLHandler handler(protein_identification, id_data, filename, peptides);


		parser->setContentHandler(&handler);
		parser->setErrorHandler(&handler);
		
		LocalFileInputSource source( Internal::StringManager().convert(filename.c_str()) );
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
		
		// Since the mascot xml can contain "peptides" without sequences the identifications 
		// without any real peptide hit are removed
  	vector<PeptideHit> peptide_hits;
  	vector<PeptideIdentification>::iterator id_it = id_data.begin();

		while(id_it != id_data.end())
		{
			peptide_hits = id_it->getHits();
			if (peptide_hits.size() == 0 || (peptide_hits.size() == 1 && peptide_hits[0].getSequence() == ""))
			{
				id_it = id_data.erase(id_it);
			}
			else
			{
				++id_it;
			}
		}         
      
  }  					 
  					 
} // namespace OpenMS
