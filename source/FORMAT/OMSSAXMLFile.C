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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/OMSSAXMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/OMSSAXMLHandler.h>
#include <OpenMS/SYSTEM/File.h>

#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>

using namespace xercesc;
using namespace std;

namespace OpenMS 
{

	OMSSAXMLFile::OMSSAXMLFile()
	{
	  	
	}
	
	OMSSAXMLFile::~OMSSAXMLFile()
	{
	}
	
  void OMSSAXMLFile::load(const String& filename, ProteinIdentification& protein_identification, vector<PeptideIdentification>& peptide_identifications) const throw (Exception::FileNotFound, Exception::ParseError)
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
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("Error during initialization: ") + Internal::StringManager().convert(toCatch.getMessage()) );
	  }

		xercesc::SAX2XMLReader* parser = xercesc::XMLReaderFactory::createXMLReader();
		parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpaces,false);
		parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpacePrefixes,false);

		Internal::OMSSAXMLHandler handler(protein_identification, peptide_identifications, filename);


		parser->setContentHandler(&handler);
		parser->setErrorHandler(&handler);
		
		xercesc::LocalFileInputSource source( Internal::StringManager().convert(filename.c_str()) );
		try 
    {
    	parser->parse(source);
    	delete(parser);
    }
    catch (const xercesc::XMLException& toCatch) 
    {
      throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("XMLException: ") + Internal::StringManager().convert(toCatch.getMessage()) );
    }
    catch (const xercesc::SAXException& toCatch) 
    {
      throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("SAXException: ") + Internal::StringManager().convert(toCatch.getMessage()) );
    }

		DateTime now;
		now.now();
		String date_string;
		now.get(date_string);
		String identifier("OMSSA_" + date_string);
	
		// post-processing
		vector<String> accessions;
		for (vector<PeptideIdentification>::iterator it = peptide_identifications.begin(); it != peptide_identifications.end(); ++it)
		{
			it->setScoreType("OMSSA");
			it->setHigherScoreBetter(false);
			it->setIdentifier(identifier);
			it->assignRanks();
			for (vector<PeptideHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
			{
				accessions.insert(accessions.end(), pit->getProteinAccessions().begin(), pit->getProteinAccessions().end());
			}
		}

		sort(accessions.begin(), accessions.end());
		vector<String>::const_iterator end_unique = unique(accessions.begin(), accessions.end());

		for (vector<String>::const_iterator it = accessions.begin(); it != end_unique; ++it)
		{
			ProteinHit hit;
			hit.setAccession(*it);
			protein_identification.insertHit(hit);
		}

		// E-values
		protein_identification.setHigherScoreBetter(false);
		protein_identification.setScoreType("OMSSA");
		
		// version of OMSSA is not available
		// Date of the search is not available -> set it to now
		protein_identification.setDateTime(now);
		protein_identification.setIdentifier(identifier);

		// search parameters are also not available
  }  					 
  					 
} // namespace OpenMS
