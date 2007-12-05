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

#include <OpenMS/FORMAT/XTandemXMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/XTandemXMLHandler.h>
#include <OpenMS/SYSTEM/File.h>

#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>

using namespace xercesc;
using namespace std;

namespace OpenMS 
{

	XTandemXMLFile::XTandemXMLFile()
	{
	  	
	}
	
	XTandemXMLFile::~XTandemXMLFile()
	{
	}
	
  void XTandemXMLFile::load(const String& filename, ProteinIdentification& protein_identification, vector<PeptideIdentification>& peptide_ids) const throw (Exception::FileNotFound, Exception::ParseError)
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

		SAX2XMLReader* parser = XMLReaderFactory::createXMLReader();
		parser->setFeature(XMLUni::fgSAX2CoreNameSpaces,false);
		parser->setFeature(XMLUni::fgSAX2CoreNameSpacePrefixes,false);

		map<UInt, vector<PeptideHit> > peptide_hits;
		map<UInt, String> descriptions;
		Internal::XTandemXMLHandler handler(protein_identification, peptide_hits, descriptions, filename);

		parser->setContentHandler(&handler);
		parser->setErrorHandler(&handler);
		
		LocalFileInputSource source(Internal::StringManager().convert(filename.c_str()));
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

		DateTime now;
		now.now();
		String date_string;
		String identifier("XTandem_" + date_string);
		//vector<String> accessions;

		// convert id -> peptide_hits into peptide hits list
		//vector<PeptideIdentification> peptide_identifications;
		PeptideIdentification().metaRegistry().registerName("spectrum_id", "the id of the spectrum counting from 1");
		for (map<UInt, vector<PeptideHit> >::const_iterator it = peptide_hits.begin(); it != peptide_hits.end(); ++it)
		{
			// reduce the hits with the same sequence to one PeptideHit
			map<String, vector<PeptideHit> > seq_to_hits;
			for (vector<PeptideHit>::const_iterator it1 = it->second.begin(); it1 != it->second.end(); ++it1)
			{
				seq_to_hits[it1->getSequence()].push_back(*it1);
			}
			
			PeptideIdentification id;
			if (descriptions.find(it->first) != descriptions.end())
			{
				id.setMetaValue("Description", descriptions[it->first]);
			}
			for (map<String, vector<PeptideHit> >::const_iterator it1 = seq_to_hits.begin(); it1 != seq_to_hits.end(); ++it1)
			{
				if (it1->second.size() > 0)
				{
					// copy the accession of all to the first hit
					PeptideHit hit = *it1->second.begin();
					vector<String> accessions;
					for (vector<PeptideHit>::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
					{
						for (vector<String>::const_iterator it3 = it2->getProteinAccessions().begin(); it3 != it2->getProteinAccessions().end(); ++it3)
						{
							accessions.push_back(*it3);
						}
					}
					hit.setProteinAccessions(accessions);
					id.insertHit(hit);
				}
			}

			id.setScoreType("XTandem");
			id.setHigherScoreBetter(true);
			id.setIdentifier(identifier);
			id.assignRanks();
			id.setMetaValue("spectrum_id", (Int)it->first);
			
			peptide_ids.push_back(id);
		}

    //sort(accessions.begin(), accessions.end());
    //vector<String>::const_iterator end_unique = unique(accessions.begin(), accessions.end());

    // E-values
    protein_identification.setHigherScoreBetter(false);
		protein_identification.sort();
    protein_identification.setScoreType("XTandem");

    // TODO version of XTandem ???? is not available from performance param section of outputfile (to be parsed)
    // TODO Date of search, dito
    protein_identification.setDateTime(now);
    protein_identification.setIdentifier(identifier);

		

    // TODO search parameters are also available



  }  					 
  					 
} // namespace OpenMS
