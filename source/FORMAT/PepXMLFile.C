// -*- Mode: C++; tab-widt: 2; -*-
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

#include <OpenMS/FORMAT/PepXMLFile.h>

#include <iostream>
#include <fstream>

using namespace std;

namespace OpenMS 
{

	PepXMLFile::PepXMLFile()
		: XMLHandler("","1.8"),
			XMLFile("/SCHEMAS/PepXML_1_8.xsd","1.8"),
			actual_title_(""),
			actual_sequence_(""),
			actual_modifications_(vector< pair<String, UInt> >()),
			peptides_(0),
			actual_aa_sequences_(vector<AASequence>()),
			modifications_(vector< pair<String, DoubleReal> >())
	{
	  	
	}

  void PepXMLFile::load(const String& filename,  map<String, vector<AASequence> >& peptides)
  {
  	//Filename for error messages in XMLHandler
  	file_ = filename;
  	
  	peptides.clear();
  	
  	peptides_ = &peptides;

		parse_(filename,this);
    
    //reset members
		actual_title_ = "";
		actual_sequence_ = "";
		actual_modifications_ = vector< pair<String, UInt> >();
		peptides_ = 0;
		modifications_ = vector< pair<String, DoubleReal> >();
  }
  					   
  void PepXMLFile::matchModification_(DoubleReal			 														mass,
  																		String&																			modification_description)
	{
		UInt i = 0;
		bool found = false;
		DoubleReal difference = 0.; 
		
		while(i < modifications_.size() && !found)
		{
			difference = modifications_[i].second - mass;
			if (difference < 0)
			{
				difference *= -1;
			}
			if (difference < 0.001)
			{
				modification_description = modifications_[i].first;
				found = true;
			}
			++i;			
		}			
	}  																		
  					   
	void PepXMLFile::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
	{		
		String element = sm_.convert(qname);
		
		//cout << "Start: " << element << endl;
		
		//SEARCH PARAMETERS
		if (element == "aminoacid_modification")
		{
			String temp_string = attributeAsString_(attributes,"variable");
			if (temp_string == "Y")
			{
				modifications_.push_back(make_pair(attributeAsString_(attributes,"description"), 
																					 attributeAsDouble_(attributes,"mass")));
			}
		}	
		//PEPTIDES
		else if (element == "spectrum_query")
		{				
			actual_title_ = attributeAsString_(attributes, "spectrum");
		}
		else if (element == "search_hit")
		{
			actual_sequence_ = attributeAsString_(attributes,"peptide");
		}
		else if (element == "mod_aminoacid_mass")
		{
			DoubleReal modification_mass = 0.;
			UInt 			 modification_position = 0;
			String 		 temp_description = "";
			
			modification_position = attributeAsInt_(attributes,"position");
			modification_mass = attributeAsDouble_(attributes,"mass");
			
			matchModification_(modification_mass, temp_description);
			
			// the modification position is 1-based
			actual_modifications_.push_back(make_pair(temp_description, modification_position));
		}
	}
	
	void PepXMLFile::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
	{
		String element = sm_.convert(qname);
				
		///SEARCH PARAMETERS
		if (element == "search_hit")
		{
			AASequence temp_aa_sequence = AASequence(actual_sequence_);
			
			// modification position is 1-based
			//actual_modifications_
			actual_aa_sequences_.push_back(temp_aa_sequence);
			
			actual_modifications_.clear();						
		}
		else if (element == "spectrum_query")
		{
			peptides_->insert(make_pair(actual_title_, actual_aa_sequences_));
			actual_aa_sequences_.clear();			
		}
	}

} // namespace OpenMS
