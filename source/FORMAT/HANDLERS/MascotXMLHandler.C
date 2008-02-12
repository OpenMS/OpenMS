// -*- Mode: C++; tab-width: 2; -*-
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

#include <OpenMS/FORMAT/HANDLERS/MascotXMLHandler.h>

#include <xercesc/sax2/Attributes.hpp>

using namespace std;
using namespace xercesc;

namespace OpenMS
{
	namespace Internal
	{
  
  MascotXMLHandler::MascotXMLHandler(ProteinIdentification& protein_identification,
								  									 vector<PeptideIdentification>& id_data, 
      								 							 const String& filename) :
    XMLHandler(filename),
    protein_identification_(protein_identification),
    id_data_(id_data),
    actual_protein_hit_(),
    actual_peptide_hit_(),
    peptide_identification_index_(0),
		tag_(),
		date_()        
  {
  	
  }
   
  MascotXMLHandler::~MascotXMLHandler()
  {
    
  }
  
  void MascotXMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const Attributes& attributes)
	{

		tag_ = String(sm_.convert(qname));
		
		if (tag_ == "protein")
		{
			String attribute_value = String(sm_.convert(attributes.getValue(0u))).trim();
 	 		actual_protein_hit_.setAccession(attribute_value);
		}
		else if (tag_ == "query")
		{
			actual_query_ = (String(sm_.convert(attributes.getValue(0u))).trim()).toInt();
		}
		else if (tag_ == "peptide" || tag_ == "u_peptide") 
		{
			if (tag_ == "peptide")
			{
				String attribute_value = String(sm_.convert(attributes.getValue(0u))).trim();
		  	peptide_identification_index_ = attribute_value.toInt() - 1;
			}
			else if (tag_ == "u_peptide")
			{
				String attribute_value = String(sm_.convert(attributes.getValue(0u))).trim();
	  		peptide_identification_index_ = attribute_value.toInt() - 1;
			}
			if (peptide_identification_index_ > id_data_.size())
			{
				throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, ".mascotXML", 
																		"No header information present: use "
																		"the show_header=1 option in the "
																		"./export_dat.pl script");  			
  		}			
		}
	}
	  
  void MascotXMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
 	{
 		tag_ = String(sm_.convert(qname)).trim();
 		 
 		if (tag_ == "protein")
 		{	
 			// since Mascot uses SwissProt IDs we set this type here
			protein_identification_.setScoreType("Mascot");
 			protein_identification_.insertHit(actual_protein_hit_);
 			actual_protein_hit_ = ProteinHit();
 		}
 		else if (tag_ == "peptide")
 		{
			bool already_stored = false;
			vector<PeptideHit>::iterator it;
 			
			vector<PeptideHit> temp_peptide_hits = id_data_[peptide_identification_index_].getHits();
				
			it = temp_peptide_hits.begin();
			while(it != temp_peptide_hits.end() && !already_stored)
			{
				if (it->getSequence() == actual_peptide_hit_.getSequence())
				{
					already_stored = true;
				}
				it++;
			}
			if (!already_stored)
			{
				id_data_[peptide_identification_index_].setScoreType("Mascot");
				actual_peptide_hit_.addProteinAccession(actual_protein_hit_.getAccession());
	 			id_data_[peptide_identification_index_].insertHit(actual_peptide_hit_); 			
			}
			else
			{
				it--;
				it->addProteinAccession(actual_protein_hit_.getAccession());
				id_data_[peptide_identification_index_].setHits(temp_peptide_hits);
			}
 			actual_peptide_hit_ = PeptideHit();
 		}
 		else if (tag_ == "u_peptide")
 		{
			id_data_[peptide_identification_index_].setScoreType("Mascot");
 			id_data_[peptide_identification_index_].insertHit(actual_peptide_hit_); 			
 			actual_peptide_hit_ = PeptideHit();
 		}
 		else if (tag_ == "mascot_search_results")
 		{
 			protein_identification_.setSearchEngine("Mascot");
 			protein_identification_.setIdentifier(identifier_);
 			protein_identification_.setSearchParameters(search_parameters_);
 		}
		tag_ = "";
 	} 

  void MascotXMLHandler::characters(const XMLCh* const chars, unsigned int /*length*/)
  {

		if (tag_ == "NumQueries")
		{
			id_data_.resize(((String) sm_.convert(chars)).trim().toInt());
			for(vector<PeptideIdentification>::iterator it = id_data_.begin();
				  it != id_data_.end();
				  ++it)
			{
				it->setIdentifier(identifier_);
			}
			tag_ = "";
		}
		else if (tag_ == "prot_score")
		{
			actual_protein_hit_.setScore(((String) sm_.convert(chars)).trim().toInt());
		}
		else if (tag_ == "pep_exp_mz")
		{
			id_data_[peptide_identification_index_].setMetaValue("MZ", ((String) sm_.convert(chars)).trim().toFloat());
			tag_ = "";
		}
		else if (tag_ == "pep_exp_z")
		{
			actual_peptide_hit_.setCharge(((String) sm_.convert(chars)).trim().toInt());
			tag_ = "";
		}
		else if (tag_ == "pep_score")
		{
			actual_peptide_hit_.setScore(((String) sm_.convert(chars)).trim().toFloat());
			tag_ = "";
		}
		else if (tag_ == "pep_expect")
		{
			actual_peptide_hit_.metaRegistry().registerName("EValue", "E-value of e.g. Mascot searches", ""); // @todo what E-value flag? (andreas)
			actual_peptide_hit_.setMetaValue("EValue", ((String)sm_.convert(chars)).trim().toFloat());
			tag_ = "";
		}
		else if (tag_ == "pep_homol")
		{			
			id_data_[peptide_identification_index_].setSignificanceThreshold(
					((String) sm_.convert(chars)).trim().toFloat());
			tag_ = "";
		}
		else if (tag_ == "pep_ident")
		{
			DoubleReal temp_homology = 0;
			DoubleReal temp_identity = 0;
			
			// According to matrixscience the homology threshold is only used if it exists and is
			// smaller than the identity threshold.
			temp_homology = 
				id_data_[peptide_identification_index_].getSignificanceThreshold();
			temp_identity = ((String) sm_.convert(chars)).trim().toFloat();
			if (temp_homology > temp_identity || temp_homology == 0)
			{
				id_data_[peptide_identification_index_].setSignificanceThreshold(
					temp_identity);				
			}
			tag_ = "";
		}
		else if (tag_ == "pep_seq")
		{
			actual_peptide_hit_.setSequence(((String) sm_.convert(chars)).trim());
			tag_ = "";
		}
		else if (tag_ == "pep_res_before")
		{
			String temp_string = ((String) sm_.convert(chars)).trim();
			if (temp_string != "")
			{
				actual_peptide_hit_.setAABefore(temp_string[0]);
			}
			tag_ = "";
		}
		else if (tag_ == "pep_res_after")
		{
			String temp_string = ((String) sm_.convert(chars)).trim();
			if (temp_string != "")
			{
				actual_peptide_hit_.setAAAfter(temp_string[0]);
			}
			tag_ = "";
		}
		else if (tag_ == "Date")
		{	
			vector< String > parts;
			
			((String) sm_.convert(chars)).trim().split('T', parts);
			if (parts.size() == 2)
			{
				date_.set(parts[0] + ' ' + parts[1].prefix('Z'));
				date_time_string_ = parts[0] + ' ' + parts[1].prefix('Z');
				identifier_ = "Mascot_" + date_time_string_;
			}
			protein_identification_.setDateTime(date_);
		}
		else if (tag_ == "StringTitle")
		{
			String title = String(sm_.convert(chars)).trim();
			vector<String> parts;
			
			title.split('_', parts);
			if (parts.size() == 2)
			{
				id_data_[actual_query_ - 1].setMetaValue("RT", parts[1].toFloat());
			}
		}
		else if (tag_ == "MascotVer")
		{
			protein_identification_.setSearchEngineVersion(((String) sm_.convert(chars)).trim());
		}
		else if (tag_ == "DB")
		{
			search_parameters_.db = (((String) sm_.convert(chars)).trim());			
		}
		else if (tag_ == "FastaVer")
		{
			search_parameters_.db_version = (((String) sm_.convert(chars)).trim());			
		}
  }

	} // namespace Internal
} // namespace OpenMS
