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
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/MascotXMLHandler.h>

#include <xercesc/sax2/Attributes.hpp>

#include <iostream>

using namespace std;
using namespace xercesc;

namespace OpenMS
{
	namespace Internal
	{
  
  MascotXMLHandler::MascotXMLHandler(ProteinIdentification& protein_identification,
								  									 vector<IdentificationData>& id_data, 
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
  
  void MascotXMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
	{

		tag_ = String(XMLString::transcode(qname));
		
		if (tag_ == "protein")
		{
			String attribute_value = String(XMLString::transcode(attributes.getValue(0u))).trim();
 	 		actual_protein_hit_.setAccession(attribute_value);
		}
		else if (tag_ == "query")
		{
			actual_query_ = (String(XMLString::transcode(attributes.getValue(0u))).trim()).toInt();
		}
		else 
		{
			if (tag_ == "peptide")
			{
				String attribute_value = String(XMLString::transcode(attributes.getValue(0u))).trim();
		  	peptide_identification_index_ = attribute_value.toInt() - 1;
			}
			else if (tag_ == "u_peptide")
			{
				String attribute_value = String(XMLString::transcode(attributes.getValue(0u))).trim();
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
 		tag_ = String(XMLString::transcode(qname)).trim();
 		 
 		if (tag_ == "protein")
 		{	
 			/// since Mascot uses SwissProt IDs we set this type here
			actual_protein_hit_.setAccessionType("SwissProt");
			actual_protein_hit_.setScoreType("Mascot");
 			protein_identification_.insertProteinHit(actual_protein_hit_);
 			actual_protein_hit_.clear();
 		}
 		else if (tag_ == "peptide")
 		{
			bool													already_stored = false;
			vector<PeptideHit>::iterator  it;
 			
			vector<PeptideHit>& temp_peptide_hits = 
				id_data_[peptide_identification_index_].id.getPeptideHits();
				
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
				actual_peptide_hit_.setScoreType("Mascot");
				actual_peptide_hit_.addProteinIndex(make_pair(date_time_string_, actual_protein_hit_.getAccession()));
	 			id_data_[peptide_identification_index_].id.insertPeptideHit(actual_peptide_hit_); 			
			}
			else
			{
				it--;
				it->addProteinIndex(make_pair(date_time_string_, actual_protein_hit_.getAccession()));
			}
 			actual_peptide_hit_.clear();
 		}
 		else if (tag_ == "u_peptide")
 		{
			actual_peptide_hit_.setScoreType("Mascot");
 			id_data_[peptide_identification_index_].id.insertPeptideHit(actual_peptide_hit_); 			
 			actual_peptide_hit_.clear();
 		}
		tag_ = "";
 	} 

  void MascotXMLHandler::characters(const XMLCh* const chars, const unsigned int /*length*/)
  {

		if (tag_ == "NumQueries")
		{
			id_data_.resize(((String) XMLString::transcode(chars)).trim().toInt());
			for(UnsignedInt i = 0; i < id_data_.size(); i++)
			{
				id_data_[i].id.setDateTime(date_);
			}
			tag_ = "";
		}
		else if (tag_ == "prot_score")
		{
			actual_protein_hit_.setScore(((String) XMLString::transcode(chars)).trim().toInt());
		}
		else if (tag_ == "pep_exp_mz")
		{
			id_data_[peptide_identification_index_].mz = ((String) XMLString::transcode(chars)).trim().toFloat();
			tag_ = "";
		}
		else if (tag_ == "pep_exp_z")
		{
			actual_peptide_hit_.setCharge(((String) XMLString::transcode(chars)).trim().toInt());
			tag_ = "";
		}
		else if (tag_ == "pep_score")
		{
			actual_peptide_hit_.setScore(((String) XMLString::transcode(chars)).trim().toFloat());
			tag_ = "";
		}
		else if (tag_ == "pep_homol")
		{			
			id_data_[peptide_identification_index_].id.setPeptideSignificanceThreshold(
					((String) XMLString::transcode(chars)).trim().toFloat());
			tag_ = "";
		}
		else if (tag_ == "pep_ident")
		{
			DoubleReal temp_homology = 0;
			DoubleReal temp_identity = 0;
			
			/// According to matrixscience the homology threshold is only used if it exists and is
			/// smaller than the identity threshold.
			temp_homology = 
				id_data_[peptide_identification_index_].id.getPeptideSignificanceThreshold();
			temp_identity = ((String) XMLString::transcode(chars)).trim().toFloat();
			if (temp_homology > temp_identity || temp_homology == 0)
			{
				id_data_[peptide_identification_index_].id.setPeptideSignificanceThreshold(
					temp_identity);				
			}
			tag_ = "";
		}
		else if (tag_ == "pep_seq")
		{
			actual_peptide_hit_.setSequence(((String) XMLString::transcode(chars)).trim());
			tag_ = "";
		}
		else if (tag_ == "Date")
		{	
			vector< String > parts;
			
			((String) XMLString::transcode(chars)).trim().split('T', parts);
			if (parts.size() == 2)
			{
				date_.set(parts[0] + ' ' + parts[1].prefix('Z'));
				date_time_string_ = parts[0] + ' ' + parts[1].prefix('Z');
			}
			protein_identification_.setDateTime(date_);
		}
		else if (tag_ == "StringTitle")
		{
			String title = String(XMLString::transcode(chars)).trim();
			vector<String> parts;
			
			title.split('_', parts);
			if (parts.size() == 2)
			{
				id_data_[actual_query_ - 1].rt = parts[1].toFloat();
			}
		}
  }

	} // namespace Internal
} // namespace OpenMS
