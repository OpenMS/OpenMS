// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
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
      								 							 const String& filename,
      								 							 map<String, vector<AASequence> >& modified_peptides) :
    XMLHandler(filename,""),
    protein_identification_(protein_identification),
    id_data_(id_data),
    actual_protein_hit_(),
    actual_peptide_hit_(),
    peptide_identification_index_(0),
		tag_(),
		date_(),
		actual_title_(""),
		modified_peptides_(modified_peptides)        
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
			String attribute_value = String(sm_.convert(attributes.getValue(XMLSize_t(0)))).trim();
 	 		actual_protein_hit_.setAccession(attribute_value);
		}
		else if (tag_ == "query")
		{
			actual_query_ = (String(sm_.convert(attributes.getValue(XMLSize_t(0)))).trim()).toInt();
		}
		else if (tag_ == "peptide" || tag_ == "u_peptide" || tag_ == "q_peptide") 
		{
			if (tag_ == "peptide")
			{
				String attribute_value = String(sm_.convert(attributes.getValue(XMLSize_t(0)))).trim();
		  	peptide_identification_index_ = attribute_value.toInt() - 1;
			}
			else if (tag_ == "u_peptide" || tag_ == "q_peptide")
			{
				String attribute_value = String(sm_.convert(attributes.getValue(XMLSize_t(0)))).trim();
	  		peptide_identification_index_ = attribute_value.toInt() - 1;
			}
			if (peptide_identification_index_ > id_data_.size())
			{
				fatalError(LOAD, "No header information present: use  the show_header=1 option in the ./export_dat.pl script");  			
  		}			
		}
	}
	  
  void MascotXMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
 	{
 		tag_ = String(sm_.convert(qname)).trim();
 		 
 		if (tag_ == "protein")
 		{	
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
 		else if (tag_ == "u_peptide" || tag_ == "q_peptide")
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

  void MascotXMLHandler::characters(const XMLCh* const chars, const XMLSize_t /*length*/)
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
			id_data_[peptide_identification_index_].setMetaValue("MZ", ((String) sm_.convert(chars)).trim().toDouble());
			tag_ = "";
		}
		else if (tag_ == "pep_scan_title") 
		{
			String title = ((String) sm_.convert(chars)).trim();
			if(title.hasSubstring("_")) 
			{
				id_data_[peptide_identification_index_].setMetaValue("RT", (title.suffix('_').toDouble()));
			}
		}
		else if (tag_ == "pep_exp_z")
		{
			actual_peptide_hit_.setCharge(((String) sm_.convert(chars)).trim().toInt());
			tag_ = "";
		}
		else if (tag_ == "pep_score")
		{
			actual_peptide_hit_.setScore(((String) sm_.convert(chars)).trim().toDouble());
			tag_ = "";
		}
		else if (tag_ == "pep_expect")
		{
			actual_peptide_hit_.metaRegistry().registerName("EValue", "E-value of e.g. Mascot searches", ""); /// @todo what E-value flag? (andreas)
			actual_peptide_hit_.setMetaValue("EValue", ((String)sm_.convert(chars)).trim().toDouble());
			tag_ = "";
		}
		else if (tag_ == "pep_homol")
		{			
			id_data_[peptide_identification_index_].setSignificanceThreshold(((String) sm_.convert(chars)).trim().toDouble());
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
			temp_identity = ((String) sm_.convert(chars)).trim().toDouble();
			actual_peptide_hit_.setMetaValue("homology_threshold", temp_homology);
			actual_peptide_hit_.setMetaValue("identity_threshold", temp_identity);
			if (temp_homology > temp_identity || temp_homology == 0)
			{
				id_data_[peptide_identification_index_].setSignificanceThreshold(temp_identity);				
			}
			tag_ = "";
		}
		else if (tag_ == "pep_seq")
		{
			AASequence temp_aa_sequence = AASequence(((String) sm_.convert(chars)).trim());
				
			// If everything is just read from the MascotXML file
			if (modified_peptides_.size() == 0)
			{			
				// fixed modifications
				for (vector<String>::const_iterator it = search_parameters_.fixed_modifications.begin(); it != search_parameters_.fixed_modifications.end(); ++it)
				{
					// e.g. Carboxymethyl (C)
					vector<String> mod_split;
					it->split(' ', mod_split);
					if (mod_split.size() == 2)
					{
						if (mod_split[1] == "(C-term)" || mod_split[1] == "(Protein C-term)")
						{
							temp_aa_sequence.setCTerminalModification(mod_split[0]);
						}
						else
						{
							if (mod_split[1] == "(N-term)" || mod_split[1] == "(Protein N-term)")
							{
								temp_aa_sequence.setNTerminalModification(mod_split[0]);
							}
							else
							{
								String origin = mod_split[1];
								origin.remove(')');
								origin.remove('(');
								for (Size i = 0; i != temp_aa_sequence.size(); ++i)
								{
									// best way we can check; because origin can be e.g. (STY)
									if (origin.hasSubstring(temp_aa_sequence[i].getOneLetterCode()))
									{
										temp_aa_sequence.setModification(i, mod_split[0]);
									}
								}
							}
						}
					}
					else
					{
						error(LOAD, String("Cannot parse fixed modification '") + *it  + "'");
					}
				}
			}
			actual_peptide_hit_.setSequence(temp_aa_sequence);
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
		else if (tag_ == "pep_var_mod_pos")
		{
			AASequence temp_aa_sequence = actual_peptide_hit_.getSequence();
			String temp_string = ((String) sm_.convert(chars)).trim();
			vector<String> parts;
			
			temp_string.split('.', parts);
			if (parts.size() == 3)
			{
				// handle internal modifications
				temp_string = parts[1];
				for (Size i = 0; i < temp_string.size(); ++i)
				{
					if (temp_string.at(i) != '0')
					{
						UInt temp_modification_index = String(temp_string[i]).toInt() - 1;
						String& temp_modification = search_parameters_.variable_modifications[temp_modification_index];

						// e.g. Carboxymethyl (C)
						vector<String> mod_split;
						temp_modification.split(' ', mod_split);

						if (mod_split.size() == 2)
						{
							// search this mod, if not directly use a general one
							temp_aa_sequence.setModification(i, mod_split[0]);
						}
						else
						{
							error(LOAD, String("Cannot parse variable modification '") + temp_modification  + "'");
						}
					}
				}

				temp_string = parts[0]; // N-term
				if (temp_string[0] != '0')
				{
					UInt temp_modification_index = String(temp_string[0]).toInt() - 1;
					String& temp_modification = search_parameters_.variable_modifications[temp_modification_index];
					vector<String> mod_split;
					temp_modification.split(' ', mod_split);

					if (mod_split.size() == 2)
					{
						temp_aa_sequence.setNTerminalModification(mod_split[0]);
					}
					else
					{
						error(LOAD, String("Cannot parse variable N-term modification '") + temp_modification  + "'");
					}
				}
				temp_string = parts[2]; // C-term
        if (temp_string[0] != '0')
        {
          UInt temp_modification_index = String(temp_string[0]).toInt() - 1;
          String& temp_modification = search_parameters_.variable_modifications[temp_modification_index];
          vector<String> mod_split;
          temp_modification.split(' ', mod_split);

          if (mod_split.size() == 2)
          {
            temp_aa_sequence.setCTerminalModification(mod_split[0]);
          }
          else
          {
            error(LOAD, String("Cannot parse variable C-term modification '") + temp_modification  + "'");
          }
        }
				
				actual_peptide_hit_.setSequence(temp_aa_sequence);
			}
			
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

			actual_title_ = title;
			if(modified_peptides_.find(title) != modified_peptides_.end())
			{
				vector<AASequence>& temp_hits = modified_peptides_[title];
				vector<PeptideHit> temp_peptide_hits = id_data_[actual_query_ - 1].getHits();
				
				if (temp_hits.size() != temp_peptide_hits.size())
				{
					warning(LOAD, "pepXML hits and Mascot hits are not the same");
				}

				// pepXML can contain more hits than MascotXML; hence we try to match all of them...
				// run-time is O(n^2) in the number of petide hits; should be a very small number
				
				for (Size i = 0; i < temp_peptide_hits.size(); ++i)
				{
					for (Size j = 0; j < temp_hits.size(); ++j)
					{
						if (temp_hits[j].isModified() && temp_hits[j].toUnmodifiedString() == temp_peptide_hits[i].getSequence().toUnmodifiedString())
						{
							temp_peptide_hits[i].setSequence(temp_hits[j]);
							break;
						}
					}
				}
				id_data_[actual_query_ - 1].setHits(temp_peptide_hits);
			}	
			title.split('_', parts);
			if (parts.size() == 2)
			{
				id_data_[actual_query_ - 1].setMetaValue("RT", parts[1].toDouble());
			}
		}
		else if (tag_ == "RTINSECONDS")
		{
			id_data_[actual_query_ - 1].setMetaValue("RT", ((String) sm_.convert(chars)).trim().toDouble());
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
		else if (tag_ == "TAXONOMY")
		{
			search_parameters_.taxonomy = (((String) sm_.convert(chars)).trim());
		}
		else if (tag_ == "CHARGE")
		{
			search_parameters_.charges = (((String) sm_.convert(chars)).trim());
		}
		else if (tag_ == "PFA")
		{
			search_parameters_.missed_cleavages = ((String) sm_.convert(chars)).trim().toInt();
		}	
		else if (tag_ == "MASS")
		{
			String temp_string = (((String) sm_.convert(chars)).trim());
			if (temp_string == "Monoisotopic")
			{
				search_parameters_.mass_type = ProteinIdentification::MONOISOTOPIC;
			}
			else if (temp_string == "Average")
			{
				search_parameters_.mass_type = ProteinIdentification::AVERAGE;
			}
		}
		else if (tag_ == "MODS")
		{
			String temp_string = (((String) sm_.convert(chars)).trim());
			temp_string.split(',', search_parameters_.fixed_modifications);
			if (search_parameters_.fixed_modifications.size() == 0 && temp_string != "")
			{
				search_parameters_.fixed_modifications.push_back(temp_string);
			}
		}
		else if (tag_ == "IT_MODS")
		{
			String temp_string = (((String) sm_.convert(chars)).trim());
			temp_string.split(',', search_parameters_.variable_modifications);
			if (search_parameters_.variable_modifications.size() == 0 && temp_string != "")
			{
				search_parameters_.variable_modifications.push_back(temp_string);
			}
		}
		else if (tag_ == "CLE")
		{
			String temp_string = (((String) sm_.convert(chars)).trim());
			if (temp_string == "Trypsin")
			{
				search_parameters_.enzyme = ProteinIdentification::TRYPSIN;
			}
			else if (temp_string == "PepsinA")
			{
				search_parameters_.enzyme = ProteinIdentification::PEPSIN_A;
			}
			else if (temp_string == "Chymotrypsin")
			{
				search_parameters_.enzyme = ProteinIdentification::CHYMOTRYPSIN;
			}
			else if (temp_string == "None")
			{
				search_parameters_.enzyme = ProteinIdentification::NO_ENZYME;
			}
			else
			{
				search_parameters_.enzyme = ProteinIdentification::UNKNOWN_ENZYME;
			}			
		}
		else if (tag_ == "TOL")
		{
			search_parameters_.precursor_tolerance = (((String) sm_.convert(chars)).trim()).toDouble();
		}
		else if (tag_ == "ITOL")
		{
			search_parameters_.peak_mass_tolerance = (((String) sm_.convert(chars)).trim()).toDouble();
		}
  }

	} // namespace Internal
} // namespace OpenMS
