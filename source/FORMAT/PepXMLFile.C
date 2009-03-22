// -*- Mode: C++; tab-widt: 2; -*-
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

#include <OpenMS/FORMAT/PepXMLFile.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

#include <iostream>
#include <fstream>

using namespace std;

namespace OpenMS 
{

	PepXMLFile::PepXMLFile()
		: XMLHandler("","1.8"),
			XMLFile("/SCHEMAS/PepXML_1_8.xsd","1.8"),
			peptides_(0)
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
		variable_modifications_ = vector< pair<String, DoubleReal> >();
		fixed_modifications_ = vector<String>();
  }

	void PepXMLFile::store(const String& filename,  std::vector<ProteinIdentification>& protein_ids, std::vector<PeptideIdentification>& peptide_ids)
	{
		ofstream f(filename.c_str());
		if (!f)
		{
			throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
		}

		f.precision(writtenDigits<DoubleReal>());


		f << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
		f << "<msms_pipeline_analysis date=\"2007-12-05T17:49:46\" xmlns=\"http://regis-web.systemsbiology.net/pepXML\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://regis-web.systemsbiology.net/pepXML http://www.matrixscience.com/xmlns/schema/pepXML_v18/pepXML_v18.xsd\" summary_xml=\".xml\">" << endl;
		f << "<msms_run_summary base_name=\"" << filename << "\" msManufacturer=\"ThermoFinnigan\" msModel=\"LCQ Classic\" msIonization=\"ESI\" msMassAnalyzer=\"Ion Trap\" msDetector=\"UNKNOWN\" raw_data_type=\"raw\" raw_data=\".mzXML\">" << endl;

		f << "<sample_enzyme name=\"trypsin\">" << endl;
		f << "<specificity cut=\"KR\" no_cut=\"P\" sense=\"C\"/>" << endl;
		f << "</sample_enzyme>" << endl;

		f << "<search_summary base_name=\"" << filename << "\" search_engine=\"SEQUEST\" precursor_mass_type=\"average\" fragment_mass_type=\"average\" out_data_type=\"out\" out_data=\".tgz\" search_id=\"1\">" << endl;
	//	f << "		<search_summary base_name=\"\" search_engine=\"Mascot\" precursor_mass_type=\"average\" fragment_mass_type=\"average\">" << endl;
		f << "		<search_database local_path=\"dbase/ipi.HUMAN.fasta.v2.31\" type=\"AA\"/>" << endl;

		// register modifications
		set<String> aa_mods;
		set<String> n_term_mods, c_term_mods;
		for (vector<PeptideIdentification>::const_iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
		{
			if (it->getHits().size() > 0)
			{
				PeptideHit h = *it->getHits().begin();
				
				if (h.getSequence().isModified())
				{
					AASequence p = h.getSequence();
					if (p.hasNTerminalModification())
					{
						n_term_mods.insert(p.getNTerminalModification());
					}
					if (p.hasCTerminalModification())
					{
						c_term_mods.insert(p.getCTerminalModification());
					}

					for (Size i = 0; i != p.size(); ++i)
					{
						if (p[i].isModified())
						{
							aa_mods.insert(p[i].getModification());
						}
					}
				}				
			}
		}

		// write modifications definitions
		// <aminoacid_modification aminoacid="C" massdiff="+58.01" mass="161.014664" variable="Y" binary="N" description="Carboxymethyl (C)"/>
		for (set<String>::const_iterator it = aa_mods.begin(); it != aa_mods.end(); ++it)
		{
			ResidueModification mod = ModificationsDB::getInstance()->getModification(*it);
			f << "<aminoacid_modification aminoacid=\"" << mod.getOrigin() << "\" massdiff=\"" << mod.getDiffMonoMass() << "\" mass=\"" << mod.getMonoMass() << "\" variable=\"Y\" binary=\"N\" description=\"" << *it << "\"/>" << endl;
		}

		for (set<String>::const_iterator it = n_term_mods.begin(); it != n_term_mods.end(); ++it)
		{
			ResidueModification mod = ModificationsDB::getInstance()->getModification(*it);
			f << "<terminal_modification terminus=\"n\" massdiff=\"" << mod.getDiffMonoMass() << "\" mass=\"" << mod.getMonoMass() << "\" variable=\"Y\" description=\"" << *it << "\" protein_terminus=\"\"/>" << endl;
		}

		for (set<String>::const_iterator it = n_term_mods.begin(); it != n_term_mods.end(); ++it)
		{
			ResidueModification mod = ModificationsDB::getInstance()->getModification(*it);
			f << "<terminal_modification terminus=\"c\" massdiff=\"" << mod.getDiffMonoMass() << "\" mass=\"" << mod.getMonoMass() << "\" variable=\"Y\" description=\"" << *it << "\" protein_terminus=\"\"/>" << endl;
		}

		f << "    </search_summary>" << endl;
		f << "    <analysis_timestamp analysis=\"peptideprophet\" time=\"2007-12-05T17:49:52\" id=\"1\"/>" << endl;
		

		UInt count(1);
		for (vector<PeptideIdentification>::const_iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it, count++)
		{
			if (it->getHits().size() > 0)
			{
				PeptideHit h = *it->getHits().begin();
				double precursor_neutral_mass(0);
				precursor_neutral_mass = h.getSequence().getAverageWeight();
				
				f << "		<spectrum_query spectrum=\"" << count << "\" start_scan=\"" << count << "\" end_scan=\"" << count << "\" precursor_neutral_mass=\"" << precursor_neutral_mass << "\" assumed_charge=\"" << h.getCharge() << "\" index=\"" << count << "\">" << endl;
				f << " 		<search_result>" << endl;
				f << "			<search_hit hit_rank=\"1\" peptide=\"" << h.getSequence().toUnmodifiedString() << "\" peptide_prev_aa=\"" << h.getAABefore() << "\" peptide_next_aa=\"" << h.getAAAfter() << "\" protein=\"Protein1\" num_tot_proteins=\"1\" num_matched_ions=\"0\" tot_num_ions=\"0\" calc_neutral_pep_mass=\"" << precursor_neutral_mass << "\" massdiff=\"\" num_tol_term=\"0\" num_missed_cleavages=\"0\" is_rejected=\"0\" protein_descr=\"Protein No. 1\">" << endl;
				if (h.getSequence().isModified())
				{
					f << "      <modification_info modified_peptide=\"" << h.getSequence() << "\"";

					if (h.getSequence().hasNTerminalModification())
					{
						f << " mod_nterm_mass=\"" << ModificationsDB::getInstance()->getModification(h.getSequence().getNTerminalModification()).getMonoMass() << "\"";
					}
					
					if (h.getSequence().hasCTerminalModification())
					{
						f << "mod_cterm_mass=\"" << ModificationsDB::getInstance()->getModification(h.getSequence().getCTerminalModification()).getMonoMass() << "\"";
					}

					f << ">" << endl;

					for (Size i = 0; i != h.getSequence().size(); ++i)
					{
						if (h.getSequence()[i].isModified())
						{
							f << "         <mod_aminoacid_mass position=\"" << i << "\" mass=\"" << ModificationsDB::getInstance()->getModification(h.getSequence()[i].getModification()).getMonoMass() << "\"/>" << endl;
						}
					}

					f << "      </modification_info>" << endl;
									
				}

				f << " 			<analysis_result analysis=\"peptideprophet\">" << endl;
				f << "			<peptideprophet_result probability=\"" << h.getScore() << "\" all_ntt_prob=\"(" << h.getScore() << "," << h.getScore() << "," << h.getScore() << ")\">" << endl;
				f << "			</peptideprophet_result>" << endl;
				f << "			</analysis_result>" << endl;
				f << "			</search_hit>" << endl;	
				f << "		</search_result>" << endl;
				f << "		</spectrum_query>" << endl;
			}
		}

		f << "	</msms_run_summary>" << endl;
		f << "</msms_pipeline_analysis>" << endl;
		
		f.close();
		
	}

  void PepXMLFile::matchModification_(DoubleReal mass,
  																		String&		 modification_description)
	{
		UInt i = 0;
		bool found = false;
		DoubleReal difference = 0.; 
		
		while(i < variable_modifications_.size() && !found)
		{
			difference = variable_modifications_[i].second - mass;
			if (difference < 0)
			{
				difference *= -1;
			}
			if (difference < 0.001)
			{
				modification_description = variable_modifications_[i].first;
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
				variable_modifications_.push_back(make_pair(attributeAsString_(attributes,"description"), 
																					 					attributeAsDouble_(attributes,"mass")));
			}
			else
			{
				fixed_modifications_.push_back(attributeAsString_(attributes,"description"));
			}
		}	

		// <terminal_modification terminus="n" massdiff="+108.05" mass="109.06" variable="N" protein_terminus="" description="dNIC (N-term)"/>
		if (element == "terminal_modification")
		{
			String temp_string = attributeAsString_(attributes, "variable");
			if (temp_string == "Y")
			{
				variable_modifications_.push_back(make_pair(attributeAsString_(attributes, "description"),
																										attributeAsDouble_(attributes, "mass")));
												
			}
			else
			{
				fixed_modifications_.push_back(attributeAsString_(attributes, "description"));
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
			for (vector<pair<String, UInt> >::const_iterator it = actual_modifications_.begin(); it != actual_modifications_.end(); ++it)
			{
				// e.g. Carboxymethyl (C)
				vector<String> mod_split;
				it->first.split(' ', mod_split);
				if (mod_split.size() == 2)
				{
					if (mod_split[1] == "(C-term)")
					{
						temp_aa_sequence.setCTerminalModification(mod_split[0]);
					}
					else
					{
						if (mod_split[1] == "(N-term)")
						{
							temp_aa_sequence.setNTerminalModification(mod_split[0]);
						}
						else
						{
							// search this mod, if not directly use a general one
							temp_aa_sequence.setModification(it->second - 1, mod_split[0]);
						}
					}
				}
				else
				{
					error(LOAD, String("Cannot parse modification '") + it->first + "@" + it->second + "'");
				}
			}

			// fixed modifications
			for (vector<String>::const_iterator it = fixed_modifications_.begin(); it != fixed_modifications_.end(); ++it)
			{
				// e.g. Carboxymethyl (C)
				vector<String> mod_split;
				it->split(' ', mod_split);
				if (mod_split.size() == 2)
				{
					if (mod_split[1] == "(C-term)")
					{
						temp_aa_sequence.setCTerminalModification(mod_split[0]);
					}
					else
					{
						if (mod_split[1] == "(N-term)")
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
					error(LOAD, String("Cannot parse fixed modification '") + *it + "'");
				}
			}

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
