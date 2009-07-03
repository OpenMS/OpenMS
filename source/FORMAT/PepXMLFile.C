// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Chris Bielow, Hendrik Weisser $
// $Authors: Chris Bielow, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/PepXMLFile.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <iostream>
#include <fstream>

using namespace std;

namespace OpenMS 
{

	PepXMLFile::PepXMLFile()
		: XMLHandler("","1.12"),
			XMLFile("/SCHEMAS/pepXML_v112.xsd","1.12"),
			protein_(0), 
			peptides_(0),
			experiment_(0),
			scan_map_(0),
			current_hit_(0),
			rt_tol_(0.5),
  		mz_tol_(0.5)
	{
	}


	PepXMLFile::~PepXMLFile()
	{
		if (scan_map_) delete scan_map_;
		if (current_hit_) delete current_hit_;
	}
	

	void PepXMLFile::load(const String& filename, ProteinIdentification& protein, vector<PeptideIdentification>& peptides, const String& experiment_name)
	{
		MSExperiment<> exp;
		load(filename, protein, peptides, experiment_name, exp);
	}

	
  void PepXMLFile::load(const String& filename, ProteinIdentification& protein, vector<PeptideIdentification>& peptides, const String& experiment_name, MSExperiment<>& experiment)
  { 
  	// initialize, load could be called several times
  	experiment_ = 0;
  	exp_name_ = "";
  	rt_tol_ = 0.5;
  	mz_tol_ = 0.5;
		// assume mass type "average" (in case element "search_summary" is missing)
		const ElementDB* db = ElementDB::getInstance();
		Element hydrogen = *db->getElement("Hydrogen");
		hydrogen_mass_ = hydrogen.getAverageWeight();	
  	
  	file_ = filename;	// filename for error messages in XMLHandler

		if (experiment_name != "")
		{
			// try and load the experiment now
			if (experiment.empty()) 
			{
				try
				{
					FileHandler().loadExperiment(experiment_name, experiment);
				}
				catch (Exception::FileNotFound ex)
				{
					warning(LOAD, String(ex.getMessage()) +	"; parsing pepXML without reference to the experiment");
				}
			}
			
			exp_name_ = File::removeExtension(experiment_name);
			
			if (!experiment.empty())
			{
				experiment_ = &experiment;
				MSExperiment<>::AreaType area = experiment_->getDataRange();
				// set tolerance to 1% of data range:
				rt_tol_ = (area.maxX() - area.minX()) * 0.01;
				mz_tol_ = (area.maxY() - area.minY()) * 0.01;
			}
		}

  	peptides.clear();
  	peptides_ = &peptides;
		protein_ = &protein;

		if (experiment_)
		{
			scan_map_ = new map<Size, Size>;
			Size scan = 0;
			for (MSExperiment<>::ConstIterator e_it = experiment_->begin(); e_it != experiment_->end(); ++e_it, ++scan)
			{
				String id = e_it->getNativeID();
				bool error = false;
				try
				{
					// expected format: "spectrum=#" (mzData) or "scan=#" (mzXML)
					Int num_id = id.suffix('=').toInt();
					if (num_id >= 0)
					{
						scan_map_->insert(scan_map_->end(), pair<Size, Size>(num_id, scan));
					}
					else
					{
						error = true;
					}
				}
				catch (Exception::ConversionError)
				{
					error = true;
				}
				if (error)
				{
					delete scan_map_;
					scan_map_ = 0;
					break;
				}
			}
		}

		wrong_experiment_ = false;
		parse_(filename, this);

		protein_->setSearchParameters(params_);
		for (set<String>::iterator a_it = accessions_.begin(); a_it != accessions_.end(); ++a_it)
		{
			ProteinHit hit;
			hit.setAccession(*a_it);
			protein_->insertHit(hit);
		}
    
    // reset members
		actual_sequence_.clear();
		actual_modifications_.clear();
		fixed_modifications_.clear();
		variable_modifications_.clear();
		exp_name_.clear();
		prot_id_.clear();
		date_.clear();
		accessions_.clear();
		params_ = ProteinIdentification::SearchParameters();
		protein_ = 0;
		peptides_ = 0;
		experiment_ = 0;
		if (scan_map_)
		{
			delete scan_map_;
			scan_map_ = 0;
		}
		if (current_hit_)
		{
			delete current_hit_;
			current_hit_ = 0;
		}
  }

	
	void PepXMLFile::store(const String& filename, std::vector<ProteinIdentification>& /* protein_ids */, std::vector<PeptideIdentification>& peptide_ids)
	{
		ofstream f(filename.c_str());
		if (!f)
		{
			throw Exception::UnableToCreateFile(__FILE__, __LINE__,
																					__PRETTY_FUNCTION__, filename);
		}

		f.precision(writtenDigits<DoubleReal>());

		f << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
		f << "<msms_pipeline_analysis date=\"2007-12-05T17:49:46\" xmlns=\"http://regis-web.systemsbiology.net/pepXML\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://regis-web.systemsbiology.net/pepXML http://www.matrixscience.com/xmlns/schema/pepXML_v18/pepXML_v18.xsd\" summary_xml=\".xml\">" << endl;
		f << "<msms_run_summary base_name=\"" << filename << "\" msManufacturer=\"ThermoFinnigan\" msModel=\"LCQ Classic\" msIonization=\"ESI\" msMassAnalyzer=\"Ion Trap\" msDetector=\"UNKNOWN\" raw_data_type=\"raw\" raw_data=\".mzXML\">"
			<< endl;

		f << "<sample_enzyme name=\"trypsin\">" << endl;
		f << "<specificity cut=\"KR\" no_cut=\"P\" sense=\"C\"/>" << endl;
		f << "</sample_enzyme>" << endl;

		f << "<search_summary base_name=\"" << filename << "\" search_engine=\"SEQUEST\" precursor_mass_type=\"average\" fragment_mass_type=\"average\" out_data_type=\"out\" out_data=\".tgz\" search_id=\"1\">" << endl;
	//	f << "		<search_summary base_name=\"\" search_engine=\"Mascot\" precursor_mass_type=\"average\" fragment_mass_type=\"average\">" << endl;
		f << "		<search_database local_path=\"dbase/ipi.HUMAN.fasta.v2.31\" type=\"AA\"/>" << endl;

		// register modifications
		set<String> aa_mods;
		set<String> n_term_mods, c_term_mods;
		for (vector<PeptideIdentification>::const_iterator it = peptide_ids.begin();
				 it != peptide_ids.end(); ++it)
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
		for (set<String>::const_iterator it = aa_mods.begin();
				 it != aa_mods.end(); ++it)
		{
			ResidueModification mod = ModificationsDB::getInstance()->
				getModification(*it);
			f << "<aminoacid_modification aminoacid=\"" << mod.getOrigin()
				<< "\" massdiff=\"" << mod.getDiffMonoMass() << "\" mass=\""
				<< mod.getMonoMass()
				<< "\" variable=\"Y\" binary=\"N\" description=\"" << *it << "\"/>"
				<< endl;
		}

		for (set<String>::const_iterator it = n_term_mods.begin();
				 it != n_term_mods.end(); ++it)
		{
			ResidueModification mod = ModificationsDB::getInstance()->
				getModification(*it);
			f << "<terminal_modification terminus=\"n\" massdiff=\""
				<< mod.getDiffMonoMass() << "\" mass=\"" << mod.getMonoMass()
				<< "\" variable=\"Y\" description=\"" << *it
				<< "\" protein_terminus=\"\"/>" << endl;
		}

		for (set<String>::const_iterator it = n_term_mods.begin();
				 it != n_term_mods.end(); ++it)
		{
			ResidueModification mod = ModificationsDB::getInstance()->
				getModification(*it);
			f << "<terminal_modification terminus=\"c\" massdiff=\""
				<< mod.getDiffMonoMass() << "\" mass=\"" << mod.getMonoMass()
				<< "\" variable=\"Y\" description=\"" << *it
				<< "\" protein_terminus=\"\"/>" << endl;
		}

		f << "    </search_summary>" << endl;
		f << "    <analysis_timestamp analysis=\"peptideprophet\" time=\"2007-12-05T17:49:52\" id=\"1\"/>" << endl;
		

		UInt count(1);
		for (vector<PeptideIdentification>::const_iterator it = peptide_ids.begin();
				 it != peptide_ids.end(); ++it, count++)
		{
			if (it->getHits().size() > 0)
			{
				PeptideHit h = *it->getHits().begin();
				double precursor_neutral_mass(0);
				precursor_neutral_mass = h.getSequence().getAverageWeight();
				
				f << "		<spectrum_query spectrum=\"" << count << "\" start_scan=\""
					<< count << "\" end_scan=\"" << count
					<< "\" precursor_neutral_mass=\"" << precursor_neutral_mass
					<< "\" assumed_charge=\"" << h.getCharge() << "\" index=\"" << count
					<< "\">" << endl;
				f << " 		<search_result>" << endl;
				f << "			<search_hit hit_rank=\"1\" peptide=\""
					<< h.getSequence().toUnmodifiedString() << "\" peptide_prev_aa=\""
					<< h.getAABefore() << "\" peptide_next_aa=\"" << h.getAAAfter()
					<< "\" protein=\"Protein1\" num_tot_proteins=\"1\" num_matched_ions=\"0\" tot_num_ions=\"0\" calc_neutral_pep_mass=\"" << precursor_neutral_mass
					<< "\" massdiff=\"\" num_tol_term=\"0\" num_missed_cleavages=\"0\" is_rejected=\"0\" protein_descr=\"Protein No. 1\">" << endl;
				if (h.getSequence().isModified())
				{
					f << "      <modification_info modified_peptide=\""
						<< h.getSequence() << "\"";

					if (h.getSequence().hasNTerminalModification())
					{
						f << " mod_nterm_mass=\"" << ModificationsDB::getInstance()->
							getModification(h.getSequence().getNTerminalModification()).
							getMonoMass() << "\"";
					}
					
					if (h.getSequence().hasCTerminalModification())
					{
						f << "mod_cterm_mass=\"" << ModificationsDB::getInstance()->
							getModification(h.getSequence().getCTerminalModification()).
							getMonoMass() << "\"";
					}

					f << ">" << endl;

					for (Size i = 0; i != h.getSequence().size(); ++i)
					{
						if (h.getSequence()[i].isModified())
						{
							f << "         <mod_aminoacid_mass position=\"" << i
								<< "\" mass=\"" << ModificationsDB::getInstance()->
								getModification(h.getSequence()[i].getModification()).
								getMonoMass() << "\"/>" << endl;
						}
					}

					f << "      </modification_info>" << endl;
									
				}

				f << " 			<analysis_result analysis=\"peptideprophet\">" << endl;
				f << "			<peptideprophet_result probability=\"" << h.getScore()
					<< "\" all_ntt_prob=\"(" << h.getScore() << "," << h.getScore()
					<< "," << h.getScore() << ")\">" << endl;
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

  void PepXMLFile::matchModification_(DoubleReal mass, String& modification_description)
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


	void PepXMLFile::startElement(const XMLCh* const /*uri*/,
																const XMLCh* const /*local_name*/,
																const XMLCh* const qname,
																const xercesc::Attributes& attributes)
	{		
		String element = sm_.convert(qname);
		
		//cout << "Start: " << element << endl;
	
		if (element == "msms_run_summary")
		{
			if (!exp_name_.empty())
			{
				String base_name = attributeAsString_(attributes, "base_name");
				if (!base_name.hasSuffix(exp_name_))
				{
					wrong_experiment_ = true;
					return;
				}
				else
				{
					wrong_experiment_ = false;
				}
			}
		}
		
		else if (wrong_experiment_)
		{
			// do nothing here (this case exists to prevent parsing of elements for
			// experiments we're not interested in)
		}
		
		// now, elements occurring more frequently are generally closer to the top
		else if (element == "search_score")
		{
			String name = attributeAsString_(attributes, "name");
			DoubleReal value = attributeAsDouble_(attributes, "value");
			// TODO: deal with different scores
			if (name == "hyperscore")
			{ // X!Tandem score
				current_hit_->setScore(value);
				current_pep_->setScoreType(name); // add "X!Tandem" to name?
				current_pep_->setHigherScoreBetter(true);
			}
			if (name == "xcorr")
			{ // Sequest score
				current_hit_->setScore(value);
				current_pep_->setScoreType(name); // add "Sequest" to name?
				current_pep_->setHigherScoreBetter(true);
			}
		}

		else if (element == "search_hit")
		{
			actual_sequence_ = attributeAsString_(attributes, "peptide");
			current_hit_ = new PeptideHit;
			current_hit_->setRank(attributeAsInt_(attributes, "hit_rank"));
			String protein = attributeAsString_(attributes, "protein");
			current_hit_->addProteinAccession(protein);
			accessions_.insert(protein);
			current_hit_->setCharge(charge_); // from "spectrum_query" tag
			String prev_aa, next_aa;
			if (optionalAttributeAsString_(prev_aa, attributes, "peptide_prev_aa")) current_hit_->setAABefore(prev_aa[0]);
			if (optionalAttributeAsString_(next_aa, attributes, "peptide_next_aa")) current_hit_->setAAAfter(next_aa[0]);
		}

		else if (element == "spectrum_query")
		{
			DoubleReal mz, rt = 0, mass = attributeAsDouble_(attributes, "precursor_neutral_mass");
			charge_ = attributeAsInt_(attributes, "assumed_charge");
			mz = (mass + hydrogen_mass_ * charge_) / charge_;
			bool rt_present = optionalAttributeAsDouble_(rt, attributes, "retention_time_sec");
			// assume only one scan, i.e. ignore "end_scan":
			Size scan = attributeAsInt_(attributes, "start_scan");
			if (scan_map_) scan = (*scan_map_)[scan];
			if (experiment_)
			{ // get precursor information
				MSSpectrum<> spec = (*experiment_)[scan];
				if ((spec.getMSLevel() == 2) && (!rt_present || Math::approximatelyEqual(spec.getRT(), rt, 0.001)))
				{
					// otherwise: wrong scan
					DoubleReal prec_mz = 0, prec_rt = 0;
					vector<Precursor> precursors = spec.getPrecursors();
					if (precursors.size()) prec_mz = precursors[0].getMZ(); // assume only one precursor
					MSExperiment<>::ConstIterator it = experiment_->getPrecursorSpectrum(experiment_->begin() + scan);
					if (it != experiment_->end())
					{
						prec_rt = it->getRT();
					}

					// check if "rt"/"mz" are similar to "prec_rt"/"prec_mz"
					// (otherwise, precursor mapping is wrong)
					if ((prec_mz > 0) && Math::approximatelyEqual(prec_mz, mz, mz_tol_)	&& (prec_rt > 0) && (!rt_present || Math::approximatelyEqual(prec_rt, rt, rt_tol_)))
					{
// 						DoubleReal diff;
// 						diff = mz - prec_mz;
// 						cout << "m/z difference: " << diff << " ("
// 								 << diff / max(mz, prec_mz) * 100 << "%)\n";
//  					diff = rt - prec_rt;
//  					cout << "RT difference: " << diff << " ("
//  							 << diff / max(rt, prec_rt) * 100 << "%)\n" << endl;
						mz = prec_mz;
						rt = prec_rt;
					}
				}
			}
			PeptideIdentification peptide;
			peptide.setMetaValue("RT", rt);
			peptide.setMetaValue("MZ", mz);
			peptide.setIdentifier(prot_id_);
			peptides_->push_back(peptide);
			current_pep_ = --(peptides_->end());
		}
	
		else if (element == "peptideprophet_result")
		{
			// PeptideProphet probability overwrites original search score!
			// maybe TODO: deal with meta data associated with PeptideProphet search
			DoubleReal value = attributeAsDouble_(attributes, "probability");
			current_hit_->setScore(value);
			current_pep_->setScoreType("PeptideProphet probability");
			current_pep_->setHigherScoreBetter(true);
		}

		else if (element == "alternative_protein")
		{
			String protein = attributeAsString_(attributes, "protein");
			current_hit_->addProteinAccession(protein);
			accessions_.insert(protein);
		}

		else if (element == "mod_aminoacid_mass")
		{
			DoubleReal modification_mass = 0.;
			UInt 			 modification_position = 0;
			String 		 temp_description = "";
			
			modification_position = attributeAsInt_(attributes, "position");
			modification_mass = attributeAsDouble_(attributes, "mass");
			
			matchModification_(modification_mass, temp_description);
			
			// the modification position is 1-based
			actual_modifications_.push_back(make_pair(temp_description, modification_position));
		}

		else if ((element == "aminoacid_modification") || (element == "terminal_modification"))
			// <terminal_modification terminus="n" massdiff="+108.05" mass="109.06" variable="N" protein_terminus="" description="dNIC (N-term)"/>
		{
			String desc, is_variable = attributeAsString_(attributes, "variable");
			bool has_desc = optionalAttributeAsString_(desc, attributes,
																								 "description");
			if (!has_desc)
			{ // generate a dummy description
				String sign, massdiff = attributeAsString_(attributes, "massdiff");
				if (!massdiff.hasPrefix("-") && !massdiff.hasPrefix("+"))
				{
					sign = "+";
				}
				desc = attributeAsString_(attributes, "aminoacid") + sign + massdiff;
			}
			if (is_variable == "Y")
			{
				variable_modifications_.push_back(make_pair(desc, attributeAsDouble_(attributes, "mass")));
				params_.variable_modifications.push_back(desc);
			}
			else
			{
				fixed_modifications_.push_back(desc);
				params_.fixed_modifications.push_back(desc);
			}
		}
		
		else if (element == "search_summary")
		{
			const ElementDB* db = ElementDB::getInstance();
			Element hydrogen = *db->getElement("Hydrogen");
			String mass_type = attributeAsString_(attributes, "precursor_mass_type");
			if (mass_type == "monoisotopic")
			{
				hydrogen_mass_ = hydrogen.getMonoWeight();
			}
			else
			{
				hydrogen_mass_ = hydrogen.getAverageWeight();
				if (mass_type != "average")
				{
					error(LOAD,	"'precursor_mass_type' attribute of 'search_summary' tag should be 'monoisotopic' or 'average', not '" + mass_type + "' (assuming 'average')");
				}
			}
			// assuming "SearchParameters::mass_type" refers to the fragment mass
			mass_type = attributeAsString_(attributes, "fragment_mass_type");
			if (mass_type == "monoisotopic") params_.mass_type = ProteinIdentification::MONOISOTOPIC;
			else if (mass_type == "average") params_.mass_type = ProteinIdentification::AVERAGE;
			else error(LOAD,	"'fragment_mass_type' attribute of 'search_summary' tag should be 'monoisotopic' or 'average', not '" + mass_type + "'");
			String search_engine = attributeAsString_(attributes, "search_engine");
			protein_->setSearchEngine(search_engine);
			// generate identifier from search engine and date:
			prot_id_ = search_engine + "_" + date_.getDate();
			protein_->setIdentifier(prot_id_);
		}

		else if (element == "sample_enzyme")
		{
			String name = attributeAsString_(attributes, "name");
			name.toLower();
			// spelling of enzyme names in pepXML?
			if (name.hasPrefix("trypsin")) params_.enzyme = ProteinIdentification::TRYPSIN;
			else if (name.hasPrefix("pepsin")) params_.enzyme = ProteinIdentification::PEPSIN_A;
			else if (name.hasPrefix("protease")) params_.enzyme = ProteinIdentification::PROTEASE_K;
			else if (name.hasPrefix("chymotrypsin")) params_.enzyme = ProteinIdentification::CHYMOTRYPSIN;
			else params_.enzyme = ProteinIdentification::UNKNOWN_ENZYME;
		}

		else if (element == "search_database")
		{
			params_.db = attributeAsString_(attributes, "local_path");
		}
		
		else if (element == "msms_pipeline_analysis")
		{
			String date = attributeAsString_(attributes, "date");
			// fix for corrupted xs:dateTime format:
			if ((date[4] == ':') && (date[7] == ':') && (date[10] == ':'))
			{
				error(LOAD, "Format of attribute 'date' in tag 'msms_pipeline_analysis' does not comply with standard 'xs:dateTime'");
				date[4] = '-';
				date[7] = '-';
				date[10] = 'T';
			}
			date_ = asDateTime_(date);
			protein_->setDateTime(date_);
			// "prot_id_" will be overwritten if element "search_summary" is present
			prot_id_ = "unknown_" + date_.getDate();
			protein_->setIdentifier(prot_id_);
		}
		
	}

	
	void PepXMLFile::endElement(const XMLCh* const /*uri*/,
															const XMLCh* const /*local_name*/,
															const XMLCh* const qname)
	{
		String element = sm_.convert(qname);

		if (wrong_experiment_)
		{
			//TODO?
		}
		else if (element == "search_hit")
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

			current_hit_->setSequence(temp_aa_sequence);
			current_pep_->insertHit(*current_hit_);
			delete current_hit_;
			current_hit_ = 0;
			actual_modifications_.clear();						
		}
	}

} // namespace OpenMS
