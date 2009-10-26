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
#include <OpenMS/CHEMISTRY/ResidueDB.h>
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
			rt_tol_(10.0),
  		mz_tol_(10.0)
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
    actual_sequence_ = "";
    actual_modifications_.clear();
    experiment_ = 0;
    exp_name_ = "";
    hydrogen_mass_ = 0;
    prot_id_ = "";
    params_ = ProteinIdentification::SearchParameters();
    charge_ = 0;
    accessions_.clear();
    rt_tol_ = 0;
		mz_tol_ = 0;
    fixed_modifications_.clear();



  	rt_tol_ = 10.0;
  	mz_tol_ = 10.0;
		// assume mass type "average" (in case element "search_summary" is missing)
		const ElementDB* db = ElementDB::getInstance();
		Element hydrogen = *db->getElement("Hydrogen");
		hydrogen_mass_ = hydrogen.getMonoWeight();	
  	
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
				// set tolerance to 1% of data range (if above a sensible minimum):
				rt_tol_ = max((area.maxX() - area.minX()) * 0.01, rt_tol_);
				mz_tol_ = max((area.maxY() - area.minY()) * 0.01, mz_tol_);
			}
		}

  	peptides.clear();
  	peptides_ = &peptides;
		protein = ProteinIdentification();
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

		if (peptides.empty())
		{
			warning(LOAD, "No data found for experiment name '" + experiment_name +
							"'");
		}
		
    // reset members
		actual_sequence_.clear();
		actual_modifications_.clear();
		fixed_modifications_.clear();
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

	
	void PepXMLFile::store(const String& filename, std::vector<ProteinIdentification>& protein_ids, std::vector<PeptideIdentification>& peptide_ids)
	{
		ofstream f(filename.c_str());
		if (!f)
		{
			throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
		}

		String search_engine_name;
		ProteinIdentification::SearchParameters search_params;
		if (protein_ids.size() != 0)
		{
			if (protein_ids.size() > 1)
			{
				warning(STORE, "More than one protein identification defined, only first one is written into pepXML more are not supported.");
			}
			search_params = protein_ids.begin()->getSearchParameters();
			search_engine_name = protein_ids.begin()->getSearchEngine();
		}

		f.precision(writtenDigits<DoubleReal>());

		f << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << "\n";
		f << "<msms_pipeline_analysis date=\"2007-12-05T17:49:46\" xmlns=\"http://regis-web.systemsbiology.net/pepXML\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://regis-web.systemsbiology.net/pepXML http://www.matrixscience.com/xmlns/schema/pepXML_v18/pepXML_v18.xsd\" summary_xml=\".xml\">" << "\n";
		f << "<msms_run_summary base_name=\"" << File::basename(filename) << "\" raw_data_type=\"raw\" raw_data=\".mzXML\" search_engine=\"" << search_engine_name << "\">" << "\n";

		f << "<sample_enzyme name=\"trypsin\">" << "\n";
		f << "<specificity cut=\"KR\" no_cut=\"P\" sense=\"C\"/>" << "\n";
		f << "</sample_enzyme>" << "\n";

		f << "<search_summary base_name=\"" << File::basename(filename);
		f << "\" search_engine=\"" << search_engine_name;	
		f << "\" precursor_mass_type=\"";
		if (search_params.mass_type == ProteinIdentification::MONOISOTOPIC)
		{
			f << "monoisotopic";
		}
		else
		{
			f << "average";
		}
		f << "\" fragment_mass_type=\"";
		if (search_params.mass_type == ProteinIdentification::MONOISOTOPIC)
		{
			f << "monoisotopic";
		}
		else
		{
			f << "average";
		}
		f << "\" out_data_type=\"\" out_data=\"\" search_id=\"1\">" << "\n";
		f << "		<search_database local_path=\"" << search_params.db << "\" type=\"AA\"/>" << "\n";


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
						n_term_mods.insert(ModificationsDB::getInstance()->getTerminalModification(p.getNTerminalModification(), ResidueModification::N_TERM).getFullId());
					}
					if (p.hasCTerminalModification())
					{
						c_term_mods.insert(ModificationsDB::getInstance()->getTerminalModification(p.getCTerminalModification(), ResidueModification::C_TERM).getFullId());
					}

					for (Size i = 0; i != p.size(); ++i)
					{
						if (p[i].isModified())
						{
							aa_mods.insert(ModificationsDB::getInstance()->getModification(p[i].getOneLetterCode(), p[i].getModification(), ResidueModification::ANYWHERE).getFullId());
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
			const ResidueModification& mod = ModificationsDB::getInstance()->getModification(*it);
			f << "<aminoacid_modification aminoacid=\"" << mod.getOrigin()
				<< "\" massdiff=\"" << precisionWrapper(mod.getDiffMonoMass()) << "\" mass=\""
				<< precisionWrapper(mod.getMonoMass())
				<< "\" variable=\"Y\" binary=\"N\" description=\"" << *it << "\"/>"
				<< "\n";
		}

		for (set<String>::const_iterator it = n_term_mods.begin(); it != n_term_mods.end(); ++it)
		{
			const ResidueModification& mod = ModificationsDB::getInstance()->
				getModification(*it);
			f << "<terminal_modification terminus=\"n\" massdiff=\""
				<< precisionWrapper(mod.getDiffMonoMass()) << "\" mass=\"" << precisionWrapper(mod.getMonoMass())
				<< "\" variable=\"Y\" description=\"" << *it
				<< "\" protein_terminus=\"\"/>" << "\n";
		}

		for (set<String>::const_iterator it = c_term_mods.begin(); it != c_term_mods.end(); ++it)
		{
			const ResidueModification& mod = ModificationsDB::getInstance()->getModification(*it);
			f << "<terminal_modification terminus=\"c\" massdiff=\""
				<< precisionWrapper(mod.getDiffMonoMass()) << "\" mass=\"" << precisionWrapper(mod.getMonoMass())
				<< "\" variable=\"Y\" description=\"" << *it
				<< "\" protein_terminus=\"\"/>" << "\n";
		}

		f << "    </search_summary>" << "\n";
		f << "    <analysis_timestamp analysis=\"peptideprophet\" time=\"2007-12-05T17:49:52\" id=\"1\"/>" << "\n";
		

		Size count(1);
		for (vector<PeptideIdentification>::const_iterator it = peptide_ids.begin();
				 it != peptide_ids.end(); ++it, count++)
		{
			if (it->getHits().size() > 0)
			{
				PeptideHit h = *it->getHits().begin();
				AASequence seq = h.getSequence();
				DoubleReal precursor_neutral_mass = seq.getMonoWeight();

				f << "		<spectrum_query spectrum=\"" << count << "\" start_scan=\""
					<< count << "\" end_scan=\"" << count
					<< "\" precursor_neutral_mass=\"" << precisionWrapper(precursor_neutral_mass)
					<< "\" assumed_charge=\"" << h.getCharge() << "\" index=\"" << count
					<< "\">" << "\n";
				f << " 		<search_result>" << "\n";
				f << "			<search_hit hit_rank=\"1\" peptide=\""
					<< seq.toUnmodifiedString() << "\" peptide_prev_aa=\""
					<< h.getAABefore() << "\" peptide_next_aa=\"" << h.getAAAfter()
					<< "\" protein=\"Protein1\" num_tot_proteins=\"1\" num_matched_ions=\"0\" tot_num_ions=\"0\" calc_neutral_pep_mass=\"" << precisionWrapper(precursor_neutral_mass)
					<< "\" massdiff=\"\" num_tol_term=\"0\" num_missed_cleavages=\"0\" is_rejected=\"0\" protein_descr=\"Protein No. 1\">" << "\n";
				if (seq.isModified())
				{
					f << "      <modification_info modified_peptide=\""
						<< seq << "\"";

					if (seq.hasNTerminalModification())
					{
						const ResidueModification& mod = ModificationsDB::getInstance()->getTerminalModification(seq.getNTerminalModification(), ResidueModification::N_TERM);
						f << " mod_nterm_mass=\"" << 
						precisionWrapper(mod.getMonoMass() + seq[(Size)0].getMonoWeight(Residue::Internal)) << "\"";
					}
					
					if (seq.hasCTerminalModification())
					{
						const ResidueModification& mod = ModificationsDB::getInstance()->getTerminalModification(seq.getCTerminalModification(), ResidueModification::C_TERM);
						f << "mod_cterm_mass=\"" << 
						precisionWrapper(mod.getMonoMass() + seq[seq.size() - 1].getMonoWeight(Residue::Internal)) << "\"";
					}

					f << ">" << "\n";

					for (Size i = 0; i != seq.size(); ++i)
					{
						if (seq[i].isModified())
						{
							const ResidueModification& mod = ModificationsDB::getInstance()->getModification(seq[i].getOneLetterCode(), seq[i].getModification(), ResidueModification::ANYWHERE);
							f << "         <mod_aminoacid_mass position=\"" << i
								<< "\" mass=\"" << 
								precisionWrapper(mod.getMonoMass() + seq[i].getMonoWeight(Residue::Internal)) << "\"/>" << "\n";
						}
					}

					f << "      </modification_info>" << "\n";
									
				}

				f << " 			<analysis_result analysis=\"peptideprophet\">" << "\n";
				f << "			<peptideprophet_result probability=\"" << h.getScore()
					<< "\" all_ntt_prob=\"(" << h.getScore() << "," << h.getScore()
					<< "," << h.getScore() << ")\">" << "\n";
				f << "			</peptideprophet_result>" << "\n";
				f << "			</analysis_result>" << "\n";
				f << "			</search_hit>" << "\n";	
				f << "		</search_result>" << "\n";
				f << "		</spectrum_query>" << "\n";
			}
		}

		f << "	</msms_run_summary>" << "\n";
		f << "</msms_pipeline_analysis>" << "\n";
		
		f.close();
		
	}

  void PepXMLFile::matchModification_(DoubleReal mass, String& modification_description, const String& origin)
	{
		DoubleReal new_mass = mass - ResidueDB::getInstance()->getResidue(origin)->getMonoWeight(Residue::Internal);
		vector<String> mods;
		ModificationsDB::getInstance()->getModificationsByDiffMonoMass(mods, origin, new_mass, 0.001);
			
		if (mods.size() == 1)
    {
			modification_description = mods[0];
    }
    else
    {
     	if (mods.size() == 0)
      {
      	error(LOAD, String("Cannot find modification '") + String(mass) + " " + String(origin) + "'");
      }
      else
      {
       	String mod_str;
        for (vector<String>::const_iterator mit = mods.begin(); mit != mods.end(); ++mit)
        {
        	mod_str += " " + *mit;
        }
        error(LOAD, String("Found more than one modification ':") + mod_str +  "'");
      }
		}
	}


	void PepXMLFile::startElement(const XMLCh* const /*uri*/,
																const XMLCh* const /*local_name*/,
																const XMLCh* const qname,
																const xercesc::Attributes& attributes)
	{		
		String element = sm_.convert(qname);
		
		//cout << "Start: " << element << "\n";
	
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
			if (name == "fval")
			{
				// SpectraST score
				current_hit_->setScore(value);
				current_pep_->setScoreType(name);
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
//  							 << diff / max(rt, prec_rt) * 100 << "%)\n" << "\n";
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
			Size 			 modification_position = 0;
			String 		 temp_description = "";
			
			modification_position = attributeAsInt_(attributes, "position");
			modification_mass = attributeAsDouble_(attributes, "mass");
			
			matchModification_(modification_mass, temp_description, String(actual_sequence_[modification_position - 1]));
			
			// the modification position is 1-based
			actual_modifications_.push_back(make_pair(temp_description, modification_position));
		}

		else if ((element == "aminoacid_modification"))
		{
			AminoAcidModification aa_mod;
			String description;
			if (optionalAttributeAsString_(description, attributes, "description"))
			{
				aa_mod.description = description;
			}
			aa_mod.massdiff = attributeAsString_(attributes, "massdiff");
			aa_mod.aminoacid = attributeAsString_(attributes, "aminoacid");
			aa_mod.mass = attributeAsDouble_(attributes, "mass");
			String is_variable = attributeAsString_(attributes, "variable");
			if (is_variable == "Y")
			{
				if (description != "")
				{
					params_.variable_modifications.push_back(description); // TODO
				}
				else
				{
					String desc = aa_mod.aminoacid;
					if (aa_mod.massdiff.toDouble() > 0)
					{
						desc += "+" + String(aa_mod.massdiff);
					}
					else
					{
						desc += String(aa_mod.massdiff);
					}
					params_.variable_modifications.push_back(desc);
				}
			}
			else
			{
				fixed_modifications_.push_back(aa_mod);
				if (description != "")
				{
					params_.fixed_modifications.push_back(description); // TODO
				}
				else
				{
          String desc = aa_mod.aminoacid;
          if (aa_mod.massdiff.toDouble() > 0)
          {
            desc += "+" + String(aa_mod.massdiff);
          }
          else
          {
            desc += String(aa_mod.massdiff);
          }
					params_.fixed_modifications.push_back(desc);
        }

			}
		}
		else if (element == "terminal_modification")
		{
			// <terminal_modification terminus="n" massdiff="+108.05" mass="109.06" variable="N" protein_terminus="" description="dNIC (N-term)"/>
			AminoAcidModification aa_mod;
			String description;
			if (optionalAttributeAsString_(description, attributes, "description"))
			{
				aa_mod.description = description;
			}
			aa_mod.massdiff = attributeAsString_(attributes, "massdiff");
      aa_mod.aminoacid = attributeAsString_(attributes, "aminoacid");
      aa_mod.mass = attributeAsDouble_(attributes, "mass");
			aa_mod.terminus = attributeAsString_(attributes, "terminus");
      String is_variable = attributeAsString_(attributes, "variable");
      if (is_variable == "Y")
      {
      	if (description != "")
        {
          params_.variable_modifications.push_back(description); // TODO
        }
        else
        {
          String desc = aa_mod.aminoacid;
          if (aa_mod.massdiff.toDouble() > 0)
          {
            desc += "+" + String(aa_mod.massdiff);
          }
          else
          {
            desc += String(aa_mod.massdiff);
          }
          params_.variable_modifications.push_back(desc);
        }
      }
      else
      {
				if (description != "")
        {
          params_.fixed_modifications.push_back(description); // TODO
        }
        else
        {
          String desc = aa_mod.aminoacid;
          if (aa_mod.massdiff.toDouble() > 0)
          {
            desc += "+" + String(aa_mod.massdiff);
          }
          else
          {
            desc += String(aa_mod.massdiff);
          }
          params_.fixed_modifications.push_back(desc);
        }

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
			if (mass_type == "monoisotopic") 
			{
				params_.mass_type = ProteinIdentification::MONOISOTOPIC;
			}
			else 
			{
				if (mass_type == "average") 
				{
					params_.mass_type = ProteinIdentification::AVERAGE;
				}
				else 
				{
					error(LOAD,	"'fragment_mass_type' attribute of 'search_summary' tag should be 'monoisotopic' or 'average', not '" + mass_type + "'");
				}
			}
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
			for (vector<pair<String, Size> >::const_iterator it = actual_modifications_.begin(); it != actual_modifications_.end(); ++it)
			{
				// e.g. Carboxymethyl (C)
				vector<String> mod_split;
				it->first.split(' ', mod_split);
				if (mod_split.size() == 2)
				{
					if (mod_split[1] == "(C-term)" || ModificationsDB::getInstance()->getModification(it->first).getTermSpecificity() == ResidueModification::C_TERM)
					{
						temp_aa_sequence.setCTerminalModification(mod_split[0]);
					}
					else
					{
						if (mod_split[1] == "(N-term)" || ModificationsDB::getInstance()->getModification(it->first).getTermSpecificity() == ResidueModification::N_TERM)
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
			for (vector<AminoAcidModification>::const_iterator it = fixed_modifications_.begin(); it != fixed_modifications_.end(); ++it)
			{
					/*if (mod_split[1] == "(C-term)")
					{
						temp_aa_sequence.setCTerminalModification(mod_split[0]);
					}
					else
					{
						if (mod_split[1] == "(N-term)")
						{
							temp_aa_sequence.setNTerminalModification(mod_split[0]);
						}*/
					DoubleReal new_mass = it->mass - ResidueDB::getInstance()->getResidue(it->aminoacid)->getMonoWeight(Residue::Internal);
      		vector<String> mods;
      		ModificationsDB::getInstance()->getModificationsByDiffMonoMass(mods, it->aminoacid, new_mass, 0.001);
					if (mods.size() == 1)
					{
						for (Size i = 0; i != temp_aa_sequence.size(); ++i)
						{
							if (it->aminoacid.hasSubstring(temp_aa_sequence[i].getOneLetterCode()))
							{
								temp_aa_sequence.setModification(i, mods[0]);
							}
						}
					}
					else
					{
						if (mods.size() == 0)
						{
							error(LOAD, String("Cannot parse modification '") + it->aminoacid + " " +  + "'");
						}
						else
						{
							String mod_str;
							for (vector<String>::const_iterator mit = mods.begin(); mit != mods.end(); ++mit)
							{
								mod_str += " " + *mit;
							}
							error(LOAD, String("Found more than one modification ':") + mod_str +  "'");
						}
					}
				//}
			}

			current_hit_->setSequence(temp_aa_sequence);
			current_pep_->insertHit(*current_hit_);
			delete current_hit_;
			current_hit_ = 0;
			actual_modifications_.clear();						
		}
	}

} // namespace OpenMS
