// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/PepXMLFile.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/SYSTEM/File.h>
#include <fstream>
#include <iostream>

using namespace std;

namespace OpenMS 
{

	PepXMLFile::PepXMLFile()
		: XMLHandler("","1.12"),
			XMLFile("/SCHEMAS/pepXML_v114.xsd","1.14"),
			proteins_(0), 
			peptides_(0),
			experiment_(0),
			scan_map_(),
			rt_tol_(10.0),
  		mz_tol_(10.0)
	{
		const ElementDB* db = ElementDB::getInstance();
		hydrogen_ = *db->getElement("Hydrogen");
	}


	PepXMLFile::~PepXMLFile()
	{
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
    if ( !protein_ids.empty() )
		{
			if (protein_ids.size() > 1)
			{
				warning(STORE, "More than one protein identification defined; only first one is written into pepXML, more are not supported.");
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

      // compute mass of modified residue
      EmpiricalFormula ef = ResidueDB::getInstance()->getResidue(mod.getOrigin())->getFormula(Residue::Internal);
      ef += mod.getDiffFormula();

      f << "      "
			  << "<aminoacid_modification aminoacid=\"" << mod.getOrigin()
        << "\" massdiff=\"" << precisionWrapper(mod.getDiffMonoMass()) << "\" mass=\""
        << precisionWrapper(ef.getMonoWeight())
				<< "\" variable=\"Y\" binary=\"N\" description=\"" << *it << "\"/>"
				<< "\n";
		}

		for (set<String>::const_iterator it = n_term_mods.begin(); it != n_term_mods.end(); ++it)
		{
			const ResidueModification& mod = ModificationsDB::getInstance()->
				getModification(*it);
      f << "      "
			  << "<terminal_modification terminus=\"n\" massdiff=\""
				<< precisionWrapper(mod.getDiffMonoMass()) << "\" mass=\"" << precisionWrapper(mod.getMonoMass())
				<< "\" variable=\"Y\" description=\"" << *it
				<< "\" protein_terminus=\"\"/>" << "\n";
		}

		for (set<String>::const_iterator it = c_term_mods.begin(); it != c_term_mods.end(); ++it)
		{
			const ResidueModification& mod = ModificationsDB::getInstance()->getModification(*it);
      f << "      "
        << "<terminal_modification terminus=\"c\" massdiff=\""
				<< precisionWrapper(mod.getDiffMonoMass()) << "\" mass=\"" << precisionWrapper(mod.getMonoMass())
				<< "\" variable=\"Y\" description=\"" << *it
				<< "\" protein_terminus=\"\"/>" << "\n";
		}

		f << "    </search_summary>" << "\n";
		f << "    <analysis_timestamp analysis=\"peptideprophet\" time=\"2007-12-05T17:49:52\" id=\"1\"/>" << "\n";
		

		Int count(1);
		for (vector<PeptideIdentification>::const_iterator it = peptide_ids.begin();
				 it != peptide_ids.end(); ++it, ++count)
		{
			if (it->getHits().size() > 0)
			{
        if (it->getHits().size() > 1)
        {
          LOG_WARN << "PepXMLFile::store() : only writing the first peptide hit of " << it->getHits().size() << " for PeptideID# " << count << "\n";
        }
				PeptideHit h = *it->getHits().begin();
				AASequence seq = h.getSequence();
				DoubleReal precursor_neutral_mass = seq.getMonoWeight();

        Int scan_index;
        if(it->metaValueExists("RT_index"))
        {
          scan_index = it->getMetaValue("RT_index");
        }
        else
        {
          scan_index = count;
        }

				f << "		<spectrum_query spectrum=\"" << count << "\""
          << " start_scan=\"" << scan_index << "\""
          << " end_scan=\"" << scan_index << "\""
					<< " precursor_neutral_mass=\"" << precisionWrapper(precursor_neutral_mass) << "\""
          << " assumed_charge=\"" << h.getCharge() << "\" index=\"" << count << "\"";

        DataValue dv = it->getMetaValue("RT");
        if (dv!=DataValue::EMPTY)
        {
          f << " retention_time_sec=\"" << dv << "\" ";
        }

        f << ">\n";
				f << " 		<search_result>" << "\n";
				f << "			<search_hit hit_rank=\"1\" peptide=\""
					<< seq.toUnmodifiedString() << "\" peptide_prev_aa=\""
					<< h.getAABefore() << "\" peptide_next_aa=\"" << h.getAAAfter()
					<< "\" protein=\"Protein1\" num_tot_proteins=\"1\" num_matched_ions=\"0\" tot_num_ions=\"0\" calc_neutral_pep_mass=\"" << precisionWrapper(precursor_neutral_mass)
          << "\" massdiff=\"0.0\" num_tol_term=\"0\" num_missed_cleavages=\"0\" is_rejected=\"0\" protein_descr=\"Protein No. 1\">" << "\n";
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
              // the modification position is 1-based
              f << "         <mod_aminoacid_mass position=\"" << (i+1)
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


  void PepXMLFile::matchModification_(const DoubleReal mass, const String& origin, String& modification_description)
	{
		DoubleReal mod_mass = mass - ResidueDB::getInstance()->getResidue(origin)->getMonoWeight(Residue::Internal);
		vector<String> mods;
		ModificationsDB::getInstance()->getModificationsByDiffMonoMass(mods, origin, mod_mass, 0.001);
			
		if (mods.size() == 1)
    {
			modification_description = mods[0];
    }
    else
    {
      if ( !mods.empty() )
      {
       	String mod_str = mods[0];
        for (vector<String>::const_iterator mit = ++mods.begin(); mit != mods.end(); ++mit)
        {
        	mod_str += ", " + *mit;
        }
				error(LOAD, "Modification '" + String(mass) + "' is not uniquely defined by the given data. Using '" + mods[0] +  "' to represent any of '" + mod_str + "'!");
				modification_description = mods[0];
      }
		}
	}
	

	void PepXMLFile::makeScanMap_()
	{
		scan_map_.clear();
		Size scan = 0;
		for (MSExperiment<>::ConstIterator e_it = experiment_->begin(); e_it != experiment_->end(); ++e_it, ++scan)
		{
			String id = e_it->getNativeID();
			bool failed = false;
			try
			{
				// expected format: "spectrum=#" (mzData) or "scan=#" (mzXML)
				Int num_id = id.suffix('=').toInt();
				if (num_id >= 0)
				{
					scan_map_.insert(scan_map_.end(), pair<Size, Size>(num_id, scan));
				}
				else
				{
					failed = true;
				}
			}
			catch (Exception::ConversionError)
			{
				failed = true;
			}
			if (failed)
			{
				scan_map_.clear();
				error(LOAD, "Could not construct mapping of native scan numbers to indexes");
			}			
		}
	}


	void PepXMLFile::readRTMZCharge_(const xercesc::Attributes& attributes)
	{
		DoubleReal mass = attributeAsDouble_(attributes, "precursor_neutral_mass");
		charge_ = attributeAsInt_(attributes, "assumed_charge");
		mz_ = (mass + hydrogen_mass_ * charge_) / charge_;
		rt_ = 0;

		bool rt_present = optionalAttributeAsDouble_(rt_, attributes, "retention_time_sec");

		if (!rt_present || use_precursor_data_) // get RT from experiment
		{
			if (!experiment_)
			{
				error(LOAD, "Cannot get precursor information - no experiment given");
				return;
			}

			// assume only one scan, i.e. ignore "end_scan":
			Size scan = attributeAsInt_(attributes, "start_scan");
			if (!scan_map_.empty()) scan = scan_map_[scan];
			MSSpectrum<> spec = (*experiment_)[scan];
			bool success = false;
			if (spec.getMSLevel() == 2)
			{
				if (!use_precursor_data_)
				{
					rt_ = spec.getRT();
					success = true;
				}
				else if (!rt_present || Math::approximatelyEqual(spec.getRT(), rt_, 0.001))
				{ 
					DoubleReal prec_mz = 0, prec_rt = 0;
					vector<Precursor> precursors = spec.getPrecursors();
          if ( !precursors.empty() )
          {
            prec_mz = precursors[0].getMZ(); // assume only one precursor
          }
					MSExperiment<>::ConstIterator it = experiment_->getPrecursorSpectrum(experiment_->begin() + scan);
					if (it != experiment_->end())
					{
						prec_rt = it->getRT();
					}

					// check if "rt"/"mz" are similar to "prec_rt"/"prec_mz"
					// (otherwise, precursor mapping is wrong)
					if ((prec_mz > 0) && Math::approximatelyEqual(prec_mz, mz_, mz_tol_)	&& (prec_rt > 0) && (!rt_present || Math::approximatelyEqual(prec_rt, rt_, rt_tol_)))
					{
						// DoubleReal diff;
						// diff = mz_ - prec_mz;
						// cout << "m/z difference: " << diff << " ("
						// 		 << diff / max(mz_, prec_mz) * 100 << "%)\n";
						// diff = rt_ - prec_rt;
						// cout << "RT difference: " << diff << " ("
						// 		 << diff / max(rt_, prec_rt) * 100 << "%)\n" << "\n";
						mz_ = prec_mz;
						rt_ = prec_rt;
						success = true;
					}
				}
			}
			if (!success)
			{
				error(LOAD, "Cannot get precursor information - scan mapping is incorrect");
			}
		}
	}


	void PepXMLFile::load(const String& filename, vector<ProteinIdentification>& 
												proteins, vector<PeptideIdentification>& peptides, 
												const String& experiment_name)
	{
		MSExperiment<> exp;
		load(filename, proteins, peptides, experiment_name, exp, false);
	}

	
  void PepXMLFile::load(const String& filename, vector<ProteinIdentification>& proteins, vector<PeptideIdentification>& peptides, const String& experiment_name, const MSExperiment<>& experiment, bool use_precursor_data)
  { 
  	// initialize, "load" could be called several times
    exp_name_ = "";
    experiment_ = 0;
		use_precursor_data_ = use_precursor_data;
    prot_id_ = "";
    charge_ = 0;
  	rt_tol_ = 10.0;
  	mz_tol_ = 10.0;
  	peptides.clear();
  	peptides_ = &peptides;
		proteins.clear();
		proteins_ = &proteins;
		// assume mass type "average" (in case element "search_summary" is missing)
		hydrogen_mass_ = hydrogen_.getAverageWeight();
  	
  	file_ = filename;	// filename for error messages in XMLHandler

		if (experiment_name != "")
		{
			exp_name_ = File::removeExtension(experiment_name);
			
			if (!experiment.empty()) // use experiment only if we know the name
			{
				experiment_ = &experiment;
				MSExperiment<>::AreaType area = experiment_->getDataRange();
				// set tolerance to 1% of data range (if above a sensible minimum):
				rt_tol_ = max((area.maxX() - area.minX()) * 0.01, rt_tol_);
				mz_tol_ = max((area.maxY() - area.minY()) * 0.01, mz_tol_);
				makeScanMap_();
			}
		}

		wrong_experiment_ = false;
		seen_experiment_ = exp_name_.empty(); // without experiment name, don't care
		parse_(filename, this);
		
		if (!seen_experiment_)
		{
			fatalError(LOAD, "Found no experiment with name '" + experiment_name + "'");
		}

		// clean up duplicate ProteinHits in ProteinIdentifications:
		// (can't use "sort" and "unique" because no "op<" defined for PeptideHit)
		for (vector<ProteinIdentification>::iterator prot_it = proteins.begin();
				 prot_it != proteins.end(); ++prot_it)
		{
			set<String> accessions;
			// modeled after "remove_if" in STL header "algorithm":
			vector<ProteinHit>::iterator first = prot_it->getHits().begin(), result = first;
			for (; first != prot_it->getHits().end(); ++first)
			{
				String accession = first->getAccession();
				bool new_element = accessions.insert(accession).second;
				if (new_element) // don't remove
				{
					*result++ = *first;
				}
			}
			prot_it->getHits().erase(result, first);
		}
		
    // reset members
		exp_name_.clear();
		prot_id_.clear();
		date_.clear();
		proteins_ = 0;
		peptides_ = 0;
		experiment_ = 0;
		scan_map_.clear();
  }


	void PepXMLFile::startElement(const XMLCh* const /*uri*/,
																const XMLCh* const /*local_name*/,
																const XMLCh* const qname,
																const xercesc::Attributes& attributes)
	{		
		String element = sm_.convert(qname);
		
		// cout << "Start: " << element << "\n";
	
		if (element == "msms_run_summary") // parent: "msms_pipeline_analysis"
		{
			if (!exp_name_.empty())
			{
				String base_name = attributeAsString_(attributes, "base_name");
				wrong_experiment_ = !base_name.hasSuffix(exp_name_);
			}
			if (wrong_experiment_) return;
			seen_experiment_ = true;

			// create a ProteinIdentification in case "search_summary" is missing:
			ProteinIdentification protein;
			protein.setDateTime(date_);
			prot_id_ = "unknown_" + date_.getDate();
			// "prot_id_" will be overwritten if elem. "search_summary" is present
			protein.setIdentifier(prot_id_);
			proteins_->push_back(protein);
			current_proteins_.clear();
			current_proteins_.push_back(--proteins_->end());
			hydrogen_mass_ = hydrogen_.getAverageWeight();
		}
		
		else if (wrong_experiment_)
		{
			// do nothing here (this case exists to prevent parsing of elements for
			// experiments we're not interested in)
		}
		
		// now, elements occurring more frequently are generally closer to the top

		else if (element == "search_score") // parent: "search_hit"
		{
			String name = attributeAsString_(attributes, "name");
			DoubleReal value;

			// TODO: deal with different scores
			if (name == "expect")
			{ // X!Tandem or Mascot E-value
				value = attributeAsDouble_(attributes, "value");
				peptide_hit_.setScore(value);
				current_peptide_.setScoreType(name);
				current_peptide_.setHigherScoreBetter(false);
			}
			// if (name == "hyperscore")
			// { // X!Tandem score
			// 	value = attributeAsDouble_(attributes, "value");
			// 	peptide_hit_.setScore(value);
			// 	current_peptide_.setScoreType(name); // add "X!Tandem" to name?
			// 	current_peptide_.setHigherScoreBetter(true);
			// }
			else if (name == "xcorr")
			{ // Sequest score
				value = attributeAsDouble_(attributes, "value");
				peptide_hit_.setScore(value);
				current_peptide_.setScoreType(name); // add "Sequest" to name?
				current_peptide_.setHigherScoreBetter(true);
			}
			else if (name == "fval")
			{ // SpectraST score
				value = attributeAsDouble_(attributes, "value");
				peptide_hit_.setScore(value);
				current_peptide_.setScoreType(name);
				current_peptide_.setHigherScoreBetter(true);
			}
		}

		else if (element == "search_hit") // parent: "search_result"
		{ // creates a new PeptideHit
			current_sequence_ = attributeAsString_(attributes, "peptide");
			current_modifications_.clear();
			peptide_hit_ = PeptideHit();
			peptide_hit_.setRank(attributeAsInt_(attributes, "hit_rank"));
			peptide_hit_.setCharge(charge_); // from parent "spectrum_query" tag
			String prev_aa, next_aa;
			if (optionalAttributeAsString_(prev_aa, attributes, "peptide_prev_aa")) peptide_hit_.setAABefore(prev_aa[0]);
			if (optionalAttributeAsString_(next_aa, attributes, "peptide_next_aa")) peptide_hit_.setAAAfter(next_aa[0]);
			String protein = attributeAsString_(attributes, "protein");
			peptide_hit_.addProteinAccession(protein);
			ProteinHit hit;
			hit.setAccession(protein);
			current_proteins_[search_id_ - 1]->insertHit(hit);
		}

		else if (element == "search_result") // parent: "spectrum_query"
		{ // creates a new PeptideIdentification
			current_peptide_ = PeptideIdentification();
			current_peptide_.setMetaValue("RT", rt_);
			current_peptide_.setMetaValue("MZ", mz_);
			search_id_ = 1; // references "search_summary"
			optionalAttributeAsUInt_(search_id_, attributes, "search_id");
			current_peptide_.setIdentifier(current_proteins_[search_id_ - 1]->getIdentifier());
		}

		else if (element == "spectrum_query") // parent: "msms_run_summary"
		{
			readRTMZCharge_(attributes); // sets "rt_", "mz_", "charge_"
		}
	
		else if (element == "peptideprophet_result") // parent: "analysis_result" (in "search_hit")
		{
			// PeptideProphet probability overwrites original search score
			// maybe TODO: deal with meta data associated with PeptideProphet search
			if (current_peptide_.getScoreType() != "InterProphet probability")
			{
				DoubleReal value = attributeAsDouble_(attributes, "probability");
				peptide_hit_.setScore(value);
				current_peptide_.setScoreType("PeptideProphet probability");
				current_peptide_.setHigherScoreBetter(true);
			}
		}

		else if (element == "interprophet_result") // parent: "analysis_result" (in "search_hit")
		{
			// InterProphet probability overwrites PeptideProphet probability and
			// original search score
			DoubleReal value = attributeAsDouble_(attributes, "probability");
			peptide_hit_.setScore(value);
			current_peptide_.setScoreType("InterProphet probability");
			current_peptide_.setHigherScoreBetter(true);
		}

		else if (element == "alternative_protein") // parent: "search_hit"
		{
			String protein = attributeAsString_(attributes, "protein");
			peptide_hit_.addProteinAccession(protein);
			ProteinHit hit;
			hit.setAccession(protein);
			current_proteins_[search_id_ - 1]->insertHit(hit);
		}

		else if (element == "mod_aminoacid_mass") // parent: "modification_info" (in "search_hit")
		{
			DoubleReal modification_mass = attributeAsDouble_(attributes, "mass");
      Size 			 modification_position = attributeAsInt_(attributes, "position");
      String     origin = String(current_sequence_[modification_position - 1]);
			String 		 temp_description = "";
			
			matchModification_(modification_mass, origin, temp_description);
			
      if (temp_description.size()>0)
      {
			  // the modification position is 1-based
			  current_modifications_.push_back(make_pair(temp_description, modification_position));
      }
      else
      {
      	error(LOAD, String("Cannot find modification '") + String(modification_mass) + " " + String(origin) + "' @" + String(modification_position));
      }
		}

		else if (element == "aminoacid_modification") // parent: "search_summary"
		{
			AminoAcidModification aa_mod;
			optionalAttributeAsString_(aa_mod.description, attributes, "description");
			aa_mod.massdiff = attributeAsString_(attributes, "massdiff");
			aa_mod.aminoacid = attributeAsString_(attributes, "aminoacid");
			aa_mod.mass = attributeAsDouble_(attributes, "mass");
			String is_variable = attributeAsString_(attributes, "variable");
			if (is_variable == "Y")
			{
				if (aa_mod.description != "")
				{
					params_.variable_modifications.push_back(aa_mod.description); // TODO
				}
				else
				{
					String desc = aa_mod.aminoacid;
					if (aa_mod.massdiff.toDouble() >= 0)
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
				if (aa_mod.description != "")
				{
					params_.fixed_modifications.push_back(aa_mod.description); // TODO
				}
				else
				{
          String desc = aa_mod.aminoacid;
          if (aa_mod.massdiff.toDouble() >= 0)
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

		else if (element == "terminal_modification") // parent: "search_summary"
		{
			// <terminal_modification terminus="n" massdiff="+108.05" mass="109.06" variable="N" protein_terminus="" description="dNIC (N-term)"/>
			AminoAcidModification aa_mod;
			optionalAttributeAsString_(aa_mod.description, attributes, "description");
			aa_mod.massdiff = attributeAsString_(attributes, "massdiff");
			optionalAttributeAsString_(aa_mod.aminoacid, attributes, "aminoacid");
      aa_mod.mass = attributeAsDouble_(attributes, "mass");
			aa_mod.terminus = attributeAsString_(attributes, "terminus");
      String is_variable = attributeAsString_(attributes, "variable");
      if (is_variable == "Y")
      {
      	if (aa_mod.description != "")
        {
          params_.variable_modifications.push_back(aa_mod.description); // TODO
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
				if (aa_mod.description != "")
        {
          params_.fixed_modifications.push_back(aa_mod.description); // TODO
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

		else if (element == "search_summary") // parent: "msms_run_summary"
		{ // creates a new ProteinIdentification (usually)
			fixed_modifications_.clear();
			variable_modifications_.clear();
			params_ = ProteinIdentification::SearchParameters();
			params_.enzyme = enzyme_;
			String mass_type = attributeAsString_(attributes, "precursor_mass_type");
			if (mass_type == "monoisotopic")
			{
				hydrogen_mass_ = hydrogen_.getMonoWeight();
			}
			else
			{
				hydrogen_mass_ = hydrogen_.getAverageWeight();
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
			// generate identifier from search engine and date:
			prot_id_ = search_engine + "_" + date_.getDate();

			search_id_ = 1;
			optionalAttributeAsUInt_(search_id_, attributes, "search_id");
			vector<ProteinIdentification>::iterator prot_it;
			if (search_id_ == 1)
			{ // ProteinIdent. was already created for "msms_run_summary" -> add to it
				prot_it = current_proteins_.front();
			}
			else
			{ // create a new ProteinIdentification
				ProteinIdentification protein;
				protein.setDateTime(date_);
				proteins_->push_back(protein);
				prot_it = --proteins_->end();
				prot_id_ = prot_id_ + "_" + search_id_; // make sure the ID is unique
			}
			prot_it->setSearchEngine(search_engine);
			prot_it->setIdentifier(prot_id_);
		}

		else if (element == "sample_enzyme") // parent: "msms_run_summary"
		{	// special case: search parameter that occurs *before* "search_summary"!
			String name = attributeAsString_(attributes, "name");
			name.toLower();
			// spelling of enzyme names in pepXML?
			if (name.hasPrefix("trypsin")) enzyme_ = ProteinIdentification::TRYPSIN;
			else if (name.hasPrefix("pepsin")) enzyme_ = ProteinIdentification::PEPSIN_A;
			else if (name.hasPrefix("protease")) enzyme_ = ProteinIdentification::PROTEASE_K;
			else if (name.hasPrefix("chymotrypsin")) enzyme_ = ProteinIdentification::CHYMOTRYPSIN;
			else enzyme_ = ProteinIdentification::UNKNOWN_ENZYME;

			ProteinIdentification::SearchParameters params = current_proteins_.front()->getSearchParameters();
			params.enzyme = enzyme_;
			current_proteins_.front()->setSearchParameters(params);
		}

		else if (element == "search_database") // parent: "search_summary"
		{
			params_.db = attributeAsString_(attributes, "local_path");
		}
		
		else if (element == "msms_pipeline_analysis") // root
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
		}
	}

	
	void PepXMLFile::endElement(const XMLCh* const /*uri*/,
															const XMLCh* const /*local_name*/,
															const XMLCh* const qname)
	{
		String element = sm_.convert(qname);

		// cout << "End: " << element << "\n";

		if (wrong_experiment_)
		{
			// do nothing here (skip all elements that belong to the wrong experiment)
		}

		else if (element == "search_hit")
		{
			AASequence temp_aa_sequence = AASequence(current_sequence_);
			
      // modification position is 1-based
			for (vector<pair<String, Size> >::const_iterator it = current_modifications_.begin(); it != current_modifications_.end(); ++it)
			{
				// e.g. Carboxymethyl (C)
				vector<String> mod_split;
				it->first.split(' ', mod_split);
				if (it->first.hasSubstring("C-term"))
				{
					temp_aa_sequence.setCTerminalModification(it->first);
				}
				else if (it->first.hasSubstring("N-term"))
				{
					temp_aa_sequence.setNTerminalModification(it->first);
				}
				else if (mod_split.size() == 2)
				{
          temp_aa_sequence.setModification(it->second - 1, mod_split[0]);
				}
				else
				{
					error(LOAD, String("Cannot parse modification '") + it->first + "@" + it->second + "'");
				}
			}

			// fixed modifications
			for (vector<AminoAcidModification>::const_iterator it = fixed_modifications_.begin(); it != fixed_modifications_.end(); ++it)
			{
				const Residue* residue = ResidueDB::getInstance()->getResidue(it->aminoacid);
				if (residue == 0)
				{
					error(LOAD, String("Cannot parse modification of unknown amino acid '") + it->aminoacid + "'");
				}
				else
				{
					DoubleReal new_mass = it->mass - residue->getMonoWeight(Residue::Internal);
      		vector<String> mods;
      		ModificationsDB::getInstance()->getModificationsByDiffMonoMass(mods, it->aminoacid, new_mass, 0.001);
					if (mods.size() == 1)
					{
						for (Size i = 0; i < temp_aa_sequence.size(); ++i)
						{
							if (it->aminoacid.hasSubstring(temp_aa_sequence[i].getOneLetterCode()))
							{
								temp_aa_sequence.setModification(i, mods[0]);
							}
						}
					}
					else if (mods.empty())
					{
						error(LOAD, String("Cannot parse modification of amino acid '") + it->aminoacid + "'");
					}
					else
					{
						String mod_str = mods[0];
						for (vector<String>::const_iterator mit = ++mods.begin(); mit != mods.end(); ++mit)
						{
							mod_str += ", " + *mit;
						}
						error(LOAD, "Modification '" + String(it->mass) + "' is not uniquely defined by the given data. Using '" + mods[0] +  "' to represent any of '" + mod_str + "'!");
						for (Size i = 0; i < temp_aa_sequence.size(); ++i)
            {
              if (it->aminoacid.hasSubstring(temp_aa_sequence[i].getOneLetterCode()))
              {
                temp_aa_sequence.setModification(i, mods[0]);
              }
            }
					}
				}
				//}
			}

			peptide_hit_.setSequence(temp_aa_sequence);
			current_peptide_.insertHit(peptide_hit_);
		}

		else if (element == "search_result")
		{
			peptides_->push_back(current_peptide_);
		}

		else if (element == "search_summary")
		{
			current_proteins_.back()->setSearchParameters(params_);
		}
	}

} // namespace OpenMS
