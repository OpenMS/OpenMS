// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/XTandemInfile.h>
#include <OpenMS/SYSTEM/File.h>

#include <OpenMS/CHEMISTRY/ModificationsDB.h>

#include <set>
#include <fstream>

using namespace xercesc;
using namespace std;

namespace OpenMS 
{

	XTandemInfile::XTandemInfile()
		: Internal::XMLFile(),
			fragment_mass_tolerance_(0.3),
			precursor_mass_tolerance_plus_(2.0),
			precursor_mass_tolerance_minus_(2.0),
      precursor_mass_type_(XTandemInfile::MONOISOTOPIC),
      precursor_mass_error_unit_(XTandemInfile::DALTONS),
      fragment_mass_error_unit_(XTandemInfile::DALTONS),
			fragment_mass_type_(XTandemInfile::MONOISOTOPIC),
      max_precursor_charge_(3),
			precursor_lower_mz_(500.0),
			fragment_lower_mz_(150.0),
			number_of_threads_(1),
			modifications_(""),
      input_filename_(""),
			output_filename_(""),
      cleavage_site_("[RK]|{P}"),
			refine_(true),
			semi_cleavage_(true),
      refine_max_valid_evalue_(1000),
      number_of_missed_cleavages_(1),
			default_parameters_file_(""),
			max_valid_evalue_(1000)

	{
	  	
	}
	
	XTandemInfile::~XTandemInfile()
	{
	}
	
  void XTandemInfile::load(const String& filename)
  {
		Internal::XTandemInfileXMLHandler handler(filename, notes_, this);
		parse_(filename, &handler);
	}
		
	void XTandemInfile::write(const String& filename)
	{
		if (!File::writable(filename))
		{
			throw (Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename));
		}
		ofstream os(filename.c_str());
		writeTo_(os);
		return;
	}

	void XTandemInfile::writeTo_(ostream& os)
	{
		set<String> used_labels; // labels which are set by OpenMS not by the default parameters file

		os << "<?xml version=\"1.0\"?>" << "\n"
			 << "<?xml-stylesheet type=\"text/xsl\" href=\"tandem-input-style.xsl\"?>" << "\n"
			 << "<bioml>" << "\n";

		//////////////// list path parameters
		writeNote_(os, "input", "list path, default parameters", default_parameters_file_);
		used_labels.insert("list path, default parameters");
		writeNote_(os, "input", "list path, taxonomy information", taxonomy_file_);
		used_labels.insert("list path, taxonomy information");
		//<note type="input" label="spectrum, path">test_spectra.mgf</note>
		writeNote_(os, "input", "spectrum, path", input_filename_);
		used_labels.insert("spectrum, path");
		////////////////////////////////////////////////////////////////////////////////


		//////////////// spectrum parameters
		//<note type="input" label="spectrum, fragment monoisotopic mass error">0.4</note>
		writeNote_(os, "input", "spectrum, fragment monoisotopic mass error", String(fragment_mass_tolerance_));
		used_labels.insert("spectrum, fragment monoisotopic mass error");
    //<note type="input" label="spectrum, parent monoisotopic mass error plus">100</note>
    writeNote_(os, "input", "spectrum, parent monoisotopic mass error plus", String(precursor_mass_tolerance_plus_));
		used_labels.insert("spectrum, parent monoisotopic mass error plus");
		//<note type="input" label="spectrum, parent monoisotopic mass error minus">100</note>
    writeNote_(os, "input", "spectrum, parent monoisotopic mass error minus", String(precursor_mass_tolerance_minus_));
		used_labels.insert("spectrum, parent monoisotopic mass error minus");
		//<note type="input" label="spectrum, parent monoisotopic mass isotope error">yes</note>
		if (precursor_mass_type_ == XTandemInfile::MONOISOTOPIC)
		{
    	writeNote_(os , "input", "spectrum, parent monoisotopic mass isotope error", "yes");
		}
		else
		{
			writeNote_(os , "input", "spectrum, parent monoisotopic mass isotope error", "no");
		}
		used_labels.insert("spectrum, parent monoisotopic mass isotope error");
		//<note type="input" label="spectrum, fragment monoisotopic mass error units">Daltons</note>
		//<note>The value for this parameter may be 'Daltons' or 'ppm': all other values are ignored</note>
    if (fragment_mass_error_unit_ == XTandemInfile::DALTONS)
		{
			writeNote_(os, "input", "spectrum, fragment monoisotopic mass error units", "Daltons");
		}
		else
		{
			writeNote_(os, "input", "spectrum, fragment monoisotopic mass error units", "ppm");
		}
		used_labels.insert("spectrum, fragment monoisotopic mass error units");
    
		//<note type="input" label="spectrum, parent monoisotopic mass error units">ppm</note>
		//<note>The value for this parameter may be 'Daltons' or 'ppm': all other values are ignored</note>
    if (precursor_mass_error_unit_ == XTandemInfile::PPM)
		{
			writeNote_(os, "input", "spectrum, parent monoisotopic mass error units", "ppm");
		}
		else
		{
			writeNote_(os, "input", "spectrum, parent monoisotopic mass error units", "Daltons");
		}
		used_labels.insert("spectrum, parent monoisotopic mass error units");
    
		//<note type="input" label="spectrum, fragment mass type">monoisotopic</note>
		//<note>values are monoisotopic|average </note>
		if (fragment_mass_type_ == XTandemInfile::MONOISOTOPIC)
		{
			writeNote_(os, "input", "spectrum, fragment mass type", "monoisotopic");
		}
		else
		{
			writeNote_(os, "input", "spectrum, fragment mass type", "average");
		}
		used_labels.insert("spectrum, fragment mass type");
		////////////////////////////////////////////////////////////////////////////////

		
		//////////////// spectrum conditioning parameters
		//<note type="input" label="spectrum, dynamic range">100.0</note>
    //<note>The peaks read in are normalized so that the most intense peak
    //is set to the dynamic range value. All peaks with values of less that
    //1, using this normalization, are not used. This normalization has the
    //overall effect of setting a threshold value for peak intensities.</note>
		//writeNote_(os, "input", "spectrum, dynamic range", String(dynamic_range_));

  	//<note type="input" label="spectrum, total peaks">50</note>
    //<note>If this value is 0, it is ignored. If it is greater than zero (lets say 50),
    //then the number of peaks in the spectrum with be limited to the 50 most intense
    //peaks in the spectrum. X! tandem does not do any peak finding: it only
    //limits the peaks used by this parameter, and the dynamic range parameter.</note>
		//writeNote_(os, "input", "spectrum, total peaks", String(total_number_peaks_));

		//<note type="input" label="spectrum, maximum parent charge">4</note>
		writeNote_(os, "input", "spectrum, maximum parent charge", String(max_precursor_charge_));
		used_labels.insert("spectrum, maximum parent charge");

		// <note type="input" label="spectrum, use noise suppression">yes</note>
		//writeNote_(os, "input", "spectrum, use noise suppression", noise_supression_);

  	//<note type="input" label="spectrum, minimum parent m+h">500.0</note>
		//writeNote_(os, "input", "spectrum, minimum parent m+h", String(precursor_lower_mz_));

		//<note type="input" label="spectrum, minimum fragment mz">150.0</note>
		//writeNote_(os, "input", "spectrum, minimum fragment mz", String(fragment_lower_mz_));
		//used_labels.insert("spectrum, minimum fragment mz");

  	//<note type="input" label="spectrum, minimum peaks">15</note>
  	//writeNote_(os, "input", "spectrum, minimum peaks", String(min_number_peaks_));
		
		//<note type="input" label="spectrum, threads">1</note>
		writeNote_(os, "input", "spectrum, threads", String(number_of_threads_));
		used_labels.insert("spectrum, threads");
  
		//<note type="input" label="spectrum, sequence batch size">1000</note>
		//writeNote_(os, "input", "spectrum, sequence batch size", String(batch_size_));
		////////////////////////////////////////////////////////////////////////////////
		

		//////////////// residue modification parameters
		//<note type="input" label="residue, modification mass">57.022@C</note>
    //<note>The format of this parameter is m@X, where m is the modfication
    //mass in Daltons and X is the appropriate residue to modify. Lists of
    //modifications are separated by commas. For example, to modify M and C
    //with the addition of 16.0 Daltons, the parameter line would be
    //+16.0@M,+16.0@C
    //Positive and negative values are allowed.
    //</note>

		String fixed_mods;
		set<ModificationDefinition> fixed_mod_defs(modifications_.getFixedModifications());
		for (set<ModificationDefinition>::const_iterator it = fixed_mod_defs.begin(); it != fixed_mod_defs.end(); ++it)
		{
			double mod_mass(ModificationsDB::getInstance()->getModification(it->getModification()).getDiffMonoMass());
			String mod_string;
			if (mod_mass >= 0)
			{
				mod_string = "+" + String(mod_mass);
			}
			else
			{
				mod_string = "-" + String(mod_mass);
			}

			ResidueModification::Term_Specificity ts = ModificationsDB::getInstance()->getModification(it->getModification()).getTermSpecificity();
			//cerr << ModificationsDB::getInstance()->getModification(it->getModification()).getTermSpecificityName(ts) << " " << ModificationsDB::getInstance()->getModification(it->getModification()).getOrigin() << "\n";
			if (ts == ResidueModification::ANYWHERE)
			{
				mod_string += "@" + ModificationsDB::getInstance()->getModification(it->getModification()).getOrigin();
			}
			else
			{
				if (ts == ResidueModification::C_TERM)
				{
					mod_string += "@]";
				}
				else
				{
					mod_string += "@[";
				}
			}

			if (fixed_mods != "")
			{
				fixed_mods += "," + mod_string;
			}
			else
			{
				fixed_mods = mod_string;
			}
		}
			
		writeNote_(os, "input", "residue, modification mass", fixed_mods);
		used_labels.insert("residue, modification mass");

  	//<note type="input" label="residue, potential modification mass"></note>
    //<note>The format of this parameter is the same as the format
    //for residue, modification mass (see above).</note>

		String var_mods;		
    set<ModificationDefinition> var_mod_defs(modifications_.getVariableModifications());
    for (set<ModificationDefinition>::const_iterator it = var_mod_defs.begin(); it != var_mod_defs.end(); ++it)
    {
      double mod_mass(ModificationsDB::getInstance()->getModification(it->getModification()).getDiffMonoMass());
      String mod_string;
      if (mod_mass >= 0)
      {
        mod_string = "+" + String(mod_mass);
      }
      else
      {
        mod_string = "-" + String(mod_mass);
      }

      mod_string += "@" + ModificationsDB::getInstance()->getModification(it->getModification()).getOrigin();

      if (var_mods != "")
      {
        var_mods += "," + mod_string;
      }
      else
      {
        var_mods = mod_string;
      }
    }
		
		writeNote_(os, "input", "residue, potential modification mass", var_mods);
		used_labels.insert("residue, potential modification mass");
		
		writeNote_(os, "input", "protein, taxon", taxon_);
		used_labels.insert("protein, taxon");

		writeNote_(os, "input", "output, path", output_filename_);
		used_labels.insert("output, path");
/*
  	//<note type="input" label="residue, potential modification motif"></note>
    //<note>The format of this parameter is similar to residue, modification mass,
    //with the addition of a modified PROSITE notation sequence motif specification.
    //For example, a value of 80@[ST!]PX[KR] indicates a modification
    //of either S or T when followed by P, and residue and the a K or an R.
    //A value of 204@N!{P}[ST]{P} indicates a modification of N by 204, if it
    //is NOT followed by a P, then either an S or a T, NOT followed by a P.
    //Positive and negative values are allowed.
    //</note>
		writeNote_(os, "input", "residue, potential modification motif", variable_modification_motif_);
		////////////////////////////////////////////////////////////////////////////////


		//////////////// protein parameters
 		//<note type="input" label="protein, taxon">other mammals</note>
    //<note>This value is interpreted using the information in taxonomy.xml.</note>
		writeNote_(os, "input", "protein, taxon", taxon_);
		used_labels.insert("protein, taxon");

  	//<note type="input" label="protein, cleavage site">[RK]|{P}</note>
    //<note>this setting corresponds to the enzyme trypsin. The first characters
    //in brackets represent residues N-terminal to the bond - the '|' pipe -
    //and the second set of characters represent residues C-terminal to the
    //bond. The characters must be in square brackets (denoting that only
    //these residues are allowed for a cleavage) or french brackets (denoting
    //that these residues cannot be in that position). Use UPPERCASE characters.
    //To denote cleavage at any residue, use [X]|[X] and reset the
    //scoring, maximum missed cleavage site parameter (see below) to something like 50.
    //</note>
		writeNote_(os, "input", "protein, cleavage site", cleavage_site_);
		*/
				//////////////// semi cleavage parameter
  	//<note type="input" label="protein, cleavage semi">yes</note>
		writeNote_(os, "input", "protein, cleavage semi", semi_cleavage_);
    used_labels.insert("protein, cleavage semi");


  	//<note type="input" label="protein, modified residue mass file"></note>
  	//writeNote_(os, "input", "protein, modified residue mass file", modified_residue_mass_file_);

		//<note type="input" label="protein, cleavage C-terminal mass change">+17.002735</note>
  	//writeNote_(os, "input", "protein, cleavage C-terminal mass change", String(cleavage_c_term_mass_change_));

		//<note type="input" label="protein, cleavage N-terminal mass change">+1.007825</note>
  	//writeNote_(os, "input", "protein, cleavage N-terminal mass change", String(cleavage_n_term_mass_change_));

		//<note type="input" label="protein, N-terminal residue modification mass">0.0</note>
  	//writeNote_(os, "input", "protein, N-terminal residue modification mass", String(protein_n_term_mod_mass_));

		//<note type="input" label="protein, C-terminal residue modification mass">0.0</note>
  	//writeNote_(os, "input", "protein, C-terminal residue modification mass", String(protein_c_term_mod_mass_));

		//<note type="input" label="protein, homolog management">no</note>
    //<note>if yes, an upper limit is set on the number of homologues kept for a particular spectrum</note>
		//writeNote_(os, "input", "protein, homolog management", protein_homolog_management_);
		////////////////////////////////////////////////////////////////////////////////



		//////////////// model refinement parameters
  	//<note type="input" label="refine">yes</note>
		writeNote_(os, "input", "refine", refine_);
    used_labels.insert("refine");


/*
  	//<note type="input" label="refine, modification mass"></note>
		writeNote_(os, "input", "refine, modification mass", String(refine_mod_mass_));
  	//<note type="input" label="refine, sequence path"></note>
		writeNote_(os, "input", "refine, sequence path", refine_sequence_path_);
  	//<note type="input" label="refine, tic percent">20</note>
		writeNote_(os, "input", "refine, tic percent", String(refine_tic_percent_));
  	//<note type="input" label="refine, spectrum synthesis">yes</note>
		writeNote_(os, "input", "refine, spectrum synthesis", refine_spectrum_sythesis_);
  	//<note type="input" label="refine, maximum valid expectation value">0.1</note>
		writeNote_(os, "input", "refine, maximum valid expectation value", String(refine_max_valid_evalue_));
  	//<note type="input" label="refine, potential N-terminus modifications">+42.010565@[</note>
		writeNote_(os, "input", "refine, potential N-terminus modifications", refine_variable_n_term_mods_);
  	//<note type="input" label="refine, potential C-terminus modifications"></note>
		writeNote_(os, "input", "refine, potential C-terminus modifications", refine_variable_c_term_mods_);
  	//<note type="input" label="refine, unanticipated cleavage">yes</note>
		writeNote_(os, "input", "refine, unanticipated cleavage", refine_unanticipated_cleavage_);
  	//<note type="input" label="refine, potential modification mass"></note>
		writeNote_(os, "input", "refine, potential modification mass", String(variable_mod_mass_));
  	//<note type="input" label="refine, point mutations">no</note>
		writeNote_(os, "input", "refine, point mutations", refine_point_mutations_);
  	//<note type="input" label="refine, use potential modifications for full refinement">no</note>
		writeNote_(os, "input", "refine, use potential modifications for full refinement", use_var_mod_for_full_refinement_);*/
  	//<note type="input" label="refine, potential modification motif"></note>
  	//<note>The format of this parameter is similar to residue, modification mass,
    //with the addition of a modified PROSITE notation sequence motif specification.
    //For example, a value of 80@[ST!]PX[KR] indicates a modification
    //of either S or T when followed by P, and residue and the a K or an R.
    //A value of 204@N!{P}[ST]{P} indicates a modification of N by 204, if it
    //is NOT followed by a P, then either an S or a T, NOT followed by a P.
    //Positive and negative values are allowed.
    //</note>

		//writeNote_(os, "input", "refine, potential modification motif", refine_var_mod_motif_);
		////////////////////////////////////////////////////////////////////////////////


		//////////////// scoring parameters
 		//<note type="input" label="scoring, minimum ion count">4</note>
		//writeNote_(os, "input", "scoring, minimum ion count", String(scoring_min_ion_count_));
 		//<note type="input" label="scoring, maximum missed cleavage sites">1</note>
		writeNote_(os, "input", "scoring, maximum missed cleavage sites", String(number_of_missed_cleavages_));
		used_labels.insert("scoring, maximum missed cleavage sites");
  	//<note type="input" label="scoring, x ions">no</note>
		//writeNote_(os, "input", "scoring, x ions", score_x_ions_);
  	//<note type="input" label="scoring, y ions">yes</note>
		//writeNote_(os, "input", "scoring, y ions", score_y_ions_);
  	//<note type="input" label="scoring, z ions">no</note>
    //writeNote_(os, "input", "scoring, z ions", score_z_ions_);
  	//<note type="input" label="scoring, a ions">no</note>
    //writeNote_(os, "input", "scoring, a ions", score_a_ions_);
  	//<note type="input" label="scoring, b ions">yes</note>
    //writeNote_(os, "input", "scoring, b ions", score_b_ions_);
  	//<note type="input" label="scoring, c ions">no</note>
    //writeNote_(os, "input", "scoring, c ions", score_c_ions_);
  	//<note type="input" label="scoring, cyclic permutation">no</note>
    //<note>if yes, cyclic peptide sequence permutation is used to pad the scoring histograms</note>
		//writeNote_(os, "input", "scoring, cyclic permutation", scoring_cyclic_permutation_);
  	//<note type="input" label="scoring, include reverse">no</note>
    //<note>if yes, then reversed sequences are searched at the same time as forward sequences</note>
		//writeNote_(os, "input", "scoring, include reverse", scoring_include_reverse_);
		////////////////////////////////////////////////////////////////////////////////


		//////////////// output parameters
  	//<note type="input" label="output, log path"></note>
  	//<note type="input" label="output, message">...</note>
		//writeNote_(os, "input", "output, message", String("..."));
  	//<note type="input" label="output, one sequence copy">no</note>
  	//<note type="input" label="output, sequence path"></note>
  	//<note type="input" label="output, path">output.xml</note>
		//writeNote_(os, "input", "output, path", output_filename_);
  	//<note type="input" label="output, sort results by">protein</note>
		writeNote_(os, "input", "output, sort results by", "spectrum");
		used_labels.insert("output, sort results by");
    //<note>values = protein|spectrum (spectrum is the default)</note>
  	//<note type="input" label="output, path hashing">yes</note>
    //<note>values = yes|no</note>
  	//<note type="input" label="output, xsl path">tandem-style.xsl</note>
  	//<note type="input" label="output, parameters">yes</note>
    //<note>values = yes|no</note>
  	//<note type="input" label="output, performance">yes</note>
    //<note>values = yes|no</note>
  	//<note type="input" label="output, spectra">yes</note>
    //<note>values = yes|no</note>
  	//<note type="input" label="output, histograms">yes</note>
    //<note>values = yes|no</note>
  	//<note type="input" label="output, proteins">yes</note>
    //<note>values = yes|no</note>
  	//<note type="input" label="output, sequences">yes</note>
    //<note>values = yes|no</note>
 		//<note type="input" label="output, one sequence copy">no</note>
    //<note>values = yes|no, set to yes to produce only one copy of each protein sequence in the output xml</note>
  	//<note type="input" label="output, results">valid</note>
		writeNote_(os, "input", "output, results", "all");
		used_labels.insert("output, results");
    //<note>values = all|valid|stochastic</note>
  	//<note type="input" label="output, maximum valid expectation value">0.1</note>
		writeNote_(os, "input", "output, maximum valid expectation value", String(max_valid_evalue_)); 
		used_labels.insert("output, maximum valid expectation value");
    //<note>value is used in the valid|stochastic setting of output, results</note>
  	//<note type="input" label="output, histogram column width">30</note>
    //<note>values any integer greater than 0. Setting this to '1' makes cutting and pasting histograms
    //into spread sheet programs easier.</note>
		//<note type="description">ADDITIONAL EXPLANATIONS</note>
  	//<note type="description">Each one of the parameters for X! tandem is entered as a labeled note
    //  node. In the current version of X!, keep those note nodes
    //  on a single line.
  	//</note>
  	//<note type="description">The presence of the type 'input' is necessary if a note is to be considered
    //  an input parameter.
  	//</note>
  	//<note type="description">Any of the parameters that are paths to files may require alteration for a
    //  particular installation. Full path names usually cause the least trouble,
    //  but there is no reason not to use relative path names, if that is the
    //  most convenient.
  	//</note>
  	//<note type="description">Any parameter values set in the 'list path, default parameters' file are
    //  reset by entries in the normal input file, if they are present. Otherwise,
    //  the default set is used.
  	//</note>
  	//<note type="description">The 'list path, taxonomy information' file must exist.
    //</note>
  	//<note type="description">The directory containing the 'output, path' file must exist: it will not be created.
    //</note>
  	//<note type="description">The 'output, xsl path' is optional: it is only of use if a good XSLT style sheet exists.
    //</note>
		////////////////////////////////////////////////////////////////////////////////

		// those of the parameters that are not set by this file adapter 
		// are just written from the default XTandem infile
		for (vector<Internal::XTandemInfileNote>::const_iterator it = notes_.begin(); it != notes_.end(); ++it)
		{
			if (it->note_type != "" && it->note_label != "" && used_labels.find(it->note_label) == used_labels.end())
			{
				writeNote_(os, it->note_type, it->note_label, it->note_value);
			}
		}

		os << "</bioml>" << "\n";

	}

	void XTandemInfile::writeNote_(ostream& os, const String& type, const String& label, const String& value)
	{
		os << "\t<note type=\"" << type << "\" label=\"" << label  << "\">" << value << "</note>" << "\n";
	}

	void XTandemInfile::writeNote_(ostream& os, const String& type, const String& label, const char* value)
	{
		String val(value);
		os << "\t<note type=\"" << type << "\" label=\"" << label  << "\">" << val << "</note>" << "\n";
	}
	
	void XTandemInfile::writeNote_(ostream& os, const String& type, const String& label, bool value)
	{
		if (value)
		{
			os << "\t<note type=\"" << type << "\" label=\"" << label  << "\">yes</note>" << "\n";
		}
		else
		{
			os << "\t<note type=\"" << type << "\" label=\"" << label  << "\">no</note>" << "\n";
		}
	}

	void XTandemInfile::setOutputFilename(const String& filename)
	{
		output_filename_ = filename;
	}

	const String& XTandemInfile::getOutputFilename() const
	{
		return output_filename_;
	}

	void XTandemInfile::setInputFilename(const String& filename)
	{
		input_filename_ = filename;
	}

	const String& XTandemInfile::getInputFilename() const
	{
		return input_filename_;
	}

	void XTandemInfile::setTaxonomyFilename(const String& filename)
	{
		taxonomy_file_ = filename;
	}

	const String& XTandemInfile::getTaxonomyFilename() const
	{
		return taxonomy_file_;
	}

	void XTandemInfile::setDefaultParametersFilename(const String& filename)
	{
		default_parameters_file_ = filename;
	}

	const String& XTandemInfile::getDefaultParametersFilename() const
	{
		return default_parameters_file_;
	}

	void XTandemInfile::setModifications(const ModificationDefinitionsSet& mods)
	{
		modifications_ = mods;
	}

	const ModificationDefinitionsSet& XTandemInfile::getModifications() const
	{
		return modifications_;
	}

	void XTandemInfile::setTaxon(const String& taxon)
	{
		taxon_ = taxon;
	}

	const String& XTandemInfile::getTaxon() const
	{
		return taxon_;
	}

	void XTandemInfile::setPrecursorMassTolerancePlus(double tolerance)
	{
		precursor_mass_tolerance_plus_ = tolerance;
	}

	double XTandemInfile::getPrecursorMassTolerancePlus() const
	{
		return precursor_mass_tolerance_plus_;
	}
	
	void XTandemInfile::setPrecursorMassToleranceMinus(double tolerance)
	{
		precursor_mass_tolerance_minus_ = tolerance;
	}

	double XTandemInfile::getPrecursorMassToleranceMinus() const
	{
		return precursor_mass_tolerance_minus_;
	}
	
	void XTandemInfile::setPrecursorMassErrorUnit(ErrorUnit unit)
	{
		precursor_mass_error_unit_ = unit;
	}

	XTandemInfile::ErrorUnit XTandemInfile::getPrecursorMassErrorUnit() const
	{
		return precursor_mass_error_unit_;
	}

	void XTandemInfile::setFragmentMassErrorUnit(ErrorUnit unit)
	{
		fragment_mass_error_unit_ = unit;
	}

	XTandemInfile::ErrorUnit XTandemInfile::getFragmentMassErrorUnit() const
	{
		return fragment_mass_error_unit_;
	}

	void XTandemInfile::setMaxPrecursorCharge(Int max_charge)
	{
		max_precursor_charge_ = max_charge;
	}

	Int XTandemInfile::getMaxPrecursorCharge() const
	{
		return max_precursor_charge_;
	}

	void XTandemInfile::setFragmentMassTolerance(double tolerance)
	{
		fragment_mass_tolerance_ = tolerance;
	}

	double XTandemInfile::getFragmentMassTolerance() const
	{
		return fragment_mass_tolerance_;
	}

	void XTandemInfile::setNumberOfThreads(UInt num_threads)
	{
		number_of_threads_ = num_threads;
	}

	UInt XTandemInfile::getNumberOfThreads() const
	{
		return number_of_threads_;
	}

	XTandemInfile::MassType XTandemInfile::getPrecursorErrorType() const
	{
		return precursor_mass_type_;
	}

	void XTandemInfile::setPrecursorErrorType(const MassType mass_type)
	{
		precursor_mass_type_ = mass_type;
	}

	void XTandemInfile::setMaxValidEValue(double value)
	{
		max_valid_evalue_ = value;
	}

	double XTandemInfile::getMaxValidEValue() const
	{
		return max_valid_evalue_;
	}

	void XTandemInfile::setNumberOfMissedCleavages(UInt missed_cleavages)
	{
		number_of_missed_cleavages_ = missed_cleavages;
	}

	UInt XTandemInfile::getNumberOfMissedCleavages() const
	{
		return number_of_missed_cleavages_;
	}
        
  bool XTandemInfile::isRefining() const
  {
    return refine_;
  }

  void XTandemInfile::setRefine(const bool refine)
  {
    refine_ = refine;
  }
  void XTandemInfile::setSemiCleavage(const bool semi_cleavage)
  {
    semi_cleavage_ = semi_cleavage;
  }
} // namespace OpenMS
