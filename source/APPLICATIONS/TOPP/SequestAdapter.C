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
// $Maintainer: Martin Langwisch $
// --------------------------------------------------------------------------


#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/PTMXMLFile.h>
#include <OpenMS/FORMAT/SequestInfile.h>
#include <OpenMS/FORMAT/SequestOutfile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/ContactPerson.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FileHandler.h>

#include <stdlib.h>
#include <vector>
#include <algorithm>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------


/**
	@page SequestAdapter SequestAdapter

	@brief Identifies peptides in MS/MS spectra via Sequest.

	This wrapper application serves for getting peptide peptide_identifications
	for MS/MS spectra. The wrapper can be executed in three different
	modes:
	<ol>
				<li>
				The whole process of ProteinIdentification via Sequest is executed.
				Inputfile is one (or more) mz file containing the MS/MS spectra
				for which peptide_identifications are to be found
				and one ore two databases in FASTA format containing
				the possible proteins.
				The results are written as an IdXML output file. This mode is selected
			 	by default.
				Note: You need a user with network access on the computer hosting sequest.
			 	</li>

				<li>
				Only the first part of the ProteinIdentification process is performed.
				This means that a Sequest input file is generated and dta files are
				created from the mz file.
				Calling an Sequest process should look like the following:

				@code sequest -P\<inputfilename\> \<path to dta files\>*.dta  @endcode

				Consult your Sequest reference manual for further details.

				This mode is selected by the <b>-sequest_in</b> option in the command line.
				</li>

				<li>
				Only the second part of the ProteinIdentification process is performed.
				This means that the output of sequest is translated into IdXML.

				This mode is selected by the <b>-sequest_out</b> option in the command line.
				</li>
	</ol>
*/

// We do not want this class to show up in the docu -> cond
// @cond

class TOPPSequestAdapter
	: public TOPPBase
{
	public:
		TOPPSequestAdapter()
			: TOPPBase("SequestAdapter", "annotates MS/MS spectra using Sequest.")
		{
		}

	protected:
		static const Int max_peptide_mass_units = 2;
		static const UInt max_dtas_per_run = 1000; // sequest has a problem when there are too many dtas, so they have to be splitted, 1000 seemed to work very good
		PointerSizeUInt dtas;

		void registerOptionsAndFlags_()
		{
			addText_("The definitions for the parameters are taken from the site:\n"
										 "http://www.grosse-coosmann.de/~florian/Parameters.html#file.");
			registerStringOption_("out", "<file>", "", "output file in IdXML format.\n"
			                                           "Note: In mode 'sequest_in' a Sequest input file is written.", false);
			registerStringOption_("in", "<file>", "", "input file(s) in mzXML or mzData format (comma-separated).\n"
					 																			"Note: In mode 'sequest_out' a directory with Sequest results files\n"
																								"(*.out) is read", false);
			registerFlag_("sequest_in", "if this flag is set the SequestAdapter will read in mzXML or mzData\n"
																								"and write an Sequest input file\n"
																								"and create dta files from the given mzXML or mzData files");
			registerFlag_("sequest_out", "if this flag is set the SequestAdapter will read in Sequest result files\n"
																									"and write IdXML");
			registerStringOption_("mzFiles", "<file>", "", "when using sequest_out the mzXML or mzData files (comma-separated)\n"
																																						"have to be given to retrieve the retention times", false);
			registerFlag_("show_enzymes", "show a list with enzymes and corresponding numbers to choose from");
			registerStringOption_("sequest_computer", "<name>", "", "the name of the computer in the network that hosts Sequest\n"
																															"(rdesktop is used to connect to this computer)", false);
			registerStringOption_("sequest_directory_win", "<dir>", "", "the windows directory in which Sequest (sequest.exe) is located", false);
			registerStringOption_("user", "<name>", "", "user name for the sequest computer (has to have access to network!)", false);
			registerStringOption_("password", "<pw>", "", "password for this user (if not given, you have to enter it at promt)", false);
			registerStringOption_("temp_data_directory", "<dir>", "", "a directory in which some temporary files can be stored", false);
			registerStringOption_("temp_data_directory_win", "<dir>", "", "windows path of the temporary data directory,\n"
																																																			"e.g. X:\\temp_data_dir", false);
			registerStringOption_("db", "<file>", "", "name of FASTA-database to search in", false);
			registerStringOption_("sequest_input", "<file>", "", "name for the input file of Sequest (may only be used in a full run)", false);
			addEmptyLine_();
			addText_("For each directory, one corresponding network drive has to be given");
			registerStringOption_("temp_data_directory_network", "<path>", "", "network path of the temporary data directory,\n"
																																																			"e.g. \\\\computername\\username\\temp_data_dir", false);
			registerStringOption_("db_directory_network", "<path>", "", "network path of the database directory", false);
			registerStringOption_("sequest_input_directory_network", "<path>", "", "network path of the sequest input file directory", false);
			addEmptyLine_();
			registerDoubleOption_("precursor_mass_tolerance", "<tol>", 2.0 , "the precursor mass tolerance", false);
			registerDoubleOption_("peak_mass_tolerance", "<tol>", 1.0, "the peak mass tolerance", false);
			registerDoubleOption_("p_value", "<prob>", 1.0, "annotations with inferior p-value are ignored", false);
		  registerStringOption_("charges", "[1>3,5]", "", "comma-seperated list of charge states (or ranges)", false);
			registerIntOption_("num_results", "<num>", 1, "the maximal number of results (peptides) to show (per scan/dta)", false);
			registerStringOption_("cleavage", "<enz>", "Trypsin", "the number of the enzyme used for digestion", false);
			registerStringOption_("enzyme_info", "<>", "", "information about the enzyme used\n"
																																							"<name>,<cut direction: N to C?>,<cuts after>,<doesn't cut before>\n"
																																							"cuts after, doesn't cut before: amino acids in 1-letter code\n"
																																							"or '-' for unspecific cleavage", false);
			registerFlag_("list_modifications", "show a list of the available modifications");
			registerStringOption_("modifications", "<mods>", "", "the colon-seperated modifications; may be\n"
																																													 "<name>,<type>, e.g.: Deamidation,opt or\n"
																																													 "<composition>,<residues>,<type>,<name>, e.g.: H(2).C(2).O,KCS,opt,Acetyl or\n"
																																													 "<mass>,<residues>,<type>,<name>, e.g.: 42.0367,KCS,opt,Acetyl or\n"
																																													 "Valid values for \"type\" are \"fix\", \"cterminal\", \"nterminal\",\n"
																																													 "and \"opt\" (the default).\n", false);
			registerFlag_("use_monoisotopic_mod_mass", "use monoisotopic masses for the modifications");
			registerStringOption_("modifications_xml_file", "<file>", "", "name of an XML file with the modifications", false);
			registerIntOption_("max_num_dif_AA_per_mod", "<num>", 0, "limits the maximum total number of\n"
																																																			 "variable modifications per amino acid", false);
			registerIntOption_("max_num_dif_mods_per_peptide", "<num>", 0, "limits the maximum total number of\n"
																																																							 "each single variable modification in one peptide", false);
			registerDoubleOption_("match_peak_tol", "", 0, "the minimal space between two peaks", false);
			registerStringOption_("neutral_loss_ABY", "[ABY]", "011", "ABY: 0 or 1 whether neutral losses of the series should be honored,\n"
																																												 "e.g.: 011", false);
			registerStringOption_("ion_series_weights", "[abcdvwxyz]", "0,1.0,0,0,0,0,0,1.0,0", "[0.0, 1.0] factor for the series,\n"
																																																			"e.g.: 0,0.5,0,0,0,0,0,1.0,0", false);
			registerDoubleOption_("ion_cutoff", "<num>", 0.0, "This value selects a cut-off below which a matching peptide is rejected.\n"
																											 "The value has to be in [0,1] and is compared with the ratio\n"
																											 "(# matching theoretical fragment peaks)/(# total theoretical fragment peaks)\n"
																											 "which means that one select a minimum coverage of matching peaks.", false);
			registerIntOption_("pep_mass_unit", "<num>", 0, "peptide mass unit: 0=amu (atomic mass unit), 1=mmu (millimass unit),\n"
																																									 "2=ppm (parts per million)", false);
			registerDoubleOption_("prot_mass", "<num>", 0, "protein mass or minimum protein mass (see below)", false);
			registerDoubleOption_("max_prot_mass_or_tol", "<num>", 0, "maximum protein mass or tolerance", false);
			registerIntOption_("max_num_int_cleav_sites", "<num>", 0, "This value is the number of cleavage positions\n"
																																																	 "that may have been ignored by the enzyme.", false);
			registerIntOption_("match_peak_count", "<num>", 0, "The highest abundant experimental peaks are checked\n"
																																													"whether they are matched by the theoretical ones.\n"
																																													"match_peak_count is the number of the top abundant peaks to check.\n"
																												"A maximum of match_peak_allowed_error may lack this test.\n", false);
			registerIntOption_("match_peak_allowed_error", "<num>", 0, "see match_peak_count", false);
			registerFlag_("show_fragment_ions", "If set the fragment peaks of the top scored peptide are listed\n"
																																"at the end of the output");
// 			registerFlag_("use_phospho_fragmentation", "?");
			registerFlag_("remove_precursor_peak", "If set the peaks near (15 amu) the precursor are removed.");
			registerFlag_("mass_type_precursor", "Set selects monoisotopic masses, not set selects average masses\n"
																																 "for calculating precursor peaks.");
			registerFlag_("mass_type_peak", "Set selects monoisotopic masses, not set selects average masses\n"
																												 "for calculating peaks.");
			registerFlag_("normalize_xcorr", "Whether to use normalized xcorr values in the out files.");
			registerFlag_("residues_in_lower_case", "Whether the residues in the FASTA database are in lower case.");
			registerStringOption_("partial_sequence", "<sequences>", "", "A comma delimited list of amino acid sequences that must occur\n"
																																																		 "in the theoretical spectra.", false);
			registerStringOption_("header_filter", "<sequences>", "", "Several elements can be splitted by commas.\n"
																																															 "Each element can be introduced by an exclamation mark (!)\n"
																																															 "meaning that this element must not appear in the header of\n"
																																															"a protein or the protein will be skipped. This test is done first.\n"
																																															"Next, all other elements are tested. The protein is processed\n"
																																															"if one filter string matches the header string.\n"
																																															"A tilde (~) in the filter string is replaced by a blank during comparison.", false);
			registerFlag_("keep_out_files", "If set the Seuest .out-files are not removed");
			registerFlag_("keep_dta_files", "If set the dta-files that were created from the mzXML or mzData files are not removed");
			registerIntOption_("nuc_reading_frame", "<num>", 0, "Format of the FASTA database:\n"
																													 "0  The FASTA file contains amino acid codes. No translation is needed.\n"
																													 "1  The DNA sequence is scanned left to right (forward direction).\n"
																													 "The amino acid code starts with the first DNA code.\n"
																													 "2  The DNA sequence is scanned left to right (forward direction).\n"
																													 "The amino acid code starts with the second DNA code.\n"
																													 "3  The DNA sequence is scanned left to right (forward direction).\n"
																													 "The amino acid code starts with the third DNA code.\n"
																													 "4  The DNA sequence is scanned right to left (backward direction\n"
																													 "for the complementary strand).\n"
																													 "The amino acid code starts with the first DNA code.\n"
																													 "5  The DNA sequence is scanned right to left (backward direction\n"
																													 "for the complementary strand).\n"
																													 "The amino acid code starts with the second DNA code.\n"
																													 "6  The DNA sequence is scanned right to left (backward direction\n"
																													 "for the complementary strand).\n"
																													 "The amino acid code starts with the third DNA code.\n"
																													 "7  Use each of the DNA translations of the codes 1, 2, 3.\n"
																													 "8  Use each of the DNA translations of the codes 4, 5, 6.\n"
																													 "9  Use each of the DNA translations of the codes 1, 2, 3, 4, 5, 6.\n", false);
			registerStringOption_("contact_name", "<name>", "unknown", "Name of the contact", false);
			registerStringOption_("contact_institution", "<name>", "unknown", "Name of the contact institution", false);
			registerStringOption_("contact_info", "<info>", "unknown", "Some information about the contact", false);
		}

		String get_composition_elements(const String& composition, vector< vector< String > >& iso_sym_occ, char seperator = ' ')
		{
			iso_sym_occ.clear();
			vector< String > substrings;
			composition.split(seperator, substrings); // get the single elements of the composition: e.g. 18O(-1) or C(3) or N
			if ( substrings.empty() ) substrings.push_back(composition);
			String::size_type pos, pos2;
			String isotope, symbol, occurences;
			// for each element, get the isotope (if used), the symbol and the occurences
			for ( vector< String >::const_iterator e_i = substrings.begin(); e_i != substrings.end(); ++e_i )
			{
				isotope.clear();
				occurences = "1";
				pos = 0;
				while ( (bool) isdigit((*e_i)[pos]) ) ++pos; // if an isotope is used, find it
				isotope = e_i->substr(0, pos);
				if ( isotope.empty() ) isotope = "0";
				pos2 = e_i->find('(', pos);
				if ( pos2 != String::npos ) // if the element occurs more than once, a bracket is found
				{
					symbol = e_i->substr(pos, pos2++ - pos);
					occurences = e_i->substr(pos2, e_i->length() - pos2 - 1 );
				}
				else
				{
					symbol = e_i->substr(pos).toLower().firstToUpper();
					occurences = "1";
				}
				// check whether this really is a chemical symbol (only characters, max length 2)
				if ( symbol.length() > 2 || (isalpha(symbol[0]) == 0) || (isalpha(symbol[symbol.length() - 1]) == 0) ) return (composition);
				// then check whether isotope and occurences are numbers
				Int i_iso, i_occ;
				try
				{
					i_iso = isotope.toInt();
					i_occ = occurences.toInt();
				}
				catch( Exception::ConversionError ce )
				{
					return composition;
				}
				if ( String(i_iso) != isotope || String(i_occ) != occurences )
				{
					return composition;
				}

				// if this is a composition, insert its elements into the vector
				iso_sym_occ.push_back(vector< String >());
				iso_sym_occ.back().push_back(isotope);
				iso_sym_occ.back().push_back(symbol);
				iso_sym_occ.back().push_back(occurences);
			}

			return String();
		}

		bool isWinFormat(const string& name)
		{
			// check for the directory and the backslash afterwards
			if ( name.length() > 1 )
			{
				if ( name[1] == ':' )
				{
					if ( name.length() > 3 )
					{
						if ( name[2] == '\\' )
						{
							// make sure there's no space within the name, as in windows 'cmd /C "command"' is used, so there's no possibility to use any more ""
							if ( name.find(" ") == String::npos ) return true;
							else return false;
						}
						else return false;
					}
					else return true;
				}
			}
			return false;
		}

		bool correctNetworkPath(String& network_path, UInt backslashes = 2)
		{
			String::size_type pos = 0;
			while ( (pos < network_path.length()) && (network_path[pos] == '\\') ) ++pos;
			if ( pos < backslashes ) network_path.insert(network_path.begin(), backslashes-pos, '\\');
			else network_path.erase(0, pos-backslashes);
			if ( network_path.length() < backslashes+1 ) return false;
			if ( network_path[network_path.length() - 1] != '\\' ) network_path.append("\\"); // if it doesn't end with a slash, append one
			return true;
		}

  UInt
		MSExperiment2DTAs(
			MSExperiment<>& msexperiment,
			const String& common_name,
			const vector< Int >& charges,
			map< String, Real >& filenames_and_precursor_retention_times,
			bool make_dtas = true)
		throw (Exception::UnableToCreateFile)
		{
			DTAFile dtafile;
			String filename;
			UInt scan_number = 0;
			UInt msms_spectra = 0;

			for ( MSExperiment<>::Iterator spec_i = msexperiment.begin(); spec_i != msexperiment.end(); ++spec_i )
			{
				++scan_number;
				if ( (spec_i->getMSLevel() == 2) && (!spec_i->empty()) )
				{
					++msms_spectra;
					if ( spec_i->getPrecursorPeak().getCharge() )
					{
						filename = common_name + "." + String(scan_number) + "." + String(spec_i->getPrecursorPeak().getCharge()) + ".dta_" + String( (PointerSizeUInt) (dtas / max_dtas_per_run) );
						if ( make_dtas )
						{
							++dtas;
							dtafile.store(filename, *spec_i);
						}
// 						filename.replace(filename.length() - 4, 4, ".out");
// 						filenames_and_precursor_retention_times[File::basename(filename)] = spec_i->getRT();
						filenames_and_precursor_retention_times[filename] = spec_i->getRT();
					}
					else
					{
						for ( vector< Int >::const_iterator i = charges.begin(); i != charges.end(); ++i )
						{
							filename = common_name + "." + String(scan_number) + "." + *i + ".dta_" + String( (PointerSizeUInt) (dtas / max_dtas_per_run) );
							if ( make_dtas )
							{
								++dtas;
								spec_i->getPrecursorPeak().setCharge(*i);
								dtafile.store(filename, *spec_i);
							}
// 							filename.replace(filename.length() - 4, 4, ".out");
// 							filenames_and_precursor_retention_times[File::basename(filename)] = spec_i->getRT();
							filenames_and_precursor_retention_times[filename] = spec_i->getRT();
						}
						spec_i->getPrecursorPeak().setCharge(0);
					}
				}
			}

			return msms_spectra;
		}

		ExitCodes main_(int , char**)
		{
			//-------------------------------------------------------------
			// (1) variables
			//-------------------------------------------------------------

			// (1.0) variables for running the program
			SequestInfile sequest_infile;
			SequestOutfile sequest_outfile;

			String
				logfile,
				output_filename,
				input_filename,
				input_file_directory_network,
				user,
				password,
				sequest_computer,
				temp_data_directory,
				temp_data_directory_win,
				temp_data_directory_network,
				sequest_directory_win,
				database,
				database_directory_network,
				out_directory,
				batch_filename,
				string_buffer,
				string_buffer2,
				modifications_filename;

			ContactPerson contact_person;

			bool
				sequest_in,
				sequest_out,
				keep_out_files,
				keep_dta_files,
				monoisotopic;

			vector< String >
				substrings,
				substrings2,
				spectra;

			vector< Int > charges;

			char char_buffer;

   Real
				Real_buffer,
				Real_buffer2;

			Int int_buffer;

			Real p_value;

			// the dta-names and their retention_times
			map< String, Real > filenames_and_precursor_retention_times;

			// filename and tag: file has to: 1 - exist  2 - be readable  4 - writable  8 - be deleted afterwards
			vector< pair< String, UInt > > files;

			//-------------------------------------------------------------
			// (2) parsing and checking parameters
			//-------------------------------------------------------------

			modifications_filename = getStringOption_("modifications_xml_file");

			if ( getFlag_("list_modifications") )
			{
				if ( modifications_filename.empty() )
				{
					writeLog_("No modifications XML file given. Aborting!");
					return INPUT_FILE_NOT_FOUND;
				}
				if ( !File::readable(modifications_filename) )
				{
					writeLog_("Modifications XML file is not readable. Aborting!");
					return INPUT_FILE_NOT_READABLE;
				}
				map< String, pair< String, String > > ptm_informations;
				try
				{
					PTMXMLFile().load(modifications_filename, ptm_informations);
				}
				catch ( Exception::ParseError pe )
				{
					writeLog_(pe.getMessage());
					return PARSE_ERROR;
				}

				// output the information
				stringstream ptm_info;
				String::size_type max_name_length, max_composition_length, max_amino_acids_length;
				max_name_length = max_composition_length = max_amino_acids_length = 0;
				for ( map< String, pair< String, String > >::const_iterator mod_i = ptm_informations.begin(); mod_i != ptm_informations.end(); ++mod_i )
				{
					max_name_length = max(max_name_length, mod_i->first.length());
					max_composition_length = max(max_composition_length, mod_i->second.first.length());
					max_amino_acids_length = max(max_amino_acids_length, mod_i->second.second.length());
				}
				ptm_info << "These modifications are taken from unimod" << endl;
				ptm_info << "name" << String(max_name_length - 4, ' ') << "\t" << "composition" << String(max_composition_length - 11, ' ') << "\t" << "amino_acids" << String(max_amino_acids_length - 11, ' ') << endl;
				for ( map< String, pair< String, String > >::const_iterator mod_i = ptm_informations.begin(); mod_i != ptm_informations.end(); ++mod_i )
				{
					ptm_info << mod_i->first << String(max_name_length - mod_i->first.length(), ' ') << "\t" << mod_i->second.first << String(max_composition_length - mod_i->second.first.length(), ' ') << "\t" << mod_i->second.second << String(max_amino_acids_length - mod_i->second.second.length(), ' ') << endl;
				}

				return EXECUTION_OK;
			}

			// only show the available enzymes, then quit
			if ( getFlag_("show_enzymes") )
			{
				writeLog_("Option show_enzymes chosen.");
				writeLog_(sequest_infile.getEnzymeInfoAsString());
				return EXECUTION_OK;
			}
			// (2.0) variables for running the program
			sequest_in = getFlag_("sequest_in");
			sequest_out = getFlag_("sequest_out");

			// a 'normal' sequest run corresponds to both sequest_in and sequest_out set
			if ( !sequest_in && !sequest_out ) sequest_in = sequest_out = true;

			logfile = getStringOption_("log");
			if ( logfile.empty() )
			{
				logfile = "temp.sequest.log";
				files.push_back(make_pair(logfile, 4+8));
			}
			files.push_back(make_pair(logfile, 4));

			string_buffer = getStringOption_("charges");
			if ( string_buffer.empty() )
			{
				writeLog_("No charge states given. Aborting!");
				return ILLEGAL_PARAMETERS;
			}
			else
			{
				Int range_start, range_end;
				string_buffer.split(',', substrings);
				if ( substrings.empty() ) substrings.push_back(string_buffer);

				for ( vector< String >::iterator s_i = substrings.begin(); s_i != substrings.end(); )
				{
					if ( s_i->empty() ) substrings.erase(s_i);
					else
					{
						s_i->split('>', substrings2);
						if ( substrings2.size() < 2 ) // only one number, no range
						{
							if ( (*s_i)[s_i->length()-1] == '-' ) charges.push_back(-1 * s_i->toInt());
							else charges.push_back(s_i->toInt());
						}
						else // range of charge states
						{
							if ( substrings2.size() > 2 )
							{
								writeLog_("Illegal range of charge states given: " + *s_i + ". Aborting!");
								return ILLEGAL_PARAMETERS;
							}

							if ( substrings2[0][substrings2[0].length()-1] == '-' ) range_start = -1 * substrings2[0].toInt();
							else range_start = substrings[0].toInt();

							if ( substrings2[1][substrings2[1].length()-1] == '-' ) range_end = -1 * substrings2[1].toInt();
							else range_end = substrings2[1].toInt();

							for ( Int i = min(range_start, range_end); i <= max(range_start, range_end); ++i )
							{
								if ( i ) charges.push_back(i);
							}
						}

						++s_i;
					}
				}

				if ( charges.empty() )
				{
					writeLog_("No charges states given. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				sort(charges.begin(), charges.end());
				for ( vector< Int >::iterator i = charges.begin(); i != --charges.end(); )
				{
					if ( (*i) == (*(i+1)) ) charges.erase(i+1);
					else ++i;
				}
			}

			temp_data_directory = getStringOption_("temp_data_directory");
			if ( temp_data_directory.empty() )
			{
				writeLog_("No directory for temporary files given. Aborting!");
				return ILLEGAL_PARAMETERS;
			}
			File::absolutePath(temp_data_directory);
			temp_data_directory.ensureLastChar('/');

			string_buffer = getStringOption_("in");
			if ( string_buffer.empty() )
			{
				writeLog_("No input file specified. Aborting!");
				return ILLEGAL_PARAMETERS;
			}
			else
			{
				if ( sequest_in ) // if sequest_in is set, in are the spectra
				{
					string_buffer.split(',', spectra);
					if ( spectra.empty() ) spectra.push_back(string_buffer);
					out_directory = temp_data_directory;
				}
				else // if only sequest_out is set, in is the out_directory
				{
					out_directory = string_buffer;
					File::absolutePath(out_directory);
					out_directory.ensureLastChar('/');

					// if only sequest_out is set, the mz files have to be given to retrieve the retention times
					string_buffer = getStringOption_("mzFiles");
					if ( string_buffer.empty() )
					{
						writeLog_("No mz files specified. Aborting!");
						return ILLEGAL_PARAMETERS;
					}
					else
					{
						string_buffer.split(',', spectra);
						if ( spectra.empty() ) spectra.push_back(string_buffer);
					}
				}
			}
			
			keep_out_files = getFlag_("keep_out_files");
			if ( sequest_out && !sequest_in ) keep_out_files = true;
			
			keep_dta_files = getFlag_("keep_dta_files");
			if ( sequest_in && !sequest_out ) keep_dta_files = true;
			
			contact_person.setName(getStringOption_("contact_name"));
			contact_person.setInstitution(getStringOption_("contact_institution"));
			contact_person.setContactInfo(getStringOption_("contact_info"));
			
			if ( sequest_in )
			{
				temp_data_directory_win = getStringOption_("temp_data_directory_win");
				temp_data_directory_win.ensureLastChar('\\');
				
				if ( !isWinFormat(temp_data_directory_win) )
				{
					writeLog_("Windows path for the directory for temporary files has wrong format: " + temp_data_directory_win + ". borting!");
					return ILLEGAL_PARAMETERS;
				}
				temp_data_directory_network = getStringOption_("temp_data_directory_network");
				if ( temp_data_directory_network.empty() )
				{
					writeLog_("No network path for the directory for temporary files given. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				if ( !correctNetworkPath(temp_data_directory_network) )
				{
					writeLog_(temp_data_directory_network + "is no network path. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				
				database = getStringOption_("db");
				if ( database.empty() )
				{
					writeLog_("No database specified. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				files.push_back(make_pair(database, 2));
				
				if ( !sequest_out )
				{
					input_filename = getStringOption_("out");
					if ( input_filename.empty() )
					{
						writeLog_("No output file specified. Aborting!");
						return ILLEGAL_PARAMETERS;
					}
					
					input_file_directory_network = getStringOption_("sequest_input_directory_network");
					if ( input_file_directory_network.empty() )
					{
						writeLog_("No network path for the directory of the Sequest input file given. Aborting!");
						return ILLEGAL_PARAMETERS;
					}
					if ( !correctNetworkPath(input_file_directory_network) )
					{
						writeLog_(input_file_directory_network + "is no network path. Aborting!");
						return ILLEGAL_PARAMETERS;
					}
				}
				else
				{
					input_filename = getStringOption_("sequest_input");
					if ( input_filename.empty() )
					{
						input_filename = temp_data_directory + "temp.sequest.in";
						files.push_back(make_pair(input_filename, 4+8));
						input_file_directory_network = temp_data_directory_network;
					}
					else
					{
						input_file_directory_network = getStringOption_("sequest_input_directory_network");
						if ( input_file_directory_network.empty() )
						{
							writeLog_("No network path for the directory of the Sequest input file given. Aborting!");
							return ILLEGAL_PARAMETERS;
						}
						files.push_back(make_pair(input_filename, 2));
					}
					if ( !correctNetworkPath(input_file_directory_network) )
					{
						writeLog_(input_file_directory_network + "is no network path. Aborting!");
						return ILLEGAL_PARAMETERS;
					}
				}
			}
			
			if ( sequest_in && sequest_out )
			{
				user = getStringOption_("user");
				
				password = getStringOption_("password");
				
				sequest_directory_win = getStringOption_("sequest_directory_win");
				if ( !sequest_directory_win.hasSuffix("sequest.exe") ) sequest_directory_win.ensureLastChar('\\');
				if ( !isWinFormat(sequest_directory_win) )
				{
					writeLog_("Windows path for the SEQUEST working directory has wrong format: " + sequest_directory_win + ". Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				else if ( sequest_directory_win.empty() )
				{
					writeLog_("No windows path for the SEQUEST working directory given. Assuming PATH variable to be set accordingly!");
					sequest_directory_win = "sequest";
				}
				
				sequest_computer = getStringOption_("sequest_computer");
				if ( sequest_computer.empty() )
				{
					writeLog_("No sequest computer name given. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
			}
			
			if ( logfile == temp_data_directory + "sequest.log")
			{
				writeLog_("The logfile must not be named " + temp_data_directory + "sequest.log. Aborting!");
				return ILLEGAL_PARAMETERS;
			}
			
			//batch_filename = getStringOption_("batchfile");
			if ( batch_filename.empty() )
			{
				batch_filename = "sequest_run.bat";
				files.push_back(make_pair(temp_data_directory + batch_filename, 4+8));
			}
			else if ( !batch_filename.hasSuffix(".bat") ) batch_filename.append(".bat");
			
			if ( sequest_in )
			{
				database_directory_network = getStringOption_("db_directory_network");
				if ( !correctNetworkPath(database_directory_network) )
				{
					writeLog_(database_directory_network + "is no network path. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				string_buffer = File::basename(database);
				if ( !database_directory_network.hasSuffix(string_buffer) ) database_directory_network.append(string_buffer);
				sequest_infile.setDatabase(database_directory_network);
				
				Real_buffer = getDoubleOption_("precursor_mass_tolerance");
				if ( Real_buffer == -1 )
				{
					writeLog_("No precursor mass tolerance specified. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				else if ( Real_buffer < 0 )
				{
					writeLog_("Precursor mass tolerance < 0. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setPrecursorMassTolerance(Real_buffer);
				
				Real_buffer = getDoubleOption_("peak_mass_tolerance");
				if ( Real_buffer == -1 )
				{
					writeLog_("No peak mass tolerance specified. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				else if ( Real_buffer < 0 )
				{
					writeLog_("peak mass tolerance < 0. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setPeakMassTolerance(Real_buffer);
				
				Real_buffer = getDoubleOption_("match_peak_tol");
				if ( Real_buffer == -1 )
				{
					writeLog_("No match peak tolerance specified. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				else if ( Real_buffer < 0 )
				{
					writeLog_("Match peak tolerance < 0. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setMatchPeakTolerance(Real_buffer);
				
				Real_buffer = getDoubleOption_("ion_cutoff");
				if ( Real_buffer < 0 || Real_buffer > 1 )
				{
					writeLog_("Ion cutoff not in [0,1]. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setIonCutoffPercentage(Real_buffer);
				
				int_buffer = getIntOption_("pep_mass_unit");
				if ( (int_buffer < 0) || (int_buffer > max_peptide_mass_units) )
				{
					writeLog_("Illegal peptide mass unit (not in [0,2]). Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setPeptideMassUnit(int_buffer);
				
				int_buffer = getIntOption_("num_results");
				if ( (int_buffer < 1) )
				{
					writeLog_("Illegal number of results (< 1). Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setOutputLines(int_buffer);
				
				string_buffer = getStringOption_("enzyme_info");
				if ( !string_buffer.empty() )
				{
					string_buffer.split(':', substrings);
					if ( substrings.empty() ) substrings.push_back(string_buffer);
					
					vector< String > enzyme_info;
					for ( vector< String >::iterator einfo_i = substrings.begin(); einfo_i != substrings.end(); ++ einfo_i )
					{
						einfo_i->split(',', enzyme_info);
						if ( (enzyme_info.size() < 3) || (enzyme_info.size() > 4) )
						{
							writeLog_("Illegal number of informations for enzyme (not in [3,4]). Aborting!");
							return ILLEGAL_PARAMETERS;
						}
						if ( !((enzyme_info[1] == "0") || (enzyme_info[1] == "1"))  )
						{
							writeLog_("Cut direction for enzyme not specified correctly (has to be 1 (N to C)) or 0 (C to N))). Aborting!");
							return ILLEGAL_PARAMETERS;
						}
						if ( enzyme_info.size() == 3 ) enzyme_info.push_back("-");
						sequest_infile.addEnzymeInfo(enzyme_info);
					}
				}
				else
				{
					substrings.clear();
					Int highest_enzyme_number = sequest_infile.setEnzyme(getStringOption_("cleavage"));
					if ( highest_enzyme_number )
					{
						writeLog_("Chosen enzym is not in list. Aborting!");
						writeLog_(sequest_infile.getEnzymeInfoAsString());
						return ILLEGAL_PARAMETERS;
					}
				}
				
				Real_buffer = getDoubleOption_("prot_mass");
				if ( Real_buffer < 0 )
				{
					writeLog_("Illegal minimum protein mass (< 0). Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				else
				{
					Real_buffer2 = getDoubleOption_("max_prot_mass_or_tol");
					if ( Real_buffer2 < 0 )
					{
						writeLog_("Illegal maximum protein mass/ tolerance (< 0). Aborting!");
						return ILLEGAL_PARAMETERS;
					}
					else if ( Real_buffer2 < Real_buffer && Real_buffer2 > 100  ) // the second value has either got to be a mass (greater than the first one), or a probability
					{
						writeLog_("Illegal tolerance (not in [0, 100]). Aborting!");
						return ILLEGAL_PARAMETERS;
					}
					else sequest_infile.setProteinMassFilter(String(Real_buffer) + " " +  String(Real_buffer2));
				}
				
				int_buffer = getIntOption_("max_num_dif_AA_per_mod");
				if ( int_buffer < 0 )
				{
					writeLog_("No maximum number of modified amino acids per different modification. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setMaxAAPerModPerPeptide(int_buffer);
				
				int_buffer = getIntOption_("max_num_dif_mods_per_peptide");
				if ( int_buffer < 0 )
				{
					writeLog_("No maximum number of differential modifications per peptide. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setMaxModsPerPeptide(int_buffer);
				
				int_buffer = getIntOption_("nuc_reading_frame");
				if ( (int_buffer < 0) || (int_buffer > 9) )
				{
					writeLog_("Illegal number for nucleotide reading frame. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setNucleotideReadingFrame(int_buffer);
				
				int_buffer = getIntOption_("max_num_int_cleav_sites");
				if ( int_buffer < 0 )
				{
					writeLog_("Illegal number of maximum internal cleavage sites. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setMaxInternalCleavageSites(int_buffer);
				
				int_buffer = getIntOption_("match_peak_count");
				if ( int_buffer < 0 )
				{
					writeLog_("Illegal number of auto-detected peaks to try matching. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setMatchPeakCount(int_buffer);
				
				int_buffer = getIntOption_("match_peak_allowed_error");
				if ( int_buffer < 0 )
				{
					writeLog_("Illegal number of allowed errors in matching auto-detected peaks. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setMatchPeakAllowedError(int_buffer);
				
				sequest_infile.setShowFragmentIons(getFlag_("show_fragment_ions"));
// 				sequest_infile.setUsePhosphoFragmentation(getFlag_("use_phospho_fragmentation"));
				sequest_infile.setRemovePrecursorNearPeaks(getFlag_("remove_precursor_peak"));
				sequest_infile.setMassTypeParent(getFlag_("mass_type_precursor"));
				sequest_infile.setMassTypeFragment(getFlag_("mass_type_peak"));
				sequest_infile.setNormalizeXcorr(getFlag_("normalize_xcorr"));
				sequest_infile.setResiduesInUpperCase(!getFlag_("residues_in_lower_case"));
				
				string_buffer = getStringOption_("neutral_loss_ABY");
				string_buffer2 = "01";
				if ( (string_buffer.size() != 3) || (string_buffer2.find(string_buffer[0], 0) == String::npos) || (string_buffer2.find(string_buffer[1], 0) == String::npos) || (string_buffer2.find(string_buffer[2], 0) == String::npos) )
				{
					writeLog_("Neutral losses for ABY-ions not given (or illegal values given). Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				else
				{
					string_buffer.insert(2, 1, ' ');
					string_buffer.insert(1, 1, ' ');
					sequest_infile.setNeutralLossesForIons(string_buffer);
				}
				
				string_buffer = getStringOption_("ion_series_weights");
				string_buffer.split(',', substrings);
				if ( substrings.size() != 9 )
				{
					writeLog_("Weights for ion series not given (or illegal values given). Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				else
				{
					for ( vector< String >::iterator s_i = substrings.begin(); s_i != substrings.end(); ++s_i )
					{
						// the values are expected to be Real, otherwise they will be seen as 0!
						Real_buffer = atof(s_i->c_str());
						if ( (Real_buffer < 0) || (Real_buffer > 1) )
						{
							writeLog_("Illegal weights for ion series given. Aborting!");
							
							return ILLEGAL_PARAMETERS;
						}
						(*s_i) = String(Real_buffer);
					}
					string_buffer.implode(substrings.begin(), substrings.end(), " ");
					sequest_infile.setIonSeriesWeights(string_buffer);
				}
				
				// modifications
				string_buffer = getStringOption_("modifications");
				monoisotopic = getFlag_("use_monoisotopic_mod_mass");
				if ( !string_buffer.empty() ) // if modifications are used get look whether whether composition and residues (and type and name) is given which needs the isotope file, the name (and type) is used (then one additionally needs the modifications file) or only the mass and residues (and type and name) is given, in which case no further file is needed
				{
					string_buffer.split(':', substrings); // get the single modifications
					
					// one vector if compositions are used (needs isotope xml file) and one vector if masses were given
					vector< vector< String > > iso_sym_occ, mass_res_type_name;
					
					// to store the informations about modifications from the ptm xml file
					map< String, pair< String, String > > ptm_informations;
					
					// to get masses from a formula
					EmpiricalFormula add_e_formula, sub_e_formula;
					
					map< char, Real > stat_mods, dyn_mods;
					map< String, Real > terminal_mods;
					
					Int comp_mass_name_given(0);
					String types = "dyn#stat#cterminal#nterminal#cterminal_dyn#nterminal_dyn#cterminal_prot#nterminal_prot#";
					
					for ( vector< String >::const_iterator mod_i = substrings.begin(); mod_i != substrings.end(); ++mod_i )
					{
						// clear the formulae
						add_e_formula = "";
						sub_e_formula = "";
						
						if ( mod_i->empty() ) continue;
						
						iso_sym_occ.clear();
						// get the components of the modification
						mod_i->split(',', substrings2);
						if ( substrings2.empty() ) substrings2.push_back(*mod_i);
						mass_res_type_name.push_back(vector< String >(4));
						
						// check whether the first component is a composition, mass or name
						// remove + signs
						if ( substrings2[0].hasPrefix("+") ) substrings2[0].erase(0, 1);
						if ( substrings2[0].hasSuffix("+") ) substrings2[0].erase(substrings2[0].length() - 1, 1);
						if ( substrings2[0].hasSuffix("-") ) // a '-' at the end will not be converted
						{
							substrings2[0].erase(substrings2[0].length() - 1, 1);
							substrings2[0].insert(0, "-");
						}
						bool go_on = false;
						try
						{
							go_on = ( String(substrings2[0].toDouble()) != substrings2[0] );
							mass_res_type_name.back()[0] = substrings2[0]; // mass
							comp_mass_name_given = 0;
						}
						catch ( Exception::ConversionError ce )
						{
							go_on = true;
						}
						if ( go_on && get_composition_elements(substrings2[0], iso_sym_occ, '.').empty() ) // if it is a composition, put it into the vector
						{
							mass_res_type_name.back()[0] = substrings2[0]; // composition
							comp_mass_name_given = 1;
							go_on = false;
						}
						if ( go_on ) // check whether it's an empirical formula
						{
							String::size_type pos = substrings2[0].find("-");
							try
							{
								if ( pos != String::npos )
								{
									add_e_formula = substrings2[0].substr(0, pos);
									sub_e_formula = substrings2[0].substr(++pos);
								}
								else
								{
									add_e_formula = substrings2[0];
								}
								// sum up the masses
								if ( monoisotopic ) mass_res_type_name.back()[0] = String(add_e_formula.getMonoWeight() - sub_e_formula.getMonoWeight());
								else mass_res_type_name.back()[0] = String(add_e_formula.getAverageWeight() - sub_e_formula.getAverageWeight());
								go_on = false;
								comp_mass_name_given = -1;
							}
							catch ( Exception::ParseError pe )
							{
								go_on = true;
							}
						}
						if ( go_on ) // if it's a name, try to find it in the ptm xml file
						{
							if ( ptm_informations.empty() ) // if the ptm xml file has not been read yet, read it
							{
								if ( modifications_filename.empty() )
								{
									writeLog_("No modifications XML file given. Aborting!");
									return INPUT_FILE_NOT_FOUND;
								}
								if ( !File::readable(modifications_filename) )
								{
									writeLog_("Modifications XML file is not readable. Aborting!");
									return INPUT_FILE_NOT_READABLE;
								}
								
								// getting all available modifications from a file
								try
								{
									PTMXMLFile().load(modifications_filename, ptm_informations);
								}
								catch ( Exception::ParseError pe )
								{
									writeLog_(pe.getMessage());
									return PARSE_ERROR;
								}
							}
							
							if ( ptm_informations.find(substrings2[0]) == ptm_informations.end() ) // if the modification cannot be found
							{
								writeLog_("The Modification " + substrings2[0] + " can not be found in file " + modifications_filename + ". Aborting!");
								return ILLEGAL_PARAMETERS;
							}
							mass_res_type_name.back()[0] = ptm_informations[substrings2[0]].first; // composition
							mass_res_type_name.back()[1] = ptm_informations[substrings2[0]].second; // residues
							mass_res_type_name.back()[3] = substrings2[0]; // name
							
							// get the type
							if ( substrings2.size() > 1 )
							{
								// if it's not a legal type
								if ( types.find(substrings2[1]) == String::npos )
								{
									writeLog_("The given type (" + substrings2[1] + ") is neither dyn, stat, cterminal, nterminal, cterminal_dyn, nterminal_dyn, cterminal_prot nor nterminal_prot. Aborting!");
									return ILLEGAL_PARAMETERS;
								}
								mass_res_type_name.back()[2] = substrings2[1];
							}
							else mass_res_type_name.back()[2] = "dyn";
							comp_mass_name_given = 2;
						}
						
						// now get the residues and, if available the type and the name
						if ( comp_mass_name_given < 2 )
						{
							if ( substrings2.size() < 2 )
							{
								writeLog_("No residues for modification given (" + *mod_i + "). Aborting!");
								return ILLEGAL_PARAMETERS;
							}
							// if the type is a terminal, there may be no residues
							if ( types.find(substrings2[1]) != String::npos ) // the second one ought be residues if it's a non-terminal mod
							{
								if ( String("dyn#stat").find(substrings2[1]) != String::npos )
								{
									writeLog_("Non-terminal modification, but no residues given. Aborting!");
									return ILLEGAL_PARAMETERS;
								}
								mass_res_type_name.back()[2] = substrings2[1];
								
								// get the name
								if ( substrings2.size() > 2 ) mass_res_type_name.back()[3] = substrings2[2];
							}
							else
							{
								// get the residues
								mass_res_type_name.back()[1] = substrings2[1];
								mass_res_type_name.back()[1].substitute('*', 'X');
								
								// get the type
								if ( substrings2.size() > 2 )
								{
									// if it's not a legal type
									if ( types.find(substrings2[2]) == String::npos )
									{
										writeLog_("The given type (" + substrings2[2] + ") is neither dyn, stat, cterminal, nterminal, cterminal_dyn, nterminal_dyn, cterminal_prot nor nterminal_prot. Aborting!");
										return ILLEGAL_PARAMETERS;
									}
									mass_res_type_name.back()[2] = substrings2[2];
									
									// get the name
									if ( substrings2.size() > 3 ) mass_res_type_name.back()[3] = substrings2[3];
								}
								else mass_res_type_name.back()[2] = "dyn";
							}
						}
						
						// if a composition is given, get the corresponding mass
						if ( comp_mass_name_given > 0 )
						{
								// get the single components of the composition, if a name was given (for not doing this work twice)
							if ( comp_mass_name_given == 2 )
							{
								if ( !get_composition_elements(mass_res_type_name.back()[0], iso_sym_occ).empty() )
								{
									writeLog_("There's something wrong with this composition: " + mass_res_type_name.back()[0] + ". Aborting!");
									return ILLEGAL_PARAMETERS;
								}
							}
							for ( vector< vector< String > >::const_iterator comp_i = iso_sym_occ.begin(); comp_i != iso_sym_occ.end(); ++comp_i )
							{
								if ( (*comp_i)[0] == "0" )
								{
									if ( (*comp_i)[2].hasPrefix("-")  ) sub_e_formula += (*comp_i)[1] + (*comp_i)[2];
									else add_e_formula += (*comp_i)[1] + (*comp_i)[2];
								}
								else // if an isotope was used, get the mass
								{
								}
							}
								// sum up the masses
							if ( monoisotopic ) mass_res_type_name.back()[0] = String(add_e_formula.getMonoWeight() - sub_e_formula.getMonoWeight());
							else mass_res_type_name.back()[0] = String(add_e_formula.getAverageWeight() - sub_e_formula.getAverageWeight());
						}
						
						// for each type, collect all masses
						if ( mass_res_type_name.back()[2] == "dyn" ) // dynamic
						{
							for ( string::const_iterator c_i = mass_res_type_name.back()[1].begin(); c_i != mass_res_type_name.back()[1].end(); ++c_i )
							{
								dyn_mods[*c_i] += mass_res_type_name.back()[0].toDouble();
							}
						}
						else if ( mass_res_type_name.back()[2] == "stat" ) // static
						{
							for ( string::const_iterator c_i = mass_res_type_name.back()[1].begin(); c_i != mass_res_type_name.back()[1].end(); ++c_i )
							{
								stat_mods[*c_i] += mass_res_type_name.back()[0].toDouble();
							}
						}
						else // terminal
						{
							terminal_mods[mass_res_type_name.back()[2]] += mass_res_type_name.back()[0].toDouble();
						}
					}
					
					// save the dynamic modifications
					map< Real, String > dyn_mods_by_mass;
					for ( map< char, Real >::const_iterator dyn_i = dyn_mods.begin(); dyn_i != dyn_mods.end(); ++dyn_i )
					{
						dyn_mods_by_mass[dyn_i->second].append(1, dyn_i->first);
					}
					if ( dyn_mods_by_mass.size() <= 6 ) // Sequest doesn't allow more than six dynamic modifications (each amino acid may only be used once as only the mass of the last occurrence of an amino acid counts: 10K 12K leads to 12K)
					{
						String dyn_mods_as_string;
						for ( map< Real, String >::const_iterator dyn_i = dyn_mods_by_mass.begin(); dyn_i != dyn_mods_by_mass.end(); ++dyn_i )
						{
							if ( dyn_i != dyn_mods_by_mass.begin() ) dyn_mods_as_string.append(" ");
							dyn_mods_as_string.append(String(dyn_i->first) + " " + dyn_i->second);
						}
						sequest_infile.setDynMods(dyn_mods_as_string);
					}
					else
					{
						writeLog_("Too many dynamic modifications used (probably at least one amino acid is used more than once. This causes some trouble to Sequest). Aborting!");
						return ILLEGAL_PARAMETERS;
					}
					
					// save the static modifications
					for ( map< char, Real >::const_iterator stat_i = stat_mods.begin(); stat_i != stat_mods.end(); ++stat_i )
					{
						char_buffer = sequest_infile.setStatMod(String(stat_i->first), stat_i->second);
						if ( char_buffer )
						{
							writeLog_("Unknown amino acid (" + String(char_buffer) + ") given. Aborting!");
							return ILLEGAL_PARAMETERS;
						}
					}
					// save the terminal modifications
					sequest_infile.setStatNTermMod(terminal_mods["nterminal"]);
					sequest_infile.setStatCTermMod(terminal_mods["cterminal"]);
					sequest_infile.setDynNTermMod(terminal_mods["nterminal_dyn"]);
					sequest_infile.setDynCTermMod(terminal_mods["cterminal_dyn"]);
					sequest_infile.setStatNTermProtMod(terminal_mods["nterminal_prot"]);
					sequest_infile.setStatCTermProtMod(terminal_mods["cterminal_prot"]);
				}
				
				string_buffer = getStringOption_("partial_sequence");
				string_buffer.split(',', substrings);
				string_buffer.implode(substrings.begin(), substrings.end(), " ");
				sequest_infile.setPartialSequence(string_buffer);
				
				string_buffer = getStringOption_("header_filter");
				string_buffer.split(',', substrings);
				string_buffer.implode(substrings.begin(), substrings.end(), " ");
				sequest_infile.setSequenceHeaderFilter(string_buffer);
			}
			
			if ( sequest_out )
			{
				output_filename = getStringOption_("out");
				if ( output_filename.empty() )
				{
					writeLog_("No output file specified. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				files.push_back(make_pair(output_filename, 4));
				
				p_value = getDoubleOption_("p_value");
				if ( (p_value <= 0) || (p_value > 1) )
				{
					writeLog_("P-value not in (0,1]. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
			}
			
			//-------------------------------------------------------------
			// running program according to parameters
			//-------------------------------------------------------------
			
			// checking accessability of files
			bool existed = false;
			UInt file_tag;
			
			for ( vector< pair< String, UInt > >::const_iterator files_i = files.begin(); files_i != files.end(); ++files_i )
			{
				string_buffer = files_i->first;
				file_tag = files_i->second;
				
				if ( (file_tag & 1) && !File::exists(string_buffer) )
				{
					throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, string_buffer);
				}
				
				if ( (file_tag & 2) && !File::readable(string_buffer) )
				{
					throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, string_buffer);
				}
				
				existed = File::exists(string_buffer);
				if ( (file_tag & 4) && !File::writable(string_buffer) )
				{
					throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, string_buffer);
				}
				else if ( !existed ) remove(string_buffer.c_str());
				existed = false;
			}
			
			// creating the input file
			if ( sequest_in )
			{
				sequest_infile.store(input_filename);
			}
			
			bool make_dtas = ( sequest_out && !sequest_in ) ? false : true; // if only sequest_out is set, just get the retention times
			// creating the dta files
			if ( make_dtas ) // if there are already .dta files in the folder, stop the adapter
			{
				vector<String> dummy;
				if (File::fileList(temp_data_directory,String("*.dta_*"),dummy))
				{
					writeLog_("There are already dta files in directory " + temp_data_directory + ". Aborting!");
					
					// deleting all temporary files
					for ( vector< pair< String, UInt > >::const_iterator files_i = files.begin(); files_i != files.end(); ++files_i )
					{
						if ( files_i->second & 8 ) remove(files_i->first.c_str());
					}
					return UNKNOWN_ERROR;
				}
			}
			
			MSExperiment<> msexperiment;
			UInt msms_spectra_in_file;
			UInt msms_spectra_altogether = 0;
			if ( make_dtas ) writeLog_("creating dta files");
			dtas = 0;
			String basename, dta_files_common_name;
			FileHandler fh;
			FileHandler::Type type;
			for ( vector< String >::const_iterator spec_i = spectra.begin(); spec_i != spectra.end(); ++spec_i )
			{
				basename = File::basename(*spec_i);
				dta_files_common_name = temp_data_directory + basename;
				
				type = fh.getTypeByContent(*spec_i);
				if ( type == FileHandler::UNKNOWN )
				{
					writeLog_("Could not determine type of the file. Aborting!");
					return PARSE_ERROR;
				}
				fh.loadExperiment(*spec_i, msexperiment, type);
				
				msms_spectra_in_file = MSExperiment2DTAs(msexperiment, dta_files_common_name, charges, filenames_and_precursor_retention_times, make_dtas);
				
				writeLog_(String(msms_spectra_in_file) + " MS/MS spectra in file " + *spec_i);
				
				msms_spectra_altogether += msms_spectra_in_file;
			}
			
			if ( !msms_spectra_altogether )
			{
				writeLog_("No MS/MS spectra found in any of the mz files. Aborting!");
				return UNKNOWN_ERROR;
			}
			
			// (3.2.3) running the program
			if ( sequest_in && sequest_out )
			{
				// creating a batch file for windows (command doesn't accept commands that are longer than 256 chars)
				String sequest_screen_output; // direct the screen-output to a file
				do
				{
					sequest_screen_output = String::random(10);
				}
				while ( File::exists(sequest_screen_output) );
				files.push_back(make_pair(temp_data_directory + sequest_screen_output, 4+8));
				
				ofstream batchfile(String(temp_data_directory + batch_filename).c_str());
				if ( !batchfile )
				{
					throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, temp_data_directory + batch_filename);
				}
				String call = "rdesktop";
				if ( !user.empty() ) call.append(" -u " + user);
				if ( !password.empty() ) call.append(" -p \"" + password + "\"");
				call.append(" -s cmd\\ /K\\ \"");
// 				call.append("echo net use " + temp_data_directory_win.substr(0,2) + " \\\\" + temp_data_directory_network + " && ");
				call.append("net use " + temp_data_directory_win.substr(0,2) + " \\\\" + temp_data_directory_network.substr(0, temp_data_directory_network.length() - 1) + " && ");
// 				call.append(" net use " + temp_data_directory_win.substr(0,2) + " " + temp_data_directory_network + " && ");
				
				batchfile << String(" cd " + temp_data_directory_win + " && " + temp_data_directory_win.substr(0,2));
				
				for ( PointerSizeUInt i = 0; i <= (PointerSizeUInt) (dtas / max_dtas_per_run); ++i )
				{
					batchfile << " && " << sequest_directory_win << "sequest.exe -P" << input_file_directory_network << File::basename(input_filename) << "  " << temp_data_directory_network << "*.dta_" << i<< " > " <<  temp_data_directory_network << sequest_screen_output << " && move sequest.log sequest.log" << i;
				}
				batchfile << " && " << sequest_directory_win.substr(0,2) << endl;
// 				batchfile << std::endl << String(" net use /delete " + temp_data_directory_win.substr(0,2));
// 				batchfile << std::endl << "logoff";
				batchfile.close();
				batchfile.clear();
				
				call.append(temp_data_directory_win + batch_filename + " && net use /delete " + temp_data_directory_win.substr(0,2) + " && logoff" +  "\" " + sequest_computer);
				writeLog_("System call: " + call);
				int status = system(call.c_str());
				
				if ( status != 0 )
				{
					writeLog_("Sequest problem. Aborting! (Details can be seen in the logfile: \"" + logfile + "\")");
					
					// deleting all temporary files
					for ( vector< pair< String, UInt > >::const_iterator files_i = files.begin(); files_i != files.end(); ++files_i )
					{
						if ( files_i->second & 8 ) remove(files_i->first.c_str());
					}
					
				// remove all dtas
					if ( !keep_dta_files )
					{
						writeLog_("removing dta files");
						for ( map< String, Real >::const_iterator dta_names_i = filenames_and_precursor_retention_times.begin(); dta_names_i != filenames_and_precursor_retention_times.end(); ++dta_names_i )
						{
							if ( !File::remove(dta_names_i->first) ) writeLog_("'" + string_buffer + "' could not be removed!");
						}
						return EXTERNAL_PROGRAM_ERROR;
					}
				}
				
				bool no_log = false;
				string_buffer.clear();
				for ( PointerSizeUInt i = 0; i <= (PointerSizeUInt) (dtas / max_dtas_per_run); ++i )
				{
					ifstream sequest_log(string(temp_data_directory + "sequest.log" + String(i)).c_str()); // write sequest log to logfile
					if ( !sequest_log )
					{
						no_log = true;
						break;
					}
					else
					{
						sequest_log.seekg (0, ios::end);
						streampos length = sequest_log.tellg();
						sequest_log.seekg (0, ios::beg);
						char * buffer = new char[length];
						sequest_log.read (buffer, length);
						sequest_log.close();
						sequest_log.clear();
						string_buffer2.assign(buffer);
						delete(buffer);
						string_buffer.append(string_buffer2.substr(string_buffer2.find("Total search time")));
						remove(string(temp_data_directory + "sequest.log" + String(i)).c_str());
					}
				}
				if ( no_log )
				{
					writeLog_("No Sequest log found!");
					
					// remove all dtas
					if ( !keep_dta_files )
					{
						writeLog_("removing dta files");
						for ( map< String, Real >::const_iterator dta_names_i = filenames_and_precursor_retention_times.begin(); dta_names_i != filenames_and_precursor_retention_times.end(); ++dta_names_i )
						{
							if ( !File::remove(dta_names_i->first) ) writeLog_("'" + string_buffer + "' could not be removed!");
						}
					}
					
					// deleting all temporary files
					for ( vector< pair< String, UInt > >::const_iterator files_i = files.begin(); files_i != files.end(); ++files_i )
					{
						if ( files_i->second & 8 ) remove(files_i->first.c_str());
					}
					return EXTERNAL_PROGRAM_ERROR;
				}
				else writeLog_(string_buffer);
			}
			
			if ( sequest_out )
			{
				// remove all dtas
				if ( !keep_dta_files )
				{
					writeLog_("removing dta files");
					for ( map< String, Real >::const_iterator dta_names_i = filenames_and_precursor_retention_times.begin(); dta_names_i != filenames_and_precursor_retention_times.end(); ++dta_names_i )
					{
						if ( !File::remove(dta_names_i->first) ) writeLog_("'" + string_buffer + "' could not be removed!");
					}
				}
				
				SequestOutfile sequest_outfile;
				vector<PeptideIdentification> peptide_identifications;
				vector<ProteinIdentification> pis;
				UInt peptide_identification_size = peptide_identifications.size();
				ProteinIdentification protein_identification;
				
				vector<String> out_files;
				if (!File::fileList(out_directory, String("*.out"), out_files))
				{
					writeLog_(String("Error: No .out files found in '") + out_directory + "'. Aborting!");
					
					// deleting all temporary files
					for ( vector< pair< String, UInt > >::const_iterator files_i = files.begin(); files_i != files.end(); ++files_i )
					{
						if ( files_i->second & 8 ) remove(files_i->first.c_str());
					}
					
					return UNKNOWN_ERROR;
				}
				
				vector< pair < String, vector< Real > > > filenames_and_pvalues;
				for ( vector< String >::iterator f_i = out_files.begin(); f_i != out_files.end(); ++f_i )
				{
					filenames_and_pvalues.push_back(make_pair(out_directory + *f_i, vector< Real >()));
				}
// 				sequest_outfile.getPValuesFromOutFiles(filenames_and_pvalues);
				
				// set the parameters
				ProteinIdentification::SearchParameters sp;
				sp.db = "Fasta";
				sp.taxonomy = sequest_infile.getSequenceHeaderFilter();
				if ( monoisotopic ) sp.mass_type = ProteinIdentification::MONOISOTOPIC;
				else sp.mass_type = ProteinIdentification::AVERAGE;
				for ( vector< Int >::const_iterator c_i = charges.begin(); c_i != charges.end(); ++c_i )
				{
					if ( *c_i > 0 ) sp.charges.append("+");
					sp.charges.append(String(*c_i));
				}
				if ( sequest_infile.getEnzyme() == "Trypsin" ) sp.enzyme = ProteinIdentification::TRYPSIN;
				else if ( sequest_infile.getEnzyme() == "No_Enzyme" ) sp.enzyme = ProteinIdentification::NO_ENZYME;
				else sp.enzyme = ProteinIdentification::UNKNOWN_ENZYME;
				sp.peak_mass_tolerance = sequest_infile.getPeakMassTolerance();
				sp.precursor_tolerance = sequest_infile.getPrecursorMassTolerance();
				protein_identification.setSearchParameters(sp);
				
				for ( vector< pair < String, vector< Real > > >::iterator fp_i = filenames_and_pvalues.begin(); fp_i != filenames_and_pvalues.end(); ++fp_i )
				{
					try
					{
						sequest_outfile.load(fp_i->first, peptide_identifications, protein_identification, p_value, fp_i->second, database);
					}
					catch( Exception::ParseError pe )
					{
						// deleting all temporary files
						for ( vector< pair< String, UInt > >::const_iterator files_i = files.begin(); files_i != files.end(); ++files_i )
						{
							if ( files_i->second & 8 ) remove(files_i->first.c_str());
						}
						writeLog_(pe.getMessage());
						return INPUT_FILE_CORRUPT;
					}
					
					// save the retention times if peptides have been identified to the p-level
					if ( peptide_identification_size != peptide_identifications.size() )
					{
						peptide_identification_size = peptide_identifications.size();
						string_buffer = fp_i->first;
						string_buffer.replace(string_buffer.length() - 3, 3, "out");
						peptide_identifications.back().setMetaValue("RT", filenames_and_precursor_retention_times[string_buffer]);
					}
				}
				
				pis.push_back(protein_identification);
				
				IdXMLFile().store(output_filename, pis, peptide_identifications);
				
				// remove all outs
				if ( !keep_out_files )
				{
					writeLog_("removing out files");
					for ( vector<String>::const_iterator i = out_files.begin(); i != out_files.end(); ++i )
					{
						if ( !File::remove(out_directory + *i) )
						{
							writeLog_(String("'") + out_directory + *i + "' could not be removed!");
						}
					}
				}
			}
			
			// deleting all temporary files
			for ( vector< pair< String, UInt > >::const_iterator files_i = files.begin(); files_i != files.end(); ++files_i )
			{
				if ( files_i->second & 8 ) remove(files_i->first.c_str());
			}
			
			return EXECUTION_OK;
		}
};

//@endcond



int main( int argc, char ** argv )
{
	TOPPSequestAdapter tool;
	
	return tool.main(argc,argv);
}
