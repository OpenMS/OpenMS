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


#include <OpenMS/FORMAT/AnalysisXMLFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/SequestInfile.h>
#include <OpenMS/FORMAT/SequestOutfile.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
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
	
	This wrapper application serves for getting peptide identifications
	for MS/MS spectra. The wrapper can be executed in three different
	modes:
	<ol>	
				<li>
				The whole process of identification via Sequest is executed. 
				Inputfile is one (or more) mzXML file containing the MS/MS spectra
				for which identifications are to be found
				and one ore two databases in FASTA format containing
				the possible proteins.
				The results are written as an analysisXML output file. This mode is selected
			 	by default.
				Note: You need a user with network access on the computer hosting sequest.
			 	</li>
				
				<li>
				Only the first part of the identification process is performed.
				This means that a Sequest input file is generated and dta files are
				created from the mzXML file.
				Calling an Sequest process should look like the following:				
				
				@code sequest -P\<inputfilename\> \<path to dta files\>*.dta  @endcode

				Consult your Sequest reference manual for further details.
				
				This mode is selected by the <b>-sequest_in</b> option in the command line.
				</li>
				
				<li>
				Only the second part of the identification process is performed.
				This means that the output of sequest is translated into analysisXML.
				
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
		static const SignedInt max_peptide_mass_units = 2;
		static const UnsignedInt max_dtas_per_run = 1000; // sequest has a problem when there are too many dtas, so they have to be splitted, 1000 seemed to work very good
		PointerSizeUInt dtas;

		void registerOptionsAndFlags_()
		{
			addText_("The definitions for the parameters are taken from the site:\n"
										 "http://www.grosse-coosmann.de/~florian/Parameters.html#file.");
			registerStringOption_("out", "<file>", "", "output file in analysisXML format.\n"
			                                           "Note: In mode 'sequest_in' a Sequest input file is written.", false);
			registerStringOption_("in", "<file>", "", "input file(s) in mzXML format (comma-separated).\n"
					 																			"Note: In mode 'sequest_out' a directory with Sequest results files\n"
																								"(*.out) is read", false);
			registerFlag_("sequest_in", "if this flag is set the SequestAdapter will read in mzXML\n"
																								"and write an Sequest input file\n"
																								"and create dta files from the given mzXML files");
			registerFlag_("sequest_out", "if this flag is set the SequestAdapter will read in Sequest result files\n"
																									"and write analysisXML");
			registerStringOption_("mzXMLs", "<file>", "", "when using sequest_out the mzXML files (comma-separated)\n"
																																						"have to be given to retrieve the retention times", false);
			registerFlag_("show_enzyme_numbers", "show a list with enzymes and corresponding numbers to choose from");
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
			registerStringOption_("peptide_prophet_directory", "<dir>", "", "the directory in which peptide prophet is located\n"
																																"peptide prophet is used to compute the p-values,\n"
																																"without it, all peptide found are accepted", false);
			addEmptyLine_();
			registerDoubleOption_("precursor_mass_tolerance", "<tol>", 2.0 , "the precursor mass tolerance", false);
			registerDoubleOption_("peak_mass_tolerance", "<tol>", 1.0, "the peak mass tolerance", false);
			registerDoubleOption_("p_value", "<prob>", 1.0, "annotations with inferior p-value are ignored", false);
		  registerStringOption_("charges", "[1>3,5]", "", "comma-seperated list of charge states (or ranges)", false);
			registerIntOption_("num_results", "<num>", 1, "the maximal number of results (peptides) to show (per scan/dta)", false);
			registerIntOption_("cleavage", "<num>", -1, "the number of the enzyme used for digestion", false);
			registerStringOption_("enzyme_info", "<>", "", "information about the enzyme used\n"
																																							"<name>,<cut direction: N to C?>,<cuts after>,<doesn't cut before>\n"
																																							"cuts after, doesn't cut before: amino acids in 1-letter code\n"
																																							"or '-' for unspecific cleavage", false);
			registerStringOption_("dyn_mods", "[44,s:80,TG]", "", "This value consists of colon-seperated pairs of variable modifications.\n"
																								"Each pair has two comma-seperated elements: A mass and a list of amino acids.\n"
																								"Sequest only applies the last modification character without warning.\n"
																								"Don't use \"44,S:80,ST\". It is interpreted as \"80,ST\"!.\n"
																								"Up to six modifications are allowed, if more are given, they are ignored.", false);
			registerDoubleOption_("dyn_N_term_mod", "", 0, "This modification mass that may be added to each N-terminus", false);
			registerDoubleOption_("dyn_C_term_mod", "", 0, "This modification mass that may be added to each C-terminus", false);
			registerDoubleOption_("stat_N_term_mod", "", 0, "This mass is added to each peptide N-terminus", false);
			registerDoubleOption_("stat_C_term_mod", "", 0, "This mass is added to each peptide C-terminus", false);
			registerDoubleOption_("stat_N_term_prot_mod", "", 0, "This mass is added to each protein N-terminus", false);
			registerDoubleOption_("stat_C_term_prot_mod", "", 0, "This mass is added to each protein C-terminus", false);
			registerStringOption_("stat_mods", "", "", "This is a colon-seperated list of amino acids in one letter code\n"
																								 "and their corrpesponding mass added: <AA_1>,<mass_1>:<AA_2>,<mass_2>:...\n"
																								 "(if several amino acids shall have the same static modification: KRLNH,10.3)", false);
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
			registerFlag_("keep_dta_files", "If set the dta-files that were created from the mzXML-files are not removed");
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
							if ( name.find(" ") == string::npos ) return true;
							else return false;
						}
						else return false;
					}
					else return true;
				}
			}
			return false;
		}

		bool correctNetworkPath(String& network_path, UnsignedInt backslashes = 2)
		{
			UnsignedInt pos = 0;
			while ( (pos < network_path.length()) && (network_path[pos] == '\\') ) ++pos;
			if ( pos < backslashes ) network_path.insert(network_path.begin(), backslashes-pos, '\\');
			else network_path.erase(0, pos-backslashes);
			if ( network_path.length() < backslashes+1 ) return false;
			return true;
		}

		UnsignedInt
		MSExperiment2DTAs(
			MSExperiment<>& msexperiment,
			const String& common_name,
			const vector< SignedInt >& charges,
			map< String, DoubleReal >& filenames_and_precursor_retention_times,
			bool make_dtas = true)
		throw (Exception::UnableToCreateFile)
		{
			DTAFile dtafile;
			String filename;
			UnsignedInt scan_number = 0;
			UnsignedInt msms_spectra = 0;
			
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
						filename.replace(filename.length() - 4, 4, ".out");
						filenames_and_precursor_retention_times[filename] = spec_i->getRetentionTime();
					}
					else
					{
						for ( vector< SignedInt >::const_iterator i = charges.begin(); i != charges.end(); ++i )
						{
							filename = common_name + "." + String(scan_number) + "." + *i + ".dta_" + String( (PointerSizeUInt) (dtas / max_dtas_per_run) );
							if ( make_dtas )
							{
								++dtas;
								spec_i->getPrecursorPeak().setCharge(*i);
								dtafile.store(filename, *spec_i);
							}
							filename.replace(filename.length() - 4, 4, ".out");
							filenames_and_precursor_retention_times[filename] = spec_i->getRetentionTime();
						}
						spec_i->getPrecursorPeak().setCharge(0);
					}
				}
			}
			
//			for (map< String, DoubleReal >::const_iterator it = filenames_and_precursor_retention_times.begin(); it!=filenames_and_precursor_retention_times.end(); ++it)
//			{
//				cout << it->first << " -> " << it->second << endl;
//			}
			
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
				peptide_prophet_directory;
			
			ContactPerson contact_person;
			
			bool
				sequest_in,
				sequest_out,
				keep_out_files,
				keep_dta_files;
			
			vector< String >
				substrings,
				substrings2,
				spectra;
			
			vector< SignedInt > charges;
			
			char char_buffer;
			
			DoubleReal
				DoubleReal_buffer,
				DoubleReal_buffer2;
			
			SignedInt int_buffer;
			
			Real p_value = 0.05;
			
			map< String, DoubleReal > filenames_and_precursor_retention_times;
			
			vector< String > tmp_names;
			
			//-------------------------------------------------------------
			// (2) parsing and checking parameters
			//-------------------------------------------------------------
			
			// only show the available enzymes, then quit
			if ( getFlag_("show_enzyme_numbers") )
			{
				writeLog_("Option show_enzyme_numbers chosen.");
				writeLog_(sequest_infile.getEnzymeInfo());
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
				tmp_names.push_back(logfile);
			}

			string_buffer = getStringOption_("charges");
			if ( string_buffer.empty() )
			{
				writeLog_("No charge states given. Aborting!");
				return ILLEGAL_PARAMETERS;
			}
			else
			{
				SignedInt range_start, range_end;
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

							for ( SignedInt i = min(range_start, range_end); i <= max(range_start, range_end); ++i )
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
				for ( vector< SignedInt >::iterator i = charges.begin(); i != --charges.end(); )
				{
					if ( (*i) == (*(i+1)) ) charges.erase(i+1);
					else ++i;
				}
			}

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
				}
				else // if only sequest_out is set, in is the out_directory
				{
					out_directory = string_buffer;
					if ( !out_directory.empty() ) out_directory.ensureLastChar('/');
					
					// if only sequest_out is set, the mzXML files have to be given to retrieve the retention times
					string_buffer = getStringOption_("mzXMLs");
					if ( string_buffer.empty() )
					{
						writeLog_("No mzXML files specified. Aborting!");
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
			
			temp_data_directory = getStringOption_("temp_data_directory");
			if ( temp_data_directory.empty() )
			{
				writeLog_("No directory for temporary files given. Aborting!");
				return ILLEGAL_PARAMETERS;
			}
			temp_data_directory.ensureLastChar('/');
			
			temp_data_directory_win = getStringOption_("temp_data_directory_win");
			temp_data_directory_win.ensureLastChar('\\');
			
			if ( !isWinFormat(temp_data_directory_win) )
			{
				writeLog_("Windows path for the directory for temporary files has wrong format: " + temp_data_directory_win + ". Aborting!");
				return ILLEGAL_PARAMETERS;
			}
			temp_data_directory_network = getStringOption_("temp_data_directory_network");
			if ( temp_data_directory_network.empty() )
			{
				writeLog_("Network path for the directory for temporary files is empty. Aborting!");
				return ILLEGAL_PARAMETERS;
			}
			if ( !correctNetworkPath(temp_data_directory_network) )
			{
				writeLog_(temp_data_directory_network + "is no network path. Aborting!");
				return ILLEGAL_PARAMETERS;
			}
			
			contact_person.setName(getStringOption_("contact_name"));
			contact_person.setInstitution(getStringOption_("contact_institution"));
			contact_person.setContactInfo(getStringOption_("contact_info"));

			if ( sequest_in )
			{
				if ( !sequest_out )
				{
					input_filename = getStringOption_("out");
					if ( input_filename.empty() )
					{
						writeLog_("No output file specified. Aborting!");
						return ILLEGAL_PARAMETERS;
					}
				}
				else
				{
					input_filename = getStringOption_("sequest_input");
					if ( input_filename.empty() )
					{
						input_filename = temp_data_directory + "temp.sequest.in";
						tmp_names.push_back(input_filename);
						input_file_directory_network = temp_data_directory_network;
					}
					else input_file_directory_network = getStringOption_("sequest_input_directory_network");
				}
			}
			
			if ( sequest_in && sequest_out )
			{
				if ( !correctNetworkPath(input_file_directory_network) )
				{
					writeLog_(input_file_directory_network + "is no network path. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				string_buffer = File::basename(input_filename);
				if ( !input_file_directory_network.hasSuffix(string_buffer) ) input_file_directory_network.append("\\" + string_buffer);
				
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
				
				peptide_prophet_directory = getStringOption_("peptide_prophet_directory");
				if ( !peptide_prophet_directory.empty() ) peptide_prophet_directory.ensureLastChar('/');
				
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
				tmp_names.push_back(batch_filename);
			}
			else if ( !batch_filename.hasSuffix(".bat") ) batch_filename.append(".bat");
			
			database = getStringOption_("db");
			if ( database.empty() )
			{
				writeLog_("No database specified. Aborting!");
				return ILLEGAL_PARAMETERS;
			}
			
			if ( sequest_in )
			{
				database_directory_network = getStringOption_("db_directory_network");
				if ( !correctNetworkPath(database_directory_network) )
				{
					writeLog_(database_directory_network + "is no network path. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				string_buffer = File::basename(database);
				if ( !database_directory_network.hasSuffix(string_buffer) ) database_directory_network.append("\\" + string_buffer);
				sequest_infile.setDatabase(database_directory_network);
				
				DoubleReal_buffer = getDoubleOption_("precursor_mass_tolerance");
				if ( DoubleReal_buffer == -1 )
				{
					writeLog_("No precursor mass tolerance specified. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				else if ( DoubleReal_buffer < 0 )
				{
					writeLog_("Precursor mass tolerance < 0. Aborting!");

					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setPeptideMassTolerance(DoubleReal_buffer);
				
				DoubleReal_buffer = getDoubleOption_("peak_mass_tolerance");
				if ( DoubleReal_buffer == -1 )
				{
					writeLog_("No peak mass tolerance specified. Aborting!");

					return ILLEGAL_PARAMETERS;
				}
				else if ( DoubleReal_buffer < 0 )
				{
					writeLog_("Fragment ion tolerance < 0. Aborting!");

					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setFragmentIonTolerance(DoubleReal_buffer);
				
				DoubleReal_buffer = getDoubleOption_("match_peak_tol");
				if ( DoubleReal_buffer == -1 )
				{
					writeLog_("No match peak tolerance specified. Aborting!");

					return ILLEGAL_PARAMETERS;
				}
				else if ( DoubleReal_buffer < 0 )
				{
					writeLog_("Match peak tolerance < 0. Aborting!");

					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setMatchPeakTolerance(DoubleReal_buffer);
				
				DoubleReal_buffer = getDoubleOption_("ion_cutoff");
				if ( DoubleReal_buffer < 0 || DoubleReal_buffer > 1 )
				{
					writeLog_("Ion cutoff not in [0,1]. Aborting!");

					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setIonCutoffPercentage(DoubleReal_buffer);
				
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
					SignedInt highest_enzyme_number = sequest_infile.setEnzymeNumber(getIntOption_("cleavage"));
					if ( highest_enzyme_number )
					{
						writeLog_("Enzyme number has to be in [0," + String(highest_enzyme_number) + "]. Aborting!");
						return ILLEGAL_PARAMETERS;
					}
				}
				
				DoubleReal_buffer = getDoubleOption_("prot_mass");
				if ( DoubleReal_buffer < 0 )
				{
					writeLog_("Illegal minimum protein mass (< 0). Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				else
				{
					DoubleReal_buffer2 = getDoubleOption_("max_prot_mass_or_tol");
					if ( DoubleReal_buffer2 < 0 )
					{
						writeLog_("Illegal maximum protein mass/ tolerance (< 0). Aborting!");
						return ILLEGAL_PARAMETERS;
					}
					else if ( DoubleReal_buffer2 < DoubleReal_buffer && DoubleReal_buffer2 > 100  ) // the second value has either got to be a mass (greater than the first one), or a probability
					{
						writeLog_("Illegal tolerance (not in [0, 100]). Aborting!");
						return ILLEGAL_PARAMETERS;
					}
					else sequest_infile.setProteinMassFilter(String(DoubleReal_buffer) + " " +  String(DoubleReal_buffer2));
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
				if ( (string_buffer.size() != 3) || (string_buffer2.find(string_buffer[0], 0) == string::npos) || (string_buffer2.find(string_buffer[1], 0) == string::npos) || (string_buffer2.find(string_buffer[2], 0) == string::npos) )
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
						DoubleReal_buffer = atof(s_i->c_str());
						if ( (DoubleReal_buffer < 0) || (DoubleReal_buffer > 1) )
						{
							writeLog_("Illegal weights for ion series given. Aborting!");

							return ILLEGAL_PARAMETERS;
						}
						(*s_i) = String(DoubleReal_buffer);
					}
					string_buffer.implode(substrings.begin(), substrings.end(), " ");
					sequest_infile.setIonSeriesWeights(string_buffer);
				}
				
				string_buffer = getStringOption_("dyn_mods");
				if ( !string_buffer.empty() )
				{
					string_buffer.split(':', substrings);
					if ( substrings.empty() ) substrings.push_back(string_buffer);
					Real f_buffer;
					char c_buffer[41]; c_buffer[40] = 0;
					for ( vector< String >::iterator s_i = substrings.begin(); s_i != substrings.end(); ++s_i )
					{
						if ( sscanf(s_i->c_str(), "%f,%40s", &f_buffer, c_buffer) != 2 )
						{
							writeLog_("Illegal number of parameters for dynamic modification given. Aborting!");

							return ILLEGAL_PARAMETERS;
						}
						(*s_i) = String(f_buffer)+" "+c_buffer;
					}
					string_buffer.implode(substrings.begin(), substrings.end(), " ");
					sequest_infile.setDynMods(string_buffer);
				}
				
				sequest_infile.setDynNTermMod(getDoubleOption_("dyn_N_term_mod"));
				sequest_infile.setDynCTermMod(getDoubleOption_("dyn_C_term_mod"));
				
				sequest_infile.setStatNTermMod(getDoubleOption_("stat_N_term_mod"));
				sequest_infile.setStatCTermMod(getDoubleOption_("stat_C_term_mod"));
				
				sequest_infile.setStatNTermProtMod(getDoubleOption_("stat_N_term_prot_mod"));
				sequest_infile.setStatCTermProtMod(getDoubleOption_("stat_C_term_prot_mod"));
				
				string_buffer = getStringOption_("stat_mods");
				if ( !string_buffer.empty() )
				{
					string_buffer.split(':', substrings);
					if ( substrings.empty() ) substrings.push_back(string_buffer);

					String::iterator ss_i;
					for ( vector< String >::iterator s_i = substrings.begin(); s_i != substrings.end(); ++s_i )
					{
						s_i->split(',', substrings2);
						if ( substrings2.size() != 2 || substrings2[0].empty() || substrings2[1].empty() )
						{
							writeLog_("Unexpected format for static modification found. Aborting!");
							return ILLEGAL_PARAMETERS;
						}
						ss_i = --substrings2[1].end();
						if ( *ss_i == '-' )
						{
							substrings2[1].erase(ss_i);
							substrings2[1].insert(0, "-");
						}
						char_buffer = sequest_infile.setStatMod(substrings2[0], substrings2[1].toFloat());
						if ( char_buffer )
						{
							writeLog_("Unknown amino acid (" + String(char_buffer) + ") given. Aborting!");
							return ILLEGAL_PARAMETERS;
						}
					}
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
				string_buffer = getStringOption_("out");
				if ( string_buffer.empty() )
				{
					writeLog_("No output file specified. Aborting!");

					return ILLEGAL_PARAMETERS;
				}
				else output_filename = string_buffer;
			
				p_value = getDoubleOption_("p_value");
				if ( (p_value <= 0) || (p_value > 1) )
				{
					writeLog_("P-value not in (0,1]. Aborting!");

					return ILLEGAL_PARAMETERS;
				}
			}
			
			//-------------------------------------------------------------
			// (3) running program according to parameters
			//-------------------------------------------------------------

			// (3.1) checking accessability of files
			if ( sequest_in )
			{
				if ( !File::writable(input_filename) )
				{
					throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, input_filename);
				}
				if ( !File::writable(temp_data_directory + batch_filename) )
				{
					throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, temp_data_directory + batch_filename);
				}
			}
			
			// (3.1.2) output file
			if ( sequest_out )
			{
				if ( !File::writable(output_filename) )
				{
					throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, output_filename);
				}
				
				// database
				if ( !File::exists(database) )
				{
					throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, database);
				}
				if ( !File::readable(database) )
				{
					throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, database);
				}
				if ( File::empty(database) )
				{
					throw Exception::FileEmpty(__FILE__, __LINE__, __PRETTY_FUNCTION__, database);
				}
			}
			
			// (3.2.1) creating the input file
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
					for ( vector< String >::const_iterator tmp_names_i = tmp_names.begin(); tmp_names_i != tmp_names.end(); ++tmp_names_i ) remove(tmp_names_i->c_str());
					return UNKNOWN_ERROR;
				}
			}
			
			MSExperiment<> msexperiment;
			UnsignedInt msms_spectra_in_file;
			UnsignedInt msms_spectra_altogether = 0;
			if ( make_dtas ) writeLog_("creating dta files");
			dtas = 0;
			String basename, dta_files_common_name;
			FileHandler fh;
			FileHandler::Type type;
			for ( vector< String >::const_iterator spec_i = spectra.begin(); spec_i != spectra.end(); ++spec_i )
			{
				basename = File::basename(*spec_i);
				dta_files_common_name = out_directory + basename;
				
				type = fh.getTypeByContent(basename);
				if ( type == FileHandler::UNKNOWN )
				{
					writeLog_("Could not determine type of the file. Aborting!");
					return PARSE_ERROR;
				}
				fh.loadExperiment(basename, msexperiment, type);
				
				msms_spectra_in_file = MSExperiment2DTAs(msexperiment, dta_files_common_name, charges, filenames_and_precursor_retention_times, make_dtas);
				
				writeLog_(String(msms_spectra_in_file) + " MS/MS spectra in file " + basename);

				msms_spectra_altogether += msms_spectra_in_file;
			}

			if ( !msms_spectra_altogether )
			{
				writeLog_("No MS/MS spectra found in any of the mzXML files. Aborting!");
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
				
				ofstream batchfile(String(temp_data_directory + batch_filename).c_str());
				if ( !batchfile )
				{
					throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, temp_data_directory + batch_filename);
				}
				String call = "rdesktop";
				if ( !user.empty() ) call.append(" -u " + user);
				if ( !password.empty() ) call.append(" -p \"" + password + "\"");
				call.append(" -s cmd\\ /C\\ \"");
				call.append("net use " + temp_data_directory_win.substr(0,2) + " \\\\" + temp_data_directory_network + " && ");
// 				call.append("net use " + temp_data_directory_win.substr(0,2) + " " + temp_data_directory_network + " && ");
				
				batchfile << String(" cd " + temp_data_directory_win + " && " + temp_data_directory_win.substr(0,2));
				
				for ( PointerSizeUInt i = 0; i <= (PointerSizeUInt) (dtas / max_dtas_per_run); ++i )
				{
					batchfile << String(" && " + sequest_directory_win + "sequest.exe -P" + input_file_directory_network + " " + temp_data_directory_network + "\\*.dta_" + String(i) + " >  " +  temp_data_directory_network +"\\" + sequest_screen_output + " && move sequest.log sequest.log" + String(i));
				}
				batchfile << String(" && " + sequest_directory_win.substr(0,2) + " &&");
				batchfile << String(" net use /delete " + temp_data_directory_win.substr(0,2));
				batchfile << " && logoff";
				batchfile.close();
				batchfile.clear();
				
				call.append(temp_data_directory_win + batch_filename + "\" " + sequest_computer);
				writeLog_("System call: " + call);
				int status = system(call.c_str());
				remove(sequest_screen_output.c_str());
				
				if ( status != 0 )
				{
					writeLog_("Sequest problem. Aborting! (Details can be seen in the logfile: \"" + logfile + "\")");
					for ( vector< String >::const_iterator tmp_names_i = tmp_names.begin(); tmp_names_i != tmp_names.end(); ++tmp_names_i ) remove(tmp_names_i->c_str());
					
					// remove all dtas
					for ( PointerSizeUInt i = 0; i <= (PointerSizeUInt) (dtas / max_dtas_per_run); ++i )
					{
						writeLog_("removing dta files");
						vector<String> to_delete;
						if (File::fileList(temp_data_directory,String("*.dta_") + i,to_delete))
						{
							for (vector<String>::const_iterator it = to_delete.begin(); it != to_delete.end(); ++it)
							{
								if ( !File::remove(temp_data_directory + *it) )
								{
									writeLog_(String("'") + temp_data_directory + *it + "' could not be removed!");
								}
							}
						}
					}
					return EXTERNAL_PROGRAM_ERROR;
				}
				
				for ( PointerSizeUInt i = 0; i <= (PointerSizeUInt) (dtas / max_dtas_per_run); ++i )
				{
					ifstream sequest_log(string(temp_data_directory + "sequest.log" + String(i)).c_str()); // write sequest log to logfile
					if ( !sequest_log )
					{
						writeLog_("No Sequest log found!");
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
						writeLog_(buffer);
						delete(buffer);
						remove(string(temp_data_directory + "sequest.log" + String(i)).c_str());
					}
				}
			}
			
			if ( sequest_out )
			{
				// remove all dtas
				if ( !keep_dta_files ) 
				{
					for ( PointerSizeUInt i = 0; i <= (PointerSizeUInt) (dtas / max_dtas_per_run); ++i )
					{
						vector<String> to_delete;
						if (File::fileList(temp_data_directory,String("*.dta_") + i,to_delete))
						{
							for (vector<String>::const_iterator it = to_delete.begin(); it != to_delete.end(); ++it)
							{
								if ( !File::remove(temp_data_directory + *it) )
								{
									writeLog_(String("'") + temp_data_directory + *it + "' could not be removed!");
								}
							}
						}
					}
				}
				
				SequestOutfile sequest_outfile;
				vector< IdentificationData > identifications;
				ProteinIdentification protein_identification;
				UnsignedInt identification_size = identifications.size();

				vector<String> out_files;
				if (!File::fileList(out_directory,String("*.out"),out_files))
				{
					writeLog_(String("Error: No .out files found in '") + out_directory + "'. Aborting!");
					for ( vector< String >::const_iterator tmp_names_i = tmp_names.begin(); tmp_names_i != tmp_names.end(); ++tmp_names_i ) remove(tmp_names_i->c_str());
					
					return UNKNOWN_ERROR;
				}
				
				String call;
				map< String, vector< Real > > filenames_and_pvalues;
				
				if ( !peptide_prophet_directory.empty() )
				{
					String summary = "tmp.summary.html";
					bool append = false;
					for ( vector<String>::const_iterator i = out_files.begin(); i != out_files.end(); ++i )
					{
						sequest_outfile.out2SummaryHtml(out_directory + *i, summary, database, append);
					}
					//sequest_outfile.finishSummaryHtml(summary);
					
					call = "cd ";
					call.append(peptide_prophet_directory);
					call.append("runPeptidProphet ");
					call.append(summary);
					SignedInt status = system(call.c_str());
					if ( status != 0 )
					{
						writeLog_("Problems with Peptide Prophet. Aborting!");
						for ( vector< String >::const_iterator tmp_names_i = tmp_names.begin(); tmp_names_i != tmp_names.end(); ++tmp_names_i ) remove(tmp_names_i->c_str());
						
						return UNKNOWN_ERROR;
					}
					filenames_and_pvalues = sequest_outfile.getPeptidePValues(temp_data_directory, summary + ".prob");
					
					remove(summary.c_str());
					remove(String(summary + ".prob").c_str());
					remove(String(summary + ".esi").c_str());
					summary.replace(summary.length() - 4, 4, "model");
					remove(summary.c_str());
				}
				String filename;
				vector< Real > pvalues;
				
				for ( vector<String>::const_iterator i = out_files.begin(); i != out_files.end(); ++i )
				{
					filename = out_directory + *i;
					if ( filenames_and_pvalues.empty() )
					{
						sequest_outfile.load(filename, identifications, protein_identification, p_value, pvalues, database);
					}
					else
					{
						sequest_outfile.load(filename, identifications, protein_identification, p_value, filenames_and_pvalues[out_directory + *i], database);
					}
					
					// save the retention times
					if ( identification_size != identifications.size() )
					{
						identification_size = identifications.size();
						identifications.back().rt = filenames_and_precursor_retention_times[filename];
						//cout << "LFF: " << filename << endl;
					}
				}
				
				vector< ProteinIdentification > pis;
				pis.push_back(protein_identification);
				
				AnalysisXMLFile().store(output_filename, pis, identifications);
				
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
			
			// (3.3) deleting all temporary files
			for ( vector< String >::const_iterator tmp_names_i = tmp_names.begin(); tmp_names_i != tmp_names.end(); ++tmp_names_i ) remove(tmp_names_i->c_str());
			
			return EXECUTION_OK;
		}
};

//@endcond



int main( int argc, char ** argv )
{
	TOPPSequestAdapter tool;

	return tool.main(argc,argv);
}
