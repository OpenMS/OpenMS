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
// $Maintainer: Andreas Bertsch $
// $Authors: Martin Langwisch $
// --------------------------------------------------------------------------


#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/SequestInfile.h>
#include <OpenMS/FORMAT/SequestOutfile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/ContactPerson.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/PTMXMLFile.h>


#include <cstdlib>
#include <vector>
#include <algorithm>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------


/**
	@page TOPP_SequestAdapter SequestAdapter

	@brief Identifies peptides in MS/MS spectra via Sequest.

	@experimental This tool has not been tested thoroughly and might behave not as expected!

	This wrapper application serves for getting peptide peptide_identifications
	for MS/MS spectra. The wrapper can be executed in three different
	modes:
	<ol>
				<li>
				The whole process of ProteinIdentification via Sequest is executed.
				Inputfile is one (or more) mz file containing the MS/MS spectra
				(Supported spectrum file formats are .mzXML, .mzData)
				for which the identifications are to be found and one database in
				FASTA format containing the possible proteins.
				The results are written as an IdXML output file. This mode is selected
			 	by default.
				Note: You need a user with network access on the computer hosting sequest.
			 	</li>

				<li>
				Only the first part of the ProteinIdentification process is performed.
				This means that a Sequest input file is generated and dta files are
				created from the mz file.
				Calling a Sequest process should look like the following:

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

	@todo Check for missing precursors (Andreas)
	
	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_SequestAdapter.cli
*/

// We do not want this class to show up in the docu -> cond
// @cond

class TOPPSequestAdapter
	: public TOPPBase
{
	public:
		TOPPSequestAdapter()
			: TOPPBase("SequestAdapter", "Annotates MS/MS spectra using Sequest.")
		{
		}

	protected:
		static const Int max_peptide_mass_units = 2;
		static const Size max_dtas_per_run = 1000; // sequest has a problem when there are too many dtas, so they have to be splitted, 1000 seemed to work very good
		Size dtas;

		void registerOptionsAndFlags_()
		{
			addText_("The definitions for the parameters are taken from the site:\n"
										 "http://www.grosse-coosmann.de/~florian/Parameters.html#file.");
			// do not change this to registerInputFile_() as it might also be a directory, which fails the property check of a file on Windows
			registerStringOption_("in", "<file>", "", "input file(s) in mzXML or mzData format (comma-separated).\n"
					 																			"Note: In mode 'sequest_out' a directory with Sequest results files\n"
																								"(*.out) is read", false);
			registerOutputFile_("out", "<file>", "", "output file in IdXML format.\n"
			                                           "Note: In mode 'sequest_in' a Sequest input file is written.", false);
			registerFlag_("sequest_in", "if this flag is set the SequestAdapter will read in mzXML or mzData\n"
																								"and write an Sequest input file\n"
																								"and create dta files from the given mzXML or mzData files");
			registerFlag_("sequest_out", "if this flag is set the SequestAdapter will read in Sequest result files\n"
																									"and write IdXML");
			registerStringOption_("mz_files", "<files>", "", "when using sequest_out the mzXML or mzData files (comma-separated)\n"
																																						"have to be given to retrieve the retention times", false);
			registerFlag_("show_enzymes", "show a list with enzymes and corresponding numbers to choose from");
			registerStringOption_("sequest_computer", "<name>", "", "the name of the computer in the network that hosts Sequest\n"
																															"(rdesktop is used to connect to this computer)", false);
			registerStringOption_("sequest_directory_win", "<dir>", "", "the windows directory in which Sequest (sequest.exe) is located", false);
			registerStringOption_("user", "<name>", "", "user name for the sequest computer (has to have access to network!)", false);
			registerStringOption_("password", "<pw>", "", "password for this user (if not given, you have to enter it at prompt)", false);
			registerStringOption_("temp_data_directory", "<dir>", "", "a directory in which some temporary files can be stored", false);
			registerStringOption_("temp_data_directory_win", "<dir>", "", "windows path of the temporary data directory,\n"
																																																			"e.g. X:\\temp_data_dir", false);
			registerStringOption_("db", "<file>", "", "name of FASTA-database to search in", false);
			registerInputFile_("sequest_input", "<file>", "", "name for the input file of Sequest (may only be used in a full run)", false);
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
																																														"<composition>,<residues>,<type>,<name>, e.g.: H2C2O,KCS,opt,Acetyl or\n"
																																														"<mass>,<residues>,<type>,<name>, e.g.: 42.0367,KCS,opt,Acetyl or\n"
																																														"Valid values for type are \"fix\" and \"opt\" (default)\n"
																																														"If you want terminal PTMs, write \"cterm\", \"nterm\", \"cterm_prot\" or \"nterm_prot\" instead of residues", false);
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
			registerFlag_("keep_out_files", "If set the Sequest .out-files are not removed");
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

		bool correctNetworkPath(String& network_path, Size backslashes = 2)
		{
			String::size_type pos(0);
			while ( (pos < network_path.length()) && (network_path[pos] == '\\') ) ++pos;
			if ( pos < backslashes ) network_path.insert(network_path.begin(), backslashes-pos, '\\');
			else network_path.erase(0, pos-backslashes);
			if ( network_path.length() < backslashes+1 ) return false;
			if ( network_path[network_path.length() - 1] != '\\' ) network_path.append("\\"); // if it doesn't end with a slash, append one
			return true;
		}

  Size
		MSExperiment2DTAs(
			MSExperiment<Peak1D>& msexperiment,
			const String& common_name,
			const vector< Int >& charges,
			map< String, DoubleReal >& outfile_names_and_precursor_retention_times,
			vector< String >& dta_filenames,
			bool make_dtas = true)
		{
			DTAFile dtafile;
			String filename;
			Size scan_number(0);
			Size msms_spectra(0);
			Size dtas(0);

			for ( MSExperiment<Peak1D>::Iterator spectra_it = msexperiment.begin(); spectra_it != msexperiment.end(); ++spectra_it )
			{
				++scan_number;
				if ( (spectra_it->getMSLevel() == 2) && (!spectra_it->empty()) )
				{
					++msms_spectra;
					if ( spectra_it->getPrecursors()[0].getCharge() )
					{
						filename = common_name + "." + String(scan_number) + "." + String(spectra_it->getPrecursors()[0].getCharge()) + ".dta_" + String( (dtas / max_dtas_per_run) );
						++dtas;
						if ( make_dtas ) dtafile.store(filename, *spectra_it);
						dta_filenames.push_back(filename);
						filename = File::basename(filename);
						filename.replace(filename.length() - 4, 4, ".out");
						outfile_names_and_precursor_retention_times[filename] = spectra_it->getRT();
					}
					else
					{
						for ( vector< Int >::const_iterator charges_it = charges.begin(); charges_it != charges.end(); ++charges_it )
						{
							filename = common_name + "." + String(scan_number) + "." + *charges_it + ".dta_" + String( (dtas / max_dtas_per_run) );
							++dtas;
							if ( make_dtas )
							{
								spectra_it->getPrecursors()[0].setCharge(*charges_it);
								dtafile.store(filename, *spectra_it);
							}
							dta_filenames.push_back(filename);
							filename = File::basename(filename);
							filename.replace(filename.length() - 4, 4, ".out");
							outfile_names_and_precursor_retention_times[filename] = spectra_it->getRT();
						}
						spectra_it->getPrecursors()[0].setCharge(0);
					}
				}
			}

			return msms_spectra;
		}

		ExitCodes main_(int , const char**)
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
				modifications_filename,
				dta_files_common_name,
				basename;

			bool
				sequest_in(false),
				sequest_out(false),
				keep_out_files(false),
				keep_dta_files(false),
				monoisotopic(false),
				make_dtas(false);

			vector< String >
				substrings,
				substrings2,
				spectra;


			Size
				msms_spectra_in_file(0),
				msms_spectra_altogether(0);

			vector< Int > charges;

			DoubleReal
				Real_buffer(0.0),
				Real_buffer2(0.0),
				p_value(1.0);

			Int int_buffer(0);

			ContactPerson contact_person;
			ExitCodes exit_code = EXECUTION_OK;
			FileHandler fh;
			FileTypes::Type type;
			MSExperiment<Peak1D> msexperiment;
			vector<PeptideIdentification> peptide_identifications;
			vector<ProteinIdentification> pis;
			ProteinIdentification protein_identification;
			StringList out_files;

			// the outfile-names and their retention_times
			map< String, DoubleReal > outfile_names_and_precursor_retention_times;

			// the names of the dta_files - used to erase them afterwards
	 		vector< String > dta_filenames;

			// filename and tag: file has to: 1 - exist  2 - be readable  4 - writable  8 - be deleted afterwards
			map< String, Size > files;
			Size const
				exist(1),
				readable(2),
				writable(4),
				delete_afterwards(8);

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
				map< String, pair< String, String > > PTM_informations;
				try
				{
					PTMXMLFile().load(modifications_filename, PTM_informations);
				}
				catch ( Exception::ParseError pe )
				{
					writeLog_(pe.getMessage());
					return PARSE_ERROR;
				}

				// output the information
				stringstream PTM_info;
				String::size_type max_name_length(4), max_composition_length(11), max_amino_acids_length(11);
				for ( map< String, pair< String, String > >::const_iterator mod_it = PTM_informations.begin(); mod_it != PTM_informations.end(); ++mod_it )
				{
					max_name_length = max(max_name_length, mod_it->first.length());
					max_composition_length = max(max_composition_length, mod_it->second.first.length());
					max_amino_acids_length = max(max_amino_acids_length, mod_it->second.second.length());
				}
				PTM_info << "name" << String(max_name_length - 4, ' ') << "\t" << "composition" << String(max_composition_length - 11, ' ') << "\t" << "amino_acids" << String(max_amino_acids_length - 11, ' ') << endl;
				for ( map< String, pair< String, String > >::const_iterator mod_it = PTM_informations.begin(); mod_it != PTM_informations.end(); ++mod_it )
				{
					PTM_info << mod_it->first << String(max_name_length - mod_it->first.length(), ' ') << "\t" << mod_it->second.first << String(max_composition_length - mod_it->second.first.length(), ' ') << "\t" << mod_it->second.second << String(max_amino_acids_length - mod_it->second.second.length(), ' ') << endl;
				}
				std::cout << PTM_info.str() << std::endl;

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
				files[logfile] = (writable | delete_afterwards);
			}
			else files[logfile] = writable;

			string_buffer = getStringOption_("charges");
			if ( string_buffer.empty() )
			{
				writeLog_("No charge states given. Aborting!");
				return ILLEGAL_PARAMETERS;
			}
			else
			{
				Int range_start(-1), range_end(-1);
				string_buffer.split(',', substrings);

				for ( vector< String >::iterator substrings_it = substrings.begin(); substrings_it != substrings.end(); )
				{
					if ( substrings_it->empty() ) substrings.erase(substrings_it);
					else
					{
						substrings_it->split('>', substrings2);
						if ( substrings2.size() < 2 ) // only one number, no range
						{
							if ( (*substrings_it)[substrings_it->length()-1] == '-' ) charges.push_back(-1 * substrings_it->toInt());
							else charges.push_back(substrings_it->toInt());
						}
						else // range of charge states
						{
							if ( substrings2.size() > 2 )
							{
								writeLog_("Illegal range of charge states given: " + *substrings_it + ". Aborting!");
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

						++substrings_it;
					}
				}

				if ( charges.empty() )
				{
					writeLog_("No charges states given. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				sort(charges.begin(), charges.end());
				for ( vector< Int >::iterator charges_it = charges.begin(); charges_it != --charges.end(); )
				{
					if ( (*charges_it) == (*(charges_it+1)) ) charges.erase(charges_it+1);
					else ++charges_it;
				}
			}

			temp_data_directory = getStringOption_("temp_data_directory");
			if ( temp_data_directory.empty() )
			{
				writeLog_("No directory for temporary files given. Aborting!");
				return ILLEGAL_PARAMETERS;
			}
			temp_data_directory = File::absolutePath(temp_data_directory);
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
					out_directory = temp_data_directory;
				}
				else // if only sequest_out is set, in is the out_directory
				{
					out_directory = string_buffer;
					out_directory = File::absolutePath(out_directory);
					out_directory.ensureLastChar('/');

					// if only sequest_out is set, the mz files have to be given to retrieve the retention times
					string_buffer = getStringOption_("mz_files");
					if ( string_buffer.empty() )
					{
						writeLog_("No mz files specified. Aborting!");
						return ILLEGAL_PARAMETERS;
					}
					else
					{
						string_buffer.split(',', spectra);
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
				files[database] = readable;

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
						files[input_filename] = (writable | delete_afterwards);
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
						files[input_filename] = readable;
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

			// the batchfile will have to be in the temporary data directory, otherwise it's hard to connect to windows and use the other files
			batch_filename = "sequest_run.bat";
			files[temp_data_directory + batch_filename] = (writable | delete_afterwards);

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
					vector< String > enzyme_info;
					for ( vector< String >::iterator einfo_it = substrings.begin(); einfo_it != substrings.end(); ++ einfo_it )
					{
						einfo_it->split(',', enzyme_info);
						if ( (enzyme_info.size() < 3) || (enzyme_info.size() > 4) )
						{
							writeLog_("Illegal number of entries for enzyme (not in [3,4]). Aborting!");
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
					SignedSize highest_enzyme_number = sequest_infile.setEnzyme(getStringOption_("cleavage"));
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
					for ( vector< String >::iterator substrings_it = substrings.begin(); substrings_it != substrings.end(); ++substrings_it )
					{
						// the values are expected to be DoubleReal, otherwise they will be seen as 0!
						Real_buffer = String(substrings_it->c_str()).toDouble();
						if ( (Real_buffer < 0) || (Real_buffer > 1) )
						{
							writeLog_("Illegal weights for ion series given. Aborting!");

							return ILLEGAL_PARAMETERS;
						}
						(*substrings_it) = String(Real_buffer);
					}
					string_buffer.concatenate(substrings.begin(), substrings.end(), " ");
					sequest_infile.setIonSeriesWeights(string_buffer);
				}
				// modifications
				string_buffer = getStringOption_("modifications");
				monoisotopic = getFlag_("use_monoisotopic_mod_mass");
				try
				{
					sequest_infile.handlePTMs(string_buffer, modifications_filename, monoisotopic);
				}
				catch ( Exception::FileNotFound fnf_e )
				{
					writeLog_("No modifications XML file given. Aborting!");
					return INPUT_FILE_NOT_FOUND;
				}
				catch ( Exception::FileNotReadable fnr_e )
				{
					writeLog_("Modifications XML file is not readable. Aborting!");
					return INPUT_FILE_NOT_READABLE;
				}
				catch ( Exception::ParseError p_e )
				{
					writeLog_(String(p_e.getMessage()) + ". Aborting!");
					return PARSE_ERROR;
				}
				string_buffer = getStringOption_("partial_sequence");
				string_buffer.substitute(',', ' ');
				sequest_infile.setPartialSequence(string_buffer);

				string_buffer = getStringOption_("header_filter");
				string_buffer.substitute(',', ' ');
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
				files[output_filename] = writable;

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
			bool existed(false);
			Size file_tag(0);

			for ( map< String, Size >::const_iterator files_it = files.begin(); files_it != files.end(); ++files_it )
			{
				string_buffer = files_it->first;
				file_tag = files_it->second;

				if ( (file_tag & exist || file_tag & readable) && !File::exists(string_buffer) )
				{
					exit_code = INPUT_FILE_NOT_FOUND;
					writeLog_(String("File ")+ string_buffer + " does not exist. Aborting!");
					break;
				}

				if ( (file_tag & readable) && !File::readable(string_buffer) )
				{
					exit_code = INPUT_FILE_NOT_READABLE;
					writeLog_(String("File ")+ string_buffer + " is not readable. Aborting!");
					break;
				}

				existed = File::exists(string_buffer);
				if ( (file_tag & writable) && !File::writable(string_buffer) )
				{
					exit_code = CANNOT_WRITE_OUTPUT_FILE;
					writeLog_(String("Cannot write file ")+ string_buffer + ". Aborting!");
					break;
				}
				else if ( !existed ) remove(string_buffer.c_str());
				existed = false;
			}
			// creating the input file
			if ( exit_code == EXECUTION_OK && sequest_in )
			{
				sequest_infile.store(input_filename);
			}

			if ( exit_code == EXECUTION_OK )
			{
				// check the Mz files, get the names for the dtas and check whether they do no already exist
				bool make_dtas = ( sequest_out && !sequest_in ) ? false : true; // if only sequest_out is set, just get the retention times
				// creating the dta files
				if ( make_dtas ) writeLog_("creating dta files");
				// first get the dta names
				for ( vector< String >::iterator spectra_it = spectra.begin(); spectra_it != spectra.end(); ++spectra_it )
				{
					*spectra_it = File::absolutePath(*spectra_it);
					type = fh.getTypeByContent(*spectra_it);
					if ( type == FileTypes::UNKNOWN )
					{
						writeLog_("Could not determine type of the file. Aborting!");
						exit_code = PARSE_ERROR;
						break;
					}
					fh.loadExperiment(*spectra_it, msexperiment, type);

					msms_spectra_in_file = MSExperiment2DTAs(msexperiment, temp_data_directory + File::basename(*spectra_it), charges, outfile_names_and_precursor_retention_times, dta_filenames, false);

					msms_spectra_altogether += msms_spectra_in_file;

					// if make_dtas is set, check whether one of them does already exist, if so, stop the adapter
					if ( make_dtas )
					{
						for ( vector< String >::const_iterator dta_names_it = dta_filenames.begin(); dta_names_it != dta_filenames.end(); ++dta_names_it )
						{
							string_buffer = temp_data_directory + *dta_names_it;
							if ( File::exists(string_buffer) )
							{
								writeLog_("The file " + string_buffer + " does already exist in directory " + temp_data_directory + ". Please remove it first. Aborting!");
								// deleting all temporary files
								for ( map< String, Size >::const_iterator files_it = files.begin(); files_it != files.end(); ++files_it )
								{
									if ( files_it->second & delete_afterwards ) remove(files_it->first.c_str());
								}
								exit_code = UNKNOWN_ERROR;
								break;
							}
						}
					}
				}
			}

			// if no msms spectra were found
			if ( exit_code == EXECUTION_OK && !msms_spectra_altogether )
			{
				writeLog_("No MS/MS spectra found in any of the mz files. Aborting!");
				exit_code = UNKNOWN_ERROR;
			}

			// if make_dtas is set and non of the dta files did already exist, create them
			if ( exit_code == EXECUTION_OK && make_dtas )
			{
				for ( vector< String >::const_iterator spectra_it = spectra.begin(); spectra_it != spectra.end(); ++spectra_it )
				{
					type = fh.getTypeByContent(*spectra_it);
					if ( type == FileTypes::UNKNOWN )
					{
						writeLog_("Could not determine type of the file. Aborting!");
						exit_code = PARSE_ERROR;
					}
					fh.loadExperiment(*spectra_it, msexperiment, type);
					basename = File::basename(*spectra_it);
					dta_files_common_name = temp_data_directory + basename;
					msms_spectra_in_file = MSExperiment2DTAs(msexperiment, dta_files_common_name, charges, outfile_names_and_precursor_retention_times, dta_filenames, make_dtas);
					writeLog_(String(msms_spectra_in_file) + " MS/MS spectra in file " + *spectra_it);
				}
			}

			// (3.2.3) running the program
			if ( exit_code == EXECUTION_OK )
			{
				if ( sequest_in && sequest_out )
				{
					// creating a batch file for windows (command doesn't accept commands that are longer than 256 chars)
					String sequest_screen_output; // direct the screen-output to a file
					do
					{
						sequest_screen_output = String::random(10);
					}
					while ( File::exists(sequest_screen_output) );
					files[temp_data_directory + sequest_screen_output] = (writable | delete_afterwards);

					ofstream batchfile(String(temp_data_directory + batch_filename).c_str());
// 					if ( !batchfile )
// 					{
// 						throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, temp_data_directory + batch_filename);
// 					}
					String call = "rdesktop";
					if ( !user.empty() ) call.append(" -u " + user);
					if ( !password.empty() ) call.append(" -p \"" + password + "\"");
					call.append(" -s cmd\\ /K\\ \"");
	// 				call.append("echo net use " + temp_data_directory_win.substr(0,2) + " \\\\" + temp_data_directory_network + " && ");
					call.append("net use " + temp_data_directory_win.substr(0,2) + " \\\\" + temp_data_directory_network.substr(0, temp_data_directory_network.length() - 1) + " && ");
	// 				call.append(" net use " + temp_data_directory_win.substr(0,2) + " " + temp_data_directory_network + " && ");

					batchfile << String(" cd " + temp_data_directory_win + " && " + temp_data_directory_win.substr(0,2));

					for ( Size i = 0; i <= Size(dtas / max_dtas_per_run); ++i )
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
						exit_code = EXTERNAL_PROGRAM_ERROR;
					}
					else
					{
						bool no_log(false);
						string_buffer.clear();
						for ( Size i = 0; i <= (Size) (dtas / max_dtas_per_run); ++i )
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
							exit_code = EXTERNAL_PROGRAM_ERROR;
						}
						else writeLog_(string_buffer);
					}
				}
			}

			if ( sequest_out )
			{
				if ( exit_code == EXECUTION_OK )
				{
					SequestOutfile sequest_outfile;
					if (!File::fileList(out_directory, String("*.out"), out_files))
					{
						writeLog_(String("Error: No .out files found in '") + out_directory + "'. Aborting!");
						exit_code = UNKNOWN_ERROR;
					}
				}
				if ( exit_code == EXECUTION_OK )
				{
					vector< pair < String, vector< DoubleReal > > > filenames_and_pvalues;
					for ( StringList::iterator out_files_it = out_files.begin(); out_files_it != out_files.end(); ++out_files_it )
					{
						filenames_and_pvalues.push_back(make_pair(out_directory + *out_files_it, vector< DoubleReal >()));
					}
	// 				sequest_outfile.getPValuesFromOutFiles(filenames_and_pvalues);

					// set the parameters
					ProteinIdentification::SearchParameters sp;
					sp.db = "Fasta";
					sp.taxonomy = sequest_infile.getSequenceHeaderFilter();
					if ( monoisotopic ) sp.mass_type = ProteinIdentification::MONOISOTOPIC;
					else sp.mass_type = ProteinIdentification::AVERAGE;
					for ( vector< Int >::const_iterator charges_it = charges.begin(); charges_it != charges.end(); ++charges_it )
					{
						if ( *charges_it > 0 ) sp.charges.append("+");
						sp.charges.append(String(*charges_it));
					}
					if ( sequest_infile.getEnzymeName() == "Trypsin" ) sp.enzyme = ProteinIdentification::TRYPSIN;
					else if ( sequest_infile.getEnzymeName() == "No_Enzyme" ) sp.enzyme = ProteinIdentification::NO_ENZYME;
					else sp.enzyme = ProteinIdentification::UNKNOWN_ENZYME;
					sp.peak_mass_tolerance = sequest_infile.getPeakMassTolerance();
					sp.precursor_tolerance = sequest_infile.getPrecursorMassTolerance();
					protein_identification.setSearchParameters(sp);

					Size peptide_identification_size = peptide_identifications.size();
					for ( vector< pair < String, vector< DoubleReal > > >::iterator filenames_and_pvalues_it = filenames_and_pvalues.begin(); filenames_and_pvalues_it != filenames_and_pvalues.end(); ++filenames_and_pvalues_it )
					{
						try
						{
							sequest_outfile.load(filenames_and_pvalues_it->first, peptide_identifications, protein_identification, p_value, filenames_and_pvalues_it->second, database);
						}
						catch( Exception::ParseError pe )
						{
							writeLog_(pe.getMessage());
							exit_code = INPUT_FILE_CORRUPT;
							break;
						}

						// save the retention times if peptides have been identified to the p-level
						if ( peptide_identification_size != peptide_identifications.size() )
						{
							peptide_identification_size = peptide_identifications.size();
							string_buffer = File::basename(filenames_and_pvalues_it->first);
							if ( outfile_names_and_precursor_retention_times.find(string_buffer) != outfile_names_and_precursor_retention_times.end() ) peptide_identifications.back().setMetaValue("RT",  outfile_names_and_precursor_retention_times[string_buffer]);
							else peptide_identifications.back().setMetaValue("RT", 0);
						}
					}
				}
				if ( exit_code == EXECUTION_OK )
				{
					pis.push_back(protein_identification);

					IdXMLFile().store(output_filename, pis, peptide_identifications);

					// remove all outs
					if ( !keep_out_files )
					{
						writeLog_("removing out files");
						for ( StringList::const_iterator out_files_it = out_files.begin(); out_files_it != out_files.end(); ++out_files_it )
						{
							if ( !File::remove(out_directory + *out_files_it) )
							{
								writeLog_(String("'") + out_directory + *out_files_it + "' could not be removed!");
							}
						}
					}
				}
			}

			if ( exit_code == EXTERNAL_PROGRAM_ERROR )
			{
				writeLog_("Sequest problem. Aborting! (Details can be seen in the logfile: \"" + logfile + "\")");
				files[logfile] = readable;
			}

			// deleting all temporary files
			writeLog_("removing temporary files");
			for ( map< String, Size >::const_iterator files_it = files.begin(); files_it != files.end(); ++files_it )
			{
				if ( files_it->second & delete_afterwards ) remove(files_it->first.c_str());
			}
			// remove all dtas
			if ( !keep_dta_files )
			{
				writeLog_("removing dta files");
				for ( vector< String >::const_iterator dta_names_it = dta_filenames.begin(); dta_names_it != dta_filenames.end(); ++dta_names_it )
				{
					if ( !File::remove(*dta_names_it) ) writeLog_("'" + string_buffer + "' could not be removed!");
				}
			}

			return exit_code;
		}
};

//@endcond



int main( int argc, const char** argv )
{
	TOPPSequestAdapter tool;

	return tool.main(argc,argv);
}
