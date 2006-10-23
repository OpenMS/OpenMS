// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: SequestAdapter.C, 2006/10/19 13:16 martinlangwisch Exp $
// $Author: martinlangwisch $
// $Maintainer: Martin Langwisch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/FORMAT/AnalysisXMLFile.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Date.h>
#include <OpenMS/METADATA/ContactPerson.h>
#include <OpenMS/FORMAT/SequestInfile.h>
#include <OpenMS/FORMAT/SequestOutfile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <stdlib.h>
#include <vector>
#include <algorithm>

#include <qfileinfo.h>
#include <qdir.h>
#include <qstringlist.h>

using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------


/**
	@page SequestAdapter SequestAdapter
	
	@brief Identifies peptides in MS/MS spectra via Sequest.
	
	This wrapper component serves for getting peptide identifications
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
				
				<ul>	
					<li>
						@code sequest -P\<inputfilename\> \<path to dta files\>*.dta  @endcode
					</li>
				</ul>
				Consult your Sequest reference manual for further details.
				
				This mode is selected by the <b>-Sequest_in</b> option in the command line.
				</li>
				
				<li>
				Only the second part of the identification process is performed.
				This means that the output of sequest is translated into analysisXML.
				
				This mode is selected by the <b>-Sequest_out</b> option in the command line.
				</li>
	</ol>
	
	@todo look for possible crash codes of sequest and catching them
	
	@ingroup TOPP
*/

// We do not want this class to show up in the docu -> @cond
// @cond 

class TOPPSequestAdapter
	: public TOPPBase
{
	public:
		TOPPSequestAdapter()
			: TOPPBase("SequestAdapter")
		{}
	
	protected:
		static const int max_peptide_mass_units = 2;
		
		struct SortRetentionTimes:
				public std::binary_function<const std::pair< String, std::vector< double > >, const std::pair< String, std::vector< double > >, bool>
		{
			bool operator()(const std::pair< String, std::vector< double > >& x, const std::pair< String, std::vector< double > >& y)
			{
				return ( x.first < y.first );
			}
		};
	
		void printToolUsage_()
		{
			std::cerr	<< std::endl
						<< tool_name_ << " -- annotates MS/MS spectra using Sequest" << std::endl
						<< std::endl
						<< "Usage:" << std::endl
						<< " " << tool_name_ << " [options]" << std::endl
						<< std::endl
						<< "Options are:" << std::endl
						<< "  [the _win parameters correspond to path of the linux directory when mounted under windows. NO SPACE ALLOWED!;" << std::endl
						<< "   the _network paramters correspond to the path network path of the directory when mounting the under windows;" << std::endl
						<< "   as Sequest runs on windows, all the files and directories used have to be mounted on the corresponding computer!" << std::endl
						<< "   rdesktop is used to connect to this computer" << std::endl << std::endl
						<< "   Sequest writes a file named 'sequest.log' into the out_dir, so you must not name the log file alike]" << std::endl
						<< "  -Sequest_in           if this flag is set the SequestAdapter will create a sequest input file" << std::endl
						<< "                        and create dta files from the given mzXML files" << std::endl
						<< "  -Sequest_out          if this flag is set the SequestAdapter will read in Sequest out files and write an analysis XML file" << std::endl
						<< "  -spectra              the names of the mzXML files" << std::endl
						<< "  -out                  the name of the analysis XML file" << std::endl
						<< "  -sequest_dir_win      the windows path where sequest.exe is located" << std::endl
						<< "  -sequest_computer     the name of the computer in the network that hosts Sequest" << std::endl
						<< "  -user                 user name for the sequest computer (has to have access to network!)" << std::endl
						<< "  -password             password for this user (if not given, you have to enter it at promt)" << std::endl
						<< "  -p_value              annotations with inferior p-value are ignored (default is 0.05)" << std::endl
						<< "  -show_enzyme_numbers  show a list with enzymes and corresponding numbers to choose from" << std::endl
						<< "  -num_results          the maximal number of results (peptides) to show" << std::endl
						<< "  -max_num_dif_AA_per_mod  limits the maximum total number of each single variable modification in one peptide" << std::endl
						<< "  -max_num_dif_mods_per_peptide  limits the maximum total number of each single variable modification in one peptide" << std::endl
						<< "  -prob_charge          the number of charge states that are used if it is unknown for a scan" << std::endl
						<< "  -pep_mass_tol         tolerance for a peptide ion" << std::endl
						<< "  -frag_ion_tol         tolerance for a fragment ion" << std::endl
						<< "  -match_peak_tol       the minimal space between two peaks" << std::endl
						<< "  -enzyme_info          <name>,<cut direction: N to C?>,<cuts after>,<doesn't cut before>;" << std::endl
						<< "                        cuts after, doesn't cut before: amino acids in 1-letter code or '-' for unspecific cleavage" << std::endl
						<< "  -enzyme_number        a number from the list (show_enzyme_numbers); if enzyme_info is used, this value is set accordingly" << std::endl
						<< "  -neutral_loss_ABY     ABY: 0 or 1 whether neutral losses of the series should be honored, eg: 011" << std::endl
						<< "  -ion_series_weights   abcdvwxyz: [0.0, 1.0] factor for the series, eg: 0,0.5,0,0,0,0,0,1.0,0" << std::endl
						<< "  -dta_dir              the directory to store the dta files" << std::endl
						<< "  -dta_dir_win          " << std::endl
						<< "  -out_dir              the directory to store the sequest output files" << std::endl
						<< "  -out_dir_win          " << std::endl
						<< "  -in                   the name of the sequest input file" << std::endl
						<< "  -in_win               " << std::endl
						<< "  -db                   the name of the database file" << std::endl
						<< "  -db_win               " << std::endl
						<< "  -snd_db               the name of the second database file" << std::endl
						<< "  -snd_db_win           " << std::endl
						<< "  -temp_data_dir        the directory to store temporary data" << std::endl
						<< "  -temp_data_dir_win    " << std::endl
						<< std::endl
						<< "  For each windows drive, one corresponding network drive has to be given, so maybe you don't need to set all the parameters below" << std::endl
						<< "  -temp_data_dir_network" << std::endl
						<< "  -out_dir_network" << std::endl
						<< "  -dta_dir_network" << std::endl
						<< "  -db_dir_network" << std::endl
						<< "  -snd_db_dir_network" << std::endl
						<< "  -in_dir_network" << std::endl;
		}


		void printToolHelpOpt_()
		{
			std::cerr	<< std::endl
			<< "  -ion_cutoff                    This value selects a cut-off below which a matching peptide is rejected." << std::endl
			<< "                                 The value compared with this value is the ratio" << std::endl
			<< "                                 (# matching theoretical fragment peaks) / (# total theoretical fragment peaks)" << std::endl
			<< "                                 which means that the user can select a minimum coverage of matching peaks." << std::endl
			<< "  -pep_mass_unit                 peptide mass unit: 0=amu (atomic mass unit), 1=mmu (millimass unit), 2=ppm (parts per million)" << std::endl
			<< "  -min_prot_mass                 minimal protein mass" << std::endl
			<< "  -max_prot_mass                 maximal protein mass" << std::endl
			<< "  -nuc_reading_frame             Format of the FASTA database:" << std::endl
			<< "                                 0  The FASTA file contains amino acid codes. No translation is needed." << std::endl
			<< "                                 1  The DNA sequence is scanned left to right (forward direction)." << std::endl
			<< "                                    The amino acid code starts with the first DNA code." << std::cout
			<< "                                 2  The DNA sequence is scanned left to right (forward direction)." << std::cout
			<< "                                    The amino acid code starts with the second DNA code." << std::endl
			<< "                                 3  The DNA sequence is scanned left to right (forward direction)." << std::cout
			<< "                                    The amino acid code starts with the third DNA code." << std::cout
			<< "                                 4  The DNA sequence is scanned right to left (backward direction for the complementary strand)." << std::cout
			<< "                                    The amino acid code starts with the first DNA code." << std::cout
			<< "                                 5  The DNA sequence is scanned right to left (backward direction for the complementary strand)." << std::cout
			<< "                                    The amino acid code starts with the second DNA code." << std::cout
			<< "                                 6  The DNA sequence is scanned right to left (backward direction for the complementary strand)." << std::cout
			<< "                                    The amino acid code starts with the third DNA code." << std::cout
			<< "                                 7  Use each of the DNA translations of the codes 1, 2, 3." << std::cout
			<< "                                 8  Use each of the DNA translations of the codes 4, 5, 6." << std::cout
			<< "                                 9  Use each of the DNA translations of the codes 1, 2, 3, 4, 5, 6." << std::cout
			<< std::cout
			<< "  -max_num_int_cleav_sites       This value is the number of cleavage positions that may have been ignored by the enzyme." << std::endl
			<< "  -match_peak_count              The highest abundant experimental peaks are checked whether they are matched by the" << std::endl
			<< "                                 theoretical ones. match_peak_count is the number of the top abundant peaks to check." << std::endl
			<< "                                 A maximum of match_peak_allowed_error may lack this test." << std::endl
			<< "  -match_peak_allowed_error      see match_peak_count" << std::endl
			<< "  -show_fragment_ions            If set to 1 the fragment peaks of the top scored peptide are listed at the end of the output" << std::endl
			<< "  -use_phospho_fragmentation     ???" << std::endl
			<< "  -remove_precursor_peak         If set to 1 the peaks near (15 amu) the precursor are removed." << std::endl
			<< "  -mass_type_parent              A value of 1 selects monoisotopic masses, 0 selects average masses for calculating precursor peaks." << std::endl
			<< "  -mass_type_fragment            A value of 1 selects monoisotopic masses, 0 selects average masses for calculating fragment peaks." << std::endl
			<< "  -normalize_xcorr               Whether to use normalized xcorr values in the out files." << std::endl
			<< "  -residues_in_upper_case        Whether the residues in the FASTA database are in upper case." << std::endl
			<< "  -dyn_mods                      This value consists of semicolon-seperated pairs of variable modifications." << std::endl
			<< "                                 Each pair has two comma-seperated elements: A mass and a list of amino acids." << std::endl
			<< "                                 Sequest only applies the last modification character without warning." << std::endl
			<< "                                 Don't use \"44 S 80 ST\". It is interpreted as \"80 ST\"!." << std::endl
			<< "                                 Sequest won't apply any modification if the first two are null." << std::endl
			<< "                                 Always put valid modifications first. Don't use \"0 X 0 X 16 M\"" << std::endl
			<< "                                 Up to six modifications are allowed, if more are given, they are ignored." << std::endl
			<< "  -dyn_N_term_mod                This is the modification (mass that may be added to each N-terminus)" << std::endl
			<< "  -dyn_C_term_mod                This is the modification (mass that may be added to each C-terminus" << std::endl
			<< "  -stat_N_term_mod               This value is the mass that is added to each peptide N-terminus" << std::endl
			<< "  -stat_C_term_mod               This value is the mass that is added to each peptide C-terminus" << std::endl
			<< "  -stat_N_term_prot_mod          This value is the mass that is added to each protein N-terminus" << std::endl
			<< "  -stat_C_term_prot_mod          This value is the mass that is added to each protein C-terminus" << std::endl
			<< "  -stat_mods                     This value consists of a semicolon-seperated list of amino acids in one letter code" << std::endl
			<< "                                 and their corrpesponding mass: <AA_1>,<mass_1>;<AA_2>,<mass_2>;..." << std::endl
			<< "  -partial_sequence              A comma delimited list of amino acid sequences that must occur in the theoretical spectra." << std::endl
			<< "  -header_filter                 Several elements can be splitted by commas. Each element can be introduced" << std::endl
			<< "                                 by an exclamation mark (!) meaning that this element must not appear" << std::endl
			<< "                                 in the header of a protein or the protein will be skipped. This test is done first." << std::endl
			<< "                                 Next, all other elements are tested. The protein is processed" << std::endl
			<< "                                 if one filter string matches the header string." << std::endl
			<< "                                 A filter string may contain a tilde (~). This is replaced by a blank during comparison." << std::endl
			<< "  -keep_out_files                If set to 1, the Seuest .out-files are not removed (default for -Sequest_out)" << std::endl
			<< "  -keep_dta_files                If set to 1, the dta-files that were created from the mzXML-files are not removed" << std::endl
			<< "                                 (default for -Sequest_in)" << std::endl
			<< "  -contactName                   " << std::endl
			<< "  -contactInstitution            " << std::endl
			<< "  -contactInfo                   " << std::endl;
		}


		void setOptionsAndFlags_()
		{
			flags_["-show_enzyme_numbers"] = "show_enzyme_numbers";
			options_["-contactName"] = "contactName";
			options_["-contactInstitution"] = "contactInstitution";
			options_["-contactInfo"] = "contactInfo";
			options_["-log"] = "log";
			options_["-sequest_dir_win"] = "sequest_dir_win";
			options_["-user"] = "user";
			options_["-password"] = "password";
			options_["-temp_data_dir"] = "temp_data_dir";
			options_["-temp_data_dir_win"] = "temp_data_dir_win";
			options_["-temp_data_dir_network"] = "temp_data_dir_network";
			options_["-dta_dir"] = "dta_dir";
			options_["-dta_dir_win"] = "dta_dir_win";
			options_["-dta_dir_network"] = "dta_dir_network";
			options_["-out_dir"] = "out_dir";
			options_["-out_dir_win"] = "out_dir_win";
			options_["-out_dir_network"] = "out_dir_network";
			options_["-in"] = "in";
			options_["-in_win"] = "in_win";
			options_["-in_dir_network"] = "in_dir_network";
			options_["-db"] = "db";
			options_["-db_win"] = "db_win";
			options_["-db_dir_network"] = "db_dir_network";
			options_["-snd_db"] = "snd_db";
			options_["-snd_db_win"] = "snd_db_win";
			options_["-snd_db_dir_network"] = "snd_db_dir_network";
			options_["-sequest_computer"] = "sequest_computer";
			flags_["-Sequest_in"] = "Sequest_in";
			flags_["-Sequest_out"] = "Sequest_out";
			options_["-spectra"] = "spectra";
			options_["-pep_mass_tol"] = "pep_mass_tol";
			options_["-frag_ion_tol"] = "frag_ion_tol";
			options_["-match_peak_tol"] = "match_peak_tol";
			options_["-ion_cutoff"] = "ion_cutoff";
			options_["-pep_mass_unit"] = "pep_mass_unit";
			options_["-num_results"] = "num_results";
			options_["-enzyme_info"] = "enzyme_info";
			options_["-enzyme_number"] = "enzyme_number";
			options_["-min_prot_mass"] = "min_prot_mass";
			options_["-max_prot_mass"] = "max_prot_mass";
			options_["-max_num_dif_AA_per_mod"] = "max_num_dif_AA_per_mod";
			options_["-max_num_dif_mods_per_peptide"] = "max_num_dif_mods_per_peptide";
			options_["-nuc_reading_frame"] = "nuc_reading_frame";
			options_["-max_num_int_cleav_sites"] = "max_num_int_cleav_sites";
			options_["-match_peak_count"] = "match_peak_count";
			options_["-match_peak_allowed_error"] = "match_peak_allowed_error";
			flags_["-show_fragment_ions"] = "show_fragment_ions";
			flags_["-use_phospho_fragmentation"] = "use_phospho_fragmentation";
			flags_["-remove_precursor_peak"] = "remove_precursor_peak";
			options_["-mass_type_parent"] = "mass_type_parent";
			options_["-mass_type_fragment"] = "mass_type_fragment";
			flags_["-normalize_xcorr"] = "normalize_xcorr";
			flags_["-residues_in_upper_case"] = "residues_in_upper_case";
			options_["-neutral_loss_ABY"] = "neutral_loss_ABY";
			options_["-ion_series_weights"] = "ion_series_weights";
			options_["-dyn_mods"] = "dyn_mods";
			options_["-dyn_N_term_mod"] = "dyn_N_term_mod";
			options_["-dyn_C_term_mod"] = "dyn_C_term_mod";
			options_["-stat_N_term_mod"] = "stat_N_term_mod";
			options_["-stat_C_term_mod"] = "stat_C_term_mod";
			options_["-stat_N_term_prot_mod"] = "stat_N_term_prot_mod";
			options_["-stat_C_term_prot_mod"] = "stat_C_term_prot_mod";
			options_["-stat_mods"] = "stat_mods";
			options_["-partial_sequence"] = "partial_sequence";
			options_["-header_filter"] = "header_filter";
			options_["-out"] = "out";
			options_["-p_value"] = "p_value";
			options_["-prob_charge"] = "prob_charge";
			flags_["-keep_out_files"] = "keep_out_files";
			flags_["-keep_dta_files"] = "keep_dta_files";
		}

		inline void ensurePathChar(std::string& path, char path_char = '/')
		{
			if ( !path.empty() && (std::string("/\\").find(path[path.length()-1], 0) == std::string::npos) ) path.append(1, path_char);
		}

		bool isWinFormat(const std::string& name)
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
							if ( name.find(" ") == std::string::npos ) return true;
							else return false;
						}
						else return false;
					}
					else return true;
				}
			}
			return false;
		}

		void correctNetworkPath(String& network_path)
		{
			if ( network_path.hasSuffix("\\") ) network_path.erase(--network_path.end());
		}
		
		long fsize(const std::string& filename)
		{
			FILE* file = fopen(filename.c_str(), "r");
			long size = 0;
			if ( file != NULL )
			{
				fseek(file, 0, SEEK_END);
				size = ftell(file);
				fclose(file);
				return size;
			}
			return -1;
		}
		
		inline bool
		emptyFile(
			const std::string& filename)
		{
			return ( fsize(filename) == 0 );
		}
		
		void
		deleteTempFiles(
			const String& input_filename,
			const String& logfile)
		{
			if ( input_filename.hasSuffix("tmp.sequest.input") ) remove(input_filename.c_str());
			if ( logfile.hasSuffix("tmp.sequest.log") ) remove(logfile.c_str());
		}

		unsigned int
		MSExperiment2DTAs(
			MSExperiment<>& msexperiment,
			const String& common_name,
			unsigned int prob_charge,
			std::vector< std::pair< String, std::vector< double > > >& filenames_and_precursor_retention_times,
			bool make_dtas = true)
		throw (Exception::UnableToCreateFile)
		{
			DTAFile dtafile;
			String filename;
			unsigned int scan_number = 0;
			unsigned int msms_spectra = 0;
			std::vector< double > retention_times;
			
			for ( MSExperiment<>::Iterator spec_i = msexperiment.begin(); spec_i != msexperiment.end(); ++spec_i )
			{
				++scan_number;
				if ( spec_i ->getMSLevel() == 2 )
				{
					++msms_spectra;
					if ( spec_i->getPrecursorPeak().getCharge() )
					{
						if ( make_dtas )
						{
							filename = common_name + "." + String(scan_number) + "." + String(spec_i->getPrecursorPeak().getCharge()) + ".dta";
							dtafile.store(filename, *spec_i);
						}
						retention_times.push_back(spec_i->getRetentionTime());
					}
					else
					{
						for ( unsigned int i = 1; i <= prob_charge; ++i )
						{
							if ( make_dtas )
							{
								spec_i->getPrecursorPeak().setCharge(i);
								filename = common_name + "." + String(scan_number) + "." + String(spec_i->getPrecursorPeak().getCharge()) + ".dta";
								dtafile.store(filename, *spec_i);
							}
							retention_times.push_back(spec_i->getRetentionTime());
						}
						spec_i->getPrecursorPeak().setCharge(0);
					}
				}
			}
			
			filenames_and_precursor_retention_times.push_back(std::make_pair(filename, retention_times));
			
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
			String logfile, output_filename, input_filename, input_filename_win, input_file_dir_network, user, password, sequest_computer;
			ContactPerson contact_person;
			bool Sequest_in, Sequest_out, keep_out_files, keep_dta_files;
			String temp_data_dir, temp_data_dir_win, temp_data_dir_network, sequest_dir_win, database, database_win, database_dir_network, snd_database, snd_database_win, snd_database_dir_network, dta_dir, dta_dir_win, dta_dir_network, out_dir, out_dir_win, out_dir_network, batch_filename;
			
			std::vector< String > substrings, spectra;
			String string_buffer, string_buffer2;
			double double_buffer;
			int	int_buffer;
			
			// (1.1) variables for the analysis
			float p_value = 0.05;
			unsigned int prob_charge = 1;
			
			std::vector< std::pair< String, std::vector< double > > > filenames_and_precursor_retention_times;

			//-------------------------------------------------------------
			// (2) parsing and checking parameters
			//-------------------------------------------------------------

			// (2.0) variables for running the program
			logfile = getParamAsString_("log", "temp.sequest.log");
			
			if ( !getParamAsString_("show_enzyme_numbers").empty() )
			{
				writeLog_("Option show_enzyme_numbers chosen. Aborting.");
				std::cout << "Enzyme numers:" << std::endl << sequest_infile.getEnzymeInfo();
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}
			
			contact_person.setName(getParamAsString_("contactName", "unknown"));
			contact_person.setInstitution(getParamAsString_("contactInstitution", "unknown"));
			contact_person.setContactInfo(getParamAsString_("contactInfo"));

			Sequest_in = getParamAsBool_("Sequest_in");
			Sequest_out = getParamAsBool_("Sequest_out");

			// a 'normal' sequest run corresponds to both Sequest_in and Sequest_out set
			if ( !Sequest_in && !Sequest_out ) Sequest_in = Sequest_out = true;
			
			if ( Sequest_in && Sequest_out )
			{
				user = getParamAsString_("user");
				if ( user.empty() )
				{
					writeLog_("No user name for Sequest computer given. Aborting!");
					std::cout << "No user name for Sequest computer given. Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				
				password = getParamAsString_("password");
				if ( password.empty() )
				{
					writeLog_("No password for user name for Sequest computer given. Aborting!");
					std::cout << "No password for user name for Sequest computer given. Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
			
				sequest_dir_win = getParamAsString_("sequest_dir_win");
				ensurePathChar(sequest_dir_win, '\\');
				if ( !isWinFormat(sequest_dir_win) && !sequest_dir_win.empty() )
				{
					writeLog_("Windows path for the SEQUEST working directory has wrong format. Aborting!");
					std::cout << "Windows path for the SEQUEST working directory has wrong format. Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else if ( sequest_dir_win.empty() )
				{
					writeLog_("No windows path for the SEQUEST working directory given. Assuming PATH variable to be set accordingly!");
					sequest_dir_win = "sequest";
				}
				
				sequest_computer = getParamAsString_("sequest_computer");
				if ( sequest_computer.empty() )
				{
					writeLog_("No sequest computer name given. Aborting!");
					std::cout << "No sequest computer name given. Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
			}
			
			keep_out_files = getParamAsBool_("keep_out_files");
			if ( Sequest_out && !Sequest_in ) keep_out_files = true;

			keep_dta_files = getParamAsBool_("keep_dta_files");
			if ( Sequest_in && !Sequest_out ) keep_dta_files = true;
			
			temp_data_dir = getParamAsString_("temp_data_dir");
			ensurePathChar(temp_data_dir);
			temp_data_dir_win = getParamAsString_("temp_data_dir_win");
			ensurePathChar(temp_data_dir_win, '\\');
			temp_data_dir_network = getParamAsString_("temp_data_dir_network");
			correctNetworkPath(temp_data_dir_network);
			
			dta_dir = getParamAsString_("dta_dir");
			ensurePathChar(dta_dir);
			dta_dir_win = getParamAsString_("dta_dir_win");
			ensurePathChar(dta_dir_win, '\\');
			dta_dir_network = getParamAsString_("dta_dir_network");
			correctNetworkPath(dta_dir_network);
			
			out_dir = getParamAsString_("out_dir");
			ensurePathChar(out_dir);
			out_dir_win = getParamAsString_("out_dir_win");
			ensurePathChar(out_dir_win, '\\');
			out_dir_network = getParamAsString_("out_dir_network");
			correctNetworkPath(out_dir_network);

			bool in_uses_temp_data_dir = false;
			
			input_filename = getParamAsString_("in");
			if ( input_filename.empty() )
			{
				if ( !Sequest_out ) // if only Sequest_in is set, a name has to be given
				{
					writeLog_("No input file specified. Aborting!");
					std::cout << "No input file specified. Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else
				{
					temp_data_dir = true;
					input_filename = temp_data_dir + "temp.sequest.in";
					input_filename_win = temp_data_dir_win + "temp.sequest.in";
				}
			}
			else if ( Sequest_out )
			{
				input_filename_win = getParamAsString_("in_win");
				if ( !isWinFormat(input_filename_win) )
				{
					writeLog_("Windows path for input file has wrong format. Aborting!");
					std::cout << "Windows path for input file has wrong format. Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
					}
			}

			
			if ( temp_data_dir.empty() && !(out_dir.empty() || dta_dir.empty()) )
			{
				writeLog_("No directory for temporary files given. Aborting!");
				std::cout << "No directory for temporary files given. Aborting!" << std::endl;
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}
			else
			{
				if ( !temp_data_dir.empty() )
				{
					if ( Sequest_in && Sequest_out )
					{
						if ( (in_uses_temp_data_dir || out_dir.empty() || dta_dir.empty()) )
						{
							if ( !isWinFormat(temp_data_dir_win) )
							{
								writeLog_("Windows path for the directory for temporary files has wrong format. Aborting!");
								std::cout << "Windows path for the directory for temporary files has wrong format. Aborting!" << std::endl;
								printUsage_();
								return ILLEGAL_PARAMETERS;
							}
							if ( temp_data_dir_network.empty() )
							{
								writeLog_("Network path for the directory for temporary files is empty. Aborting!");
								std::cout << "Network path for the directory for temporary files is empty. Aborting!" << std::endl;
								printUsage_();
								return ILLEGAL_PARAMETERS;
							}
						}
					}
					
					if ( out_dir.empty() )
					{
						out_dir = temp_data_dir;
						out_dir_win = temp_data_dir_win;
						out_dir_network = temp_data_dir_network;
					}
					
					if ( dta_dir.empty() )
					{
						dta_dir = temp_data_dir;
						dta_dir_win = temp_data_dir_win;
						dta_dir_network = temp_data_dir_network;
					}
				}
				else if ( Sequest_in && Sequest_out )
				{
					writeLog_("No directory for temporary files given. Aborting!");
					std::cout << "No directory for temporary files given. Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
			}

			
			if ( logfile == out_dir + "sequest.log")
			{
				writeLog_("The logfile must not be named " + out_dir + "sequest.log. Aborting!");
				std::cout << "The logfile must not be named " + out_dir + "sequest.log. Aborting!" << std::endl;
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}

			batch_filename = getParamAsString_("batchfile");
			if ( batch_filename.empty() ) batch_filename = "sequest_run.bat";
			else if ( !batch_filename.hasSuffix(".bat") ) batch_filename.append(".bat");
			
			database = getParamAsString_("db");
			if ( database.empty() )
			{
				writeLog_("No database specified. Aborting!");
				std::cout << "No database specified. Aborting!" << std::endl;
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}

			snd_database = getParamAsString_("snd_db");
			
			string_buffer = getParamAsString_("spectra");
			if ( string_buffer.empty() )
			{
				writeLog_("No spectrum file specified. Aborting!");
				std::cout << "No spectrum file specified. Aborting!" << std::endl;
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}
			else
			{
				string_buffer.split(',', spectra);
				if ( spectra.empty() ) spectra.push_back(string_buffer);
			}
			
			if ( Sequest_in )
			{
				database_win = getParamAsString_("db_win");
				if ( !isWinFormat(database_win) )
				{
					writeLog_("Windows path for database has wrong format. Aborting!");
					std::cout << "Windows path for database has wrong format. Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				sequest_infile.setDatabase(database_win);

				if ( !snd_database.empty() )
				{
					snd_database_win = getParamAsString_("snd_db_win");
					if ( !isWinFormat(database_win) )
					{
						writeLog_("Windows path for database has wrong format. Aborting!");
						std::cout << "Windows path for database has wrong format. Aborting!" << std::endl;
						printUsage_();
						return ILLEGAL_PARAMETERS;
					}
					sequest_infile.setSndDatabase(snd_database_win);
				}
				
				double_buffer = getParamAsDouble_("pep_mass_tol", -1);
				if ( double_buffer == -1 )
				{
					writeLog_("No peptide mass tolerance specified. Aborting!");
					std::cout << "No peptide mass tolerance specified. Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else if ( double_buffer < 0 )
				{
					writeLog_("Peptide mass tolerance < 0. Aborting!");
					std::cout << "Peptide mass tolerance < 0. Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setPeptideMassTolerance(double_buffer);
				
				double_buffer = getParamAsDouble_("frag_ion_tol", -1);
				if ( double_buffer == -1 )
				{
					writeLog_("No fragment ion tolerance specified. Aborting!");
					std::cout << "No fragment ion tolerance specified. Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else if ( double_buffer < 0 )
				{
					writeLog_("Fragment ion tolerance < 0. Aborting!");
					std::cout << "Fragment ion tolerance < 0. Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setFragmentIonTolerance(double_buffer);
				
				double_buffer = getParamAsDouble_("match_peak_tol", -1);
				if ( double_buffer == -1 )
				{
					writeLog_("No match peak tolerance specified. Aborting!");
					std::cout << "No match peak tolerance specified. Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else if ( double_buffer < 0 )
				{
					writeLog_("Match peak tolerance < 0. Aborting!");
					std::cout << "Match peak tolerance < 0. Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setMatchPeakTolerance(double_buffer);
				
				double_buffer = getParamAsDouble_("ion_cutoff", 0);
				if ( double_buffer < 0 )
				{
					writeLog_("Ion cutoff < 0. Aborting!");
					std::cout << "Ion cutoff < 0. Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setIonCutoffPercentage(double_buffer);
				
				int_buffer = getParamAsInt_("pep_mass_unit", 0);
				if ( (int_buffer < 0) || (int_buffer > max_peptide_mass_units) )
				{
					writeLog_("Illegal peptide mass unit (not in [0,2]). Aborting!");
					std::cout << "Illegal peptide mass unit (not in [0,2]). Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setPeptideMassUnits(int_buffer);
				
				int_buffer = getParamAsInt_("num_results", 0);
				if ( (int_buffer < 1) )
				{
					writeLog_("Illegal number of results (< 1). Aborting!");
					std::cout << "Illegal number of results (< 1). Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setNumOutputLines(int_buffer);
				
				string_buffer = getParamAsString_("enzyme_info");
				if ( !string_buffer.empty() )
				{
					string_buffer.split(';', substrings);
					if ( substrings.empty() ) substrings.push_back(string_buffer);
					
					std::vector< String > enzyme_info;
					for ( std::vector< String >::iterator einfo_i = substrings.begin(); einfo_i != substrings.end(); ++ einfo_i )
					{
						einfo_i->split(',', enzyme_info);
						if ( (enzyme_info.size() < 3) || (enzyme_info.size() > 4) )
						{
							writeLog_("Illegal number of informations for enzyme (not in [3,4]). Aborting!");
							std::cout << "Illegal number of informations for enzyme (not in [3,4]). Aborting!" << std::endl;
							printUsage_();
							return ILLEGAL_PARAMETERS;
						}
						if ( !((enzyme_info[1] == "0") || (enzyme_info[1] == "1"))  )
						{
							writeLog_("Cut direction for enzyme not specified correctly (has to be 1 (N to C)) or 0 (C to N))). Aborting!");
							std::cout << "Cut direction for enzyme not specified correctly (has to be 1 (N to C) or 0 (C to N)). Aborting!" << std::endl;
							printUsage_();
							return ILLEGAL_PARAMETERS;
						}
						if ( enzyme_info.size() == 3 ) enzyme_info.push_back("-");
						sequest_infile.addEnzymeInfo(enzyme_info);
					}
				}
				else
				{
					substrings.clear();
					sequest_infile.setEnzymeNumber(getParamAsInt_("enzyme_number"));
				}
				
				double_buffer = getParamAsDouble_("min_prot_mass");
				if ( double_buffer < 0 )
				{
					writeLog_("Illegal minimum protein mass (< 0). Aborting!");
					std::cout << "Illegal minimum protein mass (< 0). Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setMinimumProteinMass(double_buffer);
				
				double_buffer = getParamAsDouble_("max_prot_mass");
				if ( double_buffer < 0 )
				{
					writeLog_("Illegal maximum protein mass (< 0). Aborting!");
					std::cout << "Illegal maximum protein mass (< 0). Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setMaximumProteinMass(double_buffer);
				
				int_buffer = getParamAsInt_("max_num_dif_AA_per_mod", -1);
				if ( int_buffer < 0 )
				{
					writeLog_("No maximum number of modified amino acids per different modification. Aborting!");
					std::cout << "No maximum number of modified amino acids per different modification. Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setMaxNumDifAAPerMod(int_buffer);
				
				int_buffer = getParamAsInt_("max_num_dif_mods_per_peptide", -1);
				if ( int_buffer < 0 )
				{
					writeLog_("No maximum number of differential modifications per peptide. Aborting!");
					std::cout << "No maximum number of differential modifications per peptide. Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setMaxNumModsPerPeptide(int_buffer);
				
				int_buffer = getParamAsInt_("nuc_reading_frame", 0);
				if ( (int_buffer < 0) || (int_buffer > 9) )
				{
					writeLog_("Illegal number for nucleotide reading frame. Aborting!");
					std::cout << "Illegal number for nucleotide reading frame. Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setNucleotideReadingFrame(int_buffer);
				
				int_buffer = getParamAsInt_("max_num_int_cleav_sites", 0);
				if ( int_buffer < 0 )
				{
					writeLog_("Illegal number of maximum internal cleavage sites. Aborting!");
					std::cout << "Illegal number of maximum internal cleavage sites. Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setMaxNumInternalCleavageSites(int_buffer);
				
				int_buffer = getParamAsInt_("match_peak_count", 0);
				if ( (int_buffer < 0) && (int_buffer > 5) )
				{
					writeLog_("Illegal number of auto-detected peaks to try matching. Aborting!");
					std::cout << "Illegal number of auto-detected peaks to try matching. Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setMatchPeakCount(int_buffer);
				
				int_buffer = getParamAsInt_("match_peak_allowed_error", 0);
				if ( int_buffer < 0 )
				{
					writeLog_("Illegal number of allowed errors in matching auto-detected peaks. Aborting!");
					std::cout << "Illegal number of allowed errors in matching auto-detected peaks. Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setMatchPeakAllowedError(int_buffer);
				
				sequest_infile.setShowFragmentIons(getParamAsBool_("show_fragment_ions"));
				sequest_infile.setUsePhosphoFragmentation(getParamAsBool_("use_phospho_fragmentation"));
				sequest_infile.setRemovePrecursorPeak(getParamAsBool_("remove_precursor_peak"));
				sequest_infile.setMassTypeParent(getParamAsBool_("mass_type_parent"));
				sequest_infile.setMassTypeFragment(getParamAsBool_("mass_type_fragment"));
				sequest_infile.setNormalizeXcorr(getParamAsBool_("normalize_xcorr"));
				sequest_infile.setResiduesInUpperCase(getParamAsBool_("residues_in_upper_case", true));
				
				string_buffer = getParamAsString_("neutral_loss_ABY");
				string_buffer2 = "01";
				if ( (string_buffer.size() != 3) || (string_buffer2.find(string_buffer[0], 0) == std::string::npos) || (string_buffer2.find(string_buffer[1], 0) == std::string::npos) || (string_buffer2.find(string_buffer[2], 0) == std::string::npos) )
				{
					writeLog_("Neutral losses for ABY-ions not given (or illegal values given). Aborting!");
					std::cout << "Neutral losses for ABY-ions not given (or illegal values given). Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setNeutralLossesForIons(string_buffer);
				
				string_buffer = getParamAsString_("ion_series_weights");
				string_buffer.split(',', substrings);
				if ( substrings.size() != 9 )
				{
					writeLog_("Weights for ion series not given (or illegal values given). Aborting!");
					std::cout << "Weights for ion series not given (or illegal values given). Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else
				{
					for ( std::vector< String >::iterator s_i = substrings.begin(); s_i != substrings.end(); ++s_i )
					{
						// the values are expected to be float, otherwise they will be seen as 0!
						double_buffer = atof(s_i->c_str());
						if ( (double_buffer < 0) || (double_buffer > 1) )
						{
							writeLog_("Illegal weights for ion series given. Aborting!");
							std::cout << "Illegal weights for ion series given. Aborting!" << std::endl;
							printUsage_();
							return ILLEGAL_PARAMETERS;
						}
						(*s_i) = String(double_buffer);
					}
					string_buffer.implode(substrings.begin(), substrings.end(), " ");
					sequest_infile.setIonSeriesWeights(string_buffer);
				}
				
				string_buffer = getParamAsString_("dyn_mods");
				if ( !string_buffer.empty() )
				{
					string_buffer.split(';', substrings);
					if ( substrings.empty() ) substrings.push_back(string_buffer);
					float f_buffer;
					char c_buffer[41]; c_buffer[40] = 0;
					for ( std::vector< String >::iterator s_i = substrings.begin(); s_i != substrings.end(); ++s_i )
					{
						if ( sscanf(s_i->c_str(), "%f,%40s", &f_buffer, c_buffer) != 2 )
						{
							writeLog_("Illegal number of parameters for dynamic modification given. Aborting!");
							std::cout << "Illegal number of parameters for dynamic modification given. Aborting!" << std::endl;
							printUsage_();
							return ILLEGAL_PARAMETERS;
						}
						(*s_i) = String(f_buffer)+" "+c_buffer;
					}
					string_buffer.implode(substrings.begin(), substrings.end(), " ");
					sequest_infile.setDynMods(string_buffer);
				}
				
				sequest_infile.setDynNTermMod(getParamAsDouble_("dyn_N_term_mod"));
				sequest_infile.setDynCTermMod(getParamAsDouble_("dyn_C_term_mod"));
				
				sequest_infile.setStatNTermMod(getParamAsDouble_("stat_N_term_mod"));
				sequest_infile.setStatCTermMod(getParamAsDouble_("stat_C_term_mod"));
				
				sequest_infile.setStatNTermProtMod(getParamAsDouble_("stat_N_term_prot_mod"));
				sequest_infile.setStatCTermProtMod(getParamAsDouble_("stat_C_term_prot_mod"));
				
				string_buffer = getParamAsString_("stat_mods");
				if ( !string_buffer.empty() )
				{
					string_buffer.split(';', substrings);
					if ( substrings.empty() ) substrings.push_back(string_buffer);
					
					for ( std::vector< String >::iterator s_i = substrings.begin(); s_i != substrings.end(); ++s_i )
					{
						if ( (*s_i)[1] != ',' )
						{
							writeLog_("Unexpected format for static modification found. Aborting!");
							std::cout << "Unexpected format for static modification found. Aborting!" << std::endl;
							printUsage_();
							return ILLEGAL_PARAMETERS;
						}
						sequest_infile.setStatMod((*s_i)[0], atof(s_i->substr(2).c_str()));
					}
				}
				
				string_buffer = getParamAsString_("partial_sequence");
				string_buffer.split(',', substrings);
				string_buffer.implode(substrings.begin(), substrings.end(), " ");
				sequest_infile.setPartialSequence(string_buffer);
				
				string_buffer = getParamAsString_("header_filter");
				string_buffer.split(',', substrings);
				string_buffer.implode(substrings.begin(), substrings.end(), " ");
				sequest_infile.setSequenceHeaderFilter(string_buffer);
			}
			
			if ( Sequest_out )
			{
				int_buffer = getParamAsInt_("prob_charge");
				if ( int_buffer < 0 )
				{
					writeLog_("Maximal charge to test is less than zero. Aborting!");
					std::cout << "Maximal charge to test is less than zero. Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else if ( int_buffer ) prob_charge = int_buffer;
				else prob_charge = 1;
				
				string_buffer = getParamAsString_("out");
				if ( string_buffer.empty() )
				{
					writeLog_("No output file specified. Aborting!");
					std::cout << "No output file specified. Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else output_filename = string_buffer;
			
				double_buffer = getParamAsDouble_("p_value", -1.0);
				if ( double_buffer != -1.0 ) p_value = double_buffer;
				if ( (p_value <= 0) || (p_value > 1) )
				{
					writeLog_("P-value not in (0,1]. Aborting!");
					std::cout << "P-value not in (0,1]. Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
			}

			std::vector< String > drive_letters;
			std::vector< String > network_paths;
			if ( Sequest_in && Sequest_out )
			{
				if ( !isWinFormat(out_dir_win) )
				{
					writeLog_("Windows path for the directory for .out files has wrong format. Aborting!");
					std::cout << "Windows path for the directory for .out files has wrong format. Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				
				if ( !isWinFormat(dta_dir_win) )
				{
					writeLog_("Windows path for the directory for .dta files has wrong format. Aborting!");
					std::cout << "Windows path for the directory for .dta files has wrong format. Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}

				database_dir_network = getParamAsString_("db_dir_network");
				correctNetworkPath(database_dir_network);
				snd_database_dir_network = getParamAsString_("snd_db_dir_network");
				correctNetworkPath(snd_database_dir_network);
				input_file_dir_network = getParamAsString_("in_dir_network");
				correctNetworkPath(input_file_dir_network);

				// make a list of the directories that have to be mounted
				std::vector< String > drive_letters_all;
				std::vector< String > network_paths_all;

				drive_letters_all.push_back(out_dir_win.substr(0,2));
				network_paths_all.push_back(out_dir_network);
				drive_letters_all.push_back(dta_dir_win.substr(0,2));
				network_paths_all.push_back(dta_dir_network);
				drive_letters_all.push_back(database_win.substr(0,2));
				network_paths_all.push_back(database_dir_network);
				if ( !snd_database.empty() )
				{
					drive_letters_all.push_back(snd_database_win.substr(0,2));
					network_paths_all.push_back(snd_database_dir_network);
				}
				drive_letters_all.push_back(input_filename_win.substr(0,2));
				network_paths_all.push_back(input_file_dir_network);

				// go through the lists and search for any drive letter that has a corresponding network path
				std::vector< String >::const_iterator network_i = network_paths_all.begin();
				for ( std::vector< String >::const_iterator drive_i = drive_letters_all.begin(); drive_i != drive_letters_all.end(); ++drive_i, ++network_i )
				{
				//std::cout << *drive_i << "\t" << *network_i << std::endl;
					if ( network_i->empty() )
					{
						std::vector< String >::const_iterator network_i2 = network_paths_all.begin();
						for ( std::vector< String >::const_iterator drive_i2 = drive_letters_all.begin(); drive_i2 != drive_letters_all.end(); ++drive_i2, ++network_i2 )
						{
							if ( (*drive_i == *drive_i2) && !network_i2->empty() ) break;
						}
						if ( network_i2 == network_paths_all.end() )
						{
						//std::cout << *drive_i << "\t" << *network_i << std::endl;
							writeLog_("No network path for windows directory " + *drive_i +" given. Aborting!");
							std::cout << "No network path for windows directory " + *drive_i +" given. Aborting!" << std::endl;
							printUsage_();
							return ILLEGAL_PARAMETERS;
						}
						//std::cout << "emtpy" << std::endl;
					}
					else
					{
						drive_letters.push_back(*drive_i);
						network_paths.push_back(*network_i);
						//std::cout << "\t\t" <<*drive_i << "\t" << *network_i << std::endl;
					}
				}
				// now check whether there are drive letters with more than one network paths
				for ( std::vector< String >::iterator drive_i = drive_letters.begin(); drive_i != drive_letters.end() - 1; ++drive_i )
				{
/*std::cout << *drive_i << "\t" << *(drive_i+1) <<  std::endl;
					for ( std::vector< String >::iterator drive_i2 = drive_i+1; drive_i2 != drive_letters.end(); ++drive_i2 )
					{
					std::cout << *drive_i2 << std::endl;
						if ( *drive_i == *drive_i2 )
					}*/
						if ( std::find(drive_i+1, drive_letters.end(), *drive_i) != drive_letters.end() )
						{
							writeLog_("More than one network path for windows directory " + *drive_i +" given. Aborting!");
							std::cout << "More than one network path for windows directory " + *drive_i +" given. Aborting!" << std::endl;
							printUsage_();
							return ILLEGAL_PARAMETERS;
						}
					//}
					//std::cout << "1099" << std::endl;
				}
			}
			
			//-------------------------------------------------------------
			// (3) running program according to parameters
			//-------------------------------------------------------------

			// (3.1) checking accessability of files
			QFileInfo file_info;
			QFile file;
			
			if ( Sequest_in )
			{
				file.setName(input_filename);
				file.open(IO_WriteOnly);
				if ( !file.isWritable() )
				{
					throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, input_filename);
				}

				file.setName(out_dir + batch_filename);
				file.open(IO_WriteOnly);
				if ( !file.isWritable() )
				{
					throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, out_dir + batch_filename);
				}
			}
			
			// (3.1.2) output file
			if ( Sequest_out )
			{
				file.setName(output_filename);
				file.open(IO_WriteOnly);
				if ( !file.isWritable() )
				{
					throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, output_filename);
				}
					
				file_info.setFile(database);
				// first database
				if ( !file_info.exists() )
				{
					throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, database);
				}
				if ( !file_info.isReadable() )
				{
					throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, database);
				}
				if ( emptyFile(database) )
				{
					throw Exception::FileEmpty(__FILE__, __LINE__, __PRETTY_FUNCTION__, database);
				}
				
				// second database
				if ( !snd_database.empty() )
				{
					file_info.setFile(snd_database);
					if ( !file_info.exists() )
					{
						throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, snd_database);
					}
					if ( !file_info.isReadable() )
					{
						throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, snd_database);
					}
					if ( emptyFile(snd_database) )
					{
						throw Exception::FileEmpty(__FILE__, __LINE__, __PRETTY_FUNCTION__, snd_database);
					}
				}
			}
			
			// (3.2.1) creating the input file
			if ( Sequest_in )
			{
				sequest_infile.store(input_filename);
			}

			bool make_dtas = ( Sequest_out && !Sequest_in ) ? false : true; // if only Sequest_out is set, just get the retention times
			// creating the dta files
			MSExperiment<> msexperiment;
			unsigned int msms_spectra_in_file;
			unsigned int msms_spectra_altogether = 0;
			if ( make_dtas ) std::cout << "creating dta files" << std::endl;
			for ( std::vector< String >::const_iterator spec_i = spectra.begin(); spec_i != spectra.end(); ++spec_i )
			{
				file_info.setFile(*spec_i);
				String common_name = dta_dir + std::string(file_info.fileName().ascii());
				
				try
				{
					MzXMLFile().load(*spec_i, msexperiment);
				}
				catch ( Exception::ParseError pe )
				{
					writeLog_("Error loading mzXML file. Aborting!");
					std::cout << "Error loading mzXML file. Aborting!" << std::endl;
					printUsage_();
					return PARSE_ERROR;
				}
				
				//return UNKNOWN_ERROR;
				msms_spectra_in_file = MSExperiment2DTAs(msexperiment, common_name, prob_charge, filenames_and_precursor_retention_times, make_dtas);
				writeLog_(String(msms_spectra_in_file) + " MS/MS spectra in file " + file_info.fileName().ascii());

				msms_spectra_altogether += msms_spectra_in_file;
			}

			if ( !msms_spectra_altogether )
			{
				writeLog_("No MS/MS spectra found in any of the mzXML files. Aborting!");
				std::cout << "No MS/MS spectra found in any of the mzXML files. Aborting!" << std::endl;
				printUsage_();
				return UNKNOWN_ERROR;
			}

			// (3.2.3) running the program
			if ( Sequest_in && Sequest_out )
			{
				// creating a batch file for windows (command doesn't accept commands that are longer than 256 chars)
				std::ofstream batchfile(String(out_dir + batch_filename).c_str());
				if ( !batchfile )
				{
					throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, out_dir + batch_filename);
				}

				// get the drive and network path of out_dir
				std::vector< String >::iterator drive_i = std::find(drive_letters.begin(), drive_letters.end(), out_dir_win.substr(0,2));
				std::vector< String >::iterator network_i = network_paths.begin() + (drive_i - drive_letters.begin());
				
				String call = "rdesktop -u ";
				call.append(user);
				if ( !password.empty() ) call.append(" -p \"" + password + "\"");
				call.append(" -s cmd\\ /C\\ \"");
				call.append(" net use " + *drive_i + " \\\\" + *network_i + " && ");

				drive_letters.erase(drive_i);
				network_paths.erase(network_i);
				for ( drive_i = drive_letters.begin(), network_i = network_paths.begin(); drive_i != drive_letters.end(); ++drive_i, ++network_i )
				{
					batchfile << String(" net use " + *drive_i + " " + *network_i + " &&");
				}
				batchfile << String(" cd " + out_dir_win + " && " + out_dir_win.substr(0,2));
				batchfile << String(" && " + sequest_dir_win + "sequest.exe -P" + input_filename_win + " " + dta_dir_win + "*.dta");
				batchfile << String(" && " + sequest_dir_win.substr(0,2) + " &&");
				for ( std::vector< String >::const_iterator drive_i = drive_letters.begin(); drive_i != drive_letters.end(); ++drive_i, ++network_i )
				{
					batchfile << String(" net use /delete " + *drive_i + " &&");
				}
				batchfile << String(" net use /delete " + out_dir_win.substr(0,2));
				batchfile << " && logoff";
				batchfile.close();
				batchfile.clear();
				
				call.append(out_dir_win + batch_filename + "\" " + sequest_computer);
				std::cout << call << std::endl;
				int status = system(call.c_str());
				
				if ( status != 0 )
				{
					std::cout << "Sequest problem. Aborting! (Details can be seen " 
					<< " in the logfile: \"" << logfile << "\")" << std::endl;
					writeLog_("Sequest problem. Aborting!");
					deleteTempFiles(input_filename, logfile);
					return EXTERNAL_PROGRAM_ERROR;
				}

				std::ifstream sequest_log(std::string(out_dir + "sequest.log").c_str()); // write sequest log to logfile
				if ( !sequest_log )
				{
					std::cout << "No Sequest log found!" << std::endl;
					writeLog_("No Sequest log found!");
				}
				else
				{
					sequest_log.seekg (0, std::ios::end);
					std::streampos length = sequest_log.tellg();
					sequest_log.seekg (0, std::ios::beg);
					char * buffer = new char[length];
					sequest_log.read (buffer, length);
					sequest_log.close();
					sequest_log.clear();
					writeLog_(buffer);
					delete(buffer);
					remove(std::string(out_dir + "sequest.log").c_str());
				}

				if ( !keep_dta_files ) // remove all dtas
				{
					QDir qdir(dta_dir, "*.dta", QDir::Name, QDir::Files);
					QStringList qlist = qdir.entryList();
					QFile qfile;
					
					for ( QStringList::const_iterator i = qlist.constBegin(); i != qlist.constEnd(); ++i )
					{
						if ( !qfile.remove(QString(dta_dir.c_str() + *i)) )
						{
							std::cout << std::string((*i).ascii()) << "could not be removed!" << std::endl;
							writeLog_(std::string((*i).ascii()) + "could not be removed!");
						}
					}
				}
			}
			
			if ( Sequest_out )
			{
				if ( !keep_dta_files ) // remove all dtas
				{
					QDir qdir(dta_dir, "*.dta", QDir::Name, QDir::Files);
					QStringList qlist = qdir.entryList();
					QFile qfile;
					
					for ( QStringList::const_iterator i = qlist.constBegin(); i != qlist.constEnd(); ++i )
					{
						if ( !qfile.remove(QString(dta_dir.c_str() + *i)) )
						{
							std::cout << std::string((*i).ascii()) << "could not be removed!" << std::endl;
							writeLog_(std::string((*i).ascii()) + "could not be removed!");
						}
					}
				}
				
				AnalysisXMLFile analysisXML_file;
				
				SequestOutfile sequest_outfile;
				std::vector< Identification > identifications;
				ProteinIdentification protein_identification;
				std::vector< float > precursor_retention_times, precursor_mz_values;
				
				QDir qdir(dta_dir, "*.out", QDir::Name, QDir::Files);
				QStringList qlist = qdir.entryList();

				if ( qlist.isEmpty() )
				{
					writeLog_("No .out identified. Aborting!");
					std::cout << "No .out identified. Aborting!" << std::endl;
					qlist.clear();
					deleteTempFiles(input_filename, logfile);
					return UNKNOWN_ERROR;
				}

				for ( QStringList::const_iterator i = qlist.constBegin(); i != qlist.constEnd(); ++i )
				{
					sequest_outfile.load(out_dir + std::string((*i).ascii()), identifications, protein_identification, precursor_retention_times, precursor_mz_values,  p_value, database, snd_database);
				}

				std::sort(filenames_and_precursor_retention_times.begin(), filenames_and_precursor_retention_times.end(), SortRetentionTimes()); // sort the retention times, so they have the same order like the corresponding .out files

				// save the retention times in precursor_retention_times
				for( std::vector< std::pair< String, std::vector< double > > >::iterator pairs_i = filenames_and_precursor_retention_times.begin(); pairs_i != filenames_and_precursor_retention_times.end(); ++pairs_i )
				{
					precursor_retention_times.insert(precursor_retention_times.end(), pairs_i->second.begin(), pairs_i->second.end());
					pairs_i->second.clear();
				}
				
				std::vector< ProteinIdentification > pis;
				pis.push_back(protein_identification);

				analysisXML_file.store(output_filename, pis, identifications, precursor_retention_times, precursor_mz_values, contact_person);
				
				// remove all outs
				if ( !keep_out_files )
				{
					qdir.setPath(out_dir);
					qlist = qdir.entryList("*.out", QDir::Files, QDir::Name);
					QFile qfile;
					
					for ( QStringList::const_iterator i = qlist.constBegin(); i != qlist.constEnd(); ++i )
					{
						if ( !qfile.remove(QString(out_dir.c_str() + *i)) )
						{
							std::cout << std::string((*i).ascii()) << "could not be removed!" << std::endl;
							writeLog_(std::string((*i).ascii()) + "could not be removed!");
						}
					}
				}
				
				qlist.clear();
			}
			
			// (3.3) deleting all temporary files
			deleteTempFiles(input_filename, logfile);
			
			return EXECUTION_OK;
		}
};

//@endcond



int main( int argc, char ** argv )
{
	TOPPSequestAdapter tool;

	return tool.main(argc,argv);
}
