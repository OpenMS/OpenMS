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
// $Maintainer: Martin Langwisch $
// --------------------------------------------------------------------------


#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Date.h>
#include <OpenMS/METADATA/ContactPerson.h>
#include <OpenMS/FORMAT/SequestInfile.h>
#include <OpenMS/FORMAT/SequestOutfile.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/FORMAT/AnalysisXMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <stdlib.h>
#include <vector>
#include <algorithm>

#include <qdir.h>
#include <qstringlist.h>

using namespace OpenMS;
using namespace std;

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
	
	@todo look for possible crash codes of sequest and catching them (Martin)
	
	@ingroup TOPP
*/

// We do not want this class to show up in the docu -> cond
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
		static const unsigned int max_dtas_per_run = 1000;
		unsigned long long int dtas;
		
		struct SortRetentionTimes:
				public binary_function<const pair< String, vector< double > >, const pair< String, vector< double > >, bool>
		{
			bool operator()(const pair< String, vector< double > >& x, const pair< String, vector< double > >& y)
			{
				return ( x.first < y.first );
			}
		};
	
		void printToolUsage_() const
		{
			cerr	<< endl
						<< getToolName() << " -- annotates MS/MS spectra using Sequest" << endl
						<< "Version: " << VersionInfo::getVersion() << endl
						<< endl
						<< "Usage:" << endl
						<< " " << getToolName() << " [options]" << endl
						<< endl
						<< "Options are:" << endl
						<< "  [the _win parameters correspond to path of the linux directory when mounted under windows. NO SPACE ALLOWED!;" << endl
						<< "   the _network paramters correspond to the path network path of the directory when mounting the under windows;" << endl
						<< "   as Sequest runs on windows, all the files and directories used have to be mounted on the corresponding computer!" << endl
						<< "   rdesktop is used to connect to this computer" << endl << endl
						<< "   Sequest writes a file named 'sequest.log' into the out_dir, so you must not name the log file alike]" << endl
						<< "  -Sequest_in           if this flag is set the SequestAdapter will create a sequest input file" << endl
						<< "                        and create dta files from the given mzXML files" << endl
						<< "  -Sequest_out          if this flag is set the SequestAdapter will read in Sequest out files and write an analysis XML file" << endl
						<< "  -in                   the comma-seperated names of the mzXML files" << endl
						<< "  -out                  the name of the analysis XML file" << endl
						<< "  -show_enzyme_numbers  show a list with enzymes and corresponding numbers to choose from" << endl
						<< "  -create_dtas          creates dta files from the mzXML files" << endl
						<< "  -sequest_dir_win      the windows path where sequest.exe is located" << endl
						<< "  -sequest_computer     the name of the computer in the network that hosts Sequest" << endl
						<< "  -user                 user name for the sequest computer (has to have access to network!)" << endl
						<< "  -password             password for this user (if not given, you have to enter it at promt)" << endl
						<< "  -p_value              annotations with inferior p-value are ignored (default is 0.05)" << endl
						<< "  -num_results          the maximal number of results (peptides) to show (per scan/dta)" << endl
						<< "  -max_num_dif_AA_per_mod  limits the maximum total number of each single variable modification in one peptide" << endl
						<< "  -max_num_dif_mods_per_peptide  limits the maximum total number of each single variable modification in one peptide" << endl
						<< "  -prob_charge          the number of charge states that are used if it is unknown for a scan" << endl
						<< "  -pep_mass_tol         tolerance for a peptide ion" << endl
						<< "  -frag_ion_tol         tolerance for a fragment ion" << endl
						<< "  -match_peak_tol       the minimal space between two peaks" << endl
						<< "  -enzyme_info          <name>,<cut direction: N to C?>,<cuts after>,<doesn't cut before>;" << endl
						<< "                        cuts after, doesn't cut before: amino acids in 1-letter code or '-' for unspecific cleavage" << endl
						<< "  -enzyme_number        a number from the list (show_enzyme_numbers); if enzyme_info is used, this value is set accordingly" << endl
						<< "  -neutral_loss_ABY     ABY: 0 or 1 whether neutral losses of the series should be honored, eg: 011" << endl
						<< "  -ion_series_weights   abcdvwxyz: [0.0, 1.0] factor for the series, eg: 0,0.5,0,0,0,0,0,1.0,0" << endl
						<< "  -sequest_in           the name of the sequest input file" << endl
						<< "  -db                   the name of the database file" << endl
						<< "  -snd_db               the name of the second database file" << endl
						<< "  -temp_data_dir        the directory to store temporary data" << endl
						<< "  -temp_data_dir_win    " << endl
						<< endl
						<< "  For each windows drive, one corresponding network drive has to be given, so maybe you don't need to set all the parameters below" << endl
						<< "  -temp_data_dir_network" << endl
						<< "  -db_dir_network" << endl
						<< "  -snd_db_dir_network" << endl
						<< "  -sequest_in_dir_network" << endl;
		}


		void printToolHelpOpt_() const
		{
			cerr	<< endl
			<< "  -ion_cutoff                    This value selects a cut-off below which a matching peptide is rejected." << endl
			<< "                                 The value compared with this value is the ratio" << endl
			<< "                                 (# matching theoretical fragment peaks) / (# total theoretical fragment peaks)" << endl
			<< "                                 which means that the user can select a minimum coverage of matching peaks." << endl
			<< "  -pep_mass_unit                 peptide mass unit: 0=amu (atomic mass unit), 1=mmu (millimass unit), 2=ppm (parts per million)" << endl
			<< "  -min_prot_mass                 minimal protein mass" << endl
			<< "  -max_prot_mass                 maximal protein mass" << endl
			<< "  -nuc_reading_frame             Format of the FASTA database:" << endl
			<< "                                 0  The FASTA file contains amino acid codes. No translation is needed." << endl
			<< "                                 1  The DNA sequence is scanned left to right (forward direction)." << endl
			<< "                                    The amino acid code starts with the first DNA code." << endl
			<< "                                 2  The DNA sequence is scanned left to right (forward direction)." << endl
			<< "                                    The amino acid code starts with the second DNA code." << endl
			<< "                                 3  The DNA sequence is scanned left to right (forward direction)." << endl
			<< "                                    The amino acid code starts with the third DNA code." << endl
			<< "                                 4  The DNA sequence is scanned right to left (backward direction for the complementary strand)." << endl
			<< "                                    The amino acid code starts with the first DNA code." << endl
			<< "                                 5  The DNA sequence is scanned right to left (backward direction for the complementary strand)." << endl
			<< "                                    The amino acid code starts with the second DNA code." << endl
			<< "                                 6  The DNA sequence is scanned right to left (backward direction for the complementary strand)." << endl
			<< "                                    The amino acid code starts with the third DNA code." << endl
			<< "                                 7  Use each of the DNA translations of the codes 1, 2, 3." << endl
			<< "                                 8  Use each of the DNA translations of the codes 4, 5, 6." << endl
			<< "                                 9  Use each of the DNA translations of the codes 1, 2, 3, 4, 5, 6." << endl
			<< endl
			<< "  -max_num_int_cleav_sites       This value is the number of cleavage positions that may have been ignored by the enzyme." << endl
			<< "  -match_peak_count              The highest abundant experimental peaks are checked whether they are matched by the" << endl
			<< "                                 theoretical ones. match_peak_count is the number of the top abundant peaks to check." << endl
			<< "                                 A maximum of match_peak_allowed_error may lack this test." << endl
			<< "  -match_peak_allowed_error      see match_peak_count" << endl
			<< "  -show_fragment_ions            If set to 1 the fragment peaks of the top scored peptide are listed at the end of the output" << endl
			<< "  -use_phospho_fragmentation     ???" << endl
			<< "  -remove_precursor_peak         If set to 1 the peaks near (15 amu) the precursor are removed." << endl
			<< "  -mass_type_parent              A value of 1 selects monoisotopic masses, 0 selects average masses for calculating precursor peaks." << endl
			<< "  -mass_type_fragment            A value of 1 selects monoisotopic masses, 0 selects average masses for calculating fragment peaks." << endl
			<< "  -normalize_xcorr               Whether to use normalized xcorr values in the out files." << endl
			<< "  -residues_in_upper_case        Whether the residues in the FASTA database are in upper case." << endl
			<< "  -dyn_mods                      This value consists of semicolon-seperated pairs of variable modifications." << endl
			<< "                                 Each pair has two comma-seperated elements: A mass and a list of amino acids." << endl
			<< "                                 Sequest only applies the last modification character without warning." << endl
			<< "                                 Don't use \"44 S 80 ST\". It is interpreted as \"80 ST\"!." << endl
			<< "                                 Sequest won't apply any modification if the first two are null." << endl
			<< "                                 Always put valid modifications first. Don't use \"0 X 0 X 16 M\"" << endl
			<< "                                 Up to six modifications are allowed, if more are given, they are ignored." << endl
			<< "  -dyn_N_term_mod                This is the modification (mass that may be added to each N-terminus)" << endl
			<< "  -dyn_C_term_mod                This is the modification (mass that may be added to each C-terminus" << endl
			<< "  -stat_N_term_mod               This value is the mass that is added to each peptide N-terminus" << endl
			<< "  -stat_C_term_mod               This value is the mass that is added to each peptide C-terminus" << endl
			<< "  -stat_N_term_prot_mod          This value is the mass that is added to each protein N-terminus" << endl
			<< "  -stat_C_term_prot_mod          This value is the mass that is added to each protein C-terminus" << endl
			<< "  -stat_mods                     This value consists of a semicolon-seperated list of amino acids in one letter code" << endl
			<< "                                 and their corrpesponding mass: <AA_1>,<mass_1>;<AA_2>,<mass_2>;..." << endl
			<< "  -partial_sequence              A comma delimited list of amino acid sequences that must occur in the theoretical spectra." << endl
			<< "  -header_filter                 Several elements can be splitted by commas. Each element can be introduced" << endl
			<< "                                 by an exclamation mark (!) meaning that this element must not appear" << endl
			<< "                                 in the header of a protein or the protein will be skipped. This test is done first." << endl
			<< "                                 Next, all other elements are tested. The protein is processed" << endl
			<< "                                 if one filter string matches the header string." << endl
			<< "                                 A filter string may contain a tilde (~). This is replaced by a blank during comparison." << endl
			<< "  -keep_out_files                If set to 1, the Seuest .out-files are not removed (default for -Sequest_out)" << endl
			<< "  -keep_dta_files                If set to 1, the dta-files that were created from the mzXML-files are not removed" << endl
			<< "                                 (default for -Sequest_in)" << endl
			<< "  -contactName                   " << endl
			<< "  -contactInstitution            " << endl
			<< "  -contactInfo                   " << endl;
		}


		void setOptionsAndFlags_()
		{
			flags_["-show_enzyme_numbers"] = "show_enzyme_numbers";
			flags_["-create_dtas"] = "create_dtas";
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
			options_["-sequest_in"] = "sequest_in";
			options_["-sequest_in_dir_network"] = "sequest_in_dir_network";
			options_["-db"] = "db";
			options_["-db_dir_network"] = "db_dir_network";
			options_["-snd_db"] = "snd_db";
			options_["-snd_db_dir_network"] = "snd_db_dir_network";
			options_["-sequest_computer"] = "sequest_computer";
			flags_["-Sequest_in"] = "Sequest_in";
			flags_["-Sequest_out"] = "Sequest_out";
			options_["-in"] = "in";
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

		inline void ensurePathChar(string& path, char path_char = '/')
		{
			if ( !path.empty() && (string("/\\").find(path[path.length()-1], 0) == string::npos) ) path.append(1, path_char);
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

		bool correctNetworkPath(String& network_path, unsigned int backslashes = 2)
		{
			unsigned int pos = 0;
			while ( (pos < network_path.length()) && (network_path[pos] == '\\') ) ++pos;
			if ( pos < backslashes ) network_path.insert(network_path.begin(), backslashes-pos, '\\');
			else network_path.erase(0, pos-backslashes);
			if ( network_path.length() < backslashes+1 ) return false;
			return true;
		}
		
		long fsize(const string& filename)
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
			const string& filename)
		{
			return ( fsize(filename) == 0 );
		}

		String getRandomFilename()
		{
			String s;
			char c = 0;
			/*
			48-57 numbers
			65-90 capitals
			97-122 lower case
			*/

			QFileInfo file_info;
			do
			{
				if ( (c > 47 && c < 58 ) || (c > 64 && c < 91 ) || (c > 96 && c < 123 ) ) ++c;
				else if ( c < 48 ) c = 48;
				else if ( c > 57 && c < 65 ) c = 65;
				else if ( c > 90 && c < 97 ) c = 97;
				else if ( c > 122 )
				{
					s.append("Z");
					c = 48;
				}
				
				file_info.setFile(String(s + String(c)));
			}
			while ( file_info.exists() );

			s.append(1, c);
			return s;
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
			vector< pair< String, vector< double > > >& filenames_and_precursor_retention_times,
			bool make_dtas = true)
		throw (Exception::UnableToCreateFile)
		{
			DTAFile dtafile;
			String filename;
			unsigned int scan_number = 0;
			unsigned int msms_spectra = 0;
			vector< double > retention_times;
			
			for ( MSExperiment<>::Iterator spec_i = msexperiment.begin(); spec_i != msexperiment.end(); ++spec_i )
			{
				++scan_number;
				if ( (spec_i->getMSLevel() == 2) && (!spec_i->empty()) )
				{
					++msms_spectra;
					if ( spec_i->getPrecursorPeak().getCharge() )
					{
						if ( make_dtas )
						{
							++dtas;
							filename = common_name + "." + String(scan_number) + "." + String(spec_i->getPrecursorPeak().getCharge()) + ".dta" + String( (unsigned long long int) (dtas / max_dtas_per_run) );
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
								++dtas;
								spec_i->getPrecursorPeak().setCharge(i);
								filename = common_name + "." + String(scan_number) + "." + String(spec_i->getPrecursorPeak().getCharge()) + ".dta" + String( (unsigned long long int) (dtas / max_dtas_per_run) );
								dtafile.store(filename, *spec_i);
							}
							retention_times.push_back(spec_i->getRetentionTime());
						}
						spec_i->getPrecursorPeak().setCharge(0);
					}
				}
			}
			
			filenames_and_precursor_retention_times.push_back(make_pair(filename, retention_times));
			
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
			
			vector< String > substrings, spectra;
			String string_buffer, string_buffer2;
			double double_buffer;
			int	int_buffer;
			
			// (1.1) variables for the analysis
			float p_value = 0.05;
			unsigned int prob_charge = 1;
			
			vector< pair< String, vector< double > > > filenames_and_precursor_retention_times;

			//-------------------------------------------------------------
			// (2) parsing and checking parameters
			//-------------------------------------------------------------

			// (2.0) variables for running the program
			Sequest_in = getParamAsBool_("Sequest_in");
			Sequest_out = getParamAsBool_("Sequest_out");

			// a 'normal' sequest run corresponds to both Sequest_in and Sequest_out set
			if ( !Sequest_in && !Sequest_out ) Sequest_in = Sequest_out = true;
			
			logfile = getParamAsString_("log", "temp.sequest.log");

			int_buffer = getParamAsInt_("prob_charge");
			if ( int_buffer < 0 )
			{
				writeLog_("Maximal charge to test is less than zero. Aborting!");
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}
			else if ( int_buffer ) prob_charge = int_buffer;
			else prob_charge = 1;

			// only show the available enzymes, then quit
			if ( getParamAsBool_("show_enzyme_numbers") )
			{
				writeLog_("Option show_enzyme_numbers chosen. Aborting.");
				return EXECUTION_OK;
			}

			// the spectra
			string_buffer = getParamAsString_("in");
			if ( string_buffer.empty() )
			{
				writeLog_("No spectrum file specified. Aborting!");
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}
			else
			{
				string_buffer.split(',', spectra);
				if ( spectra.empty() ) spectra.push_back(string_buffer);
			}

			keep_out_files = getParamAsBool_("keep_out_files");
			if ( Sequest_out && !Sequest_in ) keep_out_files = true;

			keep_dta_files = getParamAsBool_("keep_dta_files");
			if ( Sequest_in && !Sequest_out ) keep_dta_files = true;
			
			temp_data_dir = getParamAsString_("temp_data_dir");
			ensurePathChar(temp_data_dir);
			if ( temp_data_dir.empty() )
			{
				writeLog_("No directory for temporary files given. Aborting!");
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}

			// only create dta files from the mzXML files, then quit
			if ( getParamAsBool_("create_dtas") )
			{
				bool make_dtas = true;
				// if there are already .dta files in the folder, stop the adapter
				QDir qdir(temp_data_dir, "*.dta", QDir::Name, QDir::Files);
				if ( !qdir.entryList().empty() )
				{
					writeLog_("There are already dta files in directory " + temp_data_dir + ". Aborting!");
					deleteTempFiles(input_filename, logfile);
					return UNKNOWN_ERROR;
				}
				
				MSExperiment<> msexperiment;
				QFileInfo file_info;
				unsigned int msms_spectra_in_file;
				unsigned int msms_spectra_altogether = 0;
				if ( make_dtas )
				{
					writeLog_("creating dta files");
				}
				dtas = 0;
				for ( vector< String >::const_iterator spec_i = spectra.begin(); spec_i != spectra.end(); ++spec_i )
				{
					file_info.setFile(*spec_i);
					String common_name = temp_data_dir + string(file_info.fileName().ascii());
					
					try
					{
						MzXMLFile().load(*spec_i, msexperiment);
					}
					catch ( Exception::ParseError pe )
					{
						writeLog_("Error loading mzXML file. Aborting!");
						printUsage_();
						return PARSE_ERROR;
					}
					
					msms_spectra_in_file = MSExperiment2DTAs(msexperiment, common_name, prob_charge, filenames_and_precursor_retention_times, make_dtas);
					writeLog_(String(msms_spectra_in_file) + " MS/MS spectra in file " + file_info.fileName().ascii());
	
					msms_spectra_altogether += msms_spectra_in_file;
				}
	
				if ( !msms_spectra_altogether )
				{
					writeLog_("No MS/MS spectra found in any of the mzXML files. Aborting!");
					printUsage_();
					return UNKNOWN_ERROR;
				}

				return EXECUTION_OK;
			}
			
			temp_data_dir_win = getParamAsString_("temp_data_dir_win");
			ensurePathChar(temp_data_dir_win, '\\');
			
			if ( !isWinFormat(temp_data_dir_win) )
			{
				writeLog_("Windows path for the directory for temporary files has wrong format. Aborting!");
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}
			temp_data_dir_network = getParamAsString_("temp_data_dir_network");
			if ( temp_data_dir_network.empty() )
			{
				writeLog_("Network path for the directory for temporary files is empty. Aborting!");
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}
				if ( !correctNetworkPath(temp_data_dir_network) )
			{
				writeLog_(temp_data_dir_network + "is no network path. Aborting!");
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}
			
			contact_person.setName(getParamAsString_("contactName", "unknown"));
			contact_person.setInstitution(getParamAsString_("contactInstitution", "unknown"));
			contact_person.setContactInfo(getParamAsString_("contactInfo"));

			if ( Sequest_in )
			{
				input_filename = getParamAsString_("sequest_in");
				if ( input_filename.empty() )
				{
					if ( !Sequest_out ) // if only Sequest_in is set, a name has to be given
					{
						writeLog_("No input file specified. Aborting!");
						printUsage_();
						return ILLEGAL_PARAMETERS;
					}
					else
					{
						input_filename = temp_data_dir + "temp.sequest.in";
						input_file_dir_network = temp_data_dir_network;
					}
				}
				else input_file_dir_network = getParamAsString_("sequest_in_dir_network");
			}
			
			if ( Sequest_in && Sequest_out )
			{
				if ( !correctNetworkPath(input_file_dir_network) )
				{
					writeLog_(input_file_dir_network + "is no network path. Aborting!");
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				QFileInfo file_info(input_filename);
				string_buffer = file_info.fileName().ascii();
				if ( !input_file_dir_network.hasSuffix(string_buffer) ) input_file_dir_network.append("\\" + string_buffer);
				
				user = getParamAsString_("user");
				if ( user.empty() )
				{
					writeLog_("No user name for Sequest computer given. Aborting!");
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				
				password = getParamAsString_("password");
			
				sequest_dir_win = getParamAsString_("sequest_dir_win");
				ensurePathChar(sequest_dir_win, '\\');
				if ( !isWinFormat(sequest_dir_win) && !sequest_dir_win.empty() )
				{
					writeLog_("Windows path for the SEQUEST working directory has wrong format. Aborting!");
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
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
			}
			
			if ( logfile == temp_data_dir + "sequest.log")
			{
				writeLog_("The logfile must not be named " + temp_data_dir + "sequest.log. Aborting!");
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
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}

			snd_database = getParamAsString_("snd_db");
			
			if ( Sequest_in )
			{
				database_dir_network = getParamAsString_("db_dir_network");
				if ( !correctNetworkPath(database_dir_network) )
				{
					writeLog_(database_dir_network + "is no network path. Aborting!");
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				QFileInfo file_info(database);
				string_buffer = file_info.fileName().ascii();
				if ( !database_dir_network.hasSuffix(string_buffer) ) database_dir_network.append("\\" + string_buffer);
				sequest_infile.setDatabase(database_dir_network);

				if ( !snd_database.empty() )
				{
					snd_database_dir_network = getParamAsString_("snd_db_dir_network");
					if ( !correctNetworkPath(snd_database_dir_network) )
					{
						writeLog_(snd_database_dir_network + "is no network path. Aborting!");
						printUsage_();
						return ILLEGAL_PARAMETERS;
					}
					QFileInfo file_info(snd_database);
					string_buffer = file_info.fileName().ascii();
// 					if ( !snd_database_dir_network.hasSuffix(string_buffer) ) snd_database_dir_network.append("\\" + string_buffer);
					sequest_infile.setSndDatabase(snd_database_dir_network);
				}
				
				double_buffer = getParamAsDouble_("pep_mass_tol", -1);
				if ( double_buffer == -1 )
				{
					writeLog_("No peptide mass tolerance specified. Aborting!");
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else if ( double_buffer < 0 )
				{
					writeLog_("Peptide mass tolerance < 0. Aborting!");
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setPeptideMassTolerance(double_buffer);
				
				double_buffer = getParamAsDouble_("frag_ion_tol", -1);
				if ( double_buffer == -1 )
				{
					writeLog_("No fragment ion tolerance specified. Aborting!");
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else if ( double_buffer < 0 )
				{
					writeLog_("Fragment ion tolerance < 0. Aborting!");
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setFragmentIonTolerance(double_buffer);
				
				double_buffer = getParamAsDouble_("match_peak_tol", -1);
				if ( double_buffer == -1 )
				{
					writeLog_("No match peak tolerance specified. Aborting!");
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else if ( double_buffer < 0 )
				{
					writeLog_("Match peak tolerance < 0. Aborting!");
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setMatchPeakTolerance(double_buffer);
				
				double_buffer = getParamAsDouble_("ion_cutoff", 0);
				if ( double_buffer < 0 )
				{
					writeLog_("Ion cutoff < 0. Aborting!");
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setIonCutoffPercentage(double_buffer);
				
				int_buffer = getParamAsInt_("pep_mass_unit", 0);
				if ( (int_buffer < 0) || (int_buffer > max_peptide_mass_units) )
				{
					writeLog_("Illegal peptide mass unit (not in [0,2]). Aborting!");
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setPeptideMassUnits(int_buffer);
				
				int_buffer = getParamAsInt_("num_results", 0);
				if ( (int_buffer < 1) )
				{
					writeLog_("Illegal number of results (< 1). Aborting!");
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setNumOutputLines(int_buffer);
				
				string_buffer = getParamAsString_("enzyme_info");
				if ( !string_buffer.empty() )
				{
					string_buffer.split(';', substrings);
					if ( substrings.empty() ) substrings.push_back(string_buffer);
					
					vector< String > enzyme_info;
					for ( vector< String >::iterator einfo_i = substrings.begin(); einfo_i != substrings.end(); ++ einfo_i )
					{
						einfo_i->split(',', enzyme_info);
						if ( (enzyme_info.size() < 3) || (enzyme_info.size() > 4) )
						{
							writeLog_("Illegal number of informations for enzyme (not in [3,4]). Aborting!");
							printUsage_();
							return ILLEGAL_PARAMETERS;
						}
						if ( !((enzyme_info[1] == "0") || (enzyme_info[1] == "1"))  )
						{
							writeLog_("Cut direction for enzyme not specified correctly (has to be 1 (N to C)) or 0 (C to N))). Aborting!");
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
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setMinimumProteinMass(double_buffer);
				
				double_buffer = getParamAsDouble_("max_prot_mass");
				if ( double_buffer < 0 )
				{
					writeLog_("Illegal maximum protein mass (< 0). Aborting!");
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setMaximumProteinMass(double_buffer);
				
				int_buffer = getParamAsInt_("max_num_dif_AA_per_mod", -1);
				if ( int_buffer < 0 )
				{
					writeLog_("No maximum number of modified amino acids per different modification. Aborting!");
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setMaxNumDifAAPerMod(int_buffer);
				
				int_buffer = getParamAsInt_("max_num_dif_mods_per_peptide", -1);
				if ( int_buffer < 0 )
				{
					writeLog_("No maximum number of differential modifications per peptide. Aborting!");
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setMaxNumModsPerPeptide(int_buffer);
				
				int_buffer = getParamAsInt_("nuc_reading_frame", 0);
				if ( (int_buffer < 0) || (int_buffer > 9) )
				{
					writeLog_("Illegal number for nucleotide reading frame. Aborting!");
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setNucleotideReadingFrame(int_buffer);
				
				int_buffer = getParamAsInt_("max_num_int_cleav_sites", 0);
				if ( int_buffer < 0 )
				{
					writeLog_("Illegal number of maximum internal cleavage sites. Aborting!");
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setMaxNumInternalCleavageSites(int_buffer);
				
				int_buffer = getParamAsInt_("match_peak_count", 0);
				if ( (int_buffer < 0) && (int_buffer > 5) )
				{
					writeLog_("Illegal number of auto-detected peaks to try matching. Aborting!");
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else sequest_infile.setMatchPeakCount(int_buffer);
				
				int_buffer = getParamAsInt_("match_peak_allowed_error", 0);
				if ( int_buffer < 0 )
				{
					writeLog_("Illegal number of allowed errors in matching auto-detected peaks. Aborting!");
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
				if ( (string_buffer.size() != 3) || (string_buffer2.find(string_buffer[0], 0) == string::npos) || (string_buffer2.find(string_buffer[1], 0) == string::npos) || (string_buffer2.find(string_buffer[2], 0) == string::npos) )
				{
					writeLog_("Neutral losses for ABY-ions not given (or illegal values given). Aborting!");
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else
				{
					string_buffer.insert(2, 1, ' ');
					string_buffer.insert(1, 1, ' ');
					sequest_infile.setNeutralLossesForIons(string_buffer);
				}
				
				string_buffer = getParamAsString_("ion_series_weights");
				string_buffer.split(',', substrings);
				if ( substrings.size() != 9 )
				{
					writeLog_("Weights for ion series not given (or illegal values given). Aborting!");
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else
				{
					for ( vector< String >::iterator s_i = substrings.begin(); s_i != substrings.end(); ++s_i )
					{
						// the values are expected to be float, otherwise they will be seen as 0!
						double_buffer = atof(s_i->c_str());
						if ( (double_buffer < 0) || (double_buffer > 1) )
						{
							writeLog_("Illegal weights for ion series given. Aborting!");
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
					for ( vector< String >::iterator s_i = substrings.begin(); s_i != substrings.end(); ++s_i )
					{
						if ( sscanf(s_i->c_str(), "%f,%40s", &f_buffer, c_buffer) != 2 )
						{
							writeLog_("Illegal number of parameters for dynamic modification given. Aborting!");
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
					
					for ( vector< String >::iterator s_i = substrings.begin(); s_i != substrings.end(); ++s_i )
					{
						if ( (*s_i)[1] != ',' )
						{
							writeLog_("Unexpected format for static modification found. Aborting!");
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
				string_buffer = getParamAsString_("out");
				if ( string_buffer.empty() )
				{
					writeLog_("No output file specified. Aborting!");
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else output_filename = string_buffer;
			
				double_buffer = getParamAsDouble_("p_value", -1.0);
				if ( double_buffer != -1.0 ) p_value = double_buffer;
				if ( (p_value <= 0) || (p_value > 1) )
				{
					writeLog_("P-value not in (0,1]. Aborting!");
					printUsage_();
					return ILLEGAL_PARAMETERS;
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

				file.setName(temp_data_dir + batch_filename);
				file.open(IO_WriteOnly);
				if ( !file.isWritable() )
				{
					throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, temp_data_dir + batch_filename);
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

			if ( make_dtas ) // if there are already .dta files in the folder, stop the adapter
			{
				QDir qdir(temp_data_dir, "*.dta", QDir::Name, QDir::Files);
				if ( !qdir.entryList().empty() )
				{
					writeLog_("There are already dta files in directory " + temp_data_dir + ". Aborting!");
					deleteTempFiles(input_filename, logfile);
					return UNKNOWN_ERROR;
				}
			}
			
			MSExperiment<> msexperiment;
			unsigned int msms_spectra_in_file;
			unsigned int msms_spectra_altogether = 0;
			if ( make_dtas ) 
			{
				writeLog_("creating dta files");
			}
			dtas = 0;
			for ( vector< String >::const_iterator spec_i = spectra.begin(); spec_i != spectra.end(); ++spec_i )
			{
				file_info.setFile(*spec_i);
				String common_name = temp_data_dir + string(file_info.fileName().ascii());
				
				try
				{
					MzXMLFile().load(*spec_i, msexperiment);
				}
				catch ( Exception::ParseError pe )
				{
					writeLog_("Error loading mzXML file. Aborting!");
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
				printUsage_();
				return UNKNOWN_ERROR;
			}

			// (3.2.3) running the program
			if ( Sequest_in && Sequest_out )
			{
				// creating a batch file for windows (command doesn't accept commands that are longer than 256 chars)
				String sequest_screen_output = getRandomFilename(); // direct the screen-output to a file
				ofstream batchfile(String(temp_data_dir + batch_filename).c_str());
				if ( !batchfile )
				{
					throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, temp_data_dir + batch_filename);
				}
				String call = "rdesktop -u ";
				call.append(user);
				if ( !password.empty() ) call.append(" -p \"" + password + "\"");
				call.append(" -s cmd\\ /C\\ \"");
				call.append(" net use " + temp_data_dir_win.substr(0,2) + " \\\\" + temp_data_dir_network + " && ");

				batchfile << String(" cd " + temp_data_dir_win + " && " + temp_data_dir_win.substr(0,2));
				
				for ( unsigned long long int i = 0; i <= (unsigned long long int) (dtas / max_dtas_per_run); ++i )
				{
					batchfile << String(" && " + sequest_dir_win + "sequest.exe -P" + input_file_dir_network + " " + temp_data_dir_network + "\\*.dta" + String(i) + " >  " +  temp_data_dir_network +"\\" + sequest_screen_output + " && move sequest.log sequest.log" + String(i));
				}
				batchfile << String(" && " + sequest_dir_win.substr(0,2) + " &&");
				batchfile << String(" net use /delete " + temp_data_dir_win.substr(0,2));
				batchfile << " && logoff";
				batchfile.close();
				batchfile.clear();
				
				call.append(temp_data_dir_win + batch_filename + "\" " + sequest_computer);
				writeLog_("System call: " + call);
				int status = system(call.c_str());
				remove(sequest_screen_output.c_str());
				
				if ( status != 0 )
				{
					writeLog_("Sequest problem. Aborting! (Details can be seen in the logfile: \"" + logfile + "\")");
					deleteTempFiles(input_filename, logfile);
					return EXTERNAL_PROGRAM_ERROR;
				}

				for ( unsigned long long int i = 0; i <= (unsigned long long int) (dtas / max_dtas_per_run); ++i )
				{
					ifstream sequest_log(string(temp_data_dir + "sequest.log" + String(i)).c_str()); // write sequest log to logfile
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
						remove(string(temp_data_dir + "sequest.log" + String(i)).c_str());
					}
				}
			}
			
			if ( Sequest_out )
			{
				if ( !keep_dta_files ) // remove all dtas
				{
					for ( unsigned long long int i = 0; i <= (unsigned long long int) (dtas / max_dtas_per_run); ++i )
					{
						writeLog_("removing dta files");
						QDir qdir(temp_data_dir, "*.dta" + String(i) , QDir::Name, QDir::Files);
						QStringList qlist = qdir.entryList();
						
						for ( QStringList::const_iterator i = qlist.constBegin(); i != qlist.constEnd(); ++i )
						{
							if ( !File::remove(temp_data_dir + *i) )
							{
								writeLog_(string((*i).ascii()) + "could not be removed!");
							}
						}
					}
				}
				
				AnalysisXMLFile analysisXML_file;
				
				SequestOutfile sequest_outfile;
				vector< Identification > identifications;
				ProteinIdentification protein_identification;
				vector< float > precursor_retention_times, precursor_mz_values;
				
				QDir qdir(temp_data_dir, "*.out", QDir::Name, QDir::Files);
				QStringList qlist = qdir.entryList();

				if ( qlist.isEmpty() )
				{
					writeLog_("No .out identified. Aborting!");
					qlist.clear();
					deleteTempFiles(input_filename, logfile);
					return UNKNOWN_ERROR;
				}

				for ( QStringList::const_iterator i = qlist.constBegin(); i != qlist.constEnd(); ++i )
				{
					sequest_outfile.load(temp_data_dir + string((*i).ascii()), identifications, protein_identification, precursor_retention_times, precursor_mz_values,  p_value, database, snd_database);
				}

				sort(filenames_and_precursor_retention_times.begin(), filenames_and_precursor_retention_times.end(), SortRetentionTimes()); // sort the retention times, so they have the same order like the corresponding .out files

				// save the retention times in precursor_retention_times
				for( vector< pair< String, vector< double > > >::iterator pairs_i = filenames_and_precursor_retention_times.begin(); pairs_i != filenames_and_precursor_retention_times.end(); ++pairs_i )
				{
					precursor_retention_times.insert(precursor_retention_times.end(), pairs_i->second.begin(), pairs_i->second.end());
					pairs_i->second.clear();
				}
				
				vector< ProteinIdentification > pis;
				pis.push_back(protein_identification);

				analysisXML_file.store(output_filename, pis, identifications, precursor_retention_times, precursor_mz_values, contact_person);
				
				// remove all outs
				if ( !keep_out_files )
				{
					qdir.setPath(temp_data_dir);
					qlist = qdir.entryList("*.out", QDir::Files, QDir::Name);
					QFile qfile;
					
					for ( QStringList::const_iterator i = qlist.constBegin(); i != qlist.constEnd(); ++i )
					{
						if ( !qfile.remove(QString(temp_data_dir.c_str() + *i)) )
						{
							writeLog_(string((*i).ascii()) + "could not be removed!");
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
