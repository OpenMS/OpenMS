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
// $Id: InspectAdapter.C,v 1.0 2006/07/12 15:58:59 martinlangwisch Exp $
// $Author: martinlangwisch $
// $Maintainer: Martin Langwisch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/FORMAT/AnalysisXMLFile.h>
#include <OpenMS/FORMAT/InspectInfile.h>
#include <OpenMS/FORMAT/InspectOutfile.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Date.h>
#include <OpenMS/METADATA/ContactPerson.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <stdlib.h>
#include <vector>

#include <qfileinfo.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

// @cond TOPPCLASSES 


/**
	@page InspectAdapter InspectAdapter
	
	@brief Identifies peptides in MS/MS spectra via Inspect.
	
	This wrapper component serves for getting peptide identifications
	for MS/MS spectra. The wrapper can be executed in three different
	modes:
	<ol>	
				<li>
				The whole process of identification via Inspect is executed. 
				Inputfile is a file (or directory with files) containing the MS/MS spectra
				(Supported spectrum file formats are .mzXML, .mzData, .ms2, .mgf, .dta,
				and .pkl. Note that multiple spectra in a single .pkl or .dta file are not supported.)
				for which identifications are to be found
				and one ore more databases in either trie, FASTA or Swissprot format containing
				the possible proteins.
				The results are written as a analysisXML output file. This mode is selected
			 	by default.
			 	</li>
				
				<li>
				Only the first part of the identification process is performed.
				This means that an Inspect input file is generated and the given databases are 
				converted (if necessary) and merged into a trie database. This file can be used
				directly with Inspect whereas the created database and the spectrum file(s)
				have to remain at the given positions.
				Calling an Inspect process should look like the following:				
				
				<ul>	
					<li>
						@code ./inspect -i  inputfilename -o outputfilename  @endcode
					</li>
				</ul>
				(Inspect may be run from anywhere adding '-r inspect_directory' but at
				the current version (20060620) this does not work properly)
				Consult your Inspect reference manual for further details.
				
				This mode is selected by the <b>-Inspect_in</b> option in the command line.
				</li>
				
				<li>
				Only the second part of the identification process is performed.
				This means that inspect is run and the output is translated into analysisXML.
				
				This mode is selected by the <b>-Inspect_out</b> option in the command line.
				</li>
	</ol>
	
	@todo look for possible crash codes of inspect and catching them; extract by-ions, read PTMs from ini file and from input, compute protein score?, catch exceptions to close files
	
	@ingroup TOPP
*/

// We do not want this class to show up in the docu -> @cond
/// @cond 

class TOPPInspectAdapter
	: public TOPPBase
{
	public:
		TOPPInspectAdapter()
			: TOPPBase("InspectAdapter")
		{}
	
	protected:
		void printToolUsage_() const
		{
			cerr	<< endl
						<< getToolName() << " -- annotates MS/MS spectra using Inspect" << endl
						<< "Version: " << VersionInfo::getVersion() << endl
						<< endl
						<< "Usage:" << endl
						<< " " << getToolName() << " [options]" << endl
						<< endl
						<< "Options are:" << endl
						<< "  -in <file>          Inspect input file" << endl
						<< "  -out <file>         output file in analysisXML" << endl
						<< "  -o <file>           direct output file from inspect" << endl
						<< "  -Inspect_in         if this flag is set the InspectAdapter will create an Inspect Input file" << endl
						<< "                      if only Inspect_in is set, a name for the trie database (see below) has to be given!" << endl
						<< "  -Inspect_out        if this flag is set the InspectAdapter will read in an Inspect Input file and write an analysisXML file." << endl
						<< "  -inspect_dir        the name of the InsPecT directory." << endl
						<< "  -temp_data_dir      the name of the directory where the temporary data will be stored." << endl
						<< "  -spectra <file>     the spectrum file OR directory to search (every file in that directory will be searched(non-recursively)" << endl
						<< "                      supported spectrum file formats are .mzXML, .mzData, .ms2, dta, and .pkl" << endl
						<< "                      multiple spectra in one .dta file are not supported" << endl
						<< "  -trie_dbs <file1>,<file2>,...      names of a databases (.trie file) to search ()" << endl
						<< "  -dbs <file>;tax1,<file2>;tax2,...  names of a other databases to search (currently FASTA and SwissProt are supported)" << endl
						<< "                                     tax - the desired taxonomy, if not given for a database, all entries are taken." << endl
						<< "  -make_trie_db       if set, the InspectAdapter will generate one trie database from all given databases." << endl
						<< "                      if you do not use this switch you may only use one FASTA database XOR one trie database" << endl
						<< "  -mods [<MASS1>,<RESIDUES1>,<TYPE1>,<NAME1>];[<MASS2>,<RESIDUES2>,<TYPE2>,<NAME2>]" << endl
						<< "                      modifications i.e. [80,STY,opt,phosphorylation] (default read from INI file)" << endl
						<< "                      MASS and RESIDUES are mandatory, make sure the modifications are seperated by a semicolon!" << endl
						<< "                      Valid values for \"type\" are \"fix\", \"cterminal\", \"nterminal\", and \"opt\" (the default)." << endl
						<< "  -blind              perform a blind search (allowing arbitrary modification masses), as this is slower than the normal search" << endl
						<< "                      A normal search is performed in advance to gain a smaller database." << endl
						<< "                      This search can only be run in full mode." << endl
						<< "  -blind_only         like blind but no prior search is performed to reduce the database size" << endl;
		}


		void printToolHelpOpt_() const
		{
			cerr	<< endl
						<< "  -instr              the instrument that was used to measure the spectra (default read from INI file)" << endl
						<< "                      (If set to QTOF, uses a QTOF-derived fragmentation model, and does not attempt to correct the parent mass.)" << endl
						<< "  -PM_tol             the precursor mass tolerance (default read from INI file)" << endl
						<< "  -ion_tol            the peak mass tolerance (default read from INI file)" << endl
						<< "  -protease           the name of a protease. \"Trypsin\", \"None\", and \"Chymotrypsin\" are the available values." << endl
						<< "                      The first four	characters of the name should be unique." << endl
						<< "  -max_mods_pp        number of PTMs permitted in a single peptide. (default: read from INI file)" << endl
						
						<< "  -p_value            annotations with inferior p-value are ignored. Default is 0.05" << endl
						<< "  -score_value        annotations with inferior score-value are ignored. Default is 1." << endl
						<< "                      (this is a workaround because sometimes inspect produces only nan as p-value;" << endl
						<< "                      a hit with score of >=1 is supposed to be good)" << endl
						<< "  -p_value_blind      used when generating the minimized database for blind search" << endl
						<< "  -score_value_blind  annotations with inferior score-value are ignored. Default is 1 (see score_value)." << endl
						<< "  -min_spp            used when generating the minimized database for blind search " << endl
						<< "                      the minimum number of spectra a protein has to annotate in order to add it to the filtered database " << endl
						<< "                      default is #spectra / #proteins * 2" << endl
						<< "  -maxptmsize         for blind search, specifies the maximum modification size (in Da) to consider (default read from INI file)" << endl
						<< "  -jumpscores <file>  file to specify PTM frequencies, for use in tag generation. This is more accurate tagging than the" << endl
						<< "                      default behavior (where tags can contain any PTM), but requires the creation of the jump frequency file" << endl
						<< "  -multicharge        attempt to guess the precursor charge and mass, and consider multiple charge states if feasible" << endl
						<< "  -twopass            use two-pass search. The first pass uses fewer tags, and produces a list of proteins" << endl
						<< "                      which are re-searched in the second pass" << endl
						<< "  -TagCountA          number of tags to generate for the first pass of a two-pass search" << endl
						<< "  -TagCountB          number of tags to generate for the second pass of a two-pass search" << endl
						<< "                      OR the number of tags to use in a one-pass search" << endl
						<< "  -cmn_conts          add the proteins from CommonContaminents.fasta (in inspect path) to the search database" << endl
						<< "  -no_tmp_dbs         no temporary databases are used" << endl
						<< "  -new_db             name of the trie database (given databases are converted and merged to one trie database)." << endl
						<< "                      This has to be set if no_tmp_dbs is set! If the name does not end with \".trie\"" << endl
						<< "                      it is extended accordingly." << endl
						<< "                      An index file with the same name but extension \".index\" will be created." << endl
						<< "  -snd_db             name of the minimized trie database generated when using blind mode." << endl
						<< "                      This has to be set if no_tmp_dbs is set!" << endl;
						//<< "  -contact		 name of the contact person" << endl
		}


		void setOptionsAndFlags_()
		{
			options_["-inspect_dir"] = "inspect_dir";
			options_["-temp_data_dir"] = "temp_data_dir";
			flags_["-Inspect_in"] = "Inspect_in";
			options_["-spectra"] = "spectra";
			options_["-trie_dbs"] = "trie_dbs";
			options_["-dbs"] = "dbs";
			options_["-new_db"] = "new_db";
			options_["-snd_db"] = "snd_db";
			options_["-tax"] = "tax";
			options_["-protease"] = "protease";
			options_["-jumpscores"] = "jumpscores";
			options_["-instrument"] = "instrument";
			options_["-mods"] = "mods";
			options_["-max_mods_pp"] = "max_mods_pp";
			options_["-PM_tol"] = "PM_tol";
			options_["-ion_tol"] = "ion_tol";
			flags_["-multicharge"] = "multicharge";
			options_["-TagCountA"] = "TagCountA";
			options_["-TagCountB"] = "TagCountB";
			flags_["-twopass"] = "twopass";
			flags_["-Inspect_out"] = "Inspect_out";
			options_["-in"] = "in";
			options_["-out"] = "out";
			options_["-o"] = "o";
			flags_["-blind_only"] = "blind_only";
			options_["-p_value"] = "p_value";
			options_["-p_value_blind"] = "p_value_blind";
			options_["-score_value"] = "score_value";
			options_["-score_value_blind"] = "score_value_blind";
			options_["-min_spp"] = "min_spp";
			options_["-maxptmsize"] = "maxptmsize";
			flags_["-blind"] = "blind";
			flags_["-cmn_conts"] = "cmn_conts";
			flags_["-no_tmp_dbs"] = "no_tmp_dbs";
			flags_["-make_trie_db"] = "make_trie_db";
			//options_["-contact"] = "contact_person";
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
		
		inline bool emptyFile(const string& filename)
		{
			return ( fsize(filename) == 0 );
		}

		string fileContent(const string& filename)
		{
			long size = fsize(filename);
			if ( size != -1 )
			{
				FILE* file = fopen(filename.c_str(), "r");
				char* buffer = new char[size+1];
				buffer[size] = 0;
				fread (buffer, size, 1, file);
				fclose(file);
				string sbuffer = buffer;
				delete(buffer);
				return sbuffer;
			}
			else return string();
		}
		
		// deleting all temporary files
		void deleteTempFiles(const String& input_filename, const String& output_filename, const String& inspect_output_filename, const String& db_filename, const String& idx_filename, const String& snd_db_filename, const String& snd_index_filename, const String& inspect_logfile)
		{
			if ( input_filename.hasSuffix("tmp.inspect.input") ) remove(input_filename.c_str());
			if ( output_filename.hasSuffix("tmp.inspect.output") ) remove(output_filename.c_str());
			if ( inspect_output_filename.hasSuffix("tmp.direct.inspect.output") ) remove(inspect_output_filename.c_str());
			if ( db_filename.hasSuffix("tmp.inspect.db.trie") ) remove(db_filename.c_str());
			if ( idx_filename.hasSuffix("tmp.inspect.db.index") ) remove(idx_filename.c_str());
			if ( snd_db_filename.hasSuffix("tmp.inspect.db.snd.trie") ) remove(snd_db_filename.c_str());
			if ( snd_index_filename.hasSuffix("tmp.inspect.db.snd.index") ) remove(snd_index_filename.c_str());
			if ( inspect_logfile.hasSuffix("tmp.inspect.log") ) remove(inspect_logfile.c_str());
		}

		ExitCodes main_(int , char**)
		{
			//-------------------------------------------------------------
			// (1) variables
			//-------------------------------------------------------------

			InspectInfile inspect_infile;

			// (1.0) general variables
			vector< String > substrings;
			String buffer, db_filename, idx_filename, snd_db_filename, snd_index_filename, common_contaminants_filename, inspect_logfile, logfile;
			ContactPerson contact_person;

			// (1.1) parameter variables
			// (1,1,0) general parameter variables
			String inspect_dir, temp_data_dir;

			// (1.1.1) Inspect_in - writing the inspect input file only and corresponding parameters
			bool Inspect_in = false;
			// (1.1.1.0) mandatory parameters
			String snd_db, snd_db_dir; // at least one of the parameters db or seq_file has to be set
			vector< String >dbs, seq_files, tax; // if several dbs are given, they are merged into one, that is then processed

			// (1.1.1.1) optional parameters
			bool make_trie_db = false;
			
			vector < vector< String > > mod; // some from ini file

			// (1.1.2) Inspect_out - executing the program only and writing xml analysis file and corresponding parameters
			double p_value_threshold = 1.0;
			double score_value_threshold = 1.0;
			bool Inspect_out = false;
			String output_filename;

			// (1.1.3) parameters corresponding to both Inspect_in and Inspect_out
			String input_filename, inspect_output_filename;

			// (1.1.4) blind_only - running inspect in blind mode only and corresponding parameters
			bool blind_only = false;

			// (1.1.5) blind - running inspect in blind mode after running a normal mode to minimize the database
			bool blind = false;
			double cutoff_p_value = 0.05;
			double cutoff_score_value = 1.0;
			int min_annotated_spectra_per_protein = -1;

			// (1.1.6) no_common_contaminants - whether to include the proteins in commonContaminants.fasta
			bool no_common_contaminants = true;

			// (1.1.7) no_tmp_dbs - whether to use temporary database files or to save them (faster if they are used more than once)
			bool no_tmp_dbs = false;


			//-------------------------------------------------------------
			// (2) parsing and checking parameters
			//-------------------------------------------------------------
			// (2.0) general variables
			Inspect_in = getParamAsBool_("Inspect_in", false);
			Inspect_out = getParamAsBool_("Inspect_out", false);
			
			// a 'normal' inspect run corresponds to both Inspect_in and Inspect_out set
			if ( !Inspect_in && !Inspect_out ) Inspect_in = Inspect_out = true;
			
			contact_person.setName(getParamAsString_("contactName", "unknown"));
			contact_person.setInstitution(getParamAsString_("contactInstitution", "unknown"));
			contact_person.setContactInfo(getParamAsString_("contactInfo"));
			
			// (2.1) parameter variables
			// (2.1,0) general parameter variables
			inspect_dir = getParamAsString_("inspect_dir");
			if ( ((Inspect_in && Inspect_out) || (Inspect_in && blind)) && inspect_dir.empty() )
			{
				writeLog_("No inspect directory file specified. Aborting!");
				cout << "No inspect directory specified. Aborting!" << endl;
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}
			inspect_infile.ensurePathChar(inspect_dir);
			
			common_contaminants_filename = inspect_dir + "CommonContaminants.fasta";
			
			temp_data_dir = getParamAsString_("temp_data_dir");
			if ( ((Inspect_in && Inspect_out) || (Inspect_in && blind)) && temp_data_dir.empty() )
			{
				writeLog_("No directory for temporary files specified. Aborting!");
				cout << "No directory for temporary files specified. Aborting!" << endl;
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}
			inspect_infile.ensurePathChar(temp_data_dir);
			
			// (2.1.3) parameters corresponding to both Inspect_in and Inspect_out
			buffer = getParamAsString_("o");
			if ( !Inspect_in && Inspect_out )
			{
				if ( buffer.empty() )
				{
					writeLog_("No InsPecT output file specified. Aborting!");
					cout << "No InsPecT output file specified. Aborting!" << endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else inspect_output_filename = buffer;
			}
			else if ( (Inspect_in && Inspect_out) || (Inspect_in && blind) )
			{
				if ( buffer.empty() ) inspect_output_filename = temp_data_dir + "tmp.direct.inspect.output";
				else inspect_output_filename = buffer;
			}
			
			buffer = getParamAsString_("in");
			// if only one mode is used, a name for the input file has to be given
			if ( Inspect_in != Inspect_out )
			{
				if ( buffer.empty() )
				{
					writeLog_("No input file specified. Aborting!");
					cout << "No input file specified. Aborting!" << endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				else input_filename = buffer;
			}
			// if both flags are set a the input file may be temporary
			else
			{
				if ( buffer.empty() ) input_filename = temp_data_dir + "tmp.inspect.input";
				else input_filename = buffer;
			}
			
			blind_only = getParamAsBool_("blind_only", false);
			
			// (2.1.1) Inspect_in - writing the inspect input file only and corresponding parameters
			if ( Inspect_in )
			{
				// (2.1.1.0) mandatory parameters
				inspect_infile.setSpectra(getParamAsString_("spectra"));
				if ( inspect_infile.getSpectra().empty() )
				{
					writeLog_("No spectrum file specified. Aborting!");
					cout << "No spectrum file specified. Aborting!" << endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}

				buffer = getParamAsString_("trie_dbs");
				if ( !buffer.empty() )
				{
					// get the single databases
					buffer.split(',', dbs);
					if ( dbs.empty() ) dbs.push_back(buffer);
				}
				
				buffer = getParamAsString_("dbs");
				if ( !buffer.empty() )
				{
					// get the single sequence files
					buffer.split(',', seq_files);
					if ( seq_files.empty() ) seq_files.push_back(buffer);
					
					// get the corresponding taxonomies
					for ( vector< String >::iterator i = seq_files.begin(); i != seq_files.end(); ++i)
					{
						substrings.clear();
						i->split(';', substrings);
						if ( !substrings.empty() )
						{
							tax.push_back(substrings[1]);
							buffer = String(tax.back());
							buffer.toUpper();
							if ( buffer == "ALL" ) tax.back() = "None";
							*i = substrings[0];
						}
						else tax.push_back("None");
					}
				}
				
				// at least one of the parameters db or seq_file has to be set
				if ( dbs.empty() && seq_files.empty() )
				{
					writeLog_("No database or sequence file specified. Aborting!");
					cout << "No database or sequence file specified. Aborting!" << endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				
				// (2.1.6) no_common_contaminants - whether to include the proteins in commonContaminants.fasta
				no_common_contaminants = !getParamAsBool_("cmn_conts", false);
				make_trie_db = getParamAsBool_("make_trie_db", false);
				if ( !make_trie_db && ((!dbs.empty()) + (!seq_files.empty()) + (!no_common_contaminants) >1) )
				{
					writeLog_("Too many databases (make_trie_db not set). Aborting!");
					cout << "Too many databases (make_trie_db not set). Aborting!" << endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				
				no_tmp_dbs = getParamAsBool_("no_tmp_dbs", false);
				if ( !make_trie_db && !dbs.empty() ) db_filename = dbs[0];
				else if ( make_trie_db )
				{
					db_filename = getParamAsString_("new_db");
					
					if ( no_tmp_dbs )
					{
						if ( db_filename.empty() )
						{
							writeLog_("No_tmp_dbs flag set but no name for database given. Aborting!");
							cout << "No_tmp_dbs flag set but no name for database given. Aborting!" << endl;
							return ILLEGAL_PARAMETERS;
						}
					}
					else
					{
						if ( db_filename.empty() )
						{
							// if only the Inspect_in flag is set, a database name has to be be given, except if inspect is run in blind mode
							if ( !Inspect_out && !blind )
							{
								writeLog_("No name for new trie database given. Aborting!");
								cout << "No name for new trie database given. Aborting!" << endl;
								return ILLEGAL_PARAMETERS;
							}
							else
							{
								db_filename = temp_data_dir + "tmp.inspect.db.trie";
								inspect_infile.setDb(db_filename);
								idx_filename = temp_data_dir + "tmp.inspect.db.index";
							}
						}
						else
						{
							if ( db_filename.hasSuffix(".trie") )
							{
								db_filename = db_filename;
								inspect_infile.setDb(db_filename);
								idx_filename = db_filename.substr(0, db_filename.size()-4) + "index";
							}
							else
							{
								db_filename = db_filename + ".trie";
								inspect_infile.setDb(db_filename);
								idx_filename = db_filename + ".index";
							}
						}
					}
				}
				// (2.1.5) blind - running inspect in blind mode after running a normal mode to minimize the database
				if ( getParamAsBool_("blind", false) )
				{
					// a blind search with prior run to minimize the database can only be run in full mode
					if ( Inspect_in && !Inspect_out )
					{
						writeLog_("A blind search with prior run to minimize the database can only be run in full mode. Aborting!");
						cout << "a blind search with prior run to minimize the database can only be run in full mode. Aborting!" << endl;
						printUsage_();
						return ILLEGAL_PARAMETERS;
					}
					blind = true;
				}
				
				if ( blind && blind_only )
				{
					writeLog_("Both blind flags set. Aborting!");
					cout << "Both blind flags set. Aborting! Only one of the two flags [-blind|-blind_only] can be set" << endl;
					return ILLEGAL_PARAMETERS;
				}
				
				snd_db = getParamAsString_("snd_db");
				
				if ( no_tmp_dbs && blind && snd_db.empty() )
				{
					writeLog_("No_tmp_dbs and blind flag set but no name for minimized database given. Aborting!");
					cout << "No_tmp_dbs and blind flag set but no name for minimized database given. Aborting!" << endl;
					return ILLEGAL_PARAMETERS;
				}
				else if ( blind && snd_db.empty() )
				{
					snd_db_filename = temp_data_dir + "tmp.inspect.db.snd.trie";
					snd_index_filename = temp_data_dir + "tmp.inspect.db.snd.index";
				}
				else if ( blind )
				{
					if ( snd_db.hasSuffix(".trie") )
					{
						snd_db_filename = snd_db;
						snd_index_filename = snd_db.substr(0, snd_db.size()-4) + "index";
					}
					else
					{
						snd_db_filename = snd_db + ".trie";
						snd_index_filename = snd_db + ".index";
					}
				}
				
				// get the single modifications
				buffer = getParamAsString_("mods");
				buffer.split(';', substrings);
				
				if ( substrings.empty() && !buffer.empty() ) substrings.push_back(buffer);
				// for each modification get the mass, residues, type (optional) and name (optional)
				for ( vector< String >::iterator i = substrings.begin(); i != substrings.end(); ++i)
				{
					mod.push_back(vector< String >());
					if ( i->hasPrefix("[") ) i->erase(0, 1);
					if ( i->hasSuffix("]") ) i->erase(i->length()-1, 1);
					i->split(',', mod.back());
				}
				if ( !blind_only && mod.empty() )
				{
					writeLog_("No modifications specified. Aborting!");
					cout << "No modifications specified. Aborting!" << endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				inspect_infile.setMod(mod);
				
				inspect_logfile = temp_data_dir + "tmp.inspect.log";

				// (2.1.1.1) optional parameters
				inspect_infile.setProtease(getParamAsString_("protease"));
				inspect_infile.setJumpscores(getParamAsString_("jumpscores"));
				inspect_infile.setInstrument(getParamAsString_("instrument"));
				
				buffer = getParamAsString_("max_mods_pp");
				if ( !buffer.empty() )
				{
					inspect_infile.setMods(getParamAsInt_("max_mods_pp"));
					if ( (inspect_infile.getMods() < 0) )
					{
						writeLog_("Illegal number of modifications (<0) given. Aborting!");
						cout << "Illegal number of modifications (<0) given. Aborting!" << endl;
						printUsage_();
						return ILLEGAL_PARAMETERS;
					}
				}

				buffer = getParamAsString_("PM_tol");
				if ( !buffer.empty() )
				{
					inspect_infile.setPMTolerance( (double) (getParam_("PM_tol")) );
					if ( (inspect_infile.getPMTolerance() < 0) )
					{
						writeLog_("Illegal parent mass tolerance (<0) given. Aborting!");
						cout << "Illegal parent mass tolerance (<0) given. Aborting!" << endl;
						printUsage_();
						return ILLEGAL_PARAMETERS;
					}
				}

				buffer = getParamAsString_("ion_tol");
				if ( !buffer.empty() )
				{
					inspect_infile.setIonTolerance( (double) (getParam_("ion_tol")) );
					if ( (inspect_infile.getIonTolerance() < 0) )
					{
						writeLog_("Illegal ion mass tolerance (<0) given. Aborting!");
						cout << "Illegal ion mass tolerance (<0) given. Aborting!" << endl;
						printUsage_();
						return ILLEGAL_PARAMETERS;
					}
				}

				if ( getParamAsBool_("multicharge", false) ) inspect_infile.setMulticharge(1);

				buffer = getParamAsString_("TagCountA");
				if ( !buffer.empty() )
				{
					inspect_infile.setTagCountA(getParamAsInt_("TagCountA"));
					if ( (inspect_infile.getTagCountA() < 0) )
					{
						writeLog_("Illegal number of tags (TagCountA <0) given. Aborting!");
						cout << "Illegal number of tags (TagCountA <0) given. Aborting!" << endl;
						printUsage_();
						return ILLEGAL_PARAMETERS;
					}
				}

				buffer = getParamAsString_("TagCountB");
				if ( !buffer.empty() ) 
				{
					inspect_infile.setTagCountB(getParamAsInt_("TagCountB"));
					if ( (inspect_infile.getTagCountB() < 0) )
					{
						writeLog_("Illegal number of tags (TagCountB <0) given. Aborting!");
						cout << "Illegal number of tags (TagCountB <0) given. Aborting!" << endl;
						printUsage_();
						return ILLEGAL_PARAMETERS;
					}
				}

				if ( getParamAsBool_("twopass", false) ) inspect_infile.setTwopass(true);
				
				buffer = getParamAsString_("maxptmsize");
				if ( !buffer.empty() )
				{
					inspect_infile.setMaxPTMsize( (double) (getParam_("maxptmsize")) );
					if ( inspect_infile.getMaxPTMsize() < 0 )
					{
						writeLog_("Illegal maximum modification size (<0). Aborting!");
						cout << "Illegal maximum modification size (<0). Aborting!" << endl;
						printUsage_();
						return ILLEGAL_PARAMETERS;
					}
				}
				
				buffer = getParamAsString_("p_value_blind");
				if ( !buffer.empty() ) cutoff_p_value = (double) (getParam_("p_value_blind"));
				if ( (cutoff_p_value < 0) || (cutoff_p_value > 1) )
				{
					writeLog_("Illegal p-value for blind search. Aborting!");
					cout << "Illegal p-value for blind search. Aborting!" << endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				
				buffer = getParamAsString_("score_value_blind");
				if ( !buffer.empty() ) cutoff_score_value = (double) (getParam_("score_value_blind"));
	
				buffer = getParamAsString_("min_spp");
				if ( !buffer.empty() ) min_annotated_spectra_per_protein = getParamAsInt_("min_spp");
			}
			
			// (2.1.2) Inspect_out - output of inspect is written xml analysis file
			if ( Inspect_out )
			{
				// get the database and sequence file name from the input file
				
				buffer = getParamAsString_("p_value");
				if ( !buffer.empty() )
				{
					p_value_threshold = (double) (getParam_("p_value"));
					if ( (p_value_threshold < 0) || (p_value_threshold > 1) )
					{
						writeLog_("Illegal p-value. Aborting!");
						cout << "Illegal p-value. Aborting!" << endl;
						printUsage_();
						return ILLEGAL_PARAMETERS;
					}
				}
				
				buffer = getParamAsString_("score_value");
				if ( !buffer.empty() )
				{
					score_value_threshold = (double) (getParam_("score_value"));
				}
				
				output_filename = getParamAsString_("out");
				if ( output_filename.empty() )
				{
					writeLog_("No output file specified. Aborting!");
					cout << "No output file specified. Aborting!" << endl;
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
			
			// (3.1.1) input file
			// if only Inspect_out is set, the file has to exist, be readable and not empty
			file_info.setFile(input_filename);
			if ( Inspect_out && !Inspect_in )
			{
				if ( !file_info.exists() )
				{
					throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, input_filename);
				}
				if ( !file_info.isReadable() )
				{
					throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, input_filename);
				}
				if ( emptyFile(input_filename) )
				{
					throw Exception::FileEmpty(__FILE__, __LINE__, __PRETTY_FUNCTION__, input_filename);
				}
			}
			// if both flags or only Inspect_in is set, the file has to be writable
			else
			{
				file.setName(input_filename);
				file.open(IO_WriteOnly);
				if ( !file.isWritable() )
				{
					throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, input_filename);
				}
				file.close();
			}
			
			// retrieve the name of the databases from the input file
			if ( Inspect_out && ! Inspect_in )
			{
				String line;
				inspect_infile.setDb("");
				ifstream get_db_names(input_filename.c_str());
				String db = "db,";
				String seq = "sequence_file,";
				while ( getline(get_db_names, line) && inspect_infile.getDb().empty() && inspect_infile.getSequenceFile().empty() )
				{
					if ( !line.empty() && (line[line.length()-1] < 33) ) line.resize(line.length()-1);
					buffer = line;
					buffer.toLower();
					if ( buffer.hasPrefix(db) )
					{
						inspect_infile.setDb(line.substr(db.length(), line.length()-db.length()));
						dbs.push_back(line.substr(db.length(), line.length()-db.length()));
					}
					else if ( buffer.hasPrefix(seq) )
					{
						inspect_infile.setSequenceFile(line.substr(seq.length(), line.length()-seq.length()));
						seq_files.push_back(line.substr(seq.length(), line.length()-seq.length()));
					}
				}
				get_db_names.close();
			}
			
			// (3.1.2.1) inspect output file
			if ( (Inspect_in && Inspect_out) || (Inspect_in && blind) )
			{
				file.setName(inspect_output_filename);
				file.open(IO_WriteOnly);
				if ( !file.isWritable() )
				{
					throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, inspect_output_filename);
				}
				file.close();
			}
			file_info.setFile(inspect_infile.getJumpscores());
			if ( (!inspect_infile.getJumpscores().empty()) && !file_info.isReadable() )
			{
				throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, inspect_infile.getJumpscores());
			}
			
			// (3.1.2) output file
			if ( Inspect_out )
			{
				file.setName(output_filename);
				file.open(IO_WriteOnly);
				if ( Inspect_out )
				{
					if ( !file.isWritable() )
					{
						throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, output_filename);
					}
				}
				file.close();
			}
			
			// (3.1.3) given databases and sequence files
			vector< String > not_accessable;
			for ( vector< String >::const_iterator i = dbs.begin(); i != dbs.end(); ++i )
			{
				file_info.setFile(i->c_str());
				if ( !file_info.exists() ) not_accessable.push_back(*i);
				else if ( !file_info.isReadable() ) not_accessable.push_back(*i);
				else if ( emptyFile(*i) ) not_accessable.push_back(*i);
			}
			
			for ( vector< String >::const_iterator i = seq_files.begin(); i != seq_files.end(); ++i )
			{
				file_info.setFile(i->c_str());
				if ( !file_info.exists() ) not_accessable.push_back(*i);
				else if ( !file_info.isReadable() ) not_accessable.push_back(*i);
				else if ( emptyFile(*i) ) not_accessable.push_back(*i);
			}
			if ( (not_accessable.size() ) == (dbs.size() + seq_files.size()) )
			{
				writeLog_("All of the given databases and sequence files are either not existent, not readable or empty. Aborting!");
				cout << "All of the given databases and sequence files are either not existent, not readable or empty. Aborting!" << endl;
				if ( dbs.empty() ) throw Exception::FileEmpty(__FILE__, __LINE__, __PRETTY_FUNCTION__, seq_files.front());
				else throw Exception::FileEmpty(__FILE__, __LINE__, __PRETTY_FUNCTION__, dbs.front());
			}
			else if ( !not_accessable.empty() )
			{
				buffer = String(not_accessable.size());
				buffer.append(" databases/sequence files are not accessable or empty. Using ");
				buffer.append(String( dbs.size()+seq_files.size()-not_accessable.size() ));
				buffer.append(" databases/sequences files only!");
				writeLog_(buffer.c_str());
				cout << buffer << endl;
			}
			
			if ( Inspect_in )
			{
				// (3.1.3.1) common contaminants
				if ( !no_common_contaminants )
				{
					file_info.setFile(common_contaminants_filename.c_str());
					if ( !file_info.exists() )
					{
						throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, common_contaminants_filename);
					}
					if ( !file_info.isReadable() )
					{
						throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, common_contaminants_filename);
					}
				}
				
				// (3.1.4) database and index
				if ( make_trie_db )
				{
					file.setName(db_filename.c_str());
					file.open(IO_WriteOnly);
					if ( !file.isWritable() )
					{
						throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, db_filename);
					}
					file.close();
					file.setName(idx_filename);
					file.open( IO_WriteOnly );
					if ( !file.isWritable() )
					{
						throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, idx_filename);
					}
					file.close();
				}
				
				// (3.1.5) second database and index
				if ( blind )
				{
					file.setName(snd_db_filename);
					file.open(IO_WriteOnly);
					if ( !file.isWritable() )
					{
						throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, snd_db_filename);
					}
					file.close();
					file.setName(snd_index_filename);
					file.open(IO_WriteOnly);
					if ( !file.isWritable() )
					{
						throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, snd_index_filename);
					}
					file.close();
				}
				
				// the on-screen output of inspect
				file.setName(inspect_logfile);
				file.open(IO_WriteOnly);
				if ( !file.isWritable() )
				{
					throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, inspect_logfile);
				}
				file.close();
			}
			
			// (3.2) running the program
			file_info.setFile(db_filename.c_str());
			String database_path = String(file_info.dirPath().ascii())+"/";
			String database_filename = file_info.fileName().ascii();
			file_info.setFile(idx_filename.c_str());
			String index_filename = file_info.fileName().ascii();
			vector< unsigned int > wanted_records;
			
			// (3.2.1) creating the input file and converting and merging the databases
			if ( Inspect_in )
			{
				if ( !no_common_contaminants )
				{
					seq_files.push_back(common_contaminants_filename);
					tax.push_back("None");
				}
				
				if ( make_trie_db )
				{
					// merging the trie databases (all but the first databases are appended)
					for ( vector< String >::const_iterator i = dbs.begin(); i != dbs.end(); ++i)
					{
						file_info.setFile(i->c_str());
						inspect_infile.compressTrieDB(file_info.fileName().ascii(), "", file_info.dirPath().ascii(), wanted_records, database_filename, index_filename, database_path, i != dbs.begin());
					}
					
					// converting and merging the other databases (all but the first databases are appended)
					vector< String >::const_iterator tax_i = tax.begin();
					for ( vector< String >::const_iterator i = seq_files.begin(); i != seq_files.end(); ++i, ++tax_i)
					{
						file_info.setFile(i->c_str());
						inspect_infile.generateTrieDB(file_info.fileName().ascii(), file_info.dirPath().ascii(), database_path, wanted_records, database_filename, index_filename, ( (i != seq_files.begin()) || (!dbs.empty()) ), *tax_i);
					}
				}
				else
				{
					if ( !dbs.empty() )
					{
						file_info.setFile(dbs[0].c_str());
						database_filename = file_info.fileName().ascii();
						database_path = String(file_info.dirPath().ascii())+"/";
					}
					else
					{
						database_filename = "";
						database_path = "";
					}
					inspect_infile.setDb(String(database_path+database_filename));
					if ( !seq_files.empty() ) inspect_infile.setSequenceFile(seq_files[0]);
				}
				
				if ( blind ) inspect_infile.setBlind(2);
				if ( blind_only ) inspect_infile.setBlind(true);
				
				inspect_infile.store(input_filename);
			}
			
			// (3.2.2) running inspect and generating a second database from the results and running inspect in blind mode on this new database
			if ( blind )
			{
				String call;
				if ( !inspect_dir.empty() )
				{
					call.append("cd ");
					call.append(inspect_dir);
					call.append(" && ./inspect");
					//call.append("inspect -r ");
					//call.append(inspect_dir);
				}
				else
				{
					writeLog_("inspect working directory not given. Aborting!");
					cout << "inspect working directory not given. Aborting!" << endl;
					return ILLEGAL_PARAMETERS;
				}
				call.append(" -i ");
				call.append(input_filename);
				call.append(" -o ");
				call.append(inspect_output_filename);
				// writing the inspect output to a temporary file
				call.append(" > ");
				call.append(inspect_logfile);
				
				int status = system(call.c_str());
				writeLog_("inspect output during running:\n");
				writeLog_(fileContent(inspect_logfile));
				
				if (status != 0)
				{
					cout << "Inspect problem. Aborting! (Details can be seen in the logfile: \"" << logfile << "\")" << endl;
					writeLog_("Inspect problem. Aborting!");
					deleteTempFiles(input_filename, output_filename, inspect_output_filename, db_filename, idx_filename, snd_db_filename, snd_index_filename, inspect_logfile);
					return EXTERNAL_PROGRAM_ERROR;
				}
				
				if ( database_filename.empty() && (!inspect_infile.getSequenceFile().empty()) )
				{
					file_info.setFile(inspect_infile.getSequenceFile().c_str());
					database_path = String(file_info.dirPath().ascii())+"/";
					database_filename = file_info.fileName().ascii();
				}
				
				file_info.setFile(snd_db_filename);
				String snd_db_path = file_info.dirPath().ascii();
				String snd_db_filename_buf = file_info.fileName().ascii();
				file_info.setFile(snd_index_filename);
				String snd_index_filename_buf = file_info.fileName().ascii();
				file_info.setFile(inspect_output_filename);
				
				inspect_infile.generateSecondDatabase(file_info.fileName().ascii(), file_info.dirPath().ascii(), database_path, database_filename, cutoff_p_value, cutoff_score_value, min_annotated_spectra_per_protein, snd_db_filename_buf, snd_index_filename_buf, snd_db_path, index_filename);
				
				if ( emptyFile(snd_db_filename) )
				{
					AnalysisXMLFile analysisXML_file;
					analysisXML_file.store(output_filename, vector< ProteinIdentification >(), vector< Identification >(), vector< float >(), vector< float >(), contact_person);
					Inspect_out = false;
					writeLog_("No proteins matching criteria for generating minimized database for blind search!");
					cout << "No proteins matching criteria for generating minimized database for blind search!" << endl;
				}
				
				// (3.2.3) setting the database name to the new database
				inspect_infile.setDb(snd_db_filename);
				inspect_infile.setSequenceFile("");
				inspect_infile.setBlind(true);
				inspect_infile.store(input_filename);
			}
			
			// (3.2.3) writing the output of inspect into an analysisXML file
			if ( Inspect_in && Inspect_out )
			{
				String call;
				if ( !inspect_dir.empty() )
				{
					call.append("cd ");
					call.append(inspect_dir);
					call.append(" && ./inspect");
					//call.append("inspect -r ");
					//call.append(inspect_dir);
				}
				//else call.append("inspect ");
				else
				{
					writeLog_("inspect working directory not given. Aborting!");
					cout << "inspect working directory not given. Aborting!" << endl;
					return ILLEGAL_PARAMETERS;
				}
				call.append(" -i ");
				call.append(input_filename);
				call.append(" -o ");
				call.append(inspect_output_filename);
				// writing the inspect output to a temporary file
				call.append(" > ");
				call.append(inspect_logfile);
				
				int status = system(call.c_str());
				writeLog_("inspect output during running:\n");
				writeLog_(fileContent(inspect_logfile));
				if (status != 0)
				{
					cout << "Inspect problem. Aborting! (Details can be seen in the logfile: \"" << logfile << "\")" << endl;
					writeLog_("Inspect problem. Aborting!");
					deleteTempFiles(input_filename, output_filename, inspect_output_filename, db_filename, idx_filename, snd_db_filename, snd_index_filename, inspect_logfile);
					return EXTERNAL_PROGRAM_ERROR;
				}
			}

			if ( Inspect_out )
			{
				AnalysisXMLFile analysisXML_file;
				
				if ( !emptyFile(inspect_output_filename) )
				{
					vector< Identification >	identifications;
					ProteinIdentification protein_identification;
					vector< float >	precursor_retention_times, precursor_mz_values;
					
					InspectOutfile inspect_outfile;
					
					file_info.setFile(inspect_infile.getDb().c_str());
					
					inspect_outfile.load(inspect_output_filename, identifications, protein_identification, precursor_retention_times, precursor_mz_values, p_value_threshold, score_value_threshold, file_info.fileName().ascii(), file_info.dirPath().ascii(), inspect_infile.getSequenceFile());
					
					vector<ProteinIdentification> protein_identifications;
					protein_identifications.push_back(protein_identification);
					
					analysisXML_file.store(output_filename, protein_identifications, identifications, precursor_retention_times, precursor_mz_values, contact_person);
				}
				else
				{
					analysisXML_file.store(output_filename, vector< ProteinIdentification >(), vector< Identification >(), vector< float >(), vector< float >(), contact_person);
					writeLog_("No proteins identified!");
					cout << "No proteins identified!" << endl;
				}
			}
			
			// (3.3) deleting all temporary files
			deleteTempFiles(input_filename, output_filename, inspect_output_filename, db_filename, idx_filename, snd_db_filename, snd_index_filename, inspect_logfile);

			return EXECUTION_OK;
		}
};

///@endcond



int main( int argc, char ** argv )
{
	TOPPInspectAdapter tool;

	return tool.main(argc,argv);
}
