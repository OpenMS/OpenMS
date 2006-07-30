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
#include "TOPPBase.h"

#include <stdlib.h>
#include <vector>

using namespace OpenMS;

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
	
	@todo look for possible crash codes of inspect and catching them; extract by-ions, read PTMs from ini file and from input, compute protein score?
	
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
		{
			
		}
	
	protected:
		void printToolUsage_()
		{
			std::cerr	<< std::endl
						<< tool_name_ << " -- annotates MS/MS spectra using Inspect" << std::endl
						<< std::endl
						<< "Usage:" << std::endl
						<< " " << tool_name_ << " [options]" << std::endl
						<< std::endl
						<< "Options are:" << std::endl
						<< "  -in <file>          Inspect input file" << std::endl
						<< "  -out <file>         output file in analysisXML" << std::endl
						<< "  -o <file>           direct output file from inspect" << std::endl
						<< "  -Inspect_in         if this flag is set the InspectAdapter will create an Inspect Input file" << std::endl
						<< "                      if only Inspect_in is set, a name for the trie database (see below) has to be given!" << std::endl
						<< "  -Inspect_out        if this flag is set the InspectAdapter will read in an Inspect Input file and write an analysisXML file." << std::endl
						<< "  -inspect_dir        the name of the InsPecT directory." << std::endl
						<< "  -temp_data_dir      the name of the directory where the temporary data will be stored." << std::endl
						<< "  -spectra <file>     the spectrum file OR directory to search (every file in that directory will be searched(non-recursively)" << std::endl
						<< "                      supported spectrum file formats are .mzXML, .mzData, .ms2, dta, and .pkl" << std::endl
						<< "                      multiple spectra in one .dta file are not supported" << std::endl
						<< "  -trie_dbs <file1>,<file2>,...      names of a databases (.trie file) to search ()" << std::endl
						<< "  -dbs <file>;tax1,<file2>;tax2,...  names of a other databases to search (currently FASTA and SwissProt are supported)" << std::endl
						<< "                                     tax - the desired taxonomy, if not given for a database, all entries are taken." << std::endl
						<< "  -make_trie_db       if set, the InspectAdapter will generate one trie database from all given databases." << std::endl
						<< "                      if you do not use this switch you may only use one FASTA database XOR one trie database" << std::endl
						<< "  -instr              the instrument that was used to measure the spectra (default read from INI file)" << std::endl
						<< "                      (If set to QTOF, uses a QTOF-derived fragmentation model, and does not attempt to correct the parent mass.)" << std::endl
						<< "  -PM_tol             the precursor mass tolerance (default read from INI file)" << std::endl
						<< "  -ion_tol            the peak mass tolerance (default read from INI file)" << std::endl
						<< "  -protease           the name of a protease. \"Trypsin\", \"None\", and \"Chymotrypsin\" are the available values." << std::endl
						<< "  -mods [<MASS1>,<RESIDUES1>,<TYPE1>,<NAME1>];[<MASS2>,<RESIDUES2>,<TYPE2>,<NAME2>]" << std::endl
						<< "                      modifications i.e. [80,STY,opt,phosphorylation] (default read from INI file)" << std::endl
						<< "                      MASS and RESIDUES are mandatory, make sure the modifications are seperated by a semicolon!" << std::endl
						<< "                      Valid values for \"type\" are \"fix\", \"cterminal\", \"nterminal\", and \"opt\" (the default)." << std::endl
						<< "                      The first four	characters of the name should be unique." << std::endl
						<< "  -max_mods_pp        number of PTMs permitted in a single peptide. (default: read from INI file)" << std::endl
						<< "  -blind              perform a blind search (allowing arbitrary modification masses), as this is slower than the normal search" << std::endl
						<< "                      A normal search is performed in advance to gain a smaller database." << std::endl
						<< "                      This search can only be run in full mode." << std::endl
						<< "  -blind_only         like blind but no prior search is performed to reduce the database size" << std::endl
						<< "  -p_value            annotations with inferior p-value are ignored. Default is 0.05" << std::endl
						<< "  -p_value_blind      used when generating the minimized database for blind search" << std::endl
						<< "  -min_spp            used when generating the minimized database for blind search " << std::endl
						<< "                      the minimum number of spectra a protein has to annotate in order to add it to the filtered database " << std::endl
						<< "                      default is #spectra / #proteins * 2" << std::endl
						<< "  -maxptmsize         for blind search, specifies the maximum modification size (in Da) to consider (default read from INI file)" << std::endl
						<< "  -jumpscores <file>  file to specify PTM frequencies, for use in tag generation. This is more accurate tagging than the" << std::endl
						<< "                      default behavior (where tags can contain any PTM), but requires the creation of the jump frequency file" << std::endl
						<< "  -multicharge        attempt to guess the precursor charge and mass, and consider multiple charge states if feasible" << std::endl
						<< "  -twopass            use two-pass search. The first pass uses fewer tags, and produces a list of proteins" << std::endl
						<< "                      which are re-searched in the second pass" << std::endl
						<< "  -TagCountA          number of tags to generate for the first pass of a two-pass search" << std::endl
						<< "  -TagCountB          number of tags to generate for the second pass of a two-pass search" << std::endl
						<< "                      OR the number of tags to use in a one-pass search" << std::endl
						<< "  -no_cmn_conts       do not add the proteins from CommonContaminents.fasta to the search database" << std::endl
						<< "  -no_tmp_dbs         no temporary databases are used" << std::endl
						<< "  -new_db             name of the trie database (given databases are converted and merged to one trie database)." << std::endl
						<< "                      This has to be set if no_tmp_dbs is set! If the name does not end with \".trie\"" << std::endl
						<< "                      it is extended accordingly." << std::endl
						<< "                      An index file with the same name but extension \".index\" will be created." << std::endl
						<< "  -snd_db             name of the minimized trie database generated when using blind mode." << std::endl
						<< "                      This has to be set if no_tmp_dbs is set!" << std::endl;
						//<< "  -contact		 name of the contact person" << std::endl
		}


		void printToolHelpOpt_()
		{
		}


		void setOptionsAndFlags_()
		{
			options_["-inspect_dir"] = "inspect_dir";
			options_["-temp_data_dir"] = "temp_data_dir";
			flags_["-Inspect_in"] = "Inspect_in";
			options_["-spectra"] = "spectra";
			options_["-trie_dbs"] = "dbs";
			options_["-dbs"] = "seq_files";
			options_["-new_db"] = "new_db";
			options_["-snd_db"] = "snd_db";
			options_["-tax"] = "tax";
			options_["-protease"] = "protease";
			options_["-jumpscores"] = "jumpscores";
			options_["-instrument"] = "instrument";
			options_["-mods"] = "mod";
			options_["-max_mods_pp"] = "mods";
			options_["-PM_tol"] = "PM_tolerance";
			options_["-ion_tol"] = "ion_tolerance";
			flags_["-multicharge"] = "multicharge";
			options_["-TagCountA"] = "TagCountA";
			options_["-TagCountB"] = "TagCountB";
			flags_["-twopass"] = "twopass";
			flags_["-Inspect_out"] = "Inspect_out";
			options_["-in"] = "in";
			options_["-out"] = "out";
			options_["-o"] = "o";
			flags_["-blind_only"] = "blind_only";
			options_["-p_value"] = "p_value_threshold";
			options_["-p_value_blind"] = "cutoff_p_value";
			options_["-min_spp"] = "min_annotated_spectra_per_protein";
			options_["-maxptmsize"] = "maxptmsize";
			flags_["-blind"] = "blind";
			flags_["-no_cmn_conts"] = "no_common_contaminants";
			flags_["-no_tmp_dbs"] = "no_tmp_dbs";
			flags_["-make_trie_db"] = "make_trie_db";
			//options_["-contact"] = "contact_person";
		}
		
		// workaround because qt string is not convertable to std::string
		bool exists(const std::string& filename)
		{
			FILE* file = fopen(filename.c_str(), "r");
			bool ret = ( file != NULL );
			if ( ret ) fclose(file);
			return ret;
		}

		bool isReadable(const std::string& filename)
		{
			bool ret;
			FILE* file = fopen(filename.c_str(), "r");
			if ( file != NULL )
			{
				fgetc(file);
				ret = !ferror(file);
				fclose(file);
				return ret;
			}
			else return false;
		}

		bool isWritable(const std::string& filename)
		{
			bool writable;
			bool existed = exists(filename);
			FILE* file = fopen(filename.c_str(), "a");
			writable = ( file != NULL );
			if ( writable ) fclose(file);
			if ( !existed ) remove(filename.c_str());
			return writable;
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
		
		bool emptyFile(const std::string& filename)
		{
			return ( fsize(filename) == 0 );
		}
		
		std::string pathDir(const std::string& filename, char slash = '/')
		{
			return filename.substr(0, filename.find_last_of(slash)+1);
		}
		
		std::string fileName(const std::string& filename, char slash = '/')
		{
			return filename.substr(filename.find_last_of(slash)+1);
		}

		std::string fileContent(const std::string& filename)
		{
			long size = fsize(filename);
			if ( size != -1 )
			{
				FILE* file = fopen(filename.c_str(), "r");
				char* buffer = new char[size+1];
				buffer[size] = 0;
				fread (buffer, size, 1, file);
				fclose(file);
				std::string sbuffer = buffer;
				delete(buffer);
				return sbuffer;
			}
			else return std::string();
		}
		
		// deleting all temporary files
		void deleteTempFiles(const String& input_filename, const String& output_filename, const String& inspect_output_filename, const String& db_filename, const String& idx_filename, const String& snd_db_filename, const String& snd_index_filename, const String& inspect_logfilename)
		{
			if ( input_filename.hasSuffix("tmp.inspect.input") ) remove(input_filename.c_str());
			if ( output_filename.hasSuffix("tmp.inspect.output") ) remove(output_filename.c_str());
			if ( inspect_output_filename.hasSuffix("tmp.direct.inspect.output") ) remove(inspect_output_filename.c_str());
			if ( db_filename.hasSuffix("tmp.inspect.db.trie") ) remove(db_filename.c_str());
			if ( idx_filename.hasSuffix("tmp.inspect.db.index") ) remove(idx_filename.c_str());
			if ( snd_db_filename.hasSuffix("tmp.inspect.db.snd.trie") ) remove(snd_db_filename.c_str());
			if ( snd_index_filename.hasSuffix("tmp.inspect.db.snd.index") ) remove(snd_index_filename.c_str());
			if ( inspect_logfilename.hasSuffix("tmp.inspect.log") ) remove(inspect_logfilename.c_str());
		}

		ExitCodes main_(int , char**)
		{
			///-------------------------------------------------------------
			// (1) variables
			///-------------------------------------------------------------

			InspectInfile inspect_infile;

			// (1.0) general variables
			std::vector< String > substrings;
			String buffer, db_filename, idx_filename, snd_db_filename, snd_index_filename, common_contaminants_filename, logfile, inspect_logfilename;
			ContactPerson contact_person;

			// (1.1) parameter variables
			// (1,1,0) general parameter variables
			String inspect_dir, temp_data_dir;

			// (1.1.1) Inspect_in - writing the inspect input file only and corresponding parameters
			bool Inspect_in = false;
			// (1.1.1.0) mandatory parameters
			String new_db, new_db_dir, snd_db, snd_db_dir; // at least one of the parameters db or seq_file has to be set
			std::vector< String >dbs, seq_files, tax; // if several dbs are given, they are merged into one, that is then processed

			// (1.1.1.1) optional parameters
			bool make_trie_db = false;
			
			std::vector < std::vector< String > > mod; // some from ini file

			// (1.1.2) Inspect_out - executing the program only and writing xml analysis file and corresponding parameters
			double p_value_threshold = 1.0;
			bool Inspect_out = false;
			String output_filename, inspect_output_filename;

			// (1.1.3) parameters corresponding to both Inspect_in and Inspect_out
			String input_filename; // in normal mode, this can be a temporary file (from ini)

			// (1.1.4) blind_only - running inspect in blind mode only and corresponding parameters
			bool blind_only = false;

			// (1.1.5) blind - running inspect in blind mode after running a normal mode to minimize the database
			bool blind = false;
			double cutoff_p_value = 0.05;
			int min_annotated_spectra_per_protein = -1;

			// (1.1.6) no_common_contaminants - whether to include the proteins in commonContaminants.fasta
			bool no_common_contaminants = false;

			// (1.1.7) no_tmp_dbs - whether to use temporary database files or to save them (faster if they are used more than once)
			bool no_tmp_dbs = false;


			///-------------------------------------------------------------
			// (2) parsing and checking parameters
			///-------------------------------------------------------------
			// (2.0) general variables
			contact_person.setName(getParamAsString_("contactName", "unknown"));
			contact_person.setInstitution(getParamAsString_("contactInstitution", "unknown"));
			contact_person.setContactInfo(getParamAsString_("contactInfo"));
			logfile = getParamAsString_("log", "Inspect.log");
			
			// (2.1) parameter variables
			// (2.1,0) general parameter variables
			inspect_dir = getParamAsString_("inspect_dir");
			
			inspect_infile.ensurePathChar(inspect_dir);
			common_contaminants_filename = inspect_dir + "CommonContaminants.fasta";
			temp_data_dir = getParamAsString_("temp_data_dir");
			if ( temp_data_dir.empty() ) temp_data_dir = "/home/bude/langwisc/inspect/temp/";
			inspect_infile.ensurePathChar(temp_data_dir);
			
			output_filename = temp_data_dir; output_filename.append("tmp.inspect.output");
			inspect_output_filename = temp_data_dir; inspect_output_filename.append("tmp.direct.inspect.output");
			input_filename = temp_data_dir; input_filename.append("tmp.inspect.input");
			
			// (2.1.1) Inspect_in - writing the inspect input file only and corresponding parameters
			if ( getParamAsString_("Inspect_in", "false") != "false" ) Inspect_in = true;
			if ( getParamAsString_("Inspect_out", "false") != "false" ) Inspect_out = true;

			// a 'normal' inspect run corresponds to both Inspect_in and Inspect_out set
			if ( !Inspect_in && !Inspect_out ) Inspect_in = Inspect_out = true;

			if ( Inspect_in )
			{
				// (2.1.1.0) mandatory parameters
				inspect_infile.setSpectra(getParamAsString_("spectra"));
				if ( inspect_infile.getSpectra().empty() )
				{
					writeLog_("No spectrum file specified. Aborting!");
					std::cout << "No spectrum file specified. Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}

				buffer = getParamAsString_("dbs");
				if ( !buffer.empty() )
				{
					// get the single databases
					buffer.split(',', dbs);
					if ( dbs.empty() ) dbs.push_back(buffer);
				}
				
				buffer = getParamAsString_("seq_files");

				if ( !buffer.empty() )
				{
					// get the single sequence files
					buffer.split(',', seq_files);
					if ( seq_files.empty() ) seq_files.push_back(buffer);
					
					// get the corresponding taxonomies
					for ( std::vector< String >::iterator i = seq_files.begin(); i != seq_files.end(); ++i)
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
					std::cout << "No database or sequence file specified. Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}

				new_db = getParamAsString_("new_db");
				snd_db = getParamAsString_("snd_db");

				// (2.1.1.1) optional parameters
				inspect_infile.setSpectra(getParamAsString_("spectra"));
				inspect_infile.setProtease(getParamAsString_("protease"));
				inspect_infile.setJumpscores(getParamAsString_("jumpscores"));
				inspect_infile.setInstrument(getParamAsString_("instrument"));

				// get the single modifications
				getParamAsString_("mod").split(';', substrings);
				
				// for each modification get the mass, residues, type (optional) and name (optional)
				for ( std::vector< String >::iterator i = substrings.begin(); i != substrings.end(); ++i)
				{
					mod.push_back(std::vector< String >());
					if ( i->hasPrefix("[") ) i->erase(0, 1);
					if ( i->hasSuffix("]") ) i->erase(i->length()-1, 1);
					i->split(',', mod.back());
				}
				inspect_infile.setMod(mod);
				
				buffer = getParamAsString_("mods");
				if ( !buffer.empty() )
				{
					inspect_infile.setMods(getParamAsInt_("mods"));
					if ( (inspect_infile.getMods() < 0) )
					{
						writeLog_("Illegal number of modifications (<0) given. Aborting!");
						std::cout << "Illegal number of modifications (<0) given. Aborting!" << std::endl;
						printUsage_();
						return ILLEGAL_PARAMETERS;
					}
				}

				buffer = getParamAsString_("PM_tolerance");
				if ( !buffer.empty() )
				{
					inspect_infile.setPMTolerance( (double) (getParam_("PM_tolerance")) );
					if ( (inspect_infile.getPMTolerance() < 0) )
					{
						writeLog_("Illegal parent mass tolerance (<0) given. Aborting!");
						std::cout << "Illegal parent mass tolerance (<0) given. Aborting!" << std::endl;
						printUsage_();
						return ILLEGAL_PARAMETERS;
					}
				}

				buffer = getParamAsString_("ion_tolerance");
				if ( !buffer.empty() )
				{
					inspect_infile.setIonTolerance( (double) (getParam_("ion_tolerance")) );
					if ( (inspect_infile.getIonTolerance() < 0) )
					{
						writeLog_("Illegal ion mass tolerance (<0) given. Aborting!");
						std::cout << "Illegal ion mass tolerance (<0) given. Aborting!" << std::endl;
						printUsage_();
						return ILLEGAL_PARAMETERS;
					}
				}

				if ( getParamAsString_("mutlicharge", "false") != "false" ) inspect_infile.setMulticharge(1);

				buffer = getParamAsString_("TagCountA");
				if ( !buffer.empty() )
				{
					inspect_infile.setTagCountA(getParamAsInt_("TagCountA"));
					if ( (inspect_infile.getTagCountA() < 0) )
					{
						writeLog_("Illegal number of tags (TagCountA <0) given. Aborting!");
						std::cout << "Illegal number of tags (TagCountA <0) given. Aborting!" << std::endl;
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
						std::cout << "Illegal number of tags (TagCountB <0) given. Aborting!" << std::endl;
						printUsage_();
						return ILLEGAL_PARAMETERS;
					}
				}

				if ( getParamAsString_("twopass", "false") != "false" ) inspect_infile.setTwopass(true);
			}

			// (2.1.2) Inspect_out - executing the program only and writing xml analysis file and corresponding parameters
			buffer = getParamAsString_("p_value_threshold");
			if ( !buffer.empty() )
			{
				p_value_threshold = (double) (getParam_("p_value_threshold"));
				if ( (p_value_threshold < 0) || (p_value_threshold > 1) )
				{
					writeLog_("Illegal p-value. Aborting!");
					std::cout << "Illegal p-value. Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
			}
			
			buffer = getParamAsString_("out");
			if ( buffer.empty() )
			{
				if ( Inspect_out || !Inspect_in )
				{
					writeLog_("No output file specified. Aborting!");
					std::cout << "No output file specified. Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
			}
			else
			{
				output_filename = buffer;
			}
			
			buffer = getParamAsString_("o");
			if ( !buffer.empty() )
			{
				inspect_output_filename = buffer;
			}

			// (2.1.3) parameters corresponding to both Inspect_in and Inspect_out
			buffer = getParamAsString_("in"); // in normal mode, this can be a temporary file (from ini)
			if ( buffer.empty() )
			{
				if ( (Inspect_in || Inspect_out) && !(Inspect_in && Inspect_out) )
				{
					writeLog_("No input file specified. Aborting!");
					std::cout << "No input file specified. Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
			}
			else
			{
				input_filename = buffer;
			}

			// (2.1.4) blind_only - running inspect in blind mode only and corresponding parameters
			if ( getParamAsString_("blind_only", "false") != "false" ) blind_only = true;

			buffer = getParamAsString_("maxptmsize");
			if ( !buffer.empty() )
			{
				inspect_infile.setMaxPTMsize( (double) (getParam_("maxptmsize")) );
				if ( inspect_infile.getMaxPTMsize() < 0 )
				{
					writeLog_("Illegal maximum modification size (<0). Aborting!");
					std::cout << "Illegal maximum modification size (<0). Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
			}

			// (2.1.5) blind - running inspect in blind mode after running a normal mode to minimize the database
			if ( getParamAsString_("blind", "false") != "false" )
			{
				// a blind search with prior run to minimize the database can only be run in full mode
				if ( Inspect_in && !Inspect_out )
				{
					writeLog_("A blind search with prior run to minimize the database can only be run in full mode. Aborting!");
					std::cout << "a blind search with prior run to minimize the database can only be run in full mode. Aborting!" << std::endl;
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				blind = true;
				inspect_infile.setBlind(true);
			}

			if ( blind && blind_only )
			{
				writeLog_("Both blind flags set. Aborting!");
				std::cout << "Both blind flags set. Aborting! Only one of the two flags [-blind|-blind_only] can be set" << std::endl;
				return ILLEGAL_PARAMETERS;
			}

			buffer = getParamAsString_("cutoff_p_value");
			if ( !buffer.empty() ) cutoff_p_value = double(getParam_("cutoff_p_value"));

			buffer = getParamAsString_("min_annotated_spectra_per_protein");
			if ( !buffer.empty() ) min_annotated_spectra_per_protein = getParamAsInt_("min_annotated_spectra_per_protein");

			// (2.1.6) no_common_contaminants - whether to include the proteins in commonContaminants.fasta
			if ( getParamAsString_("no_common_contaminants", "false") != "false" ) no_common_contaminants = true;
			
			if ( getParamAsString_("make_trie_db", "false") != "false" ) make_trie_db = true;
			if ( make_trie_db && ((!dbs.empty()) + (!seq_files.empty()) + (!no_common_contaminants) >1) )
			{
				writeLog_("Too many databases (make_trie_db not set). Aborting!");
				std::cout << "Too many databases (make_trie_db not set). Aborting!" << std::endl;
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}

			// (2.1.7) no_tmp_dbs - whether to use temporary database files or to save them (faster if they are used more than once)
			if ( getParamAsString_("no_tmp_dbs", "false") != "false" ) no_tmp_dbs = true;
			if ( no_tmp_dbs )
			{
				if ( new_db.empty() )
				{
					writeLog_("Mo_tmp_dbs flag set but no name for database given. Aborting!");
					std::cout << "No_tmp_dbs flag set but no name for database given. Aborting!" << std::endl;
					return ILLEGAL_PARAMETERS;
				}
				if ( snd_db.empty() && blind )
				{
					writeLog_("No_tmp_dbs and blind flag set but no name for minimized database given. Aborting!");
					std::cout << "No_tmp_dbs and blind flag set but no name for minimized database given. Aborting!" << std::endl;
					return ILLEGAL_PARAMETERS;
				}
			}
			else
			{
				if ( new_db.empty() )
				{
					// if only the Inspect_in flag is set, a database name has to be be given
					if ( Inspect_in && !Inspect_out )
					{
						writeLog_("No name for new trie database given. Aborting!");
						std::cout << "No name for new trie database given. Aborting!" << std::endl;
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
					if ( new_db.hasSuffix(".trie") )
					{
						db_filename = new_db;
						inspect_infile.setDb(db_filename);
						idx_filename = new_db.substr(0, new_db.size()-4) + "index";
					}
					else
					{
						db_filename = new_db + ".trie";
						inspect_infile.setDb(db_filename);
						idx_filename = new_db + ".index";
					}
				}
				
				if ( blind && snd_db.empty() )
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
			}
			
			inspect_logfilename = temp_data_dir + "tmp.inspect.log";
			
			///-------------------------------------------------------------
			// (3) running program according to parameters
			///-------------------------------------------------------------
			
			// (3.1) checking accessability of files
			// (3.1.1) input file
			// the input file has to be existent and readable if Inspect_out is set only
			if ( !Inspect_in )
			{
				if ( !exists(input_filename) )
				{
					throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, input_filename);
				}
				if ( !isReadable(input_filename) )
				{
					throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, input_filename);
				}
				if ( fsize(input_filename) == 0 )
				{
					throw Exception::FileEmpty(__FILE__, __LINE__, __PRETTY_FUNCTION__, input_filename);
				}
			}
			
			// (3.1.2) output file
			if ( Inspect_out )
			{
				if ( !isWritable(output_filename) )
				{
					throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, output_filename);
				}
			
				// (3.1.2.1) inspect output file
				if ( !isWritable(inspect_output_filename) )
				{
					throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, inspect_output_filename);
				}
				
				if ( (!inspect_infile.getJumpscores().empty()) && !isReadable(inspect_infile.getJumpscores()) )
				{
					throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, inspect_infile.getJumpscores());
				}
			}

			if ( Inspect_in )
			{
				// (3.1.3) given databases and sequence files
				std::vector< String > not_accessable;
				for ( std::vector< String >::const_iterator i = dbs.begin(); i != dbs.end(); ++i )
				{
					if ( !exists(*i) ) not_accessable.push_back(*i);
					else if ( !isReadable(*i) ) not_accessable.push_back(*i);
					else if ( fsize(*i) == 0 ) not_accessable.push_back(*i);
				}

				for ( std::vector< String >::const_iterator i = seq_files.begin(); i != seq_files.end(); ++i )
				{
					if ( !exists(*i) ) not_accessable.push_back(*i);
					else if ( !isReadable(*i) ) not_accessable.push_back(*i);
					else if ( fsize(*i) == 0 ) not_accessable.push_back(*i);
				}
				if ( (not_accessable.size() ) == (dbs.size() + seq_files.size()) )
				{
					writeLog_("All of the given databases and sequence files are either not existent, not readable or empty. Aborting!");
					std::cout << "All of the given databases and sequence files are either not existent, not readable or empty. Aborting!" << std::endl;
					if ( dbs.empty() ) throw Exception::FileEmpty(__FILE__, __LINE__, __PRETTY_FUNCTION__, seq_files.front());
					else throw Exception::FileEmpty(__FILE__, __LINE__, __PRETTY_FUNCTION__, dbs.front());
				}
				else if ( !not_accessable.empty() )
				{
					buffer = String(SignedInt(not_accessable.size()));
					buffer.append(" databases/sequence files are not accessable or empty. Using ");
					buffer.append(String( SignedInt(dbs.size()+seq_files.size()-not_accessable.size() )));
					buffer.append(" databases/sequences files only!");
					writeLog_(buffer.c_str());
					std::cout << buffer << std::endl;
				}
				
				// (3.1.3.1) common contaminants
				if ( !no_common_contaminants )
				{
					if ( !exists(common_contaminants_filename) )
					{
						throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, common_contaminants_filename);
					}
					if ( !isReadable(common_contaminants_filename) )
					{
						throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, common_contaminants_filename);
					}
				}

				// (3.1.4) database and index
				if ( make_trie_db || (!dbs.empty()) )
				{
					if ( !isWritable(db_filename) )
					{
						throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, db_filename);
					}
					
					if ( !isWritable(idx_filename) )
					{
						throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, idx_filename);
					}
				}

				// (3.1.5) second database and index
				if ( blind )
				{
					if ( !isWritable(snd_db_filename) )
					{
						throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, snd_db_filename);
					}
					
					if ( !isWritable(snd_index_filename) )
					{
						throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, snd_index_filename);
					}
				}
			}
			
			// the on-screen output of inspect
			if ( !isWritable(inspect_logfilename) )
			{
				throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, inspect_logfilename);
			}
			
			// (3.2) running the program
			
			String database_path = pathDir(db_filename);
			String database_filename = fileName(db_filename);
			String index_filename = fileName(idx_filename);
			std::vector< unsigned int > wanted_records;
			
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
					for ( std::vector< String >::const_iterator i = dbs.begin(); i != dbs.end(); ++i)
					{
						inspect_infile.compressTrieDB(fileName(*i), "", pathDir(*i), wanted_records, database_filename, index_filename, database_path, i != dbs.begin());
					}
					
					// converting and merging the other databases (all but the first databases are appended)
					std::vector< String >::const_iterator tax_i = tax.begin();
					for ( std::vector< String >::const_iterator i = seq_files.begin(); i != seq_files.end(); ++i, ++tax_i)
					{
						inspect_infile.generateTrieDB(fileName(*i), pathDir(*i), database_path, wanted_records, database_filename, index_filename, ( (i != seq_files.begin()) || (!dbs.empty()) ), *tax_i);
					}
				}
				else
				{
					if ( !dbs.empty() )
					{
						database_filename = fileName(dbs[0]);
						database_path = pathDir(dbs[0]);
						inspect_infile.ensurePathChar(database_path);
					}
					else
					{
						database_filename = "";
						database_path = "";
					}
					inspect_infile.setDb(String(database_path+database_filename));
					if ( !seq_files.empty() ) inspect_infile.setSequenceFile(seq_files[0]);
				}
				
				if ( blind ) inspect_infile.setBlind(false);
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
					std::cout << "inspect working directory not given. Aborting!" << std::endl;
					return ILLEGAL_PARAMETERS;
				}
				call.append(" -i ");
				call.append(input_filename);
				call.append(" -o ");
				call.append(inspect_output_filename);
				// writing the inspect output to a temporary file
				call.append(" > ");
				call.append(inspect_logfilename);
				
				int status = system(call.c_str());
				writeLog_("inspect output during running:\n");
				writeLog_(fileContent(inspect_logfilename));
				
				if (status != 0)
				{
					std::cout << "Inspect problem. Aborting! (Details can be seen " 
					<< " in the logfile: \"" << logfile << "\")" << std::endl;
					writeLog_("Inspect problem. Aborting!");
					deleteTempFiles(input_filename, output_filename, inspect_output_filename, db_filename, idx_filename, snd_db_filename, snd_index_filename, inspect_logfilename);
					return EXTERNAL_PROGRAM_ERROR;						
				}
				
				if ( database_filename.empty() && (!inspect_infile.getSequenceFile().empty()) )
				{
					database_filename = inspect_infile.getSequenceFile();
					database_path = pathDir(database_filename);
					database_filename = fileName(database_filename);
				}
				
				inspect_infile.generateSecondDatabase(fileName(inspect_output_filename), pathDir(inspect_output_filename), database_path, database_filename, cutoff_p_value, min_annotated_spectra_per_protein, fileName(snd_db_filename), fileName(snd_index_filename), pathDir(snd_db_filename), index_filename);
				
				
				if ( emptyFile(snd_db_filename) )
				{
					AnalysisXMLFile analysisXML_file;
					analysisXML_file.store(output_filename, std::vector< ProteinIdentification >(), std::vector< Identification >(), std::vector< float >(), std::vector< float >(), contact_person);
					Inspect_out = false;
					writeLog_("No proteins matching criteria for generating minimized database for blind search!");
					std::cout << "No proteins matching criteria for generating minimized database for blind search!" << std::endl;
				}
				
				// (3.2.3) setting the database name to the new database
				inspect_infile.setDb(snd_db_filename);
				inspect_infile.setSequenceFile("");
				inspect_infile.setBlind(true);
				inspect_infile.store(input_filename);
			}
			
			// (3.2.3) writing the output of inspect into an analysisXML file
			if ( Inspect_out )
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
					std::cout << "inspect working directory not given. Aborting!" << std::endl;
					return ILLEGAL_PARAMETERS;
				}
				call.append(" -i ");
				call.append(input_filename);
				call.append(" -o ");
				call.append(inspect_output_filename);
				// writing the inspect output to a temporary file
				call.append(" > ");
				call.append(inspect_logfilename);
				
				int status = system(call.c_str());
				writeLog_("inspect output during running:\n");
				writeLog_(fileContent(inspect_logfilename));
				if (status != 0)
				{
					std::cout << "Inspect problem. Aborting! (Details can be seen " 
					<< " in the logfile: \"" << logfile << "\")" << std::endl;
					writeLog_("Inspect problem. Aborting!");
					deleteTempFiles(input_filename, output_filename, inspect_output_filename, db_filename, idx_filename, snd_db_filename, snd_index_filename, inspect_logfilename);
					return EXTERNAL_PROGRAM_ERROR;						
				}

				// if Inspect_in is not set, retrieve the name of the database from the input file
				if ( !Inspect_in )
				{
					String line;
					inspect_infile.setDb("");
					std::ifstream get_db_name(input_filename.c_str());
					while ( getline(get_db_name, line) && inspect_infile.getDb().empty() )
					{
						buffer = line.substr(0,2);
						buffer.toUpper();
						if ( buffer == "DB" )
						{
							if ( !line.empty() ) line.resize(line.length()-1);
							inspect_infile.setDb(line.substr(3, line.length()-3));
						}
					}
					get_db_name.close();
					get_db_name.clear();
				}
				
				AnalysisXMLFile analysisXML_file;
				
				if ( !emptyFile(inspect_output_filename) )
				{
					std::vector< Identification >	identifications;
					ProteinIdentification protein_identification;
					std::vector< float >	precursor_retention_times, precursor_mz_values;
					
					InspectOutfile inspect_outfile;

					inspect_outfile.load(inspect_output_filename, identifications, protein_identification, precursor_retention_times, precursor_mz_values, p_value_threshold, fileName(inspect_infile.getDb()), pathDir(inspect_infile.getDb()), inspect_infile.getSequenceFile());

					std::vector<ProteinIdentification> protein_identifications;
					protein_identifications.push_back(protein_identification);
					
					analysisXML_file.store(output_filename, protein_identifications, identifications, precursor_retention_times, precursor_mz_values, contact_person);
				}
				else
				{
					analysisXML_file.store(output_filename, std::vector< ProteinIdentification >(), std::vector< Identification >(), std::vector< float >(), std::vector< float >(), contact_person);
					writeLog_("No proteins identified!");
					std::cout << "No proteins identified!" << std::endl;
				}
			}
			
			// (3.3) deleting all temporary files
			deleteTempFiles(input_filename, output_filename, inspect_output_filename, db_filename, idx_filename, snd_db_filename, snd_index_filename, inspect_logfilename);

			return OK;
		}
};

///@endcond



int main( int argc, char ** argv )
{
	TOPPInspectAdapter tool;

	return tool.main(argc,argv);
}
