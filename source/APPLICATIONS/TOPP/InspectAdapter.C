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

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/Date.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/AnalysisXMLFile.h>
#include <OpenMS/FORMAT/InspectInfile.h>
#include <OpenMS/FORMAT/InspectOutfile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/METADATA/ContactPerson.h>

#include <stdlib.h>
#include <vector>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

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
				(Supported spectrum file formats are .mzXML, .mzData)
				for which identifications are to be found
				and one ore more databases in either trie, FASTA or Swissprot format containing
				the possible proteins.
				The given databases are converted and merged into one trie database.
				This is done because Inspect does the conversion anyway (though with a bug) and may
				actually not use more than two databases (one of them in trie format).
				Additionally you thus can reuse the database without having Inspect done the conversion
				everytime.
				The drawback is, of course, that you need the same amount of space for the trie
				database as well, which can, in case of large and/or many databases, be a problem.
				The results are written as a analysisXML output file. This mode is selected
			 	by default.
			 	</li>
				
				<li>
				Only the first part of the identification process is performed.
				This means that an Inspect input file is generated and the given databases are 
				converted and merged into one trie database. This file can be used
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
				
				This mode is selected by the <b>-inspect_in</b> option in the command line.
				</li>
				
				<li>
				Only the second part of the identification process is performed.
				This means that the output of an Inspect run is analyzed and the result
				written to an analysisXML file.
				
				This mode is selected by the <b>-inspect_out</b> option in the command line.
				</li>
	</ol>
	
	@todo look for possible crash codes of inspect and catching them; extract by-ions, read PTMs from ini file and from input, compute protein score?, catch exceptions to close files (Martin)
	
	@ingroup TOPP
*/

// We do not want this class to show up in the docu -> cond
/// @cond TOPPCLASSES

class TOPPInspectAdapter
	: public TOPPBase
{
	public:
		TOPPInspectAdapter()
			: TOPPBase("InspectAdapter")
		{
		}
	
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
						<< "  -in <file>          the input file OR directory to search (every file in that directory will be searched (non-recursively)" << endl
						<< "                      supported file formats are .mzXML, .mzData" << endl
						<< "                      Note: In mode 'inspect_out' an Inspect result file is read" << endl
						<< "  -out <file>         output file in analysisXML" << endl
						<< "  -inspect_in         if this flag is set the InspectAdapter will write an Inspect input file and generate a trie database" << endl
						<< "  -inspect_out        if this flag is set the InspectAdapter will read an Inspect result  file and write an analysisXML file." << endl
						<< "  -inspect_dir        the Inspect directory." << endl
						<< "  -temp_data_dir      a directory in which some temporary files can be stored" << endl
						<< "  -dbs <file1>,...    names of databases(s) (FASTA and SwissProt supported)" << endl
						<< endl
						<< "  OPTIONAL PARAMETERS" << endl
						<< "  -inspect_output <file>  name for the output file of Inspect (may only be used in a full run)" << endl
						<< "  -instr              the instrument that was used to measure the spectra" << endl
						<< "                      (If set to QTOF, uses a QTOF-derived fragmentation model, and does not attempt to correct the parent mass.)" << endl
						<< "  -prcr_m_tol         the precursor mass tolerance" << endl
						<< "  -pk_m_tol           the peak mass tolerance" << endl
						<< "  -mods <MASS1>,<RESIDUES1>,<TYPE1>,<NAME1>;..." << endl
						<< "                      modifications i.e. [80,STY,opt,phosphorylation]" << endl
						<< "                      MASS and RESIDUES are mandatory" << endl
						<< "                      Valid values for \"TYPE\" are \"fix\", \"cterminal\", \"nterminal\", and \"opt\" (the default)." << endl
						<< "  -multicharge        attempt to guess the precursor charge and mass, and consider multiple charge states if feasible" << endl
						<< "  -protease           the name of a protease. (\"Trypsin\", \"None\", or \"Chymotrypsin\")" << endl
						<< "  -o <file>           direct output file from inspect" << endl
						<< "  -trie_dbs <file1>,... names of database(s) in trie format" << endl
						<< "  -max_mods_pp        number of PTMs permitted in a single peptide." << endl
//						<< "  -twopass            use two-pass search: first pass uses fewer tags, produces list of proteins" << endl
						<< "                      to be re-searched in second pass" << endl
//						<< "  -TagCountA          number of tags for the first pass" << endl
//						<< "  -TagCountB          number of tags for the second pass OR number of tags to use in a one-pass search" << endl
						<< "  -jumpscores <file>  file to specify PTM frequencies, for use in tag generation. This is more accurate tagging than the" << endl
						<< "  -no_tmp_dbs         no temporary databases are used" << endl
						<< "  -new_db             name of the merged trie database" << endl
						<< "                      an index file with extension \".index\" will be created." << endl
						<< "  -p_value            annotations with inferior p-value are ignored" << endl
						<< endl
						<< "  BLIND SEARCH" << endl
						<< "  -blind              perform a blind search (allowing arbitrary modification masses), is preceeded by a normal search to gain a smaller database." << endl
						<< "                      (can only be used in full mode)" << endl
						<< "  -blind_only         like blind but no prior search is performed to reduce the database size" << endl
						<< "  -p_value_blind      used for generating the minimized database" << endl
						<< "  -min_spp            minimum number of spectra a protein has to annotate to be added to the database" << endl
						<< "  -snd_db             name of the minimized trie database generated when using blind mode." << endl
						<< "                      (-1 is #spectra / #proteins * 2)" << endl
						<< "  -maxptmsize         maximum modification size (in Da) to consider" << endl
						<< "                      default behavior (where tags can contain any PTM), but requires the creation of the jump frequency file" << endl;
		}


		void printToolHelpOpt_() const
		{
			cerr	<< endl;
		}


		void setOptionsAndFlags_()
		{
			options_["-inspect_dir"] = "inspect_dir";
			options_["-temp_data_dir"] = "temp_data_dir";
			flags_["-inspect_in"] = "inspect_in";
			flags_["-inspect_out"] = "inspect_out";
			options_["-in"] = "in";
			options_["-trie_dbs"] = "trie_dbs";
			options_["-dbs"] = "dbs";
			options_["-new_db"] = "new_db";
			options_["-snd_db"] = "snd_db";
			options_["-protease"] = "protease";
			options_["-jumpscores"] = "jumpscores";
			options_["-instrument"] = "instrument";
			options_["-mods"] = "mods";
			options_["-max_mods_pp"] = "max_mods_pp";
			options_["-prcr_m_tol"] = "prcr_m_tol";
			options_["-pk_m_tol"] = "pk_m_tol";
			flags_["-multicharge"] = "multicharge";
// 			options_["-TagCountA"] = "TagCountA";
// 			options_["-TagCountB"] = "TagCountB";
//			flags_["-twopass"] = "twopass";
			options_["-out"] = "out";
			options_["-inspect_input"] = "inspect_input";
			options_["-inspect_output"] = "inspect_output";
			flags_["-blind_only"] = "blind_only";
			options_["-p_value"] = "p_value";
			options_["-p_value_blind"] = "p_value_blind";
			options_["-min_spp"] = "min_spp";
			options_["-maxptmsize"] = "maxptmsize";
			flags_["-blind"] = "blind";
			flags_["-cmn_conts"] = "cmn_conts";
			flags_["-no_tmp_dbs"] = "no_tmp_dbs";
		}
		
		// deleting all temporary files
		void deleteTempFiles(const String& input_filename, const String& output_filename, const String& inspect_output_filename, const String& db_filename, const String& idx_filename, const String& snd_db_filename, const String& snd_idx_filename, const String& inspect_logfile)
		{
			if ( input_filename.hasSuffix("tmp.inspect.input") ) remove(input_filename.c_str());
			if ( output_filename.hasSuffix("tmp.inspect.output") ) remove(output_filename.c_str());
			if ( inspect_output_filename.hasSuffix("tmp.direct.inspect.output") ) remove(inspect_output_filename.c_str());
			if ( db_filename.hasSuffix("tmp.inspect.db.trie") ) remove(db_filename.c_str());
			if ( idx_filename.hasSuffix("tmp.inspect.db.index") ) remove(idx_filename.c_str());
			if ( snd_db_filename.hasSuffix("tmp.inspect.db.snd.trie") ) remove(snd_db_filename.c_str());
			if ( snd_idx_filename.hasSuffix("tmp.inspect.db.snd.index") ) remove(snd_idx_filename.c_str());
 			if ( inspect_logfile.hasSuffix("tmp.inspect.log") ) remove(inspect_logfile.c_str());
		}

		ExitCodes main_(int , char**)
		{
			//-------------------------------------------------------------
			// (1) variables
			//-------------------------------------------------------------
			
			InspectInfile inspect_infile;
			InspectOutfile inspect_outfile;
			
			vector< String >
				substrings,
				dbs,
				seq_files;
			
			vector < vector< String > > mod;
			
			String
				string_buffer,
				db_filename,
				idx_filename,
				snd_db_filename,
				snd_idx_filename,
				inspect_logfile,
				logfile,
				inspect_dir,
				temp_data_dir,
				snd_db,
				snd_db_dir,
				output_filename,
				inspect_input_filename,
				inspect_output_filename;
			
			bool
				inspect_in,
				inspect_out,
				blind_only,
				blind,
				no_tmp_dbs;
			
			Real p_value_threshold = 1.0;
			Real cutoff_p_value;
			
			SignedInt
				min_annotated_spectra_per_protein;
			
			ContactPerson contact_person;
			
			
			//-------------------------------------------------------------
			// (2) parsing and checking parameters
			//-------------------------------------------------------------
			
			inspect_in = getParamAsBool_("inspect_in");
			inspect_out = getParamAsBool_("inspect_out");
			
			if ( inspect_in && inspect_out )
			{
				writeLog_("Both Inspect flags set. Aborting!\n"
				          "Only one of the two flags [-inspect_in|-inspect_out] can be set");
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}
			// a 'normal' inspect run corresponds to both inspect_in and inspect_out set
			if ( !inspect_in && !inspect_out ) inspect_in = inspect_out = true;
			
			if ( inspect_out && inspect_in )
			{
				temp_data_dir = getParamAsString_("temp_data_dir");
				if ( temp_data_dir.empty() )
				{
					writeLog_("No directory for temporary files specified. Aborting!");
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}

				File::absolutePath(temp_data_dir);
				temp_data_dir.ensureLastChar('/');
			}
			
			string_buffer = getParamAsString_("in");
			if ( string_buffer.empty() )
			{
				writeLog_("No input file specified. Aborting!");
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}
			else
			{
				File::absolutePath(string_buffer);
				if ( inspect_in )
				{
					inspect_infile.setSpectra(string_buffer);
					if ( inspect_out ) inspect_output_filename = getParamAsString_("inspect_output", temp_data_dir + "tmp.direct.inspect.output");
				}
				else inspect_output_filename = string_buffer;
			}
			
			string_buffer = getParamAsString_("out");
			if ( string_buffer.empty() )
			{
				writeLog_("No output file specified. Aborting!");
				return ILLEGAL_PARAMETERS;
			}
			else
			{
				File::absolutePath(string_buffer);
				if ( inspect_out ) output_filename = string_buffer;
				else inspect_input_filename = string_buffer;
			}
			
			if ( inspect_in && inspect_out )
			{
				inspect_input_filename = getParamAsString_("inspect_input");
				if ( inspect_input_filename.empty() )
				{
					if ( inspect_in && inspect_out ) inspect_input_filename = temp_data_dir + "tmp.inspect.input";
					else if ( inspect_in )
					{
						writeLog_("No name for the inspect input file specified. Aborting!");
						return ILLEGAL_PARAMETERS;
					}
				}
			}
			File::absolutePath(inspect_input_filename);
			
			contact_person.setName(getParamAsString_("contactName", "unknown"));
			writeDebug_(String("Contact name: ") + contact_person.getName(), 1);
			contact_person.setInstitution(getParamAsString_("contactInstitution", "unknown"));
			writeDebug_(String("Contact institution: ") + contact_person.getInstitution(), 1);
			contact_person.setContactInfo(getParamAsString_("contactInfo"));
			writeDebug_(String("Contact info: ") + contact_person.getContactInfo(), 1);
			
			inspect_dir = getParamAsString_("inspect_dir");
			if ( inspect_in && inspect_dir.empty() && inspect_out )
			{
				writeLog_("No inspect directory file specified. Aborting!");
				return ILLEGAL_PARAMETERS;
			}

			File::absolutePath(inspect_dir);
			inspect_dir.ensureLastChar('/');;

			
			blind_only = getParamAsBool_("blind_only");
			
			if ( inspect_in )
			{
				string_buffer = getParamAsString_("trie_dbs");
				if ( !string_buffer.empty() )
				{
					// get the single databases
					string_buffer.split(',', dbs);
					if ( dbs.empty() ) dbs.push_back(string_buffer);
				}
				
				string_buffer = getParamAsString_("dbs");
				if ( !string_buffer.empty() )
				{
					// get the single sequence files
					string_buffer.split(',', seq_files);
					if ( seq_files.empty() ) seq_files.push_back(string_buffer);
				}
				
				// at least one of the parameters db or seq_file has to be set
				if ( dbs.empty() && seq_files.empty() )
				{
					writeLog_("No database specified. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				
				no_tmp_dbs = getParamAsBool_("no_tmp_dbs");
				
				// blind - running inspect in blind mode after running a normal mode to minimize the database
				blind = getParamAsBool_("blind");
				if ( blind && inspect_in && !inspect_out )
				{
					blind = false;
					blind_only = true;
				}
				
				db_filename = getParamAsString_("new_db");
				if ( db_filename.empty() )
				{
					if ( !inspect_out )
					{
						if ( !blind )
						{
							writeLog_("No name for new trie database given. Aborting!");
							return ILLEGAL_PARAMETERS;
						}
					}
					else
					{
						if ( no_tmp_dbs )
						{
							writeLog_("No_tmp_dbs flag set but no name for database given. Aborting!");
							return ILLEGAL_PARAMETERS;
						}
						else
						{
							db_filename = temp_data_dir + "tmp.inspect.db.trie";
							inspect_infile.setDb(db_filename);
							idx_filename = temp_data_dir + "tmp.inspect.db.index";
						}
					}
				}
				else
				{
					File::absolutePath(db_filename);
					if ( db_filename.hasSuffix(".trie") )
					{
						inspect_infile.setDb(db_filename);
						idx_filename = db_filename.substr(0, db_filename.length()-4) + "index";
					}
					else
					{
						idx_filename = db_filename + ".index";
						db_filename = db_filename + ".trie";
						inspect_infile.setDb(db_filename);
					}
				}
				
				if ( blind && blind_only )
				{
					writeLog_("Both blind flags set. Aborting! Only one of the two flags [-blind|-blind_only] can be set");
					return ILLEGAL_PARAMETERS;
				}
				
				snd_db = getParamAsString_("snd_db");
				if ( no_tmp_dbs && blind && snd_db.empty() )
				{
					writeLog_("No_tmp_dbs and blind flag set but no name for minimized database given. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				else if ( blind && snd_db.empty() )
				{
					snd_db_filename = temp_data_dir + "tmp.inspect.db.snd.trie";
					snd_idx_filename = temp_data_dir + "tmp.inspect.db.snd.index";
				}
				else if ( blind )
				{
					File::absolutePath(snd_db_filename);
					if ( snd_db.hasSuffix(".trie") )
					{
						snd_db_filename = snd_db;
						snd_idx_filename = snd_db.substr(0, snd_db.size()-4) + "index";
					}
					else
					{
						snd_db_filename = snd_db + ".trie";
						snd_idx_filename = snd_db + ".index";
					}
				}
				
				// get the single modifications
				if ( !blind_only )
				{
					string_buffer = getParamAsString_("mods");
					string_buffer.split(';', substrings);
					
					if ( substrings.empty() && !string_buffer.empty() ) substrings.push_back(string_buffer);
					// for each modification get the mass, residues, type (optional) and name (optional)
					for ( vector< String >::iterator i = substrings.begin(); i != substrings.end(); ++i)
					{
						mod.push_back(vector< String >());
						i->split(',', mod.back());
						if ( mod.back().size() < 2 || mod.back().size() > 4 )
						{
							writeLog_("Illegal number of parameters for modification given. Aborting!");
							return ILLEGAL_PARAMETERS;
						}
						else
						{
							try
							{
								mod.back().front().toFloat();
							}
							catch ( Exception::ConversionError ce )
							{
								writeLog_("Given mass is no float. Aborting!");
								return ILLEGAL_PARAMETERS;
							}
						}
					}
					inspect_infile.setMod(mod);
				}
				
				inspect_infile.setProtease(getParamAsString_("protease"));
				inspect_infile.setJumpscores(getParamAsString_("jumpscores"));
				inspect_infile.setInstrument(getParamAsString_("instrument"));
				
				inspect_infile.setMods(getParamAsInt_("max_mods_pp", -1));
				if ( inspect_infile.getMods() < 1 && !mod.empty() )
				{
					writeLog_("Modifications specified, but max_mods_pp not set. Setting it to 1.");
					inspect_infile.setMods(1);
				}
				
				string_buffer = getParamAsString_("prcr_m_tol");
				if ( !string_buffer.empty() )
				{
					inspect_infile.setPMTolerance(getParamAsDouble_("prcr_m_tol"));
					if ( (inspect_infile.getPMTolerance() < 0) )
					{
						writeLog_("Illegal parent mass tolerance (<0) given. Aborting!");
						return ILLEGAL_PARAMETERS;
					}
				}

				string_buffer = getParamAsString_("pk_m_tol");
				if ( !string_buffer.empty() )
				{
					inspect_infile.setIonTolerance( getParamAsDouble_("pk_m_tol") );
					if ( (inspect_infile.getIonTolerance() < 0) )
					{
						writeLog_("Illegal ion mass tolerance (<0) given. Aborting!");
						return ILLEGAL_PARAMETERS;
					}
				}

				if ( getParamAsBool_("multicharge") ) inspect_infile.setMulticharge(1);

// 				string_buffer = getParamAsString_("TagCountA");
// 				if ( !string_buffer.empty() )
// 				{
// 					inspect_infile.setTagCountA(getParamAsInt_("TagCountA"));
// 					if ( (inspect_infile.getTagCountA() < 0) )
// 					{
// 						writeLog_("Illegal number of tags (TagCountA <0) given. Aborting!");
// 						return ILLEGAL_PARAMETERS;
// 					}
// 				}
// 
// 				string_buffer = getParamAsString_("TagCountB");
// 				if ( !string_buffer.empty() )
// 				{
// 					inspect_infile.setTagCountB(getParamAsInt_("TagCountB"));
// 					if ( (inspect_infile.getTagCountB() < 0) )
// 					{
// 						writeLog_("Illegal number of tags (TagCountB <0) given. Aborting!");
// 						return ILLEGAL_PARAMETERS;
// 					}
// 				}

// 				if ( getParamAsBool_("twopass") ) inspect_infile.setTwopass(true);
				
				string_buffer = getParamAsString_("maxptmsize");
				if ( !string_buffer.empty() )
				{
					inspect_infile.setMaxPTMsize(getParamAsDouble_("maxptmsize") );
					if ( inspect_infile.getMaxPTMsize() < 0 )
					{
						writeLog_("Illegal maximum modification size (<0). Aborting!");
						return ILLEGAL_PARAMETERS;
					}
				}
				
				string_buffer = getParamAsString_("min_spp");
				if ( !string_buffer.empty() ) min_annotated_spectra_per_protein = getParamAsInt_("min_spp");
			}
			
			if ( inspect_out )
			{
				p_value_threshold = getParamAsDouble_("p_value", 1.0);
				if ( (p_value_threshold < 0) || (p_value_threshold > 1) )
				{
					writeLog_("Illegal p-value. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				
				inspect_logfile = temp_data_dir + "tmp.inspect.log";
			}
			
			if ( blind && inspect_in )
			{
				cutoff_p_value = getParamAsDouble_("p_value_blind", p_value_threshold);
				if ( (cutoff_p_value < 0) || (cutoff_p_value > 1) )
				{
					writeLog_("Illegal p-value for blind search. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
			}
			
			//-------------------------------------------------------------
			// (3) running program according to parameters
			//-------------------------------------------------------------
			// checking accessability of files
			
			// the file for the inspect output
			if ( (inspect_in && inspect_out) || (inspect_in && blind) )
			{
				if ( !File::writable(inspect_output_filename) )
				{
					throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, inspect_output_filename);
				}
			}
			
			if ( !inspect_infile.getJumpscores().empty() )
			{
				if ( !File::readable(inspect_infile.getJumpscores()) )
				{
					throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, inspect_infile.getJumpscores());
				}
			}
			
			// output file
			if ( inspect_out )
			{
				if ( !File::writable(output_filename) )
				{
					throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, output_filename);
				}
			}
			
			vector< String > not_accessable, accessable_db, idx, accessable_seq;
			if ( inspect_in )
			{
				if ( !File::writable(inspect_input_filename) )
				{
					throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, inspect_input_filename);
				}
				
				// database and index
				if ( !File::writable(db_filename) )
				{
					throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, db_filename);
				}

				if ( !File::writable(idx_filename) )
				{
					throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, idx_filename);
				}
				
				// given databases and sequence files
				for ( vector< String >::const_iterator db_i = dbs.begin(); db_i != dbs.end(); ++db_i )
				{
					if ( !File::readable(*db_i) || File::empty(*db_i) ) 
					{
						not_accessable.push_back(*db_i);
					}
					else // if the file is accessable, try to find the corresponding index file and check it
					{
						if ( db_i->hasSuffix(".trie") ) 
						{
							string_buffer = db_i->substr(0, db_i->length()-4) + "index";
						}
						else 
						{
							string_buffer = *db_i + "index";
						}
						
						if ( !File::readable(string_buffer) || File::empty(string_buffer) ) 
						{
							not_accessable.push_back(*db_i);
						}
						else
						{
							accessable_db.push_back(*db_i);
							idx.push_back(string_buffer);
						}
					}
				}
				
				for ( vector< String >::const_iterator db_i = seq_files.begin(); db_i != seq_files.end(); ++db_i )
				{
					if ( !File::readable(*db_i) || File::empty(*db_i) )
					{ 
						not_accessable.push_back(*db_i);
					}
					else 
					{
						accessable_seq.push_back(*db_i);
					}
				}
				if ( (not_accessable.size() ) == (dbs.size() + seq_files.size()) )
				{
					writeLog_("All of the given databases are either not existent, not readable or empty. Aborting!");
					throw Exception::FileEmpty(__FILE__, __LINE__, __PRETTY_FUNCTION__, not_accessable.front());
				}
				else if ( !not_accessable.empty() )
				{
					string_buffer = String(not_accessable.size());
					string_buffer.append(" databases are not accessable or empty. Using ");
					string_buffer.append(String( accessable_db.size() + accessable_seq.size() ));
					string_buffer.append(" databases only!");
					writeLog_(string_buffer.c_str());
				}
				
				// second database and index
				if ( blind )
				{
					if ( !File::writable(snd_db_filename) )
					{
						throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, snd_db_filename);
					}
					if ( !File::writable(snd_idx_filename) )
					{
						throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, snd_idx_filename);
					}
				}
				
				// the on-screen output of inspect
				if ( inspect_out )
				{
					if ( !File::writable(inspect_logfile) )
					{
						writeLog_(String(" Could not write in temp data directory: ")
								+ temp_data_dir + inspect_logfile + " Aborting!");
						return ILLEGAL_PARAMETERS;
					}
				}
			}
			
			vector< unsigned int > wanted_records;
			
			// creating the input file and converting and merging the databases
			if ( inspect_in )
			{
				// merging the trie databases (all but the first databases are appended)
				vector< String >::const_iterator idx_i = idx.begin();
				for ( vector< String >::const_iterator db_i = accessable_db.begin(); db_i != accessable_db.end(); ++db_i, ++idx_i )
				{
					inspect_outfile.compressTrieDB(*db_i, *idx_i, vector< UnsignedInt >(), db_filename,  idx_filename, (db_i != accessable_db.begin()) );
				}
				
				// converting and merging the other databases (all but the first database are appended)
				for ( vector< String >::const_iterator db_i = accessable_seq.begin(); db_i != accessable_seq.end(); ++db_i )
				{
					inspect_outfile.generateTrieDB(*db_i, db_filename,  idx_filename, ( (db_i != accessable_seq.begin()) || (!accessable_db.empty()) ));
				}
				
				if ( blind_only )
				{
					inspect_infile.setMulticharge(false);
					inspect_infile.setBlind(true);
				}
				
				inspect_infile.store(inspect_input_filename);
			}
			
			// running inspect and generating a second database from the results and running inspect in blind mode on this new database
			if ( blind && inspect_in && inspect_out )
			{
				String call;
				call.append("cd ");
				call.append(inspect_dir);
				call.append(" && ./inspect");
				//call.append("inspect -r ");
				//call.append(inspect_dir);
				call.append(" -i ");
				call.append(inspect_input_filename);
				call.append(" -o ");
				call.append(inspect_output_filename);
				// writing the inspect output to a temporary file
				call.append(" > ");
				call.append(inspect_logfile);

				int status = system(call.c_str());
				
				//debug output
				writeLog_("inspect output while running:\n");
				TextFile inspect_logfile_content(inspect_logfile);
				writeLog_(inspect_logfile_content.asString());
				
				if (status != 0)
				{
					writeLog_("Inspect problem. Aborting! (Details can be seen in the logfile: \"" + logfile + "\")");
					deleteTempFiles(inspect_input_filename, output_filename, inspect_output_filename, db_filename, idx_filename, snd_db_filename, snd_idx_filename, inspect_logfile);
					return EXTERNAL_PROGRAM_ERROR;
				}
				
				vector< UnsignedInt > wanted_records = inspect_outfile.getWantedRecords(inspect_output_filename, p_value_threshold);
				
				if ( wanted_records.empty() )
				{
					AnalysisXMLFile analysisXML_file;
					analysisXML_file.store(output_filename, vector< ProteinIdentification >(), vector< Identification >(), vector< Real >(), vector< Real >(), contact_person);
					inspect_out = false;
					writeLog_("No proteins matching criteria for generating minimized database for blind search!");
					
					deleteTempFiles(inspect_input_filename, output_filename, inspect_output_filename, db_filename, idx_filename, snd_db_filename, snd_idx_filename, inspect_logfile);
				}
				inspect_outfile.compressTrieDB(db_filename, idx_filename, wanted_records, snd_db_filename, snd_idx_filename, false);
				
				// setting the database name to the new database
				inspect_infile.setDb(snd_db_filename);
				inspect_infile.setSequenceFile("");
				inspect_infile.setBlind(true);
				inspect_infile.getMod().clear();
				inspect_infile.store(inspect_input_filename);
			}
			
			// writing the output of inspect into an analysisXML file
			if ( inspect_in && inspect_out )
			{
				String call;
				call.append("cd ");
				call.append(inspect_dir);
				call.append(" && ./inspect");
				//call.append("inspect -r ");
				//call.append(inspect_dir);
				
				call.append(" -i ");
				call.append(inspect_input_filename);
				call.append(" -o ");
				call.append(inspect_output_filename);
				// writing the inspect output to a temporary file
				call.append(" > ");
				call.append(inspect_logfile);
				
				int status = system(call.c_str());
				writeLog_("inspect output while running:\n");
				TextFile inspect_logfile_content(inspect_logfile);
				writeLog_(inspect_logfile_content.asString());
				if (status != 0)
				{
					writeLog_("Inspect problem. Aborting! (Details can be seen in the logfile: \"" + logfile + "\")");
					deleteTempFiles(inspect_input_filename, output_filename, inspect_output_filename, db_filename, idx_filename, snd_db_filename, snd_idx_filename, inspect_logfile);
					return EXTERNAL_PROGRAM_ERROR;
				}
			}
			
			if ( inspect_out )
			{
				AnalysisXMLFile analysisXML_file;
				
				if ( !File::empty(inspect_output_filename) )
				{
					vector< Identification > identifications;
					ProteinIdentification protein_identification;
					vector< Real > precursor_retention_times, precursor_mz_values;
					
					try
					{
						vector< UnsignedInt > corrupted_lines = inspect_outfile.load(inspect_output_filename, identifications, protein_identification, precursor_retention_times, precursor_mz_values, p_value_threshold);
//				const std::string& database_filename)
					}
					catch( Exception::ParseError pe )
					{
						deleteTempFiles(inspect_input_filename, output_filename, inspect_output_filename, db_filename, idx_filename, snd_db_filename, snd_idx_filename, inspect_logfile);
						writeLog_(String(pe.getMessage()) + " Aborting!");
						return INPUT_FILE_CORRUPT;
					}
					vector< ProteinIdentification > protein_identifications;
					protein_identifications.push_back(protein_identification);
					
					analysisXML_file.store(output_filename, protein_identifications, identifications, precursor_retention_times, precursor_mz_values, contact_person);
				}
				else
				{
					analysisXML_file.store(output_filename, vector< ProteinIdentification >(), vector< Identification >(), vector< Real >(), vector< Real >(), contact_person);
					writeLog_("No proteins identified!");
				}
			}
			
			// (3.3) deleting all temporary files
			deleteTempFiles(inspect_input_filename, output_filename, inspect_output_filename, db_filename, idx_filename, snd_db_filename, snd_idx_filename, inspect_logfile);

			return EXECUTION_OK;
		}
};

///@endcond



int main( int argc, char ** argv )
{
	TOPPInspectAdapter tool;

	return tool.main(argc,argv);
}
