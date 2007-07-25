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
// $MaIntainer: Martin Langwisch $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/InspectInfile.h>
#include <OpenMS/FORMAT/InspectOutfile.h>
#include <OpenMS/FORMAT/PTMXMLFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/TextFile.h>

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

	This wrapper application serves for getting peptide peptide_identifications
	for MS/MS spectra. The wrapper can be executed in three different
	modes:
	<ol>
				<li>
				The whole process of ProteinIdentification via Inspect is executed.
				Inputfile is an mz file containing the MS/MS spectra
				(Supported spectrum file formats are .mzXML, .mzData)
				for which the identifications are to be found and one ore more
				databases in either trie, FASTA or Swissprot format containing
				the possible proteins.
				The given databases are converted and merged into one trie database.
				This is done because Inspect does the conversion anyway
				(though with a bug) and may actually not use more than two
				databases (one of them in trie format). Additionally you thus can
				reuse the database without having Inspect done the conversion
				everytime.
				The drawback is, of course, that you need the same amount of space
				for the trie database as well, which can, in case of large and/or many
				databases, be a problem.
				The results are written as a IdXML output file. This mode is selected
			 	by default.
			 	</li>

				<li>
				Only the first part of the ProteinIdentification process is performed.
				This means that an Inspect input file is generated and the given databases are
				converted and merged Into one trie database. This file can be used
				directly with Inspect whereas the created database and the spectrum file(s)
				have to remain at the given positions.
				Calling an Inspect process should look like the following:

				@code ./inspect -i  inputfilename -o outputfilename  @endcode

				Consult your Inspect reference manual for further details.

				This mode is selected by the <b>-inspect_in</b> option in the command line.
				</li>

				<li>
				Only the second part of the ProteinIdentification process is performed.
				This means that the output of an Inspect run is analyzed and the result
				written to an IdXML file.

				This mode is selected by the <b>-inspect_out</b> option in the command line.
				</li>
	</ol>

	@todo extract by-ions, compute protein score? (Martin)
*/

// We do not want this class to show up in the docu -> cond
/// @cond TOPPCLASSES

class TOPPInspectAdapter
	: public TOPPBase
{
	public:
		TOPPInspectAdapter()
			: TOPPBase("InspectAdapter", "annotates MS/MS spectra using Inspect.")
		{
		}

	protected:
		void registerOptionsAndFlags_()
		{
			registerStringOption_("out", "<file>", "", "output file in IdXML format.\n"
			                                           "Note: In mode 'inspect_in' an Inspect input file is written.");
			registerStringOption_("in", "<file>", "", "input file in mzXML or mzData format.\n"
					 																			"Note: In mode 'inspect_out' an Inspect results file is read");
			registerFlag_("inspect_in", "if this flag is set the InspectAdapter will read in mzXML,\n"
																							 "write an Inspect input file and generate a trie database");
			registerFlag_("inspect_out", "if this flag is set the InspectAdapter will read in a Inspect results file\n"
																								 "and write IdXML");
			registerStringOption_("inspect_directory", "<dir>", "", "the directory in which Inspect is located", false);
			registerStringOption_("temp_data_directory", "<dir>", "", "a directory in which some temporary files can be stored", false);
			registerStringOption_("dbs", "<file>", "", "name(s) of database(s) to search in (FASTA and SwissProt supported)", false);
			registerStringOption_("trie_dbs", "<file>", "", "name(s) of databases(s) to search in (trie-format)", false);
			registerStringOption_("new_db", "<file>", "", "name of the merged trie database", false);
			registerStringOption_("instrument", "<i>", "", "the instrument that was used to measure the spectra\n"
																										 "(If set to QTOF, uses a QTOF-derived fragmentation model,\n"
																										 "and does not attempt to correct the parent mass.)", false);
			registerDoubleOption_("precursor_mass_tolerance", "<tol>", 2.0 , "the precursor mass tolerance", false);
			registerDoubleOption_("peak_mass_tolerance", "<tol>", 1.0, "the peak mass tolerance", false);
			registerFlag_("list_modifications", "show a list of the available modifications");
			registerStringOption_("modifications", "<mods>", "", "the colon-seperated modifications; may be\n"
																																														"<name>,<type>, e.g.: Deamidation,opt or\n"
																																														"<composition>,<residues>,<type>,<name>, e.g.: H2C2O,KCS,opt,Acetyl or\n"
																																														"<mass>,<residues>,<type>,<name>, e.g.: 42.0367,KCS,opt,Acetyl or\n"
																																														"Valid values for type are \"fix\" and \"opt\" (default)\n"
																																														"If you want terminal PTMs, write \"cterm\" or \"nterm\" instead of residues", false);
			registerFlag_("use_monoisotopic_mod_mass", "use monoisotopic masses for the modifications");
			registerStringOption_("modifications_xml_file", "<file>", "", "name of an XML file with the modifications", false);
			registerStringOption_("cleavage", "<enz>", "Trypsin", "the enzyme used for digestion", false);
			registerStringOption_("inspect_output", "<file>", "", "name for the output file of Inspect (may only be used in a full run)", false);
			registerStringOption_("inspect_input", "<file>", "", "name for the input file of Inspect (may only be used in a full run)", false);
			registerFlag_("multicharge", "attempt to guess the precursor charge and mass,\n"
																								 "and consider multiple charge states if feasible");
			registerIntOption_("max_modifications_pp", "<num>", -1 ,"number of PTMs permitted in a single peptide.", false);
			registerIntOption_("tag_count", "<num>", -1, "number of tags to generate", false);
			registerFlag_("no_tmp_dbs", "no temporary databases are used");
			registerDoubleOption_("p_value", "<prob>", 1.0, "annotations with inferior p-value are ignored", false);
			addEmptyLine_();
			addText_("Options for blind search");
			registerFlag_("blind", "perform a blind search (allowing arbitrary modification masses),\n"
																			 "is preceeded by a normal search to gain a smaller database.\n"
														 "(in full mode only)");
			registerFlag_("blind_only", "like blind but no prior search is performed to reduce the database size");
			registerDoubleOption_("p_value_blind", "<prob>", 1.0, "used for generating the minimized database", false);
			registerIntOption_("min_spp", "<num>", -1, "minimum number of spectra a protein has to annotate\n"
																																					 "to be added to the database", false);
			registerStringOption_("snd_db", "<file>", "", "name of the minimized trie database generated when using blind mode.", false);
			registerDoubleOption_("maxPTMsize", "<num>", 250.0, "maximum modification size (in Da) to consider", false);
			registerStringOption_("contact_name", "<name>", "unknown", "Name of the contact", false);
			registerStringOption_("contact_institution", "<name>", "unknown", "Name of the contact institution", false);
			registerStringOption_("contact_info", "<info>", "unknown", "Some information about the contact", false);
		}

		ExitCodes main_(Int , char**)
		{
			//-------------------------------------------------------------
			// (1) variables
			//-------------------------------------------------------------

			InspectInfile inspect_infile;
			InspectOutfile inspect_outfile;

			vector< String >
				substrings,
				substrings2,
				dbs,
				seq_files,
				not_accessable,
				accessable_db,
				idx, accessable_seq;

			String
				string_buffer,
				db_filename,
				idx_filename,
				snd_db_filename,
				snd_idx_filename,
				inspect_logfile,
				logfile,
				inspect_directory,
				temp_data_directory,
				snd_db,
				snd_db_directory,
				output_filename,
				inspect_input_filename,
				inspect_output_filename,
				modifications_filename,
				isotope_filename;

			bool
				inspect_in(false),
				inspect_out(false),
				blind_only(false),
				blind(false),
				no_tmp_dbs(false),
				monoisotopic(false);

			Real p_value_threshold(1.0);
			Real cutoff_p_value(1.0);

			char separator = '/';

			ContactPerson contact_person;

			// filename and tag: file has to: 1 - exist  2 - be readable  4 - writable  8 - be deleted afterwards
			vector< pair< String, UInt > > files;
			UInt const
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
				for ( map< String, pair< String, String > >::const_iterator mod_i = PTM_informations.begin(); mod_i != PTM_informations.end(); ++mod_i )
				{
					max_name_length = max(max_name_length, mod_i->first.length());
					max_composition_length = max(max_composition_length, mod_i->second.first.length());
					max_amino_acids_length = max(max_amino_acids_length, mod_i->second.second.length());
				}
				PTM_info << "name" << String(max_name_length - 4, ' ') << "\t" << "composition" << String(max_composition_length - 11, ' ') << "\t" << "amino_acids" << String(max_amino_acids_length - 11, ' ') << endl;
				for ( map< String, pair< String, String > >::const_iterator mod_i = PTM_informations.begin(); mod_i != PTM_informations.end(); ++mod_i )
				{
					PTM_info << mod_i->first << String(max_name_length - mod_i->first.length(), ' ') << "\t" << mod_i->second.first << String(max_composition_length - mod_i->second.first.length(), ' ') << "\t" << mod_i->second.second << String(max_amino_acids_length - mod_i->second.second.length(), ' ') << endl;
				}
				std::cout << PTM_info.str() << std::endl;

				return EXECUTION_OK;
			}

			inspect_in = getFlag_("inspect_in");
			inspect_out = getFlag_("inspect_out");

			if ( inspect_in && inspect_out )
			{
				writeLog_("Both Inspect flags set. Only one of the two flags [-inspect_in|-inspect_out] can be set. Aborting!");
				return ILLEGAL_PARAMETERS;
			}

			if ( inspect_in ) writeDebug_("Inspect flag: mascot_in (reads in MzXML/MzData, writes Inspect generic format)", 1);
			else if ( inspect_out ) writeDebug_("Inspect flag: mascot_in (reads in Inspect result file, writes IdXML file)", 1);
			else writeDebug_("No Inspect flag set: reads in MzXML/MzData, writes IdXML file", 1);

			// a 'normal' inspect run corresponds to both inspect_in and inspect_out set
			if ( !inspect_in && !inspect_out ) inspect_in = inspect_out = true;

			if ( inspect_out && inspect_in )
			{
				temp_data_directory = getStringOption_("temp_data_directory");
				if ( temp_data_directory.empty() )
				{
					writeLog_("No directory for temporary files specified. Aborting!");
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
				File::absolutePath(temp_data_directory);
				temp_data_directory.ensureLastChar(separator);
			}

			string_buffer = getStringOption_("in");
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
					MSExperiment<> experiment;
					String type;
					try
					{
						inspect_outfile.getExperiment(experiment, type, string_buffer); // may throw an exception if the filetype could not be determined
					}
					catch(Exception::ParseError pe )
					{
						writeLog_(pe.getMessage());
						return PARSE_ERROR;
					}
					if ( type != "mzXML" )
					{
						string_buffer.append(".mzXML");
						MzXMLFile().store(string_buffer, experiment);
						files.push_back(make_pair(string_buffer, writable | delete_afterwards));
					}
					inspect_infile.setSpectra(string_buffer);

					if ( inspect_out )
					{
						inspect_output_filename = getStringOption_("inspect_output");
						if ( inspect_output_filename.empty() )
						{
							inspect_output_filename = temp_data_directory + "tmp.direct.inspect.output";
							files.push_back(make_pair(inspect_output_filename, writable | delete_afterwards));
						}
						else
						{
							File::absolutePath(inspect_output_filename);
							files.push_back(make_pair(inspect_output_filename, writable));
						}
					}
				}
				else
				{
					inspect_output_filename = string_buffer;
					File::absolutePath(inspect_output_filename);
					files.push_back(make_pair(inspect_output_filename, readable));
				}
			}

			string_buffer = getStringOption_("out");
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
				files.push_back(make_pair(string_buffer, writable));
			}

			if ( inspect_in && inspect_out )
			{
				inspect_input_filename = getStringOption_("inspect_input");
				if ( inspect_input_filename.empty() )
				{
					inspect_input_filename = temp_data_directory + "tmp.inspect.input";
					files.push_back(make_pair(inspect_input_filename, writable | delete_afterwards));
				}
				else
				{
					File::absolutePath(inspect_input_filename);
					files.push_back(make_pair(inspect_input_filename, writable));
				}
			}

			inspect_directory = getStringOption_("inspect_directory");
			if ( inspect_in && inspect_directory.empty() && inspect_out )
			{
				writeLog_("No inspect directory file specified. Aborting!");
				return ILLEGAL_PARAMETERS;
			}
			File::absolutePath(inspect_directory);
			inspect_directory.ensureLastChar(separator);

			blind_only = getFlag_("blind_only");

			contact_person.setName(getStringOption_("contact_name"));
			contact_person.setInstitution(getStringOption_("contact_institution"));
			contact_person.setContactInfo(getStringOption_("contact_info"));

			if ( inspect_in )
			{
				string_buffer = getStringOption_("trie_dbs");
				if ( !string_buffer.empty() )
				{
					// get the single databases
					string_buffer.split(',', dbs);
					if ( dbs.empty() ) dbs.push_back(string_buffer);
				}

				string_buffer = getStringOption_("dbs");
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

				no_tmp_dbs = getFlag_("no_tmp_dbs");

				// blind - running inspect in blind mode after running a normal mode to minimize the database
				blind = getFlag_("blind");
				if ( blind && inspect_in && !inspect_out )
				{
					blind = false;
					blind_only = true;
				}

				db_filename = getStringOption_("new_db");
				if ( db_filename.empty() && (!seq_files.empty() || dbs.size() != 1) )
				{
					if ( !inspect_out )
					{
								writeLog_("No name for new trie database given. Aborting!");
								return ILLEGAL_PARAMETERS;
					}
					else
					{
						if ( no_tmp_dbs )
						{
							writeLog_("no_tmp_dbs flag set but no name for database given. Aborting!");
							return ILLEGAL_PARAMETERS;
						}
						else
						{
							db_filename = temp_data_directory + "tmp.inspect.db.trie";
							files.push_back(make_pair(db_filename, writable | delete_afterwards));
							inspect_infile.setDb(db_filename);
							idx_filename = temp_data_directory + "tmp.inspect.db.index";
							files.push_back(make_pair(idx_filename, writable | delete_afterwards));
						}
					}
				}
				else
				{
					// if only one trie database is given, this one is used
					if ( db_filename.empty() ) db_filename = dbs.front();

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
					writeLog_("Both blind flags set. Only one of the two flags [-blind|-blind_only] can be set. Aborting!");
					return ILLEGAL_PARAMETERS;
				}

				snd_db = getStringOption_("snd_db");
				if ( no_tmp_dbs && blind && snd_db.empty() )
				{
					writeLog_("No_tmp_dbs and blind flag set but no name for minimized database given. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				else if ( blind && snd_db.empty() )
				{
					snd_db_filename = temp_data_directory + "tmp.inspect.db.snd.trie";
					snd_idx_filename = temp_data_directory + "tmp.inspect.db.snd.index";
					files.push_back(make_pair(snd_db_filename, writable | delete_afterwards));
					files.push_back(make_pair(snd_idx_filename, writable | delete_afterwards));
				}
				else if ( blind )
				{
					File::absolutePath(snd_db);
					if ( snd_db.hasSuffix(".trie") )
					{
						snd_db_filename = snd_db;
						snd_idx_filename = snd_db.substr(0, snd_db.size()-4) + "index";
						files.push_back(make_pair(snd_db_filename, writable));
						files.push_back(make_pair(snd_idx_filename, writable));
					}
					else
					{
						snd_db_filename = snd_db + ".trie";
						snd_idx_filename = snd_db + ".index";
						files.push_back(make_pair(snd_db_filename, writable));
						files.push_back(make_pair(snd_idx_filename, writable));
					}
				}

				// get the known modifications
				monoisotopic = getFlag_("use_monoisotopic_mod_mass");
				if ( !blind_only )
				{
					// modifications
					string_buffer = getStringOption_("modifications");
					try
					{
						inspect_infile.handlePTMs(string_buffer, modifications_filename, monoisotopic);
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
				}

				inspect_infile.setEnzyme(getStringOption_("cleavage"));
				inspect_infile.setInstrument(getStringOption_("instrument"));

				inspect_infile.setModificationsPerPeptide(getIntOption_("max_modifications_pp"));
				if ( inspect_infile.getModificationsPerPeptide() < 1 && !inspect_infile.getModifications().empty() )
				{
					writeLog_("Modifications specified, but max_modifications_pp not set. Setting it to 1.");
					inspect_infile.setModificationsPerPeptide(1);
				}

				inspect_infile.setPrecursorMassTolerance(getDoubleOption_("precursor_mass_tolerance"));
				if ( (inspect_infile.getPrecursorMassTolerance() < 0 && inspect_infile.getPrecursorMassTolerance() != -1) )
				{
					writeLog_("Illegal precursor mass tolerance (<0) given. Aborting!");
					return ILLEGAL_PARAMETERS;
				}

				inspect_infile.setPeakMassTolerance( getDoubleOption_("peak_mass_tolerance") );
				if ( (inspect_infile.getPeakMassTolerance() < 0 && inspect_infile.getPeakMassTolerance() != -1) )
				{
					writeLog_("Illegal peak mass tolerance (<0) given. Aborting!");
					return ILLEGAL_PARAMETERS;
				}

				if ( getFlag_("multicharge") ) inspect_infile.setMulticharge(1);

				inspect_infile.setTagCount(getIntOption_("tag_count"));
				if ( (inspect_infile.getTagCount() < 0 && inspect_infile.getTagCount() != -1) )
				{
					writeLog_("Illegal number of tags (tag_count <0) given. Aborting!");
					return ILLEGAL_PARAMETERS;
				}

				inspect_infile.setMaxPTMsize(getDoubleOption_("maxPTMsize") );
				if ( (inspect_infile.getMaxPTMsize() < 10 || inspect_infile.getMaxPTMsize() > 2000) && inspect_infile.getMaxPTMsize() != -1)
				{
					writeLog_("Illegal maximum modification size (not in [10,2000]). Aborting!");
					return ILLEGAL_PARAMETERS;
				}
			}

			if ( inspect_out )
			{
				p_value_threshold = getDoubleOption_("p_value");
				if ( (p_value_threshold < 0) || (p_value_threshold > 1) )
				{
					writeLog_("Illegal p-value. Aborting!");
					return ILLEGAL_PARAMETERS;
				}

				inspect_logfile = temp_data_directory + "tmp.inspect.log";
				files.push_back(make_pair(inspect_logfile, writable | delete_afterwards));
			}

			if ( blind && inspect_in )
			{
				cutoff_p_value = getDoubleOption_("p_value_blind");
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

			bool existed(false);
			UInt file_tag(0);

			for ( vector< pair< String, UInt > >::const_iterator files_i = files.begin(); files_i != files.end(); ++files_i )
			{
				string_buffer = files_i->first;
				file_tag = files_i->second;

				if ( (file_tag & exist || file_tag & readable) && !File::exists(string_buffer) )
				{
					throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, string_buffer);
				}

				if ( (file_tag & readable) && !File::readable(string_buffer) )
				{
					throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, string_buffer);
				}

				existed = File::exists(string_buffer);
				if ( (file_tag & writable) && !File::writable(string_buffer) )
				{
					throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, string_buffer);
				}
				else if ( !existed ) remove(string_buffer.c_str());
				existed = false;
			}

			if ( inspect_in )
			{
				// given databases and sequence files
				for ( vector< String >::iterator db_i = dbs.begin(); db_i != dbs.end(); ++db_i )
				{
					File::absolutePath(*db_i);

					if ( !File::exists(*db_i) ) not_accessable.push_back(*db_i);
					else if ( !File::readable(*db_i) ) not_accessable.push_back(*db_i);
					else if ( File::empty(*db_i) ) not_accessable.push_back(*db_i);
					else // if the file is accessable, try to find the corresponding index file and check it
					{
						if ( db_i->hasSuffix(".trie") ) string_buffer = db_i->substr(0, db_i->length()-4) + "index";
						else string_buffer = *db_i + "index";

						if ( !File::exists(string_buffer) ) not_accessable.push_back(*db_i);
						else if ( !File::readable(string_buffer) ) not_accessable.push_back(*db_i);
						else if ( File::empty(string_buffer) ) not_accessable.push_back(*db_i);
						else
						{
							accessable_db.push_back(*db_i);
							idx.push_back(string_buffer);
						}
					}
				}
			}

			vector< UInt > wanted_records;

			// creating the input file and converting and merging the databases
			if ( inspect_in )
			{
				if ( !seq_files.empty() || dbs.size() != 1 ) // don't do it, if only one trie database is given
				{
					// merging the trie databases (all but the first databases are appended)
					vector< String >::const_iterator idx_i = idx.begin();
					vector< UInt > v;
					for ( vector< String >::const_iterator db_i = accessable_db.begin(); db_i != accessable_db.end(); ++db_i, ++idx_i )
					{
						inspect_outfile.compressTrieDB(*db_i, *idx_i, v, db_filename,  idx_filename, (db_i != accessable_db.begin()) );
					}

					// converting and merging the other databases (all but the first database are appended)
					for ( vector< String >::const_iterator db_i = accessable_seq.begin(); db_i != accessable_seq.end(); ++db_i )
					{
						inspect_outfile.generateTrieDB(*db_i, db_filename,  idx_filename, ( (db_i != accessable_seq.begin()) || (!accessable_db.empty()) ));
					}
				}

				if ( blind_only ) inspect_infile.setBlind(true);

				inspect_infile.store(inspect_input_filename);
			}

			// running inspect and generating a second database from the results and running inspect in blind mode on this new database
			if ( blind && inspect_in && inspect_out )
			{
				writeLog_("Searching and generating minimised database for blind mode ...");
				writeDebug_("The Inspect process created the following output:", 1);
				String call;
				call.append(inspect_directory);
				call.append("inspect -r ");
				call.append(inspect_directory);
				call.append(" -i ");
				call.append(inspect_input_filename);
				call.append(" -o ");
				call.append(inspect_output_filename);
				// writing the inspect output to a temporary file
				call.append(" -e ");
				call.append(inspect_logfile);

				Int status = system(call.c_str());

				if (status != 0)
				{
					string_buffer = TextFile(inspect_logfile).asString();
					writeLog_("Inspect problem: " + string_buffer + " Aborting!");

					// deleting all temporary files
					for ( vector< pair< String, UInt > >::const_iterator files_i = files.begin(); files_i != files.end(); ++files_i )
					{
						if ( files_i->second & delete_afterwards ) remove(files_i->first.c_str());
					}
					return EXTERNAL_PROGRAM_ERROR;
				}

				vector< UInt > wanted_records = inspect_outfile.getWantedRecords(inspect_output_filename, p_value_threshold);

				if ( wanted_records.empty() )
				{
					IdXMLFile IdXML_file;
					IdXML_file.store(output_filename, vector<ProteinIdentification>(), vector<PeptideIdentification>());
					inspect_out = false;
					writeLog_("No proteins matching criteria for generating minimized database for blind search!");

					// deleting all temporary files
					for ( vector< pair< String, UInt > >::const_iterator files_i = files.begin(); files_i != files.end(); ++files_i )
					{
						if ( files_i->second & delete_afterwards ) remove(files_i->first.c_str());
					}
				}
				inspect_outfile.compressTrieDB(db_filename, idx_filename, wanted_records, snd_db_filename, snd_idx_filename, false);

				// setting the database name to the new database
				inspect_infile.setDb(snd_db_filename);
				inspect_infile.setBlind(true);
				inspect_infile.store(inspect_input_filename);
			}

			// writing the output of inspect Into an IdXML file
			if ( inspect_in && inspect_out )
			{
				String call;
				call.append(inspect_directory);
				call.append("inspect -r ");
				call.append(inspect_directory);
				call.append(" -i ");
				call.append(inspect_input_filename);
				call.append(" -o ");
				call.append(inspect_output_filename);
				// writing the inspect output to a temporary file
				call.append(" -e ");
				call.append(inspect_logfile);

				writeLog_("Searching ...");
				writeDebug_("The Inspect process created the following output:", 1);

				Int status = system(call.c_str());

				if (status != 0)
				{
					string_buffer = TextFile(inspect_logfile).asString();
					writeLog_("Inspect problem: " + string_buffer + " Aborting!");

					// deleting all temporary files
					for ( vector< pair< String, UInt > >::const_iterator files_i = files.begin(); files_i != files.end(); ++files_i )
					{
						if ( files_i->second & delete_afterwards ) remove(files_i->first.c_str());
					}
					return EXTERNAL_PROGRAM_ERROR;
				}
			}

			if ( inspect_out )
			{
				vector<PeptideIdentification> peptide_identifications;
				ProteinIdentification protein_identification;
				IdXMLFile IdXML_file;

				if ( inspect_in ) // the version can only be retrieved by running inspect without parameters
				{
					// first get the InsPecT version
					String call = inspect_directory;
					call.append("inspect > ");
					// writing the inspect output to a temporary file
					call.append(inspect_logfile);
					Int status = system(call.c_str());

					if (status != 0)
					{
						string_buffer = TextFile(inspect_logfile).asString();
						writeLog_("Inspect problem: " + string_buffer + " Aborting!");

						// deleting all temporary files
						for ( vector< pair< String, UInt > >::const_iterator files_i = files.begin(); files_i != files.end(); ++files_i )
						{
							if ( files_i->second & delete_afterwards ) remove(files_i->first.c_str());
						}
						return EXTERNAL_PROGRAM_ERROR;
					}

					// set the search engine and its version and the score type
					inspect_outfile.getSearchEngineAndVersion(inspect_logfile, protein_identification);
				}
				else protein_identification.setSearchEngine("InsPecT");

				if ( !File::empty(inspect_output_filename) )
				{
					// set the parameters
					ProteinIdentification::SearchParameters sp;
					if ( monoisotopic ) sp.mass_type = ProteinIdentification::MONOISOTOPIC;
					else sp.mass_type = ProteinIdentification::AVERAGE;
					if ( inspect_infile.getEnzyme() == "Trypsin" ) sp.enzyme = ProteinIdentification::TRYPSIN;
					else if ( inspect_infile.getEnzyme() == "No_Enzyme" ) sp.enzyme = ProteinIdentification::NO_ENZYME;
					else sp.enzyme = ProteinIdentification::UNKNOWN_ENZYME;
					sp.peak_mass_tolerance = inspect_infile.getPeakMassTolerance();
					sp.precursor_tolerance = inspect_infile.getPrecursorMassTolerance();
					protein_identification.setSearchParameters(sp);

					try
					{
						vector< UInt > corrupted_lines = inspect_outfile.load(inspect_output_filename, peptide_identifications, protein_identification, p_value_threshold, inspect_infile.getDb());
					}
					catch( Exception::ParseError pe )
					{
						// deleting all temporary files
						for ( vector< pair< String, UInt > >::const_iterator files_i = files.begin(); files_i != files.end(); ++files_i )
						{
							if ( files_i->second & delete_afterwards ) remove(files_i->first.c_str());
						}
						writeLog_(pe.getMessage());
						return INPUT_FILE_CORRUPT;
					}

					vector< ProteinIdentification > protein_identifications(1, protein_identification);
					IdXML_file.store(output_filename, protein_identifications, peptide_identifications);
				}
				else
				{
					IdXML_file.store(output_filename, vector<ProteinIdentification>(), vector<PeptideIdentification>());
					writeLog_("No proteins identified!");
				}
			}

			// deleting all temporary files
			for ( vector< pair< String, UInt > >::const_iterator files_i = files.begin(); files_i != files.end(); ++files_i )
			{
				if ( files_i->second & delete_afterwards ) remove(files_i->first.c_str());
			}

			return EXECUTION_OK;
		}
};

///@endcond



Int main( Int argc, char ** argv )
{
	TOPPInspectAdapter tool;

	return tool.main(argc,argv);
}
