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

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/InspectInfile.h>
#include <OpenMS/FORMAT/InspectOutfile.h>
#include <OpenMS/FORMAT/PTMXMLFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <cstdlib>
#include <vector>

#include <QtCore/QProcess>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_InspectAdapter InspectAdapter

	@brief Identifies peptides in MS/MS spectra via Inspect.

<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ InspectAdapter \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any signal-/preprocessing tool @n (in mzML format)</td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
		</tr>
	</table>
</CENTER>

	@experimental This tool has not been tested thoroughly and might behave not as expected!

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

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_InspectAdapter.cli
*/

// We do not want this class to show up in the docu -> cond
/// @cond TOPPCLASSES

class TOPPInspectAdapter
	: public TOPPBase
{
	public:
		TOPPInspectAdapter()
			: TOPPBase("InspectAdapter", "Annotates MS/MS spectra using Inspect.")
		{
		}

	protected:
		void registerOptionsAndFlags_()
		{
			registerInputFile_("in", "<file>", "", "input file in mzXML or mzData format.\n"
					 																			"Note: In mode 'inspect_out' an Inspect results file is read.");
			registerOutputFile_("out", "<file>", "", "output file in IdXML format.\n"
			                                           "Note: In mode 'inspect_in' an Inspect input file is written.");
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
      setValidStrings_("instrument",StringList::create("ESI-ION-TRAP,QTOF,FT-Hybrid"));
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
			registerOutputFile_("inspect_output", "<file>", "", "name for the output file of Inspect (may only be used in a full run)", false);
			registerInputFile_("inspect_input", "<file>", "", "name for the input file of Inspect (may only be used in a full run)", false);
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
			registerDoubleOption_("max_ptm_size", "<num>", 250.0, "maximum modification size (in Da) to consider", false);
			registerStringOption_("contact_name", "<name>", "unknown", "Name of the contact", false);
			registerStringOption_("contact_institution", "<name>", "unknown", "Name of the contact institution", false);
			registerStringOption_("contact_info", "<info>", "unknown", "Some information about the contact", false);
		}

		ExitCodes main_(Int , const char**)
		{
			//-------------------------------------------------------------
			// (1) variables
			//-------------------------------------------------------------

			InspectInfile inspect_infile;
			InspectOutfile inspect_outfile;

			vector< String >
				substrings,
				substrings2,
				trie_database_filenames,
				sequence_database_filenames,
				index_filenames;

			String
				string_buffer,
				trie_database_filename,
				index_filename,
				snd_trie_database_filename,
				snd_index_filename,
				inspect_logfile,
				logfile,
				inspect_directory,
				temp_data_directory,
				snd_trie_database,
				snd_trie_database_directory,
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

			DoubleReal p_value_threshold(1.0);
			DoubleReal cutoff_p_value(1.0);

			char separator = '/';

			ContactPerson contact_person;

			ExitCodes exit_code = EXECUTION_OK;

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
				catch ( Exception::ParseError& pe )
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

			logfile = getStringOption_("log");
			if ( logfile.empty() )
			{
				logfile = "temp.inspect.log";
				files[logfile] = (writable | delete_afterwards);
			}
			else files[logfile] = writable;

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
				temp_data_directory = File::absolutePath(temp_data_directory);
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
				string_buffer = File::absolutePath(string_buffer);
				if ( inspect_in )
				{
					MSExperiment<Peak1D> experiment;
					String type;
					try
					{
						inspect_outfile.getExperiment(experiment, type, string_buffer); // may throw an exception if the filetype could not be determined
					}
					catch(Exception::ParseError& pe )
					{
						writeLog_(pe.getMessage());
						return PARSE_ERROR;
					}
					if ( type != "mzXML" )
					{
						string_buffer.append(".mzXML");
						MzXMLFile().store(string_buffer, experiment);
						files[string_buffer] = (writable | delete_afterwards);
					}
					inspect_infile.setSpectra(string_buffer);

					if ( inspect_out )
					{
						inspect_output_filename = getStringOption_("inspect_output");
						if ( inspect_output_filename.empty() )
						{
							inspect_output_filename = temp_data_directory + "tmp.direct.inspect.output";
							files[inspect_output_filename] = (writable | delete_afterwards);
						}
						else
						{
							inspect_output_filename = File::absolutePath(inspect_output_filename);
							files[inspect_output_filename] = writable;
						}
					}
				}
				else
				{
					inspect_output_filename = string_buffer;
					inspect_output_filename = File::absolutePath(inspect_output_filename);
					files[inspect_output_filename] = readable;
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
				string_buffer = File::absolutePath(string_buffer);
				if ( inspect_out ) output_filename = string_buffer;
				else inspect_input_filename = string_buffer;
				files[string_buffer] = writable;
			}

			if ( inspect_in && inspect_out )
			{
				inspect_input_filename = getStringOption_("inspect_input");
				if ( inspect_input_filename.empty() )
				{
					inspect_input_filename = temp_data_directory + "tmp.inspect.input";
					files[inspect_input_filename] = (writable | delete_afterwards);
				}
				else
				{
					inspect_input_filename = File::absolutePath(inspect_input_filename);
					files[inspect_input_filename] = writable;
				}
			}

			inspect_directory = getStringOption_("inspect_directory");
			if ( inspect_in && inspect_directory.empty() && inspect_out )
			{
				writeLog_("No inspect directory file specified. Aborting!");
				return ILLEGAL_PARAMETERS;
			}
			inspect_directory = File::absolutePath(inspect_directory);
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
					string_buffer.split(',', trie_database_filenames);
					if ( trie_database_filenames.empty() ) trie_database_filenames.push_back(string_buffer);

					// the database files have to be readable, (by the way changing the names using the absolute path)
					for ( vector< String >::iterator trie_database_filenames_it = trie_database_filenames.begin(); trie_database_filenames_it != trie_database_filenames.end(); ++trie_database_filenames_it )
					{
						*trie_database_filenames_it = File::absolutePath(*trie_database_filenames_it);
						files[*trie_database_filenames_it] = readable;

						// get the according index file
						if ( trie_database_filenames_it->hasSuffix(".trie") ) string_buffer = trie_database_filenames_it->substr(0, trie_database_filenames_it->length()-4) + "index";
						else string_buffer = *trie_database_filenames_it + "index";
						index_filenames.push_back(string_buffer);
						files[string_buffer] = readable;
					}
				}

				string_buffer = getStringOption_("dbs");
				if ( !string_buffer.empty() )
				{
					// get the single sequence files
					string_buffer.split(',', sequence_database_filenames);
					if ( sequence_database_filenames.empty() ) sequence_database_filenames.push_back(string_buffer);
					// the sequence files have to be readable, (by the way changing the names using the absolute path)
					for ( vector< String >::iterator sequence_database_filenames_it = sequence_database_filenames.begin(); sequence_database_filenames_it != sequence_database_filenames.end(); ++sequence_database_filenames_it )
					{
						*sequence_database_filenames_it = File::absolutePath(*sequence_database_filenames_it);
						files[*sequence_database_filenames_it] = readable;
					}
				}

				// at least one of the parameters db or sequence_file has to be set
				if ( trie_database_filenames.empty() && sequence_database_filenames.empty() )
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

				trie_database_filename = getStringOption_("new_db");
				if ( trie_database_filename.empty() && (!sequence_database_filenames.empty() || trie_database_filenames.size() != 1) )
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
							trie_database_filename = temp_data_directory + "tmp.inspect.db.trie";
							files[trie_database_filename] = (writable | delete_afterwards);
							inspect_infile.setDb(trie_database_filename);
							index_filename = temp_data_directory + "tmp.inspect.db.index";
							files[index_filename] = (writable | delete_afterwards);
						}
					}
				}
				else
				{
					// if only one trie database is given, this one is used
					if ( trie_database_filename.empty() ) trie_database_filename = trie_database_filenames.front();

					trie_database_filename = File::absolutePath(trie_database_filename);
					if ( trie_database_filename.hasSuffix(".trie") )
					{
						inspect_infile.setDb(trie_database_filename);
						index_filename = trie_database_filename.substr(0, trie_database_filename.length()-4) + "index";
					}
					else
					{
						index_filename = trie_database_filename + ".index";
						trie_database_filename = trie_database_filename + ".trie";
						inspect_infile.setDb(trie_database_filename);
					}
					files[trie_database_filename] = writable;
					files[index_filename] = writable;
				}

				if ( blind && blind_only )
				{
					writeLog_("Both blind flags set. Only one of the two flags [-blind|-blind_only] can be set. Aborting!");
					return ILLEGAL_PARAMETERS;
				}

				snd_trie_database = getStringOption_("snd_db");
				if ( no_tmp_dbs && blind && snd_trie_database.empty() )
				{
					writeLog_("No_tmp_dbs and blind flag set but no name for minimized database given. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				else if ( blind && snd_trie_database.empty() )
				{
					snd_trie_database_filename = temp_data_directory + "tmp.inspect.db.snd.trie";
					snd_index_filename = temp_data_directory + "tmp.inspect.db.snd.index";
					files[snd_trie_database_filename] = (writable | delete_afterwards);
					files[snd_index_filename] = (writable | delete_afterwards);
				}
				else if ( blind )
				{
					snd_trie_database = File::absolutePath(snd_trie_database);
					if ( snd_trie_database.hasSuffix(".trie") )
					{
						snd_trie_database_filename = snd_trie_database;
						snd_index_filename = snd_trie_database.substr(0, snd_trie_database.size()-4) + "index";
						files[snd_trie_database_filename] = writable;
						files[snd_index_filename] = writable;
					}
					else
					{
						snd_trie_database_filename = snd_trie_database + ".trie";
						snd_index_filename = snd_trie_database + ".index";
						files[snd_trie_database_filename] = writable;
						files[snd_index_filename] = writable;
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
					catch ( Exception::FileNotFound& /*fnf_e*/ )
					{
						writeLog_("No modifications XML file given. Aborting!");
						return INPUT_FILE_NOT_FOUND;
					}
					catch ( Exception::FileNotReadable& /*fnr_e*/ )
					{
						writeLog_("Modifications XML file is not readable. Aborting!");
						return INPUT_FILE_NOT_READABLE;
					}
					catch ( Exception::ParseError& p_e )
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

				inspect_infile.setMaxPTMsize(getDoubleOption_("max_ptm_size") );
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
				files[inspect_logfile] = (writable | delete_afterwards);
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

			vector< Size > wanted_records;

			// creating the input file and converting and merging the databases
			if ( exit_code == EXECUTION_OK && inspect_in )
			{
				if ( !sequence_database_filenames.empty() || trie_database_filenames.size() != 1 ) // don't do it, if only one trie database is given
				{
					// merging the trie databases (all but the first databases are appended)
					vector< String >::const_iterator index_filenames_itt = index_filenames.begin();
					for ( vector< String >::const_iterator trie_database_filenames_it = trie_database_filenames.begin(); trie_database_filenames_it != trie_database_filenames.end(); ++trie_database_filenames_it, ++index_filenames_itt )
					{
						inspect_outfile.compressTrieDB(*trie_database_filenames_it, *index_filenames_itt, wanted_records, trie_database_filename,  index_filename, (trie_database_filenames_it != trie_database_filenames.begin()) );
					}

					// converting and merging the other databases (all but the first database are appended)
					for ( vector< String >::const_iterator sequence_database_filenames_it = sequence_database_filenames.begin(); sequence_database_filenames_it != sequence_database_filenames.end(); ++sequence_database_filenames_it )
					{
						inspect_outfile.generateTrieDB(*sequence_database_filenames_it, trie_database_filename,  index_filename, ( (sequence_database_filenames_it != sequence_database_filenames.begin()) || (!sequence_database_filenames.empty()) ));
					}
				}

				if ( blind_only ) inspect_infile.setBlind(true);

				inspect_infile.store(inspect_input_filename);
			}

			// running inspect and generating a second database from the results and running inspect in blind mode on this new database
			if ( exit_code == EXECUTION_OK && blind && inspect_in && inspect_out )
			{
				writeLog_("Searching and generating minimised database for blind mode ...");
				writeDebug_("The Inspect process created the following output:", 1);
				String call;
				call.append(" -r ");
				call.append(inspect_directory);
				call.append(" -i ");
				call.append(inspect_input_filename);
				call.append(" -o ");
				call.append(inspect_output_filename);
				// writing the inspect output to a temporary file
				call.append(" -e ");
				call.append(inspect_logfile);

        Int status = QProcess::execute((inspect_directory+"inspect").toQString(), QStringList(call.toQString().split(" ", QString::SkipEmptyParts))); // does automatic escaping etc...
				if (status != 0)
				{
					string_buffer = TextFile(inspect_logfile).concatenate();
					writeLog_("Inspect problem: " + string_buffer + " Aborting!");

					exit_code = EXTERNAL_PROGRAM_ERROR;
				}

				wanted_records = inspect_outfile.getWantedRecords(inspect_output_filename, p_value_threshold);

				if ( wanted_records.empty() )
				{
					IdXMLFile IdXML_file;
					IdXML_file.store(output_filename, vector<ProteinIdentification>(), vector<PeptideIdentification>());
					inspect_out = false;
					writeLog_("No proteins matching criteria for generating minimized database for blind search. Aborting!");
					exit_code = UNKNOWN_ERROR;
				}
				else
				{
					inspect_outfile.compressTrieDB(trie_database_filename, index_filename, wanted_records, snd_trie_database_filename, snd_index_filename, false);

					// setting the database name to the new database
					inspect_infile.setDb(snd_trie_database_filename);
					inspect_infile.setBlind(true);
					inspect_infile.store(inspect_input_filename);
				}
			}

			// writing the output of inspect Into an IdXML file
			if ( exit_code == EXECUTION_OK && inspect_in && inspect_out )
			{
				String call;
				call.append(" -r ");
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

        Int status = QProcess::execute((inspect_directory+"inspect").toQString(), QStringList(call.toQString().split(" ", QString::SkipEmptyParts))); // does automatic escaping etc...
				if (status != 0)
				{
					string_buffer = TextFile(inspect_logfile).concatenate();
					writeLog_("Inspect problem: " + string_buffer + ". Aborting!");
					exit_code =  EXTERNAL_PROGRAM_ERROR;
				}
			}

			if ( exit_code == EXECUTION_OK && inspect_out )
			{
				vector<PeptideIdentification> peptide_identifications;
				ProteinIdentification protein_identification;
				IdXMLFile IdXML_file;

				if ( inspect_in ) // the version can only be retrieved by running inspect without parameters
				{
					// first get the InsPecT version
          QProcess builder;
          builder.start((inspect_directory+"inspect").toQString(), QStringList()); // does automatic escaping etc...

          if (!builder.waitForFinished(-1))
          {
						writeLog_("Inspect problem: " + String(QString(builder.readAll())) + ". Aborting!");
						exit_code = EXTERNAL_PROGRAM_ERROR;
					}
					else
					{
            QString output = builder.readAll();
            // set the search engine and its version and the score type
            if (!inspect_outfile.getSearchEngineAndVersion(output, protein_identification)) LOG_WARN << "Could not read version of InsPecT from:\n" << String(output) << "\n\n";
					}
				}
				else protein_identification.setSearchEngine("InsPecT");

				if ( exit_code == EXECUTION_OK )
				{
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
							vector< Size > corrupted_lines = inspect_outfile.load(inspect_output_filename, peptide_identifications, protein_identification, p_value_threshold, inspect_infile.getDb());
						}
						catch( Exception::ParseError& pe )
						{
							writeLog_(pe.getMessage());
							exit_code = INPUT_FILE_CORRUPT;
						}
						if ( exit_code == EXECUTION_OK )
						{
							vector< ProteinIdentification > protein_identifications(1, protein_identification);
							IdXML_file.store(output_filename, protein_identifications, peptide_identifications);
						}
					}
					else
					{
						IdXML_file.store(output_filename, vector<ProteinIdentification>(), vector<PeptideIdentification>());
						writeLog_("No proteins identified!");
					}
				}
			}

			// if an external program error occured, the logfile must not be deleted
			if ( exit_code == EXTERNAL_PROGRAM_ERROR )
			{
				writeLog_("PepNovo problem. Aborting! (Details can be seen in the logfile: \"" + logfile + "\")");
				files[logfile] = readable;
			}
			// deleting all temporary files
			for ( map< String, Size >::const_iterator files_it = files.begin(); files_it != files.end(); ++files_it )
			{
				if ( files_it->second & delete_afterwards ) remove(files_it->first.c_str());
			}

			return exit_code;
		}
};

///@endcond



Int main( Int argc, const char** argv )
{
	TOPPInspectAdapter tool;

	return tool.main(argc,argv);
}
