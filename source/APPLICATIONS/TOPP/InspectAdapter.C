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

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/AnalysisXMLFile.h>
#include <OpenMS/FORMAT/InspectInfile.h>
#include <OpenMS/FORMAT/InspectOutfile.h>
#include <OpenMS/FORMAT/IsotopeXMLFile.h>
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
	
	This wrapper application serves for getting peptide identifications
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
				
				@code ./inspect -i  inputfilename -o outputfilename  @endcode

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
			registerStringOption_("out", "<file>", "", "output file in analysisXML format.\n"
			                                           "Note: In mode 'inspect_in' an Inspect input file is written.");
			registerStringOption_("in", "<file>", "", "input file in mzXML format OR directory to search in.\n"
					 																			"Note: In mode 'inspect_out' an Inspect results file is read");
			registerFlag_("inspect_in", "if this flag is set the InspectAdapter will read in mzXML,\n"
																							 "write an Inspect input file and generate a trie database");
			registerFlag_("inspect_out", "if this flag is set the InspectAdapter will read in a Inspect results file\n"
																								 "and write analysisXML");
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
																																													 "<composition>,<residues>,<type>,<name>, e.g.: H(2).C(2).O,KCS,opt,Acetyl or\n"
																																													 "<mass>,<residues>,<type>,<name>, e.g.: 42.0367,KCS,opt,Acetyl or\n"
																																													 "Valid values for \"type\" are \"fix\", \"cterminal\", \"nterminal\",\n"
																																													 "and \"opt\" (the default).\n", false);
			registerFlag_("use_monoisotopic_mod_mass", "use monoisotopic masses for the modifications");
			registerStringOption_("modifications_xml_file", "<file>", "", "name of an XML file with the modifications", false);
			registerStringOption_("isotopes_xml_file", "<file>", "", "name of an XML file with the masses and probabilities of isotopes", false);
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
			registerDoubleOption_("maxptmsize", "<num>", 250.0, "maximum modification size (in Da) to consider", false);
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
			SignedInt pos, pos2;
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
				if ( pos2 != String::NPOS ) // if the element occurs more than once, a bracket is found
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
				SignedInt i_iso, i_occ;
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

		ExitCodes main_(int , char**)
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
				no_tmp_dbs(false);
			
			Real p_value_threshold = 1.0;
			Real cutoff_p_value;
			
			char separator = '/';
			
			ContactPerson contact_person;
			
			vector< String > tmp_names;
			
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
			}
			
			inspect_in = getFlag_("inspect_in");
			inspect_out = getFlag_("inspect_out");
			
			if ( inspect_in && inspect_out )
			{
				writeLog_("Both Inspect flags set. Only one of the two flags [-inspect_in|-inspect_out] can be set. Aborting!");
				return ILLEGAL_PARAMETERS;
			}
			
			if ( inspect_in ) writeDebug_("Inspect flag: mascot_in (reads in MzXML/MzData, writes Inspect generic format)", 1);
			else if ( inspect_out ) writeDebug_("Inspect flag: mascot_in (reads in Inspect result file, writes analysisXML file)", 1);
			else writeDebug_("No Inspect flag set: reads in MzXML/MzData, writes analysisXML file", 1);
			
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
						//tmp_names.push_back(string_buffer);
					}
					inspect_infile.setSpectra(string_buffer);
					
					if ( inspect_out )
					{
						inspect_output_filename = getStringOption_("inspect_output");
						if ( inspect_output_filename.empty() )
						{
							inspect_output_filename = temp_data_directory + "tmp.direct.inspect.output";
							tmp_names.push_back(inspect_output_filename);
						}
					}
				}
				else inspect_output_filename = string_buffer;
				File::absolutePath(inspect_output_filename);
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
			}
			
			if ( inspect_in && inspect_out )
			{
				inspect_input_filename = getStringOption_("inspect_input");
				if ( inspect_input_filename.empty() )
				{
					inspect_input_filename = temp_data_directory + "tmp.inspect.input";
					tmp_names.push_back(inspect_input_filename);
				}
				
				File::absolutePath(inspect_input_filename);
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
							tmp_names.push_back(db_filename);
							inspect_infile.setDb(db_filename);
							idx_filename = temp_data_directory + "tmp.inspect.db.index";
							tmp_names.push_back(idx_filename);
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
					tmp_names.push_back(snd_db_filename);
					tmp_names.push_back(snd_idx_filename);
				}
				else if ( blind )
				{
					File::absolutePath(snd_db);
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
				
				// get the known modifications
				if ( !blind_only )
				{
					string_buffer = getStringOption_("modifications");
					bool monoisotopic = getFlag_("use_monoisotopic_mod_mass");
					if ( !string_buffer.empty() ) // if modifications are used get look whether whether composition and residues (and type and name) is given which needs the isotope file, the name (and type) is used (then one additionally needs the modifications file) or only the mass and residues (and type and name) is given, in which case no further file is needed
					{
						string_buffer.split(':', substrings); // get the single modifications
						
						// one vector if compositions are used (needs isotope xml file) and one vector if masses were given
						vector< vector< String > > iso_sym_occ, mass_res_type_name;
						
						// to store the informations about modifications from the ptm xml file
						map< String, pair< String, String > > ptm_informations;
						
						// to store the informations about isotopes from the isotopes xml file
						map< String, vector< pair< DoubleReal, DoubleReal > > > isotopes_mass_and_probability;
						map< String, DoubleReal > isotope_masses;
						
						UnsignedInt comp_mass_name_given;
						String types = "opt#fix#cterminal#nterminal";
						
						for ( vector< String >::const_iterator mod_i = substrings.begin(); mod_i != substrings.end(); ++mod_i )
						{
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
							bool is_mass = false;
							try
							{
								is_mass = ( String(substrings2[0].toDouble()) == substrings2[0] );
							}
							catch ( Exception::ConversionError ce ) {}
							if ( is_mass ) // if it's a mass
							{
								mass_res_type_name.back()[0] = substrings2[0]; // mass
								comp_mass_name_given = 0;
							}
							else if ( get_composition_elements(substrings2[0], iso_sym_occ, '.').empty() ) // if it is a composition, put it into the vector
							{
								mass_res_type_name.back()[0] = substrings2[0]; // composition
								comp_mass_name_given = 1;
							}
							else // if it's a name, try to find it in the ptm xml file
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
									if ( types.find(substrings2[1]) == string::npos )
									{
										writeLog_("The given type (" + substrings2[1] + ") is neither opt, fix, cterminal nor nterminal. Aborting!");
										return ILLEGAL_PARAMETERS;
									}
									mass_res_type_name.back()[2] = substrings2[1];
								}
								else mass_res_type_name.back()[2] = "opt";
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
								// get the residues
								mass_res_type_name.back()[1] = substrings2[1];
								
								// get the type
								if ( substrings2.size() > 2 )
								{
									// if it's not a legal type
									if ( types.find(substrings2[2]) == string::npos )
									{
										writeLog_("The given type (" + substrings2[2] + ") is neither opt, fix, cterminal nor nterminal. Aborting!");
										return ILLEGAL_PARAMETERS;
									}
									mass_res_type_name.back()[2] = substrings2[2];
									
									// get the name
									if ( substrings2.size() > 3 ) mass_res_type_name.back()[3] = substrings2[3];
								}
								else mass_res_type_name.back()[2] = "opt";
							}
							
							// if a composition is given, get the corresponding mass
							if ( comp_mass_name_given )
							{
								if ( isotope_masses.empty() ) // if the isotopes xml file has not been read yet, read it
								{
									string isotopes_filename = getStringOption_("isotopes_xml_file");
									if ( isotopes_filename.empty() )
									{
										writeLog_("No isotopes XML file given. Aborting!");
										return INPUT_FILE_NOT_FOUND;
									}
									if ( !File::readable(isotopes_filename) )
									{
										writeLog_("Isotopes XML file is not readable. Aborting!");
										return INPUT_FILE_NOT_READABLE;
									}
									
									try
									{
										IsotopeXMLFile().load(isotopes_filename, isotopes_mass_and_probability);
									}
									catch ( Exception::ParseError pe )
									{
										writeLog_(pe.getMessage());
										return PARSE_ERROR;
									}
									
									DoubleReal mass, probability; // compute the monoisotopic or average mass
									for ( map< String, vector< pair< DoubleReal, DoubleReal > > >::const_iterator iso_i = isotopes_mass_and_probability.begin(); iso_i != isotopes_mass_and_probability.end(); ++iso_i )
									{
										mass = probability = 0;
										for ( vector< pair< DoubleReal, DoubleReal > >::const_iterator mp_i = iso_i->second.begin(); mp_i != iso_i->second.end(); ++mp_i )
										{
											if ( monoisotopic )
											{
												if ( probability < mp_i->second )
												{
													mass = mp_i->first;
													probability = mp_i->second;
												}
											}
											else mass += mp_i->first * mp_i->second;
										}
										isotope_masses[iso_i->first] = mass;
									}
								}
								
								// compute the mass
								DoubleReal mass = 0;
								// get the single components of the composition
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
										mass += isotope_masses[(*comp_i)[1]] * (*comp_i)[2].toDouble();
									}
									else // if an isotope was used, get the mass
									{
										for ( vector< pair< DoubleReal, DoubleReal > >::const_iterator iso_i = isotopes_mass_and_probability[(*comp_i)[1]].begin(); iso_i != isotopes_mass_and_probability[(*comp_i)[1]].end(); ++iso_i )
										{
											if ( ((SignedInt) (iso_i->first + 0.5)) == (*comp_i)[0].toDouble() ) // round the mass
											{
												mass += iso_i->first * (*comp_i)[2].toDouble();
												break;
											}
										}
									}
								}
								mass_res_type_name.back()[0] = String(mass);
							}
						}
						
						inspect_infile.setMod(mass_res_type_name);
					}
				}
				
				inspect_infile.setProtease(getStringOption_("cleavage"));
				inspect_infile.setInstrument(getStringOption_("instrument"));
				
				inspect_infile.setMods(getIntOption_("max_modifications_pp"));
				if ( inspect_infile.getMods() < 1 && !mod.empty() )
				{
					writeLog_("Modifications specified, but max_modifications_pp not set. Setting it to 1.");
					inspect_infile.setMods(1);
				}
				
				inspect_infile.setPMTolerance(getDoubleOption_("precursor_mass_tolerance"));
				if ( (inspect_infile.getPMTolerance() < 0 && inspect_infile.getPMTolerance() != -1) )
				{
					writeLog_("Illegal precursor mass tolerance (<0) given. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				
				inspect_infile.setIonTolerance( getDoubleOption_("peak_mass_tolerance") );
				if ( (inspect_infile.getIonTolerance() < 0 && inspect_infile.getIonTolerance() != -1) )
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
				
				inspect_infile.setMaxPTMsize(getDoubleOption_("maxptmsize") );
				if ( inspect_infile.getMaxPTMsize() < 0 && inspect_infile.getMaxPTMsize() != -1)
				{
					writeLog_("Illegal maximum modification size (<0). Aborting!");
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
				tmp_names.push_back(inspect_logfile);
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
			
			// the file for the inspect output
			if ( (inspect_in && inspect_out) || (inspect_in && blind) )
			{
				if ( !File::writable(inspect_output_filename) )
				{
					throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, inspect_output_filename);
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
				if ( !seq_files.empty() || dbs.size() != 1 )
				{
					if ( !File::writable(db_filename) )
					{
						throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, db_filename);
					}
					if ( !File::writable(idx_filename) )
					{
						throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, idx_filename);
					}
				}
				
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
				
				for ( vector< String >::iterator db_i = seq_files.begin(); db_i != seq_files.end(); ++db_i )
				{
					File::absolutePath(*db_i);
					
					if ( !File::exists(*db_i) ) not_accessable.push_back(*db_i);
					else if ( !File::readable(*db_i) ) not_accessable.push_back(*db_i);
					else if ( File::empty(*db_i) ) not_accessable.push_back(*db_i);
					else accessable_seq.push_back(*db_i);
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
								+ temp_data_directory + inspect_logfile
								+ String(" Aborting!"));
						return ILLEGAL_PARAMETERS;
					}
				}
			}
			
			vector< unsigned int > wanted_records;
			
			// creating the input file and converting and merging the databases
			if ( inspect_in )
			{
				if ( !seq_files.empty() || dbs.size() != 1 ) // don't do it, if only one trie database is given
				{
					// merging the trie databases (all but the first databases are appended)
					vector< String >::const_iterator idx_i = idx.begin();
					vector< UnsignedInt > v;
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
				writeDebug_("Searching and generating minimised database for blind mode ...", 1);
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

				int status = system(call.c_str());

				if (status != 0)
				{
					string_buffer = TextFile(inspect_logfile).asString();
					writeLog_("Inspect problem: " + string_buffer + " Aborting!");
					for ( vector< String >::const_iterator tmp_names_i = tmp_names.begin(); tmp_names_i != tmp_names.end(); ++tmp_names_i ) remove(tmp_names_i->c_str());
					return EXTERNAL_PROGRAM_ERROR;
				}
				
				vector< UnsignedInt > wanted_records = inspect_outfile.getWantedRecords(inspect_output_filename, p_value_threshold);
				
				if ( wanted_records.empty() )
				{
					AnalysisXMLFile analysisXML_file;
					analysisXML_file.store(output_filename, vector< ProteinIdentification >(), vector< IdentificationData >());
					inspect_out = false;
					writeLog_("No proteins matching criteria for generating minimized database for blind search!");
					
					for ( vector< String >::const_iterator tmp_names_i = tmp_names.begin(); tmp_names_i != tmp_names.end(); ++tmp_names_i ) remove(tmp_names_i->c_str());
				}
				inspect_outfile.compressTrieDB(db_filename, idx_filename, wanted_records, snd_db_filename, snd_idx_filename, false);
				
				// setting the database name to the new database
				inspect_infile.setDb(snd_db_filename);
				inspect_infile.setBlind(true);
				inspect_infile.getMod().clear();
				inspect_infile.store(inspect_input_filename);
			}
			
			// writing the output of inspect into an analysisXML file
			if ( inspect_in && inspect_out )
			{
				writeDebug_("Searching ...", 1);
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
				
				int status = system(call.c_str());
				
				if (status != 0)
				{
					string_buffer = TextFile(inspect_logfile).asString();
					writeLog_("Inspect problem: " + string_buffer + " Aborting!");
					for ( vector< String >::const_iterator tmp_names_i = tmp_names.begin(); tmp_names_i != tmp_names.end(); ++tmp_names_i ) remove(tmp_names_i->c_str());
					return EXTERNAL_PROGRAM_ERROR;
				}
			}
			
			if ( inspect_out )
			{
				AnalysisXMLFile analysisXML_file;
				
				if ( !File::empty(inspect_output_filename) )
				{
					vector< IdentificationData > identifications;
					ProteinIdentification protein_identification;
					
					try
					{
						vector< UnsignedInt > corrupted_lines = inspect_outfile.load(inspect_output_filename, identifications, protein_identification, p_value_threshold);
					}
					catch( Exception::ParseError pe )
					{
						for ( vector< String >::const_iterator tmp_names_i = tmp_names.begin(); tmp_names_i != tmp_names.end(); ++tmp_names_i ) remove(tmp_names_i->c_str());
						writeLog_(pe.getMessage());
						return INPUT_FILE_CORRUPT;
					}
					
					vector< ProteinIdentification > protein_identifications;
					protein_identifications.push_back(protein_identification);
					
					analysisXML_file.store(output_filename, protein_identifications, identifications);
				}
				else
				{
					analysisXML_file.store(output_filename, vector< ProteinIdentification >(), vector< IdentificationData >());
					writeLog_("No proteins identified!");
				}
			}
			
			// (3.3) deleting all temporary files
			for ( vector< String >::const_iterator tmp_names_i = tmp_names.begin(); tmp_names_i != tmp_names.end(); ++tmp_names_i ) remove(tmp_names_i->c_str());
			
			return EXECUTION_OK;
		}
};

///@endcond



int main( int argc, char ** argv )
{
	TOPPInspectAdapter tool;

	return tool.main(argc,argv);
}
