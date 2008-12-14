// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// --------------------------------------------------------------------------


#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/PepNovoInfile.h>
#include <OpenMS/FORMAT/PepNovoOutfile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/ContactPerson.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/PTMXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>

#include <cstdlib>
#include <vector>
#include <algorithm>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------


/**
	@page TOPP_PepNovoAdapter PepNovoAdapter

	@brief Identifies peptides in MS/MS spectra via PepNovo.

	This wrapper application serves for getting peptide identifications
	for MS/MS spectra. The wrapper can be executed in three different
	modes:
	<ol>
				<li>
				The whole process of identification via PepNovo is executed.
				Inputfile is one (or more) mz file containing the MS/MS spectra
				(Supported spectrum file formats are .mzXML, .mzData)
				for which the identifications are to be found. The results are written
				as an idXML output file. This mode is selected by default.
			 	</li>

				<li>
				Only the first part of the ProteinIdentification process is performed.
				This means that a PepNovo input file is generated and dta files are
				created from the mz file.
				The call for the corresponding DeNovo process is written to standard
				output.

				Consult your PepNovo reference manual for further details.

				This mode is selected by the <b>-pepnovo_in</b> option in the command line.
				</li>

				<li>
				Only the second part of the ProteinIdentification process is performed.
				This means that the output of pepnovo is translated into idXML.

				This mode is selected by the <b>-pepnovo_out</b> option in the command line.
				</li>
	</ol>
	
	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_PepNovoAdapter.cli
*/

// We do not want this class to show up in the docu -> cond
// @cond

class TOPPPepNovoAdapter
	: public TOPPBase
{
	public:
		TOPPPepNovoAdapter()
			: TOPPBase("PepNovoAdapter", "Annotates MS/MS spectra using PepNovo.")
		{
		}

	protected:

		void registerOptionsAndFlags_()
		{
			addText_("The definitions for the parameters are taken from the site:\n"
										 "http://www.grosse-coosmann.de/~florian/Parameters.html#file.");
			registerInputFile_("in", "<file>", "", "input file(s) in mzXML or mzData format (comma-separated).\n"
					 																			"Note: In mode 'pepnovo_out' a directory with PepNovo results files\n"
																								"(*.out) is read", false);
			registerOutputFile_("out", "<file>", "", "output file in idXML format.\n"
			                                           "Note: In mode 'pepnovo_in' a PepNovo input file is written.", false);
			registerFlag_("pepnovo_in", "if this flag is set the PepNovoAdapter will read in mzXML or mzData\n"
																								"and write an PepNovo input file\n"
																								"and create dta files from the given mzXML or mzData files");
			registerFlag_("pepnovo_out", "if this flag is set the PepNovoAdapter will read in PepNovo result files\n"
																									"and write idXML");
			registerStringOption_("mz_files", "<file>", "", "when using pepnovo_out the mzXML or mzData files (comma-separated)\n"
																																							"have to be given to retrieve the retention times", false);
			registerStringOption_("pepnovo_directory", "<dir>", "", "the PepNovo working directory", false);
			registerStringOption_("temp_data_directory", "<dir>", "", "a directory in which some temporary files can be stored", false);
			registerStringOption_("charges", "[1,3,5]", "", "comma-seperated list of charge states (or ranges).", false);
			registerStringOption_("model_directory", "<file>", "", "name of the directory where the model files are kept.");
			registerFlag_("list_models", "show a list of the available models");
			registerStringOption_("model", "<file>", "", "name of the model that should be used (e.g. tryp_model.txt).");
			registerStringOption_("cleavage", "<enz>", "Trypsin", "the name of the enzyme used for digestion (currently there's only distinction\nbetween Trypsin and everything else)", false);
			registerIntOption_("max_number_of_tags", "<num>", -1, "maximal number of tags used (zero means not set).", false);
			registerStringOption_("dta_list", "<file>", "", "name of the file that holds the names of the dta files (created from the input) to be\nsearched. This name has to be given, if pepnovo_in is used only!", false);
			registerDoubleOption_("precursor_mass_tolerance", "<tol>", -1 , "the precursor mass tolerance", false);
			registerDoubleOption_("peak_mass_tolerance", "<tol>", -1, "the peak mass tolerance", false);
			registerFlag_("list_modifications", "show a list of the available modifications");
			registerStringOption_("modifications", "<mods>", "", "the colon-seperated modifications; may be\n"
																																														"<name>,<type>, e.g.: Deamidation,opt or\n"
																																														"<composition>,<residues>,<type>,<name>, e.g.: H2C2O,KCS,opt,Acetyl or\n"
																																														"<mass>,<residues>,<type>,<name>, e.g.: 42.0367,KCS,opt,Acetyl or\n"
																																														"Valid values for type are \"fix\" and \"opt\" (default)\n"
																																														"If you want terminal PTMs, write \"cterm\" or \"nterm\" instead of residues", false);
			registerFlag_("use_monoisotopic_mod_mass", "use monoisotopic masses for the modifications");
			registerStringOption_("modifications_xml_file", "<file>", "", "name of an XML file with the modifications", false);
			registerDoubleOption_("p_value", "<prob>", 1.0, "annotations with inferior p-value are ignored", false);
			registerIntOption_("min_sequence_length", "<min>", 3, "minimal number of amino acids in predicted sequence", false);
			registerIntOption_("max_sequence_length", "<max>", 40, "maximal number of amino acids in predicted sequence", false);
			registerIntOption_("num_results", "<num>", 20, "the number of possible peptides per scan", false);
			registerFlag_("keep_dta_files", "If set, the dta-files that were created from the mzXML or mzData files are not removed");
			registerOutputFile_("pepnovo_output", "<file>", "", "name for the output file of PepNovo (may only be used in a full run)", false);
			registerInputFile_("pepnovo_input", "<file>", "", "name for the input file of PepNovo (may only be used in a full run)", false);
			registerStringOption_("contact_name", "<name>", "unknown", "Name of the contact", false);
			registerStringOption_("contact_institution", "<name>", "unknown", "Name of the contact institution", false);
			registerStringOption_("contact_info", "<info>", "unknown", "Some information about the contact", false);
		}

  UInt MSExperiment2DTAs(MSExperiment<Peak1D>& msexperiment, const String& common_name,	const vector< Int >& charges,	map< String, Real >& dta_filenames_and_precursor_retention_times, bool make_dtas = true)
		throw (Exception::UnableToCreateFile)
		{
			DTAFile dtafile;
			String filename;
			UInt scan_number(0);
			UInt msms_spectra(0);

			for ( MSExperiment<Peak1D>::Iterator spec_it = msexperiment.begin(); spec_it != msexperiment.end(); ++spec_it )
			{
				++scan_number;
				if ( (spec_it->getMSLevel() == 2) && (!spec_it->empty()) )
				{
					++msms_spectra;
					if ( spec_it->getPrecursorPeak().getCharge() )
					{
						filename = common_name + "." + String(scan_number) + "." + String(spec_it->getPrecursorPeak().getCharge()) + ".dta";
						if ( make_dtas ) dtafile.store(filename, *spec_it);
						dta_filenames_and_precursor_retention_times[File::basename(filename)] = spec_it->getRT();
					}
					else
					{
						for ( vector< Int >::const_iterator i = charges.begin(); i != charges.end(); ++i )
						{
							filename = common_name + "." + String(scan_number) + "." + *i + ".dta";
							// for PepNovo the precursor mass may not be less than the highest peak mass
							if ( spec_it->back().getPosition()[0] < ((spec_it->getPrecursorPeak().getPosition()[0] - 1.0) * (*i) +1.0) )
							{
								if ( make_dtas )
								{
									spec_it->getPrecursorPeak().setCharge(*i);
									dtafile.store(filename, *spec_it);
								}
								dta_filenames_and_precursor_retention_times[File::basename(filename)] = spec_it->getRT();
							}
						}
						spec_it->getPrecursorPeak().setCharge(0);
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
			PepNovoInfile pepnovo_infile;
			PepNovoOutfile pepnovo_outfile;

			String
				logfile,
				output_filename,
				input_filename,
				pepnovo_output_filename,
				temp_data_directory,
				string_buffer,
				pepnovo_directory,
				dta_list,
				model,
				model_directory,
				modifications_filename,
				default_model,
				cleavage,
				basename,
				dta_files_common_name,
				pepnovo_modifications_filename,
				call,
				abbreviation_string;

			Int
				max_number_of_tags(0),
				min_sequence_length(0),
				max_sequence_length(0),
				num_results(0);
			
			UInt
				msms_spectra_altogether(0),
				msms_spectra_in_file(0);

			Real
				p_value(1.0),
				precursor_mass_tolerance(0.0),
				peak_mass_tolerance(0.0);

			bool
				pepnovo_in(false),
				pepnovo_out(false),
				keep_dta_files(false),
				monoisotopic(false),
				make_dtas(false);

			vector<String>
				substrings,
				substrings2,
				spectra,
				models;
			
			FileHandler fh;
			FileHandler::Type type;
			MSExperiment<Peak1D> msexperiment;
			vector< PeptideIdentification > peptide_identifications;
			ProteinIdentification protein_identification;
			ContactPerson contact_person;
			ExitCodes exit_code = EXECUTION_OK;

			// filename and tag: file has to: 1 - exist  2 - be readable  4 - writable  8 - be deleted afterwards
			map< String, UInt > files;
			UInt const
				exist(1),
				readable(2),
				writable(4),
				delete_afterwards(8);

			vector<Int> charges;

			map<String, Real> dta_filenames_and_precursor_retention_times;

			/*
				LTQ - linear quadrupole ion trap
				FT - fourier transformation
				ORBI - Orbitrap
			*/

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

			if ( getFlag_("list_models") )
			{
				model_directory = getStringOption_("model_directory");
				if ( model_directory.empty() )
				{
					writeLog_("No model directory given. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				File::absolutePath(model_directory);
				model_directory.ensureLastChar('/');
				if ( File::fileList(model_directory, String("*_config.txt"), models) )
				{
					for ( vector< String >::iterator model_it = models.begin(); model_it != models.end(); ++model_it )
					{
						model_it->erase(model_it->length() - strlen("_config.txt"));
					}
				}
				if ( models.empty() )
				{
					writeLog_("No models found in the model directory (" + model_directory + "). Aborting!");
				}
				else
				{
					cout << "Available Models:" << endl;
					for ( vector< String >::iterator model_it = models.begin(); model_it != models.end(); ++model_it )
					{
						cout << *model_it << endl;
					}
				}
				return EXECUTION_OK;
			}

			pepnovo_in = getFlag_("pepnovo_in");
			pepnovo_out = getFlag_("pepnovo_out");

			// a 'normal' pepnovo run corresponds to both pepnovo_in and pepnovo_out set
			if ( !pepnovo_in && !pepnovo_out ) pepnovo_in = pepnovo_out = true;

			logfile = getStringOption_("log");
			if ( logfile.empty() )
			{
				logfile = "temp.pepnovo.log";
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
				Int range_start(0), range_end(0);
				string_buffer.split(',', substrings);
				if ( substrings.empty() ) substrings.push_back(string_buffer);
				
				for ( vector< String >::iterator substrings_it = substrings.begin(); substrings_it != substrings.end(); )
				{
					if ( substrings_it->empty() ) substrings.erase(substrings_it);
					else
					{
						substrings_it->split('}', substrings2);
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
				for ( vector< Int >::iterator i = charges.begin(); i != --charges.end(); )
				{
					if ( (*i) == (*(i+1)) ) charges.erase(i+1);
					else ++i;
				}
				for ( vector< Int >::iterator i = charges.begin(); i != charges.end(); ++i )
				{
					if ( (*i) < 1 || (*i) > 3 )
					{
						writeLog_("Charges states allowed in [1,3] only. Aborting!");
						return ILLEGAL_PARAMETERS;
					}
				}
			}

			temp_data_directory = getStringOption_("temp_data_directory");
			if ( temp_data_directory.empty() )
			{
				writeLog_("No directory for temporary files given. Aborting!");
				return ILLEGAL_PARAMETERS;
			}
			File::absolutePath(temp_data_directory);
			temp_data_directory.ensureLastChar('/');

			string_buffer = getStringOption_("in");
			if ( string_buffer.empty() )
			{
				writeLog_("No input file specified. Aborting!");
				return ILLEGAL_PARAMETERS;
			}
			else
			{
				if ( pepnovo_in ) // if pepnovo_in is set, in are the spectra
				{
					string_buffer.split(',', spectra);
					if ( spectra.empty() ) spectra.push_back(string_buffer);
					for ( vector< String >::const_iterator spectra_it = spectra.begin(); spectra_it != spectra.end(); ++spectra_it )
					{
						files[*spectra_it] = readable;
					}
				}
				else // otherwise the pepnovo output is the input
				{
					pepnovo_output_filename = string_buffer;

					// if only pepnovo_out is set, the mz files have to be given to retrieve the retention times
					string_buffer = getStringOption_("mz_files");
					if ( string_buffer.empty() )
					{
						writeLog_("No mz files specified. Aborting!");
						return ILLEGAL_PARAMETERS;
					}
					else
					{
						string_buffer.split(',', spectra);
						if ( spectra.empty() ) spectra.push_back(string_buffer);
						for ( vector< String >::const_iterator spectra_it = spectra.begin(); spectra_it != spectra.end(); ++spectra_it )
						{
							files[*spectra_it] = readable;
						}
					}
				}
			}
			
			keep_dta_files = getFlag_("keep_dta_files");
			if ( pepnovo_in && !pepnovo_out ) keep_dta_files = true;

			contact_person.setName(getStringOption_("contact_name"));
			contact_person.setInstitution(getStringOption_("contact_institution"));
			contact_person.setContactInfo(getStringOption_("contact_info"));

			min_sequence_length = getIntOption_("min_sequence_length");
			if ( min_sequence_length < 3 || min_sequence_length > 40 )
			{
				writeLog_("min_sequence_length not in [3, 40]. Aborting!");
				return ILLEGAL_PARAMETERS;
			}

			max_sequence_length = getIntOption_("max_sequence_length");
			if ( max_sequence_length < 3 || max_sequence_length > 40 )
			{
				writeLog_("max_sequence_length not in [3, 40]. Aborting!");
				return ILLEGAL_PARAMETERS;
			}
			if ( max_sequence_length < min_sequence_length )
			{
				writeLog_("max_sequence_length is less than min_sequence_length. Aborting!");
				return ILLEGAL_PARAMETERS;
			}

			if ( pepnovo_in )
			{
				// if pepnovo_in is set (independ whether pepnovo_out is set)
				precursor_mass_tolerance = getDoubleOption_("precursor_mass_tolerance");
				if ( precursor_mass_tolerance != -1 && precursor_mass_tolerance < 0 )
				{
					writeLog_("Precursor mass tolerance < 0. Aborting!");
					return ILLEGAL_PARAMETERS;
				}

				peak_mass_tolerance = getDoubleOption_("peak_mass_tolerance");
				if ( peak_mass_tolerance != -1 && peak_mass_tolerance < 0 )
				{
					writeLog_("peak mass tolerance < 0. Aborting!");
					return ILLEGAL_PARAMETERS;
				}

				num_results = getIntOption_("num_results");
				if ( (num_results < 1) )
				{
					writeLog_("Illegal number of results (< 1). Aborting!");
					return ILLEGAL_PARAMETERS;
				}

				pepnovo_directory = getStringOption_("pepnovo_directory");
				if ( pepnovo_directory.empty() ) writeLog_("PepNovo working directory not given. Assuming PATH variable to be set accordingly.");
				else
				{
					File::absolutePath(pepnovo_directory);
					pepnovo_directory.ensureLastChar('/');
				}

				// set the protease (trypsin or not trypsin)
				cleavage = getStringOption_("cleavage");

				// maximal number of tags to use for identification
				max_number_of_tags = getIntOption_("max_number_of_tags");
				if ( max_number_of_tags != -1 && (max_number_of_tags < 0 || max_number_of_tags > 200) )
				{
					writeLog_("Maximal number of tags not in [1,200]. Aborting!");
					return ILLEGAL_PARAMETERS;
				}

				// model directory and model used if tryptic model (default) is used, pepnovo v1 is used, otherwise pepnovo v2 is used)
				model_directory = getStringOption_("model_directory");
				if ( model_directory.empty() )
				{
					writeLog_("No model directory given. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				File::absolutePath(model_directory);
				model_directory.ensureLastChar('/');
				if ( File::fileList(model_directory, String("*_config.txt"), models) )
				{
					for ( vector< String >::iterator models_it = models.begin(); models_it != models.end(); ++models_it )
					{
						models_it->erase(models_it->length() - strlen("_config.txt"));
					}
				}
				if ( models.empty() )
				{
					writeLog_("No models found in the model directory (" + model_directory + "). Aborting!");
					return INPUT_FILE_EMPTY;
				}
				else
				{
					model = getStringOption_("model");
					if ( model.empty() )
					{
						writeLog_("No model file given. Aborting!");
						return ILLEGAL_PARAMETERS;
					}
					else if ( find(models.begin(), models.end(), model) == models.end() ) // if a model was given, that's not in the model directory, abort
					{
						writeLog_("No model file given. Aborting!");
						writeLog_("Available Models:");
						for ( vector< String >::iterator models_it = models.begin(); models_it != models.end(); ++models_it )
						{
							writeLog_(*models_it);
						}
						return ILLEGAL_PARAMETERS;
					}
					else // if a correct model was given, check what maximal charge may be used
					{
						if ( !File::readable(model_directory + model + "_break_score.txt") )
						{
							return INPUT_FILE_NOT_READABLE;
						}
						else
						{
							String model_filename = model_directory + model + "_break_score.txt";
							ifstream model_file( model_filename.c_str() );
							while ( getline(model_file, string_buffer) )
							{
								if ( string_buffer.hasPrefix("#MAX_CHARGE ") )
								{
									if ( !string_buffer.empty() && (string_buffer[string_buffer.length()-1] < 33) ) string_buffer.resize(string_buffer.length()-1);
									string_buffer.trim();
									while ( charges.back() > string_buffer.substr(strlen("#MAX_CHARGE ")).toInt() ) charges.pop_back();
									break;
								}
							}
							model_file.close();
						}
					}
				}

				// the list with the names of the dta files to be analyzed
				dta_list = getStringOption_("dta_list");
				if ( dta_list.empty() )
				{
					if ( !pepnovo_out ) // if only pepnovo_in is given, the dta_list name has to be given
					{
						writeLog_("No name for dta list given (has to be given when only pepnovo_in is set). Aborting!");
						return ILLEGAL_PARAMETERS;
					}
					dta_list = temp_data_directory + "tmp.dta.list";
					files[dta_list] = (writable | delete_afterwards);
				}
				else
				{
					File::absolutePath(dta_list);
					files[dta_list] = writable;
				}

				// modifications
				string_buffer = getStringOption_("modifications");
				monoisotopic = getFlag_("use_monoisotopic_mod_mass");
				try
				{
					pepnovo_infile.handlePTMs(string_buffer, modifications_filename, monoisotopic);
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
					writeLog_(p_e.getMessage());
					return PARSE_ERROR;
				}

				if ( !pepnovo_infile.getModifications().empty() )
				{
					pepnovo_modifications_filename = model_directory + "PepNovo_PTMs.txt";
					files[pepnovo_modifications_filename] = writable;
				}
			}

			if ( pepnovo_out )
			{
				output_filename = getStringOption_("out");
				if ( output_filename.empty() )
				{
					writeLog_("No output file specified. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				File::absolutePath(output_filename);
				files[output_filename] = writable;

				// if only pepnovo out is set, -in gives the pepnovo_output_filename
				if ( pepnovo_output_filename.empty() ) pepnovo_output_filename = getStringOption_("pepnovo_output");
				if ( pepnovo_in )
				{
					if ( pepnovo_output_filename.empty() )
					{
						pepnovo_output_filename = temp_data_directory + "tmp.pepnovo.output";
						files[pepnovo_output_filename] = (writable | delete_afterwards);
					}
					else
					{
						File::absolutePath(pepnovo_output_filename);
						files[pepnovo_output_filename] = writable;
					}
				}
				else
				{
					File::absolutePath(pepnovo_output_filename);
					files[pepnovo_output_filename] = readable;
				}

				p_value = getDoubleOption_("p_value");
				if ( (p_value <= 0) || (p_value > 1) )
				{
					writeLog_("P-value not in (0, 1]. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
			}

			//-------------------------------------------------------------
			// (3) running program according to parameters
			//-------------------------------------------------------------

			// (3.1) checking accessability of files
			bool existed(false);
			UInt file_tag(0);

			for ( map< String, UInt >::const_iterator files_it = files.begin(); files_it != files.end(); ++files_it )
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

			if ( exit_code == EXECUTION_OK )
			{
				// check the Mz files, get the names for the dtas and check whether they do no already exist
				make_dtas = ( pepnovo_out && !pepnovo_in ) ? false : true; // if only pepnovo_out is set, just get the retention times
				if ( make_dtas ) writeLog_("creating dta files");
				// first get the dta names
				for ( vector< String >::iterator spectra_it = spectra.begin(); spectra_it != spectra.end(); ++spectra_it )
				{
					File::absolutePath(*spectra_it);
					type = fh.getTypeByContent(*spectra_it);
					if ( type == FileHandler::UNKNOWN )
					{
						writeLog_("Could not determine type of the file. Aborting!");
						exit_code = PARSE_ERROR;
						break;
					}
					fh.loadExperiment(*spectra_it, msexperiment, type);

					msms_spectra_in_file = MSExperiment2DTAs(msexperiment, temp_data_directory + File::basename(*spectra_it), charges, dta_filenames_and_precursor_retention_times, false);

					msms_spectra_altogether += msms_spectra_in_file;

					// if make_dtas is set, check whether one of them does already exist, if so, stop the adapter
					if ( make_dtas )
					{
						for ( map< String, Real >::const_iterator dta_names_it = dta_filenames_and_precursor_retention_times.begin(); dta_names_it != dta_filenames_and_precursor_retention_times.end(); ++dta_names_it )
						{
							string_buffer = temp_data_directory + dta_names_it->first;
							if ( File::exists(string_buffer) )
							{
								writeLog_("The file " + string_buffer + " does already exist in directory " + temp_data_directory + ". Please remove it first. Aborting!");
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
					if ( type == FileHandler::UNKNOWN )
					{
						writeLog_("Could not determine type of the file. Aborting!");
						exit_code = PARSE_ERROR;
						break;
					}
					fh.loadExperiment(*spectra_it, msexperiment, type);
					basename = File::basename(*spectra_it);
					dta_files_common_name = temp_data_directory + basename;
					msms_spectra_in_file = MSExperiment2DTAs(msexperiment, dta_files_common_name, charges, dta_filenames_and_precursor_retention_times, make_dtas);
					writeLog_(String(msms_spectra_in_file) + " MS/MS spectra in file " + *spectra_it);
				}

				if ( exit_code == EXECUTION_OK )
				{
					// make a list of all dtas
					ofstream dta_list_file(dta_list.c_str());
					if ( !dta_list_file )
					{
						exit_code = CANNOT_WRITE_OUTPUT_FILE;
						writeLog_(String("Cannot write file ")+ dta_list + ". Aborting!");
					}
					else
					{
						for ( map< String, Real >::const_iterator filenames_it = dta_filenames_and_precursor_retention_times.begin(); filenames_it != dta_filenames_and_precursor_retention_times.end(); ++filenames_it )
						{
							string_buffer = temp_data_directory + filenames_it->first;
							dta_list_file << string_buffer << endl;
						}
						dta_list_file.close();
						dta_list_file.clear();
					}
				}
			}

			if ( exit_code == EXECUTION_OK )
			{
				if ( pepnovo_in && !pepnovo_infile.getModifications().empty() )
				{
					try
					{
						abbreviation_string = pepnovo_infile.store(pepnovo_modifications_filename);
					}
					catch( Exception::UnableToCreateFile utcr_e )
					{
						writeLog_("Cannot write file " + pepnovo_modifications_filename + ". Aborting!");
						exit_code = CANNOT_WRITE_OUTPUT_FILE;
						keep_dta_files = false;
					}
				}
			}

			if ( exit_code == EXECUTION_OK )
			{
				if ( pepnovo_out ) // try to get the program version by starting the program without parameters and reading the output
				{
					// use output_filename as a temporary file
					call = pepnovo_directory;
					call.append("PepNovo_bin > " + output_filename);
					int status = system(call.c_str());
					if ( status != 256 )
					{
						pepnovo_directory.append("src/");
						call = pepnovo_directory;
						call.append("PepNovo_bin > " + output_filename);
						status = system(call.c_str());
					}

					if ( status == 256 ) pepnovo_outfile.getSearchEngineAndVersion(output_filename, protein_identification);
				}

				// how to call the program (if only pepnovo_in is set, this is returned to the user, if both flags are set, this is exectued)
				call = pepnovo_directory;
				call.append("PepNovo_bin -list " + dta_list);
				call.append(" -model " + model);
				if ( peak_mass_tolerance != -1 ) call.append(" -fragment_tolerance " + String(peak_mass_tolerance));
				if ( precursor_mass_tolerance != -1 ) call.append(" -pm_tolerance " + String(precursor_mass_tolerance));
				if ( !pepnovo_infile.getModifications().empty() ) call.append(" -PTMs " + abbreviation_string);
				if ( cleavage != "Trypsin" ) call.append(" -digest NON_SPECIFIC ");
				call.append(" -num_solutions " + String(num_results));
				call.append(" -min_length " + String(min_sequence_length));
				call.append(" -max_length " + String(max_sequence_length));
				call.append(" -model_dir " + model_directory);
				call.append(" -denovo_mode ");
				call.append(" > " + pepnovo_output_filename);

				// if only pepnovo_in is set, output the call of pepnovo
				if ( pepnovo_in )
				{
					if ( pepnovo_out )
					{
						// running the program
						writeLog_("System call: " + call);
						int status = system(call.c_str());

						if ( status != 0 ) exit_code = EXTERNAL_PROGRAM_ERROR;
					}
					else
					{
						writeLog_("Use this line to call PepNovo: ");
						writeLog_(call);
					}
				}
			}

			if ( exit_code == EXECUTION_OK && pepnovo_out )
			{
				// set the parameters
				ProteinIdentification::SearchParameters sp;
				if ( monoisotopic ) sp.mass_type = ProteinIdentification::MONOISOTOPIC;
				else sp.mass_type = ProteinIdentification::AVERAGE;
				for ( vector< Int >::const_iterator charges_it = charges.begin(); charges_it != charges.end(); ++charges_it )
				{
					if ( *charges_it > 0 ) sp.charges.append("+");
					sp.charges.append(String(*charges_it));
				}
				if ( cleavage == "Trypsin" ) sp.enzyme = ProteinIdentification::TRYPSIN;
				else if ( cleavage == "No_Enzyme" ) sp.enzyme = ProteinIdentification::NO_ENZYME;
				else sp.enzyme = ProteinIdentification::UNKNOWN_ENZYME;
				sp.peak_mass_tolerance = peak_mass_tolerance;
				sp.precursor_tolerance = precursor_mass_tolerance;
				protein_identification.setSearchParameters(sp);

				pepnovo_outfile.load(pepnovo_output_filename, peptide_identifications, protein_identification, p_value, dta_filenames_and_precursor_retention_times);

				vector< ProteinIdentification > identifications;
				identifications.push_back(protein_identification);
				
				IdXMLFile().store(output_filename, identifications, peptide_identifications);
			}
			
			if ( exit_code == EXTERNAL_PROGRAM_ERROR )
			{
				writeLog_("PepNovo problem. Aborting! (Details can be seen in the logfile: \"" + logfile + "\")");
				files[logfile] = readable;
			}
			
			// deleting all temporary files
			writeLog_("removing temporary files");
			for ( map< String, UInt >::const_iterator files_it = files.begin(); files_it != files.end(); ++files_it )
			{
				if ( files_it->second & delete_afterwards ) remove(files_it->first.c_str());
			}
			// remove all dtas
			if ( !keep_dta_files )
			{
				writeLog_("removing dta files");
				for ( map< String, Real >::const_iterator dta_names_it = dta_filenames_and_precursor_retention_times.begin(); dta_names_it != dta_filenames_and_precursor_retention_times.end(); ++dta_names_it )
				{
					string_buffer = temp_data_directory + dta_names_it->first;
					if ( !File::remove(string_buffer) ) writeLog_("'" + string_buffer + "' could not be removed!");
				}
			}
			
			return exit_code;
		}
};

//@endcond



int main( int argc, const char** argv )
{
	TOPPPepNovoAdapter tool;

	return tool.main(argc,argv);
}
