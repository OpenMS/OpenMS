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
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/AnalysisXMLFile.h>
#include <OpenMS/FORMAT/MascotXMLFile.h>
#include <OpenMS/FORMAT/MascotInfile.h>
#include <OpenMS/FORMAT/MascotOutfile.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Date.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/ANALYSIS/ID/MSExperimentAnnotator.h>
#include <OpenMS/FILTERING/ID/IDFilter.h>
#include <OpenMS/METADATA/ContactPerson.h>
#include "TOPPBase.h"

#include <qfileinfo.h>
#include <qfile.h>

#include <map>
#include <iostream>
#include <fstream>
#include <string>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

// @cond TOPPCLASSES 

/**
	@page MascotAdapter MascotAdapter
	
	@brief Identifies peptides in MS/MS spectra via Mascot.
	
	This wrapper component serves for getting peptide identifications
	for MS/MS spectra. The wrapper can be executed in three different
	modes:
	<ol>	
				<li>
				The whole process of identification via Mascot is executed. 
				Inputfile is a mzData file containing the MS/MS spectra
			 	for which the identifications are to be found. The results
			 	are written as a analysisXML output file. This mode is selected
			 	by default.
			 	</li>
				
				<li>
				Only the first part of the identification process is performed.
				This means that the MS/MS data is transformed into Mascot
				Generic Format (mgf) which can be used directly with Mascot.
				Being in the cgi directory of the Mascot directory calling a Mascot
				process should look like the following:				
				
				<ul>	
					<li>
						@code ./nph-mascot.exe 1 -commandline -f outputfilename < inputfilename @endcode
					</li>
				</ul>
				Consult your Mascot reference manual for further details.
				
				This mode is selected by the <b>-mascot_in</b> option in the command line.
				</li>
				
				<li>
				Only the second part of the identification process is performed.
				This means that the outputfile of the Mascot server is
				translated into analysisXML.
				
				This mode is selected by the <b>-mascot_out</b> option in the command line.
				</li>
	</ol>

	<br>			
	If your Mascot server is installed on the same computer as the 
	TOPP components the MascotAdapter can be executed in mode 1. 
	Otherwise the Mascot engine has to be executed manually assisted
	by mode 2 and mode 3. The identification steps then look like:
	
	<ul>
		<li>
			# execute MascotAdapter in mode 2
			@code ./MascotAdapter -in mzDataFile -out mascotGenericFormatFile -mascot_in @endcode	
		</li>
		<li>
			copy mascotGenericFormatFile to your Mascot server
		</li>
		<li>
			# call your Mascot server process:
			@code ./nph-mascot.exe 1 -commandline -f mascotOutFile < mascotGenericFormatFile @endcode
		</li>
		<li>
			copy mascotOutFile to the server on which the TOPP components are installed
		</li>
		<li>
			# execute MascotAdapter in mode 3			
			@code ./MascotAdapter -in mascotOutFile -out analysisXMLFile -mascot_out @endcode
		</li>
	</ul>

	<p>
	For mode 1 you have to specify the directory in which the Mascot
	server is installed. This is done by setting the option <b>mascot_dir</b> 
	in the ini file. Furthermore you have to specify a folder in which
	the user has write permissions. This is done by setting the option 
	<b>temp_data_directory</b> in the ini file. 
	Two temporary files will be created in this directory during execution 
	but deleted at the end of execution.
	<br>
	
	You can specify the Mascot parameters <b>precursor_mass_tolerance</b> 
	(the peptide mass tolerance), <b>peak_mass_tolerance</b> (the MS/MS tolerance), 
	<b>taxonomy</b> (restriction to a certain subset of the database), <b>modifications</b>, 
	<b>variable_modifications</b>, <b>charges</b> (the possible charge variants), 
	<b>db</b> (database where the peptides are searched in), <b>hits</b> (number of hits), 
	<b>cleavage</b> (the cleavage enzyme), <b>missed_cleavages</b> (number of missed cleavages) 
	and <b>mass_type</b> (Monoisotopic or Average) via the ini file.
	
	<br>			
	Known problems with Mascot server execution:
	<ul>
		<li>	
		getting error message:
		"FATAL_ERROR: M00327
		 The ms-monitor daemon/service is not running, please start it."
		</li>
		
		<li>
		Possible explanations:
		</li>
		<ul>
			<li>
			Your ms-monitor is really not running => consult your Mascot
																							 reference manual for
																							 details about starting 
																							 the Mascot server.
			</li>
			<li>
			(Suppose you have Mascot installed in directory mascot.)
			mascot/data/mascot.control is not writable for the current user.
			This has to be changed. Otherwise you will not be able to 
			use the Mascot server via the shell and receive the above error
			message.<br>			
			=> Change write permissions of the file mascot/data/mascot.control
				 such that the current user has write permissions to it.
			</li>
		</ul>
	</ul>		
	
	@todo Fix --help-opt output (Nico)
	
	@ingroup TOPP
*/
class TOPPMascotAdapter
	: public TOPPBase
{
	public:
		TOPPMascotAdapter()
			: TOPPBase("MascotAdapter")
		{
			
		}
	
	protected:
		void printToolUsage_()
		{
			cerr << endl
		       << tool_name_ << " -- annotates MS/MS spectra using Mascot" << endl
		       << endl
		       << "Usage:" << endl
					 << " " << tool_name_ << " [options]" << endl
					 << endl
					 << "Options are:" << endl
					 << "  -in <file>           input file in mzData/Mascot resultsfile " 
					 << "(default read from INI file)" << endl
					 << "  -out <file>          output file in analysisXML/Mascot generic format "
					 << "(default read from INI file)" << endl
					 << "  -mascot_in           if this flag is set the MascotAdapter will read in "
					 << "mzData and write Mascot generic format" << endl
					 << "  -mascot_out          if this flag is set the MascotAdapter will read in "
					 << "a Mascot resultsfile and write analysisXML." << endl
					 << "  -instr               the instrument that was used to measure the spectra (default read from INI file)" << endl
					 <<	"  -prcr_m_tol          the precursor mass tolerance (default read from INI file)" << endl
					 << "  -pk_m_tol            the peak mass tolerance (default read from INI file)" << endl
					 << "  -tax                 the taxonomy (default read from INI file)" << endl
					 << "  -mods                the modifications i.e. Carboxymethyl (C) (default read from INI file)" << endl
					 << "  -vmods               the variable modifications i.e. Carboxymethyl (C) (default read from INI file)" << endl
					 << "  -charges             the different charge states separated by comma ( (default read from INI file)" << endl
					 << "  -mascot_directory    the directory in which mascot is located" << endl
					 << "  -temp_data_directory a directory in which some temporary files can be stored" << endl
					 << endl ;
		}

		void printToolHelpOpt_()
		{
			cerr << endl
		       << tool_name_ << endl
		       << endl
		       << "INI options:" << endl
					 << "  in                        input file" << endl
					 << "  out                       output file" << endl
					 << "  mascot_in                 if this flag is set the MascotAdapter will read in "
					 << "mzData and write Mascot generic format" << endl
					 << "  mascot_out                if this flag is set the MascotAdapter will read in "
					 << "a Mascot resultsfile and write analysisXML." << endl
					 << "  instrument                the instrument that was used to measure the spectra" << endl
					 <<	"  precursor_mass_tolerance  the precursor mass tolerance" << endl
					 << "  peak_mass_tolerance       the peak mass tolerance" << endl
					 << "  taxonomy                  the taxonomy" << endl
					 << "  modifications             the modifications i.e. Carboxymethyl (C)" << endl
					 << "  variable_modifications    the variable modifications i.e. Carboxymethyl (C)" << endl
					 << "  charges                   the different charge states separated by comma" << endl
					 << "  mascot_directory          the directory in which mascot is located" << endl
					 << "  temp_data_directory       a directory in which some temporary files can be stored" << endl
					 << endl
					 << "INI File example section:" << endl
					 << "  <ITEM name=\"in\" value=\"input.mzData\" type=\"string\"/>" << endl
					 << "  <ITEM name=\"out\" value=\"output.analysisXML\" type=\"string\"/>" << endl
					 << "  <ITEM name=\"instrument\" value=\"ESI-TRAP\" type=\"string\"/>" << endl
           << "  <ITEM name=\"precursor_mass_tolerance\" value=\"1.3\" type=\"float\"/>" << endl
           << "  <ITEM name=\"peak_mass_tolerance\" value=\"0.3\" type=\"float\"/>" << endl
           << "  <ITEM name=\"taxonomy\" value=\". . . . . . Chordata (vertebrates and relatives)\" type=\"string\"/>" << endl
           << "  <ITEM name=\"modifications\" value=\"Carboxymethyl (C)\" type=\"string\"/>" << endl
           << "  <ITEM name=\"charges\" value=\"1+,2+,3+\" type=\"string\"/>" << endl
           << "  <ITEM name=\"db\" value=\"MSDB\" type=\"string\"/>" << endl
           << "  <ITEM name=\"hits\" value=\"AUTO\" type=\"string\"/>" << endl
           << "  <ITEM name=\"cleavage\" value=\"Trypsin\" type=\"string\"/>" << endl
           << "  <ITEM name=\"missed_cleavages\" value=\"1\" type=\"UnsignedInt\"/>" << endl
           << "  <ITEM name=\"mass_type\" value=\"Monoisotopic\" type=\"string\"/>" << endl
		       << "  <ITEM name=\"mascot_directory\" value=\"/local/mascot/\" type=\"string\"/>" << endl
		       << "  <ITEM name=\"temp_data_directory\" value=\"/local/mascot/tmp/\" type=\"string\"/>" << endl;
 
		}		

		void setOptionsAndFlags_()
		{
			options_["-out"] = "out";
			options_["-in"] = "in";
			options_["-additional_in"] = "additional_in";
			options_["-instr"] = "instrument";
			options_["-prcr_m_tol"] = "precursor_mass_tolerance";
			options_["-pk_m_tol"] = "peak_mass_tolerance";
			options_["-tax"] = "taxonomy";
			options_["-mods"] = "modifications";
			options_["-vmods"] = "variable_modifications";
			options_["-charges"] = "charges";
			options_["-db"] = "db";
			options_["-hits"] = "hits";
			options_["-cleavage"] = "cleavage";
			options_["-missed_cleavages"] = "missed_cleavages";
			options_["-mass_type"] = "mass_type";
			options_["-mascot_in"] = "mascot_in";
			options_["-mascot_out"] = "mascot_out";
			options_["-mascot_directory"] = "mascot_directory";
			options_["-temp_data_directory"] = "temp_data_directory";
		}

		ExitCodes main_(int , char**)
		{
			// instance specific location of settings in INI file (e.g. 'TOPP_Skeleton:1:')
			String ini_location;
			// path to the log file
			String logfile = "mascot.log";
			// log filestream (as long as the real logfile is not determined yet)
			ofstream log;		
			String inputfile_name;
			String outputfile_name;
			String mascot_infile_name = "tmp.mascot_in";
			String mascot_outfile_name = "tmp_mascot_in.out";
			String mascot_output_name = "tmp_mascot.output";
			String mascot_cgi_dir;
			String mascot_data_dir;
			String call;
			String instrument;
			String taxonomy;
			String temp_string;
			String mascotXML_file_name = "";
			MzDataFile mzdata_infile;
			MSExperiment< DPeak<1> > experiment;
			MSExperimentAnnotator annotator;
			IDFilter filter;
			AnalysisXMLFile analysisXML_file;
			MascotXMLFile mascotXML_file;
			MascotInfile* mascot_infile;
			MascotOutfile* mascot_outfile;
			ContactPerson contact_person;
			vector<String> mods;
			vector<String> variable_mods;
			ProteinIdentification protein_identification;
			vector<Identification> identifications;
			vector<Real> precursor_retention_times;
			vector<Real> precursor_mz_values;
			vector<SignedInt> charges;
			vector<String> parts;
			DoubleReal precursor_mass_tolerance = 2;
			DoubleReal peak_mass_tolerance = 1;
			String temp_charge;
			QFileInfo file_info;
			QFile file;
			string db = "MSDB";
			string hits = "20";
			string cleavage = "Trypsin";
			UnsignedInt missed_cleavages = 1;
			string mass_type = "Monoisotopic";
			int status = 0;
			bool mascot_in = false;
			bool mascot_out = false;
			DateTime date_time;
			String date_time_string;
			String time_string;
			
			date_time.now();
			date_time.get(date_time_string);
			date_time_string.split(' ', parts);
			
			mascot_infile_name = parts[0] + "_" + parts[1] + "_" + mascot_infile_name;
			mascot_outfile_name	= parts[0] + "_" + parts[1] + "_" + mascot_outfile_name;
			mascot_output_name = parts[0] + "_" + parts[1] + "_" + mascot_output_name;
			parts.clear();
				
			//-------------------------------------------------------------
			// parsing parameters
			//-------------------------------------------------------------
						
			inputfile_name = getParamAsString_("in");			
			writeDebug_(String("Input file: ") + inputfile_name, 1);
			if (inputfile_name == "")
			{
				writeLog_("No input file specified. Aborting!");
				cout << "No input file specified. Aborting!" << endl;
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}
	
			outputfile_name = getParamAsString_("out");
			writeDebug_(String("Output file: ") + outputfile_name, 1);
			if (outputfile_name == "")
			{
				writeLog_("No output file specified. Aborting!");
				cout << "No output file specified. Aborting!" << endl;
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}				

			mascotXML_file_name = getParamAsString_("additional_in");
			writeDebug_(String("Additional input file: ") + mascotXML_file_name, 1);

			if (!getParam_("log").isEmpty())
			{
				logfile = getParamAsString_("log");
			}

			if (getParamAsString_("mascot_in", "false") != "false")
			{
				mascot_in = true;
			}
			
			if (getParamAsString_("mascot_out", "false") != "false")
			{
				mascot_out = true;
				if (mascot_in)
				{
					writeLog_("Both Mascot flags set. Aborting!");
					cout << "Both Mascot flags set. Aborting! Only one of the two "
						<< "flags [-mascot_in|-mascot_out] can be set" << endl;
					return ILLEGAL_PARAMETERS;
				}				
			}
			else
			{		
				instrument = getParamAsString_("instrument", "Default");
				writeDebug_(String("Instrument: ") + instrument, 1);
				precursor_mass_tolerance = getParamAsString_("precursor_mass_tolerance", "2.0f").toDouble();
				writeDebug_(String("Precursor mass tolerance: ") + 
					String(precursor_mass_tolerance), 1);


				peak_mass_tolerance = getParamAsString_("peak_mass_tolerance", "1.0f").toDouble();
				writeDebug_(String("Peak mass tolerance: ") + String(peak_mass_tolerance), 1);

				taxonomy = getParamAsString_("taxonomy", "All entries");
				writeDebug_(String("Taxonomy: ") + taxonomy, 1);

				/// fixed modifications
				temp_string = getParamAsString_("modifications");
				temp_string.split(',', mods);
				if (mods.size() == 0 && temp_string != "")
				{
					mods.push_back(temp_string);
				}
				writeDebug_(String("Modifications: ") + temp_string, 1);

				/// variable modifications
				temp_string = getParamAsString_("variable_modifications");
				temp_string.split(',', variable_mods);
				if (variable_mods.size() == 0 && temp_string != "")
				{
					mods.push_back(temp_string);
				}				
				writeDebug_(String("Variable modifications: ") + temp_string, 1);							

				///charges
				temp_string = getParamAsString_("charges");
				temp_string.split(',', parts);
				if (parts.size() == 0 && temp_string != "")
				{
					temp_charge = temp_string;
					if (temp_charge[temp_charge.size() - 1] == '-'
							|| temp_charge[0] == '-')
					{
						charges.push_back(-1 * (temp_charge.toInt()));
					}
					else
					{
						charges.push_back(temp_charge.toInt());						
					}
				}									
				else if (temp_string != "")
				{
					for(UnsignedInt i = 0; i < parts.size(); i++)
					{
						temp_charge = parts[i];
						if (temp_charge[temp_charge.size() - 1] == '-'
								|| temp_charge[0] == '-')
						{
							charges.push_back(-1 * (temp_charge.toInt()));
						}
						else
						{
							charges.push_back(temp_charge.toInt());						
						}
					}
				}
				if (charges.size() == 0)
				{
					writeLog_(String("No charge states specified for ") + 
										String("Mascot search. Aborting!"));
					cout << "No charge states specified for "
						<< "Mascot search. Aborting!" << endl;
					return ILLEGAL_PARAMETERS;			
				}
				writeDebug_(String("Charges: ") + temp_string, 1);
				  				
				db = getParamAsString_("db", "MSDB");
				writeDebug_(String("Database: ") + db, 1);
				
				hits = getParamAsString_("hits", "AUTO");
				writeDebug_(String("Hits: ") + hits, 1);
								
				cleavage = getParamAsString_("cleavage", "Trypsin");
				writeDebug_(String("Cleavage: ") + cleavage, 1);
				
				missed_cleavages = (UnsignedInt) getParamAsString_("missed_cleavages", "0").toInt();
				writeDebug_(String("Number of allowed missed cleavages: ") + String(missed_cleavages), 1);
				
				mass_type = getParamAsString_("mass_type", "Monoisotopic");
				writeDebug_(String("Precursor mass type: ") + mass_type, 1);
	
			}
			if (mascot_in)
			{
				mascot_infile_name = outputfile_name;
				writeDebug_(String("Mascot flag: ") + 
					String("mascot_in (reads in MzData writes Mascot generic format)"), 1);
			}
			else if (mascot_out)
			{
				mascot_outfile_name = inputfile_name;
				writeDebug_(String("Mascot flag: ") + 
					String("mascot_out (reads in Mascot results file writes analysisXML file)"), 1);
			}
			else
			{
				writeDebug_(String("No Mascot flag set: ") + 
					String("reads in MzData writes analysisXML file"), 1);
			}
			if (!mascot_in && !mascot_out)
			{
				
				mascot_cgi_dir = getParamAsString_("mascot_directory");
				if (mascot_cgi_dir == "")
				{
					writeLog_(String("No Mascot directory specified. Aborting!"));
					cout << "No Mascot directory specified."
						<< " Aborting!" << endl;
					return ILLEGAL_PARAMETERS;
				}
				writeDebug_(String("Mascot directory: ") + mascot_cgi_dir, 1);
				mascot_cgi_dir += "/cgi/";

				mascot_data_dir = getParamAsString_("temp_data_directory");

				if (mascot_data_dir == "")
				{
					writeLog_("No temp directory specified. Aborting!");
					cout << "No temp directory specified."
						<< " Aborting!" << endl;
					return ILLEGAL_PARAMETERS;
				}
				
				writeDebug_(String("Temp directory: ") + mascot_data_dir, 1);

				file.setName(mascot_data_dir + "/" + mascot_outfile_name);
				file.open( IO_WriteOnly );
				if (!file.isWritable())
				{
					writeLog_(String(" Could not write in temp data directory: ")
						+ mascot_data_dir + "/" + mascot_outfile_name
						+ String(" Aborting!"));
					cout << "Could not write in temp data directory: "
						<< mascot_data_dir
						<< " Aborting!" 
						<< endl;
					file.close();				
					return ILLEGAL_PARAMETERS;
				}
				file.close();
				mascotXML_file_name = mascot_data_dir + "/" + mascot_outfile_name + ".mascotXML";				
			}

			contact_person.setName(getParamAsString_("contactName", "unknown"));
			writeDebug_(String("Contact name: ") + contact_person.getName(), 1);

			contact_person.setInstitution(getParamAsString_("contactInstitution", "unknown"));
			writeDebug_(String("Contact institution: ") + contact_person.getInstitution(), 1);
			
			contact_person.setContactInfo(getParamAsString_("contactInfo"));
			writeDebug_(String("Contact info: ") + contact_person.getContactInfo(), 1);
					
			//-------------------------------------------------------------
			// testing whether input and output files are accessible
			//-------------------------------------------------------------
	
			file_info.setFile(inputfile_name.c_str());
			if (!file_info.exists())
			{
				throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, inputfile_name);
			}
			if (!file_info.isReadable())
			{
				throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, inputfile_name);			
			}
	    if (file_info.size() == 0)
	    {
	      throw Exception::FileEmpty(__FILE__, __LINE__, __PRETTY_FUNCTION__, inputfile_name);
	    }		
			file.setName(outputfile_name.c_str());
			file.open( IO_WriteOnly );
			if (!file.isWritable())
			{
				throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, outputfile_name);
			}
			file.close();				
	
			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------
	
			if(!mascot_out)
			{
				mzdata_infile.load(inputfile_name, experiment);
					
			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------
			
			
				mascot_infile = new MascotInfile();
			
				mascot_infile->setInstrument(instrument);
				mascot_infile->setPrecursorMassTolerance(precursor_mass_tolerance);
				mascot_infile->setPeakMassTolerance(peak_mass_tolerance);
				if (mods.size() > 0)
				{
					mascot_infile->setModifications(mods);
				}
				if (variable_mods.size() > 0)
				{
					mascot_infile->setVariableModifications(variable_mods);
				}
				mascot_infile->setTaxonomy(taxonomy);
				mascot_infile->setDB(db);
				mascot_infile->setHits(hits);
				mascot_infile->setCleavage(cleavage);
				mascot_infile->setMissedCleavages(missed_cleavages);
				mascot_infile->setMassType(mass_type);
				mascot_infile->setCharges(charges);
				if (!mascot_in)
				{
					mascot_infile->store(mascot_data_dir + "/" + mascot_infile_name, 
															 experiment, 
															 "OpenMS search");
					file_info.setFile(logfile.c_str());
					writeDebug_("The Mascot process created the following output:", 0);
					/// calling the Mascot process
					call = "cd " + mascot_cgi_dir + "; ./nph-mascot.exe 1 -commandline -f " + 
						mascot_data_dir + "/" + mascot_outfile_name + " < " + 
						mascot_data_dir + "/" + mascot_infile_name + 
						" >> " + String(file_info.absFilePath().ascii()) + ";"
						+ "./export_dat.pl do_export=1 export_format=XML file=" + mascot_data_dir + 
						"/" + mascot_outfile_name + " _showsubset=1 show_same_sets=1 show_unassigned=1 " + 
						"prot_score=1 pep_exp_z=1 pep_score=1 pep_homol=1 pep_ident=1 pep_seq=1 " + 
						"show_header=1 > " + mascotXML_file_name + ";";
					status = system(call.c_str());
					if (status != 0)
					{
						cout << "Mascot server problem. Aborting! (Details can be seen " 
						<< " in the logfile: \"" << logfile << "\")" << endl;
						writeLog_("Mascot server problem. Aborting!");
						call = "rm " + mascot_data_dir + "/" 
										+ mascot_infile_name + "; rm " + mascotXML_file_name + ";";
						system(call.c_str());
						return EXTERNAL_PROGRAM_ERROR;						
					}
					
				} // from if(!mascot_in)
				else
				{
					mascot_infile->store(mascot_infile_name, 
															 experiment, 
															 "OpenMS search");		
				}
		
			} // from if(!mascot_out)
			if (!mascot_in)
			{
	
				/// Reading in the Mascot outfile
				/// Since Mascot does not store the retention times in the mascot xml format
				/// we need also the old out-file to retrieve them. All other data is loaded 
				/// from the mascot xml file (if possible).
				mascot_outfile = new MascotOutfile();

				vector<Real> wrong_retention_times;				
				if (!mascot_out)
				{
//					mascot_outfile->load(mascot_data_dir + "/" + mascot_outfile_name,
//															identifications,
//															precursor_retention_times,
//															precursor_mz_values);
					mascotXML_file.load(mascotXML_file_name,
															&protein_identification,
															&identifications,
															&precursor_retention_times,
															&precursor_mz_values);			
				}
				else
				{
					if (mascotXML_file_name != "")
					{
						mascotXML_file.load(mascotXML_file_name,
															&protein_identification,
															&identifications,
															&precursor_retention_times,
															&precursor_mz_values);			
						
					}
					else
					{
						mascot_outfile->load(mascot_outfile_name,
																identifications,
																precursor_retention_times,
																precursor_mz_values);
						
					}
							
				}
				wrong_retention_times.clear();
			
			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------
				vector<ProteinIdentification> protein_identifications;
				protein_identifications.push_back(protein_identification);
				analysisXML_file.store(outputfile_name,
															 protein_identifications, 
													 		 identifications, 
													 		 precursor_retention_times, 
													 		 precursor_mz_values,
													 		 contact_person);
													 		 												 		 
				/// Deletion of temporary Mascot files
				if (!mascot_out)
				{
					call = "rm " + mascot_data_dir + "/" + mascot_infile_name + ";"
						+ "rm " + mascot_data_dir + "/" + mascot_outfile_name + ";"
						+ "rm " + mascotXML_file_name + ";";
					system(call.c_str());
				}
			
			} // from if(!mascot_in)
			return OK;	
		}
};


int main( int argc, char ** argv )
{
	TOPPMascotAdapter tool;

	return tool.main(argc,argv);
}

// @endcond
