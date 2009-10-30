// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

#include <QtCore/QFile>
#include <QtCore/QDir>

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

	@experimental This tool has not been tested thoroughly and might behave not as expected!

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
			registerInputFile_("in", "<file>", "", "input file ");
			setValidFormats_("in",StringList::create("mzXML"));

			registerOutputFile_("out", "<file>", "", "output file ");
			setValidFormats_("out",StringList::create("idXML"));

			registerInputFile_("pepnovo_executable","<file>", "", "The \"PepNovo\" executable of the PepNovo installation", true);
			registerStringOption_("temp_data_directory", "<dir>", "", "Directory were temporary data can be stored. If not set the directory were startet is used.", true);
			registerFlag_("correct_pm", "find optimal precursor mass and charge values.");
			registerFlag_("use_spectrum_charge", "do not correct charge");
			registerFlag_("use_spectrum_mz", "do not correct the precursor m/z value that appears in the file.");
			registerFlag_("no_quality_filter", "do not remove low quality spectra.");
			registerDoubleOption_("fragment_tolerance", "<Float>", -1.0, "the fragment tolerance (between 0 and 0.75 Da. Set to -1.0 to use model's default setting)", false, false);
			registerDoubleOption_("pm_tolerance", "<Float>", -1.0, "the precursor mass tolerance (between 0 and 5.0 Da. Set to -1.0 to use model's default setting)", false, false);
			registerStringOption_("model_directory", "<file>", " ", "name of the directory where the model files are kept.",true);
			registerStringOption_("model", "<file>", "CID_IT_TRYP", "name of the model that should be used", false);

			registerStringOption_("digest", "", "TRYPSIN", "enzyme used for digestion (default TRYPSIN)", false);
			setValidStrings_("digest", StringList::create("TRYPSIN,NON_SPECIFIC"));

			registerIntOption_("tag_length", "<num>", -1, "returns peptide sequence of the specified length (only lengths 3-6 are allowed)", false);

			registerIntOption_("num_solutions", "<num>", 20, "number of solutions to be computed", false);
			setMinInt_("num_solutions",1);
			setMaxInt_("num_solutions",2000);

			std::vector<String>all_possible_modifications;
			ModificationsDB::getInstance()->getAllSearchModifications(all_possible_modifications);
			registerStringList_("fixed_modifications", "<mod1,mod2,...>", StringList::create(""), "list of fixed modifications", false);
			setValidStrings_("fixed_modifications", all_possible_modifications);
			registerStringList_("variable_modifications", "<mod1,mod2,...>", StringList::create(""), "list of fixed modifications", false);
			setValidStrings_("variable_modifications", all_possible_modifications);
		}


		ExitCodes main_(int , const char**)
		{

			// path to the log file
			String logfile(getStringOption_("log"));
			String pepnovo_executable(getStringOption_("pepnovo_executable"));

			//ofstream log;
			String inputfile_name, outputfile_name, model_directory;
			PeakMap exp;

			inputfile_name = getStringOption_("in");
			writeDebug_(String("Input file: ") + inputfile_name, 1);
			if (inputfile_name == "")
			{
				writeLog_("No input file specified. Aborting!");
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}

			outputfile_name = getStringOption_("out");
			writeDebug_(String("Output file: ") + outputfile_name, 1);
			if (outputfile_name == "")
			{
				writeLog_("No output file specified. Aborting!");
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}

			model_directory = getStringOption_("model_directory");
			writeDebug_(String("model directory: ") + model_directory, 1);
			if (model_directory == "")
			{
				writeLog_("No model directory specified. Aborting!");
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}

			String model_name = getStringOption_("model");
			writeDebug_(String("model directory: ") + model_name, 1);
			if (model_name == "")
			{
				writeLog_("No model specified. Aborting!");
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}

			DoubleReal fragment_tolerance = getDoubleOption_("fragment_tolerance");
			if(fragment_tolerance!=-1.0 && (fragment_tolerance<0 || fragment_tolerance>0.75))
			{
				writeLog_("Invalid fragment tolerance");
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}

			DoubleReal pm_tolerance = getDoubleOption_("pm_tolerance");
			if(pm_tolerance!=-1.0 && (pm_tolerance<0.0 || pm_tolerance>5.0))
			{
				writeLog_("Invalid fragment tolerance");
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}

			Int tag_length = getIntOption_("tag_length");
			if( tag_length!=-1 && (tag_length<3 || tag_length>6))
			{
				writeLog_("Invalid fragment tolerance");
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}

			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------

			// only load msLevel 2
			MzXMLFile mzdata_infile;
			mzdata_infile.getOptions().addMSLevel(2);
			mzdata_infile.setLogType(log_type_);
			mzdata_infile.load(inputfile_name, exp);

			// we need to replace the native id with a simple numbering schema, to be able to
			// map the IDs back to the spectra (RT, and MZ infomration)
			std::map<String, Real>id_to_rt;
			Size native_id(1);
			for (PeakMap::Iterator it = exp.begin(); it != exp.end(); ++it)
			{
				id_to_rt[native_id]=it->getRT();
				it->setNativeID(native_id++);
			}

			logfile = getStringOption_("log");
			
			QString temp_data_directory = getStringOption_("temp_data_directory").c_str();
			if ( temp_data_directory=="")
			{
				writeLog_("No directory for temporary files given. Aborting!");
				return ILLEGAL_PARAMETERS;
			}

			QDir qdir_temp(temp_data_directory);
			QDir qdir_models_source(model_directory.c_str());

			if(!qdir_temp.exists())
			{
				writeLog_("The temporary directory does not exist");
				return INPUT_FILE_NOT_FOUND;
			}
			if(!qdir_temp.exists())
			{
				writeLog_("The model directory does not exist");
				return INPUT_FILE_NOT_FOUND;
			}

			try{

			  //temporary File to store PepNovo output
			  String temp_pepnovo_outfile = qdir_temp.absoluteFilePath("tmp_pepnovo_out.txt");
			  String tmp_models_dir=qdir_temp.absoluteFilePath("Models");

				if(qdir_temp.cd("Models"))
				{
					writeLog_("The temporary directory already contains \"Model\" Folder. Please delete it and re-run. Aborting!");
					return CANNOT_WRITE_OUTPUT_FILE;
				}
				else
				{
					qdir_temp.mkdir("Models");
					qdir_temp.cd("Models");
				}

				//copy the Models folder of OpenMS into the temp_data_directory
				QStringList pepnovo_files = qdir_models_source.entryList(QDir::Dirs | QDir::Files|QDir::NoDotAndDotDot);
				if(pepnovo_files.empty())
				{
					writeLog_("The \"Model\" directory does not contain model files. Aborting!");
					return INPUT_FILE_NOT_FOUND;
				}

				for(QStringList::ConstIterator file_it=pepnovo_files.begin(); file_it!=pepnovo_files.end(); ++file_it)
				{
				  if(qdir_models_source.cd(*file_it))
				  {
				    //std::cout<<"if: "<<file_it->toStdString()<<std::endl;
				    qdir_temp.mkdir(*file_it);
				    qdir_temp.cd(*file_it);
				    QStringList subdir_files = qdir_models_source.entryList(QDir::Dirs | QDir::Files|QDir::NoDotAndDotDot);
				    for(QStringList::ConstIterator subdir_file_it=subdir_files.begin(); subdir_file_it!=subdir_files.end(); ++subdir_file_it)
            {
				      //std::cout<<subdir_file_it->toStdString()<<std::endl;
				      QFile::copy(qdir_models_source.filePath(*subdir_file_it), qdir_temp.filePath(*subdir_file_it));
            }
				    qdir_temp.cdUp();
				    qdir_models_source.cdUp();
				  }
          else
          {
            //std::cout<<"else: "<<file_it->toStdString()<<std::endl;
            QFile::copy(qdir_models_source.filePath(*file_it), qdir_temp.filePath(*file_it));
          }
				}

				//generate PTM File and store in temp directory
				PepNovoInfile p_novo_infile;
				String ptm_command;
				if(!getStringList_("fixed_modifications").empty() || !getStringList_("variable_modifications").empty())
				{
					p_novo_infile.setModifications(getStringList_("fixed_modifications"), getStringList_("variable_modifications"));
					p_novo_infile.store(qdir_temp.filePath("PepNovo_PTMs.txt"));
					pepnovo_files.append("PepNovo_PTMs.txt");
					std::map<String, String>mods_and_keys;
					p_novo_infile.getModifications(mods_and_keys);

					for(std::map<String, String>::const_iterator key_it=mods_and_keys.begin(); key_it!=mods_and_keys.end();++key_it)
					{
						if(ptm_command!="")
						{
							ptm_command+=":";
						}
						ptm_command+= key_it->second;
					}
				}

				//-------------------------------------------------------------
				// (3) running program according to parameters
				//-------------------------------------------------------------

				String call;
				call = pepnovo_executable;
				//call.append("PepNovo_bin");
				call.append(" -file " + inputfile_name);
				call.append(" -model " + model_name);
				if (pm_tolerance != -1 ) call.append(" -pm_tolerance " + String(pm_tolerance));
				if (fragment_tolerance != -1 ) call.append(" -fragment_tolerance " + String(fragment_tolerance));
				if (!ptm_command.empty()) call.append(" -PTMs " + ptm_command);
				call.append(" -digest "+ getStringOption_("digest"));
				call.append(" -num_solutions " + String(getIntOption_("num_solutions")));
				if(tag_length!=-1)call.append(" -tag_length " + String(tag_length));
				call.append(" -model_dir " + tmp_models_dir);
				call.append(String(" > ") + temp_pepnovo_outfile);

				writeLog_("Use this line to call PepNovo: ");
				writeLog_(call);

				Int status=system(call.c_str());

				if (status != 0)
				{
					writeLog_("PepNovo problem. Aborting! (Details can be seen in the logfile: \"" + logfile + "\")");
					// clean temporary files
					/*
					for(QStringList::ConstIterator file_it=pepnovo_files.begin(); file_it!=pepnovo_files.end(); ++file_it)
					{
						qdir_temp.remove(*file_it);
					}
					qdir_temp.cdUp();
					//qdir_temp.remove("tmp_pepnovo_out.txt");
					qdir_temp.rmdir("Models");
					*/
					return EXTERNAL_PROGRAM_ERROR;
				}

				//if PepNovo finished succesfully use PepNovoOutfile to parse the results and generate idxml
				std::vector< PeptideIdentification > peptide_identifications;
				ProteinIdentification protein_identification;

				PepNovoOutfile p_novo_outfile;
				p_novo_outfile.load(temp_pepnovo_outfile, peptide_identifications, protein_identification, std::numeric_limits<Real>::max(), id_to_rt);
				IdXMLFile().store(outputfile_name,std::vector<ProteinIdentification>(1,protein_identification),peptide_identifications);

				for(QStringList::ConstIterator file_it=pepnovo_files.begin(); file_it!=pepnovo_files.end(); ++file_it)
        {
          //std::cout<<file_it->toStdString()<<std::endl;
          if(qdir_temp.cd(*file_it))
          {
            QStringList subdir_files = qdir_temp.entryList(QDir::Dirs | QDir::Files|QDir::NoDotAndDotDot);
            for(QStringList::ConstIterator subdir_file_it=subdir_files.begin(); subdir_file_it!=subdir_files.end(); ++subdir_file_it)
            {
              qdir_temp.remove(*subdir_file_it);
            }
            qdir_temp.cdUp();
            qdir_temp.rmdir(*file_it);
          }
          else
          {
            qdir_temp.remove(*file_it);
          }
        }
				qdir_temp.cdUp();
				qdir_temp.remove("tmp_pepnovo_out.txt");
				qdir_temp.rmdir("Models");

				return EXECUTION_OK;

			}
			catch(...)
			{
				//remove all possibly created files and folders ion case of uexpected behavior
				qdir_temp.setPath(temp_data_directory);
				if(qdir_temp.cd("Models"))
				{
					QStringList pepnovo_files = qdir_temp.entryList(QDir::Dirs | QDir::Files|QDir::NoDotAndDotDot);
					for(QStringList::ConstIterator file_it=pepnovo_files.begin(); file_it!=pepnovo_files.end(); ++file_it)
          {
            //std::cout<<file_it->toStdString()<<std::endl;
            if(qdir_temp.cd(*file_it))
            {
              QStringList subdir_files = qdir_temp.entryList();
              for(QStringList::ConstIterator subdir_file_it=subdir_files.begin(); subdir_file_it!=subdir_files.end(); ++subdir_file_it)
              {
                qdir_temp.remove(*subdir_file_it);
              }
              qdir_temp.cdUp();
              qdir_temp.rmdir(*file_it);
            }
            else
            {
              qdir_temp.remove(*file_it);
            }
					qdir_temp.cdUp();
					qdir_temp.remove("tmp_pepnovo_out.txt");
					qdir_temp.rmdir("Models");
				}
			}
      return EXTERNAL_PROGRAM_ERROR;
		}
		}

/*
			c
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
					if ( type == FileTypes::UNKNOWN )
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
					if ( type == FileTypes::UNKNOWN )
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
		*/
};

//@endcond



int main( int argc, const char** argv )
{
	TOPPPepNovoAdapter tool;

	return tool.main(argc,argv);
}
