// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Sandro Andreotti $
// $Authors: Sandro Andreotti $
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

#include <OpenMS/ANALYSIS/ID/IDMapper.h>

#include <QtCore/QFile>
#include <QtCore/QDir>
#include <QtCore/QProcess>

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

<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ PepNovoAdapter \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any signal-/preprocessing tool @n (in mzXML format)</td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
		</tr>
	</table>
</CENTER>

	This wrapper application serves for getting peptide identifications
	for MS/MS spectra.

	The whole process of identification via PepNovo is executed.
	Inputfile is one mzXML file containing the MS/MS spectra
	for which the identifications are to be found. The results are written
	as an idXML output file.

	The resulting idXML file can then be directly mapped to the spectra using the
	IDMapper class.

	Consult your PepNovo reference manual for further details about parameter meanings.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_PepNovoAdapter.cli
	<B>INI file documentation of this tool:</B>
	@htmlinclude TOPP_PepNovoAdapter.html
*/

// We do not want this class to show up in the docu -> cond
// @cond

class TOPPPepNovoAdapter
	: public TOPPBase
{
	public:
		TOPPPepNovoAdapter()
			: TOPPBase("PepNovoAdapter", "Adapter to PepNovo supporting all PepNovo command line parameters. The results are converted from the PepNovo text outfile format into the idXML format.")
		{
		}

	protected:

		void registerOptionsAndFlags_()
		{
			registerInputFile_("in", "<file>", "", "input file ");
      setValidFormats_("in",StringList::create("mzXML"));

			registerOutputFile_("out", "<file>", "", "output file ");
			setValidFormats_("out",StringList::create("idXML"));

			registerInputFile_("pepnovo_executable","<file>", "", "The \"PepNovo\" executable of the PepNovo installation", true, false, StringList::create("skipexists"));
			registerStringOption_("temp_data_directory", "<dir>", "", "Directory were temporary data can be stored. If not set the directory were startet is used.", true);
			registerStringOption_("model_directory", "<file>", "", "Name of the directory where the model files are kept.",true);
      addEmptyLine_ ();
      addText_("PepNovo Parameters");
      registerFlag_("correct_pm", "Find optimal precursor mass and charge values.");
      registerFlag_("use_spectrum_charge", "Do not correct charge");
      registerFlag_("use_spectrum_mz", "Do not correct the precursor m/z value that appears in the file.");
      registerFlag_("no_quality_filter", "Do not remove low quality spectra.");
      registerDoubleOption_("fragment_tolerance", "<Float>", -1.0, "The fragment tolerance (between 0 and 0.75 Da. Set to -1.0 to use model's default setting)", false, false);
      registerDoubleOption_("pm_tolerance", "<Float>", -1.0, "The precursor mass tolerance (between 0 and 5.0 Da. Set to -1.0 to use model's default setting)", false, false);
      registerStringOption_("model", "<file>", "CID_IT_TRYP", "Name of the model that should be used", false);

			registerStringOption_("digest", "", "TRYPSIN", "Enzyme used for digestion (default TRYPSIN)", false);
			setValidStrings_("digest", StringList::create("TRYPSIN,NON_SPECIFIC"));

			registerIntOption_("tag_length", "<num>", -1, "Returns peptide sequence of the specified length (only lengths 3-6 are allowed)", false);

			registerIntOption_("num_solutions", "<num>", 20, "Number of solutions to be computed", false);
			setMinInt_("num_solutions",1);
			setMaxInt_("num_solutions",2000);

			std::vector<String>all_possible_modifications;
			ModificationsDB::getInstance()->getAllSearchModifications(all_possible_modifications);
			registerStringList_("fixed_modifications", "<mod1,mod2,...>", StringList::create(""), "List of fixed modifications", false);
			setValidStrings_("fixed_modifications", all_possible_modifications);
			registerStringList_("variable_modifications", "<mod1,mod2,...>", StringList::create(""), "List of variable modifications", false);
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

			outputfile_name = getStringOption_("out");
			writeDebug_(String("Output file: ") + outputfile_name, 1);

			model_directory = getStringOption_("model_directory");
			writeDebug_(String("model directory: ") + model_directory, 1);

			String model_name = getStringOption_("model");
			writeDebug_(String("model directory: ") + model_name, 1);

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
      String digest = getStringOption_("digest");
      Size num_solutions=getIntOption_("num_solutions");

			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------

			// only load msLevel 2
      MzXMLFile mzdata_infile;
			mzdata_infile.getOptions().addMSLevel(2);
			mzdata_infile.setLogType(log_type_);
			mzdata_infile.load(inputfile_name, exp);

			// we map the native id to the MZ and RT to be able to
			// map the IDs back to the spectra (RT, and MZ Meta Information)
			std::map<String, pair<DoubleReal, DoubleReal> >id_to_rt;
			for (PeakMap::Iterator it = exp.begin(); it != exp.end(); ++it)
			{
			  Int valid_id;
			  Size num_pos=0;
			  String native_id=it->getNativeID();

			  while(!isdigit(native_id[num_pos]) && num_pos<native_id.length())
			  {
			    ++num_pos;
			  }
			  if(num_pos==native_id.length())
			  {
			    writeLog_("No valid NativeId for spectrum. Aborting!");
          return INPUT_FILE_CORRUPT;
			  }
			  else
			  {
			    valid_id=native_id.substr(num_pos).toInt();
			  }
			  id_to_rt[valid_id]=make_pair(it->getRT(), it->getPrecursors()[0].getPosition()[0]); //set entry <RT, MZ>
				//std::cout<<"stored id: "<<valid_id<<std::endl;
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

			try
			{
			  //temporary File to store PepNovo output
			  String temp_pepnovo_outfile = qdir_temp.absoluteFilePath("tmp_pepnovo_out.txt");
			  String tmp_models_dir=qdir_temp.absoluteFilePath("Models");

        std::map<String, String>mods_and_keys; //, key_to_id;

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
				    qdir_temp.mkdir(*file_it);
				    qdir_temp.cd(*file_it);
				    QStringList subdir_files = qdir_models_source.entryList(QDir::Dirs | QDir::Files|QDir::NoDotAndDotDot);
				    for(QStringList::ConstIterator subdir_file_it=subdir_files.begin(); subdir_file_it!=subdir_files.end(); ++subdir_file_it)
            {
				      QFile::copy(qdir_models_source.filePath(*subdir_file_it), qdir_temp.filePath(*subdir_file_it));
            }
				    qdir_temp.cdUp();
				    qdir_models_source.cdUp();
				  }
          else
          {
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
					p_novo_infile.getModifications(mods_and_keys);

					for(std::map<String, String>::const_iterator key_it=mods_and_keys.begin(); key_it!=mods_and_keys.end();++key_it)
					{
						if(ptm_command!="")
						{
							ptm_command+=":";
						}
						ptm_command+= key_it->first;
						//key_to_id[key_it->second]=key_it->first;
					}
				}

				//-------------------------------------------------------------
				// (3) running program according to parameters
				//-------------------------------------------------------------
        QStringList arguments;

        arguments<<"-file" << inputfile_name.toQString();
        arguments<<"-model" << model_name.toQString();
        if (pm_tolerance != -1 ) arguments<<"-pm_tolerance"<<String(pm_tolerance).toQString();
        if (fragment_tolerance != -1 ) arguments<<"-fragment_tolerance" <<String(fragment_tolerance).toQString();
        if (!ptm_command.empty()) arguments<<"-PTMs" <<ptm_command.toQString();
        if(getFlag_("correct_pm")) arguments<<"-correct_pm";
        if(getFlag_("use_spectrum_charge")) arguments<<"-use_spectrum_charge";
        if(getFlag_("use_spectrum_mz")) arguments<<"-use_spectrum_mz";
        if(getFlag_("no_quality_filter")) arguments<<"-no_quality_filter";
        arguments<<"-digest" << digest.toQString();
        arguments<<"-num_solutions" << String(num_solutions).toQString();
        if(tag_length!=-1)arguments<<"-tag_length" << String(tag_length).toQString();
        arguments<<"-model_dir" << tmp_models_dir.toQString();
        //arguments<<">" << temp_pepnovo_outfile.toQString();

				writeLog_("Use this line to call PepNovo: ");
        writeLog_(arguments.join(" "));
        QProcess process;
        process.setStandardOutputFile(temp_pepnovo_outfile.toQString());
        process.setStandardErrorFile(temp_pepnovo_outfile.toQString());
        process.start(pepnovo_executable.toQString(), arguments); // does automatic escaping etc...
        if (process.waitForFinished(-1))
				{
          //if PepNovo finished succesfully use PepNovoOutfile to parse the results and generate idxml
          std::vector< PeptideIdentification > peptide_identifications;
          ProteinIdentification protein_identification;

          PepNovoOutfile p_novo_outfile;

          //resolve PTMs (match them back to the OpenMs Identifier String)
          std::vector<ProteinIdentification>prot_ids;
          p_novo_outfile.load(temp_pepnovo_outfile, peptide_identifications, protein_identification, (-1)*std::numeric_limits<DoubleReal>::max(), id_to_rt, mods_and_keys);
          prot_ids.push_back(protein_identification);
          IdXMLFile().store(outputfile_name,prot_ids, peptide_identifications);
        }

        //remove the temporary files
				for(QStringList::ConstIterator file_it=pepnovo_files.begin(); file_it!=pepnovo_files.end(); ++file_it)
        {
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
				qdir_temp.rmdir("Models");

        if(process.exitStatus() == 0)
				{
					qdir_temp.remove("tmp_pepnovo_out.txt");
				  return EXECUTION_OK;
				}
				else
				{
				  writeLog_("PepNovo problem. Aborting! (Details can be seen in outfile: \"" + qdir_temp.absoluteFilePath("tmp_pepnovo_out.txt") + "\")");
				  return EXTERNAL_PROGRAM_ERROR;
				}
      }
			catch(Exception::BaseException &exc)
			{
				//remove all possibly created files and folders ion case of unexpected error
				qdir_temp.setPath(temp_data_directory);
				if(qdir_temp.cd("Models"))
				{
					QStringList pepnovo_files = qdir_temp.entryList(QDir::Dirs | QDir::Files|QDir::NoDotAndDotDot);
					for(QStringList::ConstIterator file_it=pepnovo_files.begin(); file_it!=pepnovo_files.end(); ++file_it)
          {
					  std::cout<<qdir_temp.absolutePath().toStdString()<<std::endl;
            std::cout<<file_it->toStdString()<<std::endl;
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
        }
        writeLog_(exc.what());
        return EXTERNAL_PROGRAM_ERROR;
      }
    }
};

//@endcond

int main( int argc, const char** argv )
{
	TOPPPepNovoAdapter tool;

	return tool.main(argc,argv);
}
