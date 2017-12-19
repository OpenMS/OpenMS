// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Sandro Andreotti, Chris Bielow $
// --------------------------------------------------------------------------


#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/PepNovoInfile.h>
#include <OpenMS/FORMAT/PepNovoOutfile.h>
#include <OpenMS/FORMAT/MascotGenericFile.h>
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
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any signal-/preprocessing tool @n (in mzML format)</td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
    </tr>
  </table>
</CENTER>

  This wrapper application serves for getting peptide identifications
  for MS/MS spectra.

  The whole process of identification via PepNovo is executed.
  Input file is one mzML file containing the MS/MS spectra
  for which the identifications are to be found. The results are written
  as an idXML output file.

  The resulting idXML file can then be directly mapped to the spectra using the
  IDMapper class.

  Consult your PepNovo reference manual for further details about parameter meanings.

 @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_PepNovoAdapter.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_PepNovoAdapter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPPepNovoAdapter :
  public TOPPBase
{
  public:
  TOPPPepNovoAdapter() :
    TOPPBase("PepNovoAdapter", "Adapter to PepNovo supporting all PepNovo command line parameters. The results are converted from the PepNovo text outfile format into the idXML format.")
    {
    }

  protected:

    void registerOptionsAndFlags_() override
    {
      registerInputFile_("in", "<file>", "", "input file ");
      setValidFormats_("in", ListUtils::create<String>("mzML"));

      registerOutputFile_("out", "<file>", "", "output file ");
      setValidFormats_("out",ListUtils::create<String>("idXML"));

      registerInputFile_("pepnovo_executable","<file>", "", "The \"PepNovo\" executable of the PepNovo installation", true, false, ListUtils::create<String>("skipexists"));
      registerStringOption_("model_directory", "<file>", "", "Name of the directory where the model files are kept.",true);

      addEmptyLine_ ();
      registerFlag_("correct_pm", "Find optimal precursor mass and charge values.");
      registerFlag_("use_spectrum_charge", "Do not correct charge");
      registerFlag_("use_spectrum_mz", "Do not correct the precursor m/z value that appears in the file.");
      registerFlag_("no_quality_filter", "Do not remove low quality spectra.");
      registerDoubleOption_("fragment_tolerance", "<Float>", -1.0, "The fragment tolerance (between 0 and 0.75 Da. Set to -1.0 to use model's default setting)", false, false);
      registerDoubleOption_("pm_tolerance", "<Float>", -1.0, "The precursor mass tolerance (between 0 and 5.0 Da. Set to -1.0 to use model's default setting)", false, false);
      registerStringOption_("model", "<file>", "CID_IT_TRYP", "Name of the model that should be used", false);

      registerStringOption_("digest", "", "TRYPSIN", "Enzyme used for digestion (default TRYPSIN)", false);
      setValidStrings_("digest", ListUtils::create<String>("TRYPSIN,NON_SPECIFIC"));

      registerIntOption_("tag_length", "<num>", -1, "Returns peptide sequence of the specified length (only lengths 3-6 are allowed)", false);

      registerIntOption_("num_solutions", "<num>", 20, "Number of solutions to be computed", false);
      setMinInt_("num_solutions",1);
      setMaxInt_("num_solutions", 2000);

      std::vector<String>all_possible_modifications;
      ModificationsDB::getInstance()->getAllSearchModifications(all_possible_modifications);
      registerStringList_("fixed_modifications", "<mod1,mod2,...>", ListUtils::create<String>(""), "List of fixed modifications", false);
      setValidStrings_("fixed_modifications", all_possible_modifications);
      registerStringList_("variable_modifications", "<mod1,mod2,...>", ListUtils::create<String>(""), "List of variable modifications", false);
      setValidStrings_("variable_modifications", all_possible_modifications);
    }

    ExitCodes main_(int , const char**) override
    {

      // path to the log file
      String logfile(getStringOption_("log"));
      String pepnovo_executable(getStringOption_("pepnovo_executable"));

      PeakMap exp;

      String inputfile_name = getStringOption_("in");
      writeDebug_(String("Input file: ") + inputfile_name, 1);

      String outputfile_name = getStringOption_("out");
      writeDebug_(String("Output file: ") + outputfile_name, 1);

      String model_directory = getStringOption_("model_directory");
      writeDebug_(String("model directory: ") + model_directory, 1);

      String model_name = getStringOption_("model");
      writeDebug_(String("model directory: ") + model_name, 1);

      double fragment_tolerance = getDoubleOption_("fragment_tolerance");
      if (fragment_tolerance!=-1.0 && (fragment_tolerance<0 || fragment_tolerance>0.75))
      {
        writeLog_("Invalid fragment tolerance");
        printUsage_();
        return ILLEGAL_PARAMETERS;
      }

      double pm_tolerance = getDoubleOption_("pm_tolerance");
      if (pm_tolerance!=-1.0 && (pm_tolerance<0.0 || pm_tolerance>5.0))
      {
        writeLog_("Invalid fragment tolerance");
        printUsage_();
        return ILLEGAL_PARAMETERS;
      }

      Int tag_length = getIntOption_("tag_length");
      if ( tag_length!=-1 && (tag_length<3 || tag_length>6))
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
      MzMLFile mzml_infile;
      mzml_infile.getOptions().addMSLevel(2);
      mzml_infile.setLogType(log_type_);
      mzml_infile.load(inputfile_name, exp);

      // we map the native id to the MZ and RT to be able to
      // map the IDs back to the spectra (RT, and MZ Meta Information)
      PepNovoOutfile::IndexPosMappingType index_to_precursor;
      for (Size i = 0; i < exp.size(); ++i)
      {
        index_to_precursor[i]= make_pair(exp[i].getRT(), exp[i].getPrecursors()[0].getPosition()[0]); //set entry <RT, MZ>
      }

      logfile = getStringOption_("log");
      
      QDir qdir_models_source(model_directory.c_str());
      if (!qdir_models_source.exists())
      {
        writeLog_("The model directory does not exist");
        return INPUT_FILE_NOT_FOUND;
      }
      
      // create temp directory
      QDir qdir_temp(File::getTempDirectory().toQString());
      String temp_data_directory = File::getUniqueName();
      qdir_temp.mkdir(temp_data_directory.toQString());
      qdir_temp.cd(temp_data_directory.toQString());
      temp_data_directory  = File::getTempDirectory() + "/" + temp_data_directory; // delete later

      String mgf_file = temp_data_directory + "/" + File::getUniqueName() + ".mgf"; // the mzXML parser of PepNovo is somewhat broken.. don't use mzXML
      MascotGenericFile().store(mgf_file, exp);

      bool error(false);

      try
      {
        //temporary File to store PepNovo output
        String temp_pepnovo_outfile = qdir_temp.absoluteFilePath("tmp_pepnovo_out.txt");
        String tmp_models_dir = qdir_temp.absoluteFilePath("Models");

        std::map<String, String>mods_and_keys; //, key_to_id;

        if (qdir_temp.cd("Models"))
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
        if (pepnovo_files.empty())
        {
          writeLog_("The \"Model\" directory does not contain model files. Aborting!");
          return INPUT_FILE_NOT_FOUND;
        }

        for (QStringList::ConstIterator file_it=pepnovo_files.begin(); file_it!=pepnovo_files.end(); ++file_it)
        {
          if (qdir_models_source.cd(*file_it))
          {
            qdir_temp.mkdir(*file_it);
            qdir_temp.cd(*file_it);
            QStringList subdir_files = qdir_models_source.entryList(QDir::Dirs | QDir::Files|QDir::NoDotAndDotDot);
            for (QStringList::ConstIterator subdir_file_it=subdir_files.begin(); subdir_file_it!=subdir_files.end(); ++subdir_file_it)
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
        if (!getStringList_("fixed_modifications").empty() || !getStringList_("variable_modifications").empty())
        {
          p_novo_infile.setModifications(getStringList_("fixed_modifications"), getStringList_("variable_modifications"));
          p_novo_infile.store(qdir_temp.filePath("PepNovo_PTMs.txt"));
          pepnovo_files.append("PepNovo_PTMs.txt");
          p_novo_infile.getModifications(mods_and_keys);

          for (std::map<String, String>::const_iterator key_it=mods_and_keys.begin(); key_it!=mods_and_keys.end();++key_it)
          {
            if (ptm_command!="")
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

        arguments << "-file" << mgf_file.toQString();
        arguments << "-model" << model_name.toQString();
        if (pm_tolerance != -1 ) arguments << "-pm_tolerance"<<String(pm_tolerance).toQString();
        if (fragment_tolerance != -1 ) arguments << "-fragment_tolerance" <<String(fragment_tolerance).toQString();
        if (!ptm_command.empty()) arguments <<"-PTMs" <<ptm_command.toQString();
        if (getFlag_("correct_pm")) arguments << "-correct_pm";
        if (getFlag_("use_spectrum_charge")) arguments << "-use_spectrum_charge";
        if (getFlag_("use_spectrum_mz")) arguments << "-use_spectrum_mz";
        if (getFlag_("no_quality_filter")) arguments << "-no_quality_filter";
        arguments << "-digest" << digest.toQString();
        arguments << "-num_solutions" << String(num_solutions).toQString();
        if (tag_length!=-1) arguments<<"-tag_length" << String(tag_length).toQString();
        arguments<<"-model_dir" << tmp_models_dir.toQString();
        //arguments<<">" << temp_pepnovo_outfile.toQString();

        writeDebug_("Use this line to call PepNovo: ", 1);
        writeDebug_(pepnovo_executable + " " + String(arguments.join(" ")), 1);
        QProcess process;
        process.setStandardOutputFile(temp_pepnovo_outfile.toQString());
        process.setStandardErrorFile(temp_pepnovo_outfile.toQString());
        process.start(pepnovo_executable.toQString(), arguments); // does automatic escaping etc...
        if (process.waitForFinished(-1))
        {
          //if PepNovo finished successfully use PepNovoOutfile to parse the results and generate idXML
          std::vector< PeptideIdentification > peptide_identifications;
          ProteinIdentification protein_identification;
          StringList ms_runs;
          exp.getPrimaryMSRunPath(ms_runs);
          protein_identification.setPrimaryMSRunPath(ms_runs);

          PepNovoOutfile p_novo_outfile;

          //resolve PTMs (match them back to the OpenMs Identifier String)
          std::vector<ProteinIdentification>prot_ids;
          p_novo_outfile.load(temp_pepnovo_outfile, peptide_identifications, protein_identification, -1e5, index_to_precursor, mods_and_keys);
          prot_ids.push_back(protein_identification);
          IdXMLFile().store(outputfile_name, prot_ids, peptide_identifications);
        }

        if (process.exitStatus() != 0)  error = true;
       
      }
      catch(Exception::BaseException &exc)
      {
        writeLog_(exc.what());
        LOG_ERROR << "Error occurred: " << exc.what() << std::endl;
        error = true;
      }
      
      if (!error)
      {
        File::removeDirRecursively(temp_data_directory);
        return EXECUTION_OK;
      }
      else
      {
        writeLog_("PepNovo problem. Aborting! (Details can be seen in outfiles: '" + temp_data_directory + "')");
        return EXTERNAL_PROGRAM_ERROR;
      }

    }

};


int main( int argc, const char** argv )
{
  TOPPPepNovoAdapter tool;
  return tool.main(argc,argv);
}

/// @endcond

