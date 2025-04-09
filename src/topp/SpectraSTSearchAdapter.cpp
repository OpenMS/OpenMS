// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// -----------------------------------------
// $Maintainer: Lukas Zimmermann $
// $Authors: Lukas Zimmermann $
// --------------------------------------------------------------------------
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/SYSTEM/File.h>

#include <sstream> 

#include <QDir>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_SpectraSTSearchAdapter SpectraSTSearchAdapter

@brief This util provides an interface to the 'SEARCH' mode of the SpectraST program.
       All non-advanced parameters of the executable of SpectraST were translated into
       parameters of this util.

SpectraST: Version: 5

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_SpectraSTSearchAdapter.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_SpectraSTSearchAdapter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPSpectraSTSearchAdapter :
  public TOPPBase
{
 public:
  // Define parameter name
  static const String param_executable;
  static const String param_spectra_files;
  static const String param_library_file;
  static const String param_sequence_database_file;
  static const String param_sequence_database_type;
  static const String param_search_file;
  static const String param_params_file;
  static const String param_precursor_mz_tolerance;
  static const String param_use_isotopically_averaged_mass;
  static const String param_use_all_charge_states;
  static const String param_output_files;
  static const vector<String> param_output_file_formats;
  static const vector<String> param_input_file_formats;
  static const String param_user_mod_file;

  TOPPSpectraSTSearchAdapter() :
    TOPPBase("SpectraSTSearchAdapter", "Interface to the SEARCH Mode of the SpectraST executable")
  {
  }

protected:

  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_() override
  {
      StringList empty;

      // Handle executable
      registerInputFile_(TOPPSpectraSTSearchAdapter::param_executable, "<path>", "spectrast", "Path to the SpectraST executable to use; may be empty if the executable is globally available.", true, false, ListUtils::create<String>("skipexists"));

      // register spectra input files
      registerInputFileList_(TOPPSpectraSTSearchAdapter::param_spectra_files, "<SearchFileName1> [ <SearchFileName2> ... <SearchFileNameN> ]", empty, "File names(s) of spectra to be searched.", true, false);
      setValidFormats_(TOPPSpectraSTSearchAdapter::param_spectra_files, ListUtils::create<String>(TOPPSpectraSTSearchAdapter::param_input_file_formats), false);

      // register Output files
      registerOutputFileList_(TOPPSpectraSTSearchAdapter::param_output_files, "<OutputFile1> [ <OutputFileName2> ... <OutputFileNameN> ]", empty, "Output files. Make sure to specify one output file for each input file", true, false);
      setValidFormats_(TOPPSpectraSTSearchAdapter::param_output_files, TOPPSpectraSTSearchAdapter::param_output_file_formats, false);

      // Require library file to be searched
      registerInputFile_(TOPPSpectraSTSearchAdapter::param_library_file, "<lib_file>.splib", "", "Specify library file.", true, false);
      setValidFormats_(TOPPSpectraSTSearchAdapter::param_library_file, ListUtils::create<String>("splib"), false);

      // Sequence database file
      registerInputFile_(TOPPSpectraSTSearchAdapter::param_sequence_database_file, "<sequencedb_file>.fasta", "", "The sequence database.", false, false);
      setValidFormats_(TOPPSpectraSTSearchAdapter::param_sequence_database_file, ListUtils::create<String>("fasta"), false);

      registerStringOption_(TOPPSpectraSTSearchAdapter::param_sequence_database_type, "<sequencedb_type>", "AA", "Specify type of sequence database", false, false);
      setValidStrings_(TOPPSpectraSTSearchAdapter::param_sequence_database_type, ListUtils::create<String>("DNA,AA"));

      // Search file
      registerInputFile_(TOPPSpectraSTSearchAdapter::param_search_file, "<search_file>", "", "Only search a subset of the query spectra in the search file", false, false);
      setValidFormats_(TOPPSpectraSTSearchAdapter::param_search_file, ListUtils::create<String>("txt, dat"), false);

      // params file
      registerInputFile_(TOPPSpectraSTSearchAdapter::param_params_file, "<params_file>", "", "Read search options from file. All options set in the file will be overridden by command-line options, if specified.", false, false);
      setValidFormats_(TOPPSpectraSTSearchAdapter::param_params_file, ListUtils::create<String>("params"), false);

      // Precursor m/z tolerance
      registerDoubleOption_(TOPPSpectraSTSearchAdapter::param_precursor_mz_tolerance, "<precursor_mz_tolerance>", 3, "m/z (in Th) tolerance within which candidate entries are compared to the query. Monoisotopic mass is assumed.", false, false);
      setMinFloat_(TOPPSpectraSTSearchAdapter::param_precursor_mz_tolerance, 0);

      //Whether to use isotope average instead of monoisotopic mass
      registerFlag_(TOPPSpectraSTSearchAdapter::param_use_isotopically_averaged_mass, "Use isotopically averaged mass instead of monoisotopic mass", true);

      // Whether to use all charge states
      registerFlag_(TOPPSpectraSTSearchAdapter::param_use_all_charge_states, "Search library spectra of all charge states, i.e., ignore specified charge state (if any) of the query spectrum", true);

      // User defined modifications file
      registerInputFile_(TOPPSpectraSTSearchAdapter::param_user_mod_file, "<user_mod_file>", "", "Specify name of user-defined modifications file. Default is \"spectrast.usermods\".", false, true);
  }

  // the main_ function is called after all parameters are read
  ExitCodes main_(int, const char **) override
  {
     // Assemble command line for SpectraST
     QStringList arguments;

     // Executable
     String executable = getStringOption_(TOPPSpectraSTSearchAdapter::param_executable);
     if (executable.empty())
     {
         executable = "spectrast";
     }

     // Make Search Mode explicit
     arguments << "-s";

     double precursor_mz_tolerance = getDoubleOption_(TOPPSpectraSTSearchAdapter::param_precursor_mz_tolerance);
     arguments << QString::number(precursor_mz_tolerance).prepend("-sM");

     // Set the parameter file if present
     String params_file = getStringOption_(TOPPSpectraSTSearchAdapter::param_params_file);
     if (! params_file.empty())
     {
         arguments << params_file.toQString().prepend("-sF");
     }

     // Add library file argument, terminate if the corresponding spidx is not present
     String library_file = getStringOption_(TOPPSpectraSTSearchAdapter::param_library_file);
     String index_file = FileHandler::stripExtension(library_file).append(".spidx");
     if (! File::exists(index_file))
     {
         OPENMS_LOG_ERROR << "ERROR: Index file required by spectrast not found:\n" << index_file << endl;
         return INPUT_FILE_NOT_FOUND;
     }
     arguments << library_file.toQString().prepend("-sL");

     // Add Sequence Database file if exists
     String sequence_database_file  = getStringOption_(TOPPSpectraSTSearchAdapter::param_sequence_database_file);
     if (! sequence_database_file.empty())
     {
         String sequence_database_type = getStringOption_(TOPPSpectraSTSearchAdapter::param_sequence_database_type);

         // Check empty or invalid sequence database type
         if (sequence_database_type.empty())
         {
            OPENMS_LOG_ERROR << "ERROR: Sequence database type invalid or not provided" << endl;
            return MISSING_PARAMETERS;
         }
         arguments << sequence_database_type.toQString().prepend("-sT");
         arguments << sequence_database_file.toQString().prepend("-sD");
     }

     // Set the number of threads in SpectraST
     Int threads = getIntOption_("threads");
     arguments << (threads > 1 ?  QString::number(threads).prepend("-sP") : "-sP!");

     // Set the search file
     String search_file = getStringOption_(TOPPSpectraSTSearchAdapter::param_search_file);
     if (! search_file.empty())
     {
         arguments << search_file.toQString().prepend("-sS");
     }

     // Flags
     arguments << (getFlag_(TOPPSpectraSTSearchAdapter::param_use_isotopically_averaged_mass) ? "-sA" : "-sA!");
     arguments << (getFlag_(TOPPSpectraSTSearchAdapter::param_use_all_charge_states) ? "-sz" : "-sz!");

     // User mod file
     String user_mod_file = getStringOption_(TOPPSpectraSTSearchAdapter::param_user_mod_file);
     if (! user_mod_file.empty())
     {
        arguments << user_mod_file.toQString().prepend("-M");
     }

     // Input and output files, errors if lists are not equally long
     StringList spectra_files = getStringList_(TOPPSpectraSTSearchAdapter::param_spectra_files);
     StringList output_files = getStringList_(TOPPSpectraSTSearchAdapter::param_output_files);
     if (spectra_files.size() != output_files.size())
     {
        OPENMS_LOG_ERROR << "ERROR: Number of output files does not match number of input files." << endl;
        return ILLEGAL_PARAMETERS;
     }
     if (spectra_files.empty())
     {
         OPENMS_LOG_ERROR << "ERROR: At least one file containing spectra to be searched must be provided." << endl;
         return ILLEGAL_PARAMETERS;
     }
     String first_output_file = output_files[0];
     String outputFormat;
     // Determine the output format with the first file
     for (StringList::const_iterator it = TOPPSpectraSTSearchAdapter::param_output_file_formats.begin();
         it != TOPPSpectraSTSearchAdapter::param_output_file_formats.end(); ++it)
     {
          String format = *it;
          if (first_output_file.hasSuffix(format))
          {
              outputFormat = format;
          }
     }
     if (outputFormat.empty())
     {
         OPENMS_LOG_ERROR << "ERROR: Unrecognized output format from file: " << first_output_file << endl;
         return ILLEGAL_PARAMETERS;
     }
     // Output files must agree on format
     for (StringList::const_iterator it = output_files.begin(); it != output_files.end(); ++it)
     {
         String output_file = *it;
         if (! output_file.hasSuffix(outputFormat))
         {
             OPENMS_LOG_ERROR << "ERROR: Output filename does not agree in format: "
                       << output_file << " is not " << outputFormat << endl;
             return ILLEGAL_PARAMETERS;
         }
     }

     String temp_dir = File::getTempDirectory();
     arguments << outputFormat.toQString().prepend("-sE");
     arguments << temp_dir.toQString().prepend("-sO");

     // Check whether input files agree in format
     String first_input_file = spectra_files[0];
     String inputFormat;

     // Determine the output format with the first file
     for (StringList::const_iterator it = TOPPSpectraSTSearchAdapter::param_input_file_formats.begin();
         it != TOPPSpectraSTSearchAdapter::param_input_file_formats.end(); ++it)
     {
          String format = *it;
          if (first_input_file.hasSuffix(format))
          {
              inputFormat = format;
          }
     }

     // Exit if the input file format is invalid
     if (inputFormat.empty())
     {
         OPENMS_LOG_ERROR << "ERROR: Unrecognized input format from file: " << first_input_file << endl;
         return ILLEGAL_PARAMETERS;
     }

     for (vector<String>::const_iterator it = spectra_files.begin(); it != spectra_files.end(); ++it)
     {
      String input_file = *it;
      if (! input_file.hasSuffix(inputFormat))
      {
        OPENMS_LOG_ERROR << "ERROR: Input filename does not agree in format: "
                       << input_file << " is not " << inputFormat << endl;
        return ILLEGAL_PARAMETERS;
      }
      arguments << input_file.toQString();
     }

     // Writing the final SpectraST command to the DEBUG LOG
     std::stringstream ss;
     ss << "COMMAND: " << executable;
     for (QStringList::const_iterator it = arguments.begin(); it != arguments.end(); ++it)
     {
         ss << " " << it->toStdString();
     }
     OPENMS_LOG_DEBUG << ss.str() << endl;

     // Run SpectraST
     TOPPBase::ExitCodes exit_code = runExternalProcess_(executable.toQString(), arguments);
     if (exit_code != EXECUTION_OK)
     {
       return exit_code;
     }

     // Copy the output files to the specified location
     QDir temp_dir_qt = QDir(temp_dir.toQString());
     for (size_t i = 0; i < spectra_files.size(); i++)
     {
        String spectra_file = spectra_files[i];
        QString actual_path = temp_dir_qt.filePath(FileHandler::stripExtension(File::basename(spectra_file)).toQString().append(".").append(outputFormat.toQString()));

        std::ifstream ifs(actual_path.toStdString().c_str(), std::ios::in | std::ios::binary);
        std::ofstream ofs(output_files[i].c_str(), std::ios::out | std::ios::binary);
        ofs << ifs.rdbuf();
     }

    // Exit the tool
    return EXECUTION_OK;
  }

};
// End of Tool definition

// Definition of static members
const String TOPPSpectraSTSearchAdapter::param_executable = "executable";
const String TOPPSpectraSTSearchAdapter::param_spectra_files = "spectra_files";
const String TOPPSpectraSTSearchAdapter::param_library_file = "library_file";
const String TOPPSpectraSTSearchAdapter::param_sequence_database_file = "sequence_database_file";
const String TOPPSpectraSTSearchAdapter::param_sequence_database_type = "sequence_database_type";
const String TOPPSpectraSTSearchAdapter::param_search_file = "search_file";
const String TOPPSpectraSTSearchAdapter::param_params_file = "params_file";
const String TOPPSpectraSTSearchAdapter::param_precursor_mz_tolerance = "precursor_mz_tolerance";
const String TOPPSpectraSTSearchAdapter::param_use_isotopically_averaged_mass = "use_isotopically_averaged_mass";
const String TOPPSpectraSTSearchAdapter::param_use_all_charge_states = "use_all_charge_states";
const String TOPPSpectraSTSearchAdapter::param_output_files = "output_files";
const String TOPPSpectraSTSearchAdapter::param_user_mod_file = "user_mod_file";
const StringList TOPPSpectraSTSearchAdapter::param_output_file_formats = ListUtils::create<String>("txt,tsv,xml,pepXML,html");
const StringList TOPPSpectraSTSearchAdapter::param_input_file_formats = ListUtils::create<String>("mzML,mzXML,mzData,mgf,dta,msp");

// the actual main function needed to create an executable
int main(int argc, const char ** argv)
{
  TOPPSpectraSTSearchAdapter tool;
  return tool.main(argc, argv);
}

///@endcond
