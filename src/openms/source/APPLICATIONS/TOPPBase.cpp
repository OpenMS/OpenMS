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
// $Authors: Marc Sturm, Clemens Groepl, Johannes Junker, Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/SYSTEM/SysInfo.h>
#include <OpenMS/SYSTEM/UpdateCheck.h>

#include <OpenMS/DATASTRUCTURES/Date.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>

#include <OpenMS/KERNEL/ConsensusMap.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/FORMAT/VALIDATORS/XMLValidator.h>

#include <OpenMS/APPLICATIONS/ConsoleUtils.h>

#include <iostream>

#if  defined(__APPLE__)
  #include <QCoreApplication.h> // needed to disable plugin loading on Mac OSX
#endif

#include <QDir>
#include <QFile>

#include <boost/math/special_functions/fpclassify.hpp>

#include <time.h>
#include <stdio.h>
#include <stdlib.h>

// OpenMP support
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef OPENMS_WINDOWSPLATFORM
#undef min
#undef max
#endif

#include <cmath>

using namespace std;
  
namespace OpenMS
{

  using namespace Exception;

  String TOPPBase::topp_ini_file_ = String(QDir::homePath()) + "/.TOPP.ini";
  const Citation TOPPBase::cite_openms_ = { "Rost HL, Sachsenberg T, Aiche S, Bielow C et al.",
      "OpenMS: a flexible open-source software platform for mass spectrometry data analysis",
      "Nat Meth. 2016; 13, 9: 741-748",
      "10.1038/nmeth.3959" };

  void TOPPBase::setMaxNumberOfThreads(int
#ifdef _OPENMP
                                       num_threads // to avoid the unused warning we enable this
                                                   // argument only if openmp is available
#endif
                                       )
  {
#ifdef _OPENMP
    omp_set_num_threads(num_threads);
#endif
  }

  TOPPBase::TOPPBase(const String& tool_name, const String& tool_description, bool official, const std::vector<Citation>& citations) :
    tool_name_(tool_name),
    tool_description_(tool_description),
    instance_number_(-1),
    version_(""),
    verboseVersion_(""),
    official_(official),
    citations_(citations),
    log_type_(ProgressLogger::NONE),
    test_mode_(false),
    debug_level_(-1)
  {
    version_ = VersionInfo::getVersion();
    verboseVersion_ = version_ + " " + VersionInfo::getTime();

    // if the revision info is meaningful, show it as well
    if (!VersionInfo::getRevision().empty() && VersionInfo::getRevision() != "exported")
    {
      verboseVersion_ += String(", Revision: ") + VersionInfo::getRevision() + "";
    }

    //check if tool is in official tools list
    if (official_ && tool_name_ != "GenericWrapper" && !ToolHandler::getTOPPToolList().count(tool_name_))
    {
      writeLog_(String("Warning: Message to maintainer - If '") + tool_name_ + "' is an official TOPP tool, add it to the tools list in ToolHandler. If it is not, set the 'official' flag of the TOPPBase constructor to false.");
    }

#if  defined(__APPLE__)
    // we do not want to load plugins as this leads to serious problems
    // when shipping on mac os x
    QCoreApplication::setLibraryPaths(QStringList());
#endif
  }

  TOPPBase::~TOPPBase()
  {
    //delete log file if empty
    StringList log_files;
    if (!getParam_("log").isEmpty())
      log_files.push_back((String)(getParam_("log")));
    for (Size i = 0; i < log_files.size(); ++i)
    {
      if (File::empty(log_files[i]))
      {
        File::remove(log_files[i]);
      }
    }
  }

  TOPPBase::ExitCodes TOPPBase::main(int argc, const char** argv)
  {
    //----------------------------------------------------------
    //parse command line
    //----------------------------------------------------------

    //register values from derived TOPP tool
    registerOptionsAndFlags_();
    addEmptyLine_();
    //common section for all tools
    if (ToolHandler::getTOPPToolList().count(tool_name_))
      addText_("Common TOPP options:");
    else
      addText_("Common UTIL options:");
    registerStringOption_("ini", "<file>", "", "Use the given TOPP INI file", false);
    registerStringOption_("log", "<file>", "", "Name of log file (created only when specified)", false, true);
    registerIntOption_("instance", "<n>", 1, "Instance number for the TOPP INI file", false, true);
    registerIntOption_("debug", "<n>", 0, "Sets the debug level", false, true);
    registerIntOption_("threads", "<n>", 1, "Sets the number of threads allowed to be used by the TOPP tool", false);
    registerStringOption_("write_ini", "<file>", "", "Writes the default configuration file", false);
    registerStringOption_("write_ctd", "<out_dir>", "", "Writes the common tool description file(s) (Toolname(s).ctd) to <out_dir>", false, true);
    registerFlag_("no_progress", "Disables progress logging to command line", true);
    registerFlag_("force", "Overwrite tool specific checks.", true);
    registerFlag_("test", "Enables the test mode (needed for internal use only)", true);
    registerFlag_("-help", "Shows options");
    registerFlag_("-helphelp", "Shows all options (including advanced)", false);

    // parse command line parameters:
    try
    {
      param_cmdline_ = parseCommandLine_(argc, argv);
    }
    catch (Exception::BaseException& e)
    {
      writeLog_("Invalid parameter values (" + String(e.getName()) + "): " + String(e.getMessage()) + ". Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    // for now command line is all we have, final assembly will follow below
    param_ = param_cmdline_;

    // assign instance number
    *const_cast<int*>(&instance_number_) = getParamAsInt_("instance", 1);
    writeDebug_(String("Instance: ") + String(instance_number_), 1);

    // assign ini location
    *const_cast<String*>(&ini_location_) = tool_name_ + ':' + String(instance_number_) + ':';
    writeDebug_(String("Ini_location: ") + getIniLocation_(), 1);

    // set debug level
    debug_level_ = getParamAsInt_("debug", 0);
    writeDebug_(String("Debug level: ") + String(debug_level_), 1);

    // print command line to console
    StringList args;
    for (int i = 0; i < argc; ++i)
    {
      if (String(argv[i]).has(' ')) args.push_back(String("\"") + argv[i] + String("\"")); // surround with quotes if argument contains a space
      else args.push_back(argv[i]);
    }
    writeDebug_(String(" >> ") + ListUtils::concatenate(args, " "), 1);


    // test if no options were given
    if (argc == 1)
    {
      writeLog_("No options given. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    // '--help' given
    if (param_cmdline_.exists("-help") || param_cmdline_.exists("-helphelp"))
    {
      printUsage_();
      return EXECUTION_OK;
    }

    // test if unknown options were given
    if (param_cmdline_.exists("unknown"))
    {
      writeLog_(String("Unknown option(s) '") + getParamAsString_("unknown") + "' given. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    // test if unknown text argument were given (we do not use them)
    if (param_cmdline_.exists("misc"))
    {
      writeLog_(String("Trailing text argument(s) '") + getParamAsString_("misc") + "' given. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    ExitCodes result;
#ifndef DEBUG_TOPP
    try
    {
#endif
    // '-write_ini' given
    if (param_cmdline_.exists("write_ini"))
    {
      String write_ini_file = param_cmdline_.getValue("write_ini");
      outputFileWritable_(write_ini_file, "write_ini");
      Param default_params = getDefaultParameters_();

      // check if augmentation with -ini param is needed
      DataValue in_ini;
      if (param_cmdline_.exists("ini"))
      {
        in_ini = param_cmdline_.getValue("ini");
      }
      if (!in_ini.isEmpty())
      {
        Param ini_params;
        ParamXMLFile paramFile;
        paramFile.load((String)in_ini, ini_params);

        // check if ini parameters are applicable to this tool
        checkIfIniParametersAreApplicable_(ini_params);
        // update default params with outdated params given in -ini and be verbose
        default_params.update(ini_params, false);
      }
      ParamXMLFile paramFile;
      paramFile.store(write_ini_file, default_params);
      return EXECUTION_OK;
    }

    // '-write_ctd' given
    if (param_cmdline_.exists("write_ctd"))
    {
      if (!writeCTD_())
      {
        writeLog_("Error: Could not write CTD file!");
        return INTERNAL_ERROR;
      }
      return EXECUTION_OK;
    }

    //-------------------------------------------------------------
    // load INI file
    //-------------------------------------------------------------
    {
      DataValue value_ini;

      if (param_cmdline_.exists("ini"))
      {
        value_ini = param_cmdline_.getValue("ini");
      }
      if (!value_ini.isEmpty())
      {
        writeDebug_("INI file: " + (String)value_ini, 1);
        writeDebug_("INI location: " + getIniLocation_(), 1);
        ParamXMLFile paramFile;

        paramFile.load((String)value_ini, param_inifile_);
        checkIfIniParametersAreApplicable_(param_inifile_);

        // dissect loaded INI parameters
        param_instance_ = param_inifile_.copy(getIniLocation_(), true);
        writeDebug_("Parameters from instance section:", param_instance_, 2);
        param_common_tool_ = param_inifile_.copy("common:" + tool_name_ + ":", true);
        writeDebug_("Parameters from common section with tool name:", param_common_tool_, 2);
        param_common_ = param_inifile_.copy("common:", true);
        writeDebug_("Parameters from common section without tool name:", param_common_, 2);

        checkIfIniParametersAreApplicable_(param_inifile_);

        // set type on command line if given in .ini file
        if (param_inifile_.exists(getIniLocation_() + "type") && !param_cmdline_.exists("type"))
          param_cmdline_.setValue("type", param_inifile_.getValue(getIniLocation_() + "type"));
      }

      // construct the set of final parameters as they will be available in the main_ method
      Param finalParam;

      // 1. the CMD parameters
      writeDebug_("Initialize final param with cmd line:", param_cmdline_, 2);
      finalParam = param_cmdline_;

      // 2. the instance values from the ini-file
      writeDebug_("Merging instance section into param:", param_instance_, 2);
      finalParam.merge(param_instance_);

      // 3. the tools data from the common section
      writeDebug_("Merging common section with tool name into param:", param_common_tool_, 2);
      finalParam.merge(param_common_tool_);

      // 4. everything else from the common section
      writeDebug_("Merging common section without tool name into param:", param_common_, 2);
      finalParam.merge(param_common_);


      finalParam.remove("ini"); // not contained in default params; remove to avoid "unknown param" in update()

      // finally: augment default values with INI/CLI values
      // note the copy(getIniLocation_(),..) as we want the param tree without instance
      // information
      param_ = this->getDefaultParameters_().copy(getIniLocation_(), true);
      if (!param_.update(finalParam, false, false, true, true, LOG_WARN))
      {
        LOG_ERROR << "Parameters passed to '" << this->tool_name_ << "' are invalid. To prevent usage of wrong defaults, please update/fix the parameters!" << std::endl;
        return ILLEGAL_PARAMETERS;
      }

      if (finalParam.exists("type"))
      {
        param_.setValue("type", finalParam.getValue("type"));
      }

      // check if all parameters are registered and have the correct type
      checkParam_(param_instance_, (String)value_ini, getIniLocation_());
      checkParam_(param_common_tool_, (String)value_ini, "common:" + tool_name_ + "::");
      checkParam_(param_common_, (String)value_ini, "common:");

      // check if the version of the parameters file matches the version of this tool
      // the parameters and values are all ok, but there might be more valid values now or new parameters which are currently not visible in the outdated INI
      String file_version = "";
      if (param_inifile_.exists(tool_name_ + ":version"))
      {
        file_version = param_inifile_.getValue(tool_name_ + ":version");
        if (file_version != version_)
        {
          writeLog_(String("Warning: Parameters file version (") + file_version + ") does not match the version of this tool (" + version_ + ").\n"
                    "Your current parameters are still valid, but there might be new valid values or even new parameters. Upgrading the INI might be useful.");
        }
      }
    }

    // 'test' flag is set
    if (getFlag_("test"))
    {
      test_mode_ = true;

      // initialize the random generator as early as possible!
      UniqueIdGenerator::setSeed(19991231235959);
    }

    // enable / disable collection of usage statistics by build variable
#ifdef ENABLE_UPDATE_CHECK
    // disable collection of usage statistics if environment variable is present
    char* disable_usage = getenv("OPENMS_DISABLE_UPDATE_CHECK");
 
    // only perform check if variable is not set or explicitly enabled by setting it to "OFF"  
    if (!test_mode_ && (disable_usage == NULL || strcmp(disable_usage, "OFF") == 0))
    {
      UpdateCheck::run(tool_name_, version_, debug_level_);
    }
#endif

    //-------------------------------------------------------------
    // determine and open the real log file
    //-------------------------------------------------------------

    {
      DataValue const& value_log = getParam_("log");
      if (!value_log.isEmpty())
      {
        writeDebug_("Log file: " + (String)value_log, 1);
        log_.close();
        log_.open(((String)value_log).c_str(), ofstream::out | ofstream::app);
        writeDebug_("Writing to '" + (String)value_log + '\'', 1);
      }
    }

    //-------------------------------------------------------------
    // debug level
    //-------------------------------------------------------------
    debug_level_ = getParamAsInt_("debug", 0);
    writeDebug_(String("Debug level (after ini file): ") + String(debug_level_), 1);
    if (debug_level_ > 0) Log_debug.insert(cout); // allows to use LOG_DEBUG << "something" << std::endl;

    //-------------------------------------------------------------
    //progress logging
    //-------------------------------------------------------------
    if (!getFlag_("no_progress"))
    {
      log_type_ = ProgressLogger::CMD;
    }

    //----------------------------------------------------------
    //threads
    //----------------------------------------------------------
    Int threads = getParamAsInt_("threads", 1);
    TOPPBase::setMaxNumberOfThreads(threads);

    //----------------------------------------------------------
    //main
    //----------------------------------------------------------
    StopWatch sw;
    sw.start();
    result = main_(argc, argv);
    sw.stop();
    LOG_INFO << this->tool_name_ << " took " << sw.toString() << "." << std::endl;

    // useful for benchmarking
    if (debug_level_ >= 1)
    {
      size_t mem_virtual(0);
      writeLog_(String("Peak Memory Usage: ") + (SysInfo::getProcessPeakMemoryConsumption(mem_virtual) ? String(mem_virtual / 1024) + " MB" : "<unknown>"));
    }


#ifndef DEBUG_TOPP
  }

  //----------------------------------------------------------
  //error handling
  //----------------------------------------------------------
  // Errors caused by the user
  catch (UnableToCreateFile& e)
  {
    writeLog_(String("Error: Unable to write file (") + e.what() + ")");
    writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ")!", 1);
    return CANNOT_WRITE_OUTPUT_FILE;
  }
  catch (FileNotFound& e)
  {
    writeLog_(String("Error: File not found (") + e.what() + ")");
    writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !", 1);
    return INPUT_FILE_NOT_FOUND;
  }
  catch (FileNotReadable& e)
  {
    writeLog_(String("Error: File not readable (") + e.what() + ")");
    writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !", 1);
    return INPUT_FILE_NOT_READABLE;
  }
  catch (FileEmpty& e)
  {
    writeLog_(String("Error: File empty (") + e.what() + ")");
    writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !", 1);
    return INPUT_FILE_EMPTY;
  }
  catch (ParseError& e)
  {
    writeLog_(String("Error: Unable to read file (") + e.what() + ")");
    writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !", 1);
    return INPUT_FILE_CORRUPT;
  }
  catch (RequiredParameterNotGiven& e)
  {
    String what = e.what();
    if (!what.hasPrefix("'"))
      what = "'" + what + "'";
    writeLog_(String("Error: The required parameter ") + what + " was not given or is empty!");
    writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !", 1);
    return MISSING_PARAMETERS;
  }
  catch (InvalidParameter& e)
  {
    writeLog_(String("Invalid parameter: ") + e.what());
    writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !", 1);
    return ILLEGAL_PARAMETERS;
  }
  // Internal errors because of wrong use of this class
  catch (UnregisteredParameter& e)
  {
    writeLog_(String("Internal error: Request for unregistered parameter '") + e.what() + "'");
    writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !", 1);
    return INTERNAL_ERROR;
  }
  catch (WrongParameterType& e)
  {
    writeLog_(String("Internal error: Request for parameter with wrong type '") + e.what() + "'");
    writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !", 1);
    return INTERNAL_ERROR;
  }
  // All other errors
  catch (BaseException& e)
  {
    writeLog_(String("Error: Unexpected internal error (") + e.what() + ")");
    writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !", 1);
    return UNKNOWN_ERROR;
  }
#endif

    log_.close();

    return result;
  }

  void TOPPBase::printUsage_()
  {
    // show advanced options?
    bool verbose = getFlag_("-helphelp");

    // common output
    cerr << "\n"
         << ConsoleUtils::breakString(tool_name_ + " -- " + tool_description_, 0, 10) << "\n"
         << "Version: " << verboseVersion_ << "\n"
         << "To cite OpenMS:\n  " << cite_openms_.toString() << "\n";
    if (!citations_.empty())
    {
      cerr << "To cite " << tool_name_ << ":\n";
      for (const Citation& c : citations_) cerr << "  " << c.toString() << "\n";
    }
    cerr << "\n";
    cerr << "Usage:" << "\n"
         << "  " << tool_name_ << " <options>" << "\n"
         << "\n";

    // print warning regarding not shown parameters
    if (!subsections_.empty() && !verbose)
      cerr << ConsoleUtils::breakString("This tool has algorithm parameters that are not shown here! Please check the ini file for a detailed description or use the --helphelp option.", 0, 10) + "\n\n";

    if (verbose)
    {
      // add all subsection parameters to the command line
      try
      {
        Param p = getSubsectionDefaults_();
        registerFullParam_(p);
      }
      catch (BaseException& /*e*/)
      {
        writeDebug_("Failed to add subsection parameters", 1);
      }
    }

    cerr << "Options (mandatory options marked with '*'):" << "\n";

    //determine max length of parameters (including argument) for indentation
    UInt max_size = 0;
    for (vector<ParameterInformation>::const_iterator it = parameters_.begin(); it != parameters_.end(); ++it)
    {
      if ((!it->advanced) || (it->advanced && verbose))
      {
        max_size = max((UInt)max_size, (UInt)(it->name.size() + it->argument.size() + it->required));
      }
    }

    //offset of the descriptions
    UInt offset = 6 + max_size;
    //keep track of the current subsection we are in, to display the subsection help when a new section starts
    String current_TOPP_subsection("");

    // PRINT parameters && description, restrictions and default
    for (vector<ParameterInformation>::const_iterator it = parameters_.begin(); it != parameters_.end(); ++it)
    {
      if (!((!it->advanced) || (it->advanced && verbose)))
        continue;

      //new subsection?
      String subsection = getSubsection_(it->name);
      if (!subsection.empty() && current_TOPP_subsection != subsection)
      {
        current_TOPP_subsection = subsection;
        map<String, String>::const_iterator subsec_it = subsections_TOPP_.find(current_TOPP_subsection);
        if (subsec_it == subsections_TOPP_.end())
        {
          throw ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "'" + current_TOPP_subsection + "' (TOPP subsection not registered)");
        }
        cerr << "\n"; // print newline for new subsection

        String subsection_description = subsec_it->second;
        if (subsection_description.length() == 0)
        {
          subsection_description = current_TOPP_subsection;
        }

        cerr << ConsoleUtils::breakString(subsection_description, 0, 10) << ":\n"; // print subsection description
      }
      else if (subsection.empty() && !current_TOPP_subsection.empty()) // subsection ended and normal parameters start again
      {
        current_TOPP_subsection = "";
        cerr << "\n"; // print newline to separate ending subsection
      }

      //NAME + ARGUMENT
      String str_tmp = "  -";
      str_tmp += it->name + " " + it->argument;
      if (it->required)
        str_tmp += '*';
      if (it->type == ParameterInformation::NEWLINE)
        str_tmp = "";

      //OFFSET
      str_tmp.fillRight(' ', offset);
      if (it->type == ParameterInformation::TEXT)
        str_tmp = "";

      //DESCRIPTION
      String desc_tmp = it->description;
      desc_tmp.firstToUpper();

      //DEFAULT
      StringList addons;
      switch (it->type)
      {
      case ParameterInformation::STRING:
      case ParameterInformation::DOUBLE:
      case ParameterInformation::INT:
      case ParameterInformation::STRINGLIST:
      case ParameterInformation::INTLIST:
      case ParameterInformation::DOUBLELIST:
      {
        String tmp_s = it->default_value.toString().substitute(", ", " ");
        if (tmp_s != "" && tmp_s != "[]")
        {
          addons.push_back(String("default: '") + tmp_s + "'");
        }
      }
      break;

      default:
        break;
      }

      //RESTRICTIONS
      switch (it->type)
      {
      case ParameterInformation::STRING:
      case ParameterInformation::INPUT_FILE:
      case ParameterInformation::OUTPUT_FILE:
      case ParameterInformation::STRINGLIST:
      case ParameterInformation::INPUT_FILE_LIST:
      case ParameterInformation::OUTPUT_FILE_LIST:
        if (it->valid_strings.size() != 0)
        {
          StringList copy = it->valid_strings;
          for (StringList::iterator str_it = copy.begin();
               str_it != copy.end(); ++str_it)
          {
            str_it->quote('\'');
          }

          String add = "";
          if (it->type == ParameterInformation::INPUT_FILE || it->type == ParameterInformation::OUTPUT_FILE ||
              it->type == ParameterInformation::INPUT_FILE_LIST || it->type == ParameterInformation::OUTPUT_FILE_LIST)
            add = " formats";

          addons.push_back(String("valid") + add + ": " + ListUtils::concatenate(copy, ", ")); // concatenate restrictions by comma
        }
        break;

      case ParameterInformation::INT:
      case ParameterInformation::INTLIST:
        if (it->min_int != -std::numeric_limits<Int>::max())
        {
          addons.push_back(String("min: '") + it->min_int + "'");
        }
        if (it->max_int != std::numeric_limits<Int>::max())
        {
          addons.push_back(String("max: '") + it->max_int + "'");
        }
        break;

      case ParameterInformation::DOUBLE:
      case ParameterInformation::DOUBLELIST:
        if (it->min_float != -std::numeric_limits<double>::max())
        {
          addons.push_back(String("min: '") + it->min_float + "'");
        }
        if (it->max_float != std::numeric_limits<double>::max())
        {
          addons.push_back(String("max: '") + it->max_float + "'");
        }
        break;

      default:
        break;
      }

      //add DEFAULT and RESTRICTIONS
      if (addons.size() != 0)
      {
        desc_tmp += String(" (") + ListUtils::concatenate(addons, " ") + ")";
      }

      if (it->type == ParameterInformation::TEXT)
        cerr << ConsoleUtils::breakString(str_tmp + desc_tmp, 0, 10); // no indentation for text
      else
        cerr << ConsoleUtils::breakString(str_tmp + desc_tmp, offset, 10);
      cerr << "\n";
    }

    // SUBSECTION's at the end
    if (subsections_.size() != 0 && !verbose)
    {
      //determine indentation of description
      UInt indent = 0;
      for (map<String, String>::const_iterator it = subsections_.begin(); it != subsections_.end(); ++it)
      {
        indent = max((UInt)it->first.size(), indent);
      }
      indent += 6;

      //output
      cerr << "\n"
           << "The following configuration subsections are valid:" << "\n";
      for (map<String, String>::const_iterator it = subsections_.begin(); it != subsections_.end(); ++it)
      {
        String tmp = String(" - ") + it->first;
        tmp.fillRight(' ', indent);
        cerr << ConsoleUtils::breakString(tmp  + it->second, indent, 10);
        cerr << "\n";
      }
      cerr << "\n"
           << ConsoleUtils::breakString("You can write an example INI file using the '-write_ini' option.", 0, 10) << "\n"
           << ConsoleUtils::breakString("Documentation of subsection parameters can be found in the doxygen documentation or the INIFileEditor.", 0, 10) << "\n"
           << ConsoleUtils::breakString("Have a look at the OpenMS documentation for more information.", 0, 10) << "\n";
    }
    cerr << endl;
  }

  ParameterInformation TOPPBase::paramEntryToParameterInformation_(const Param::ParamEntry& entry, const String& argument, const String& full_name) const
  {
    String name = full_name.empty() ? entry.name : full_name;
    bool advanced = entry.tags.count("advanced");
    // special case for flags:
    if ((entry.value.valueType() == DataValue::STRING_VALUE) &&
        (entry.value == "false") && (entry.valid_strings.size() == 2) &&
        (entry.valid_strings[0] == "true") && (entry.valid_strings[1] == "false"))
    {
      return ParameterInformation(name, ParameterInformation::FLAG, "", "", entry.description, false, advanced);
    }

    bool input_file = entry.tags.count("input file");
    bool output_file = entry.tags.count("output file");
    if (input_file && output_file)
    {
      throw InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Parameter '" + full_name + "' marked as both input and output file");
    }
    enum ParameterInformation::ParameterTypes type = ParameterInformation::NONE;
    switch (entry.value.valueType())
    {
    case DataValue::STRING_VALUE:
      if (input_file)
        type = ParameterInformation::INPUT_FILE;
      else if (output_file)
        type = ParameterInformation::OUTPUT_FILE;
      else
        type = ParameterInformation::STRING;
      break;

    case DataValue::INT_VALUE:
      type = ParameterInformation::INT;
      break;

    case DataValue::DOUBLE_VALUE:
      type = ParameterInformation::DOUBLE;
      break;

    case DataValue::STRING_LIST:
      if (input_file)
        type = ParameterInformation::INPUT_FILE_LIST;
      else if (output_file)
        type = ParameterInformation::OUTPUT_FILE_LIST;
      else
        type = ParameterInformation::STRINGLIST;
      break;

    case DataValue::INT_LIST:
      type = ParameterInformation::INTLIST;
      break;

    case DataValue::DOUBLE_LIST:
      type = ParameterInformation::DOUBLELIST;
      break;

    case DataValue::EMPTY_VALUE:
      type = ParameterInformation::NONE;
      break;
    }
    bool required = entry.tags.count("required");
    ParameterInformation param(name, type, argument, entry.value, entry.description, required, advanced);
    param.valid_strings = entry.valid_strings;
    // here, we rely on the fact that defaults (meaning "not set") are the same for both:
    param.min_int = entry.min_int;
    param.max_int = entry.max_int;
    param.min_float = entry.min_float;
    param.max_float = entry.max_float;
    return param;
  }

  String TOPPBase::getParamArgument_(const Param::ParamEntry& entry) const
  {
    String argument = "";
    switch (entry.value.valueType())
    {
    case DataValue::STRING_VALUE:
      if (entry.valid_strings.empty())
        argument = "<text>"; // name?
      else
        argument = "<choice>";
      break;

    case DataValue::INT_VALUE:
      argument = "<number>"; // integer?
      break;

    case DataValue::DOUBLE_VALUE:
      argument = "<value>"; // float?
      break;

    case DataValue::STRING_LIST:
      argument = "<list>";
      break;

    case DataValue::INT_LIST:
      argument = "<numbers>";
      break;

    case DataValue::DOUBLE_LIST:
      argument = "<values>";
      break;

    case DataValue::EMPTY_VALUE:
      argument = "";
      break;
    }
    return argument;
  }

  std::vector<ParameterInformation> TOPPBase::paramToParameterInformation_(const Param& param) const
  {
    std::vector<ParameterInformation> parameter_information;
    for (Param::ParamIterator it = param.begin(); it != param.end(); ++it)
    {
      String full_name = it.getName();
      // make up a value for "argument":
      String argument = getParamArgument_(*it);
      // transform to ParameterInformation and register
      parameter_information.push_back(paramEntryToParameterInformation_(*it, argument, full_name));
    }
    return parameter_information;
  }

  void TOPPBase::registerParamSubsectionsAsTOPPSubsections_(const Param& param)
  {
    for (Param::ParamIterator it = param.begin(); it != param.end(); ++it)
    {
      String full_name = it.getName();
      String subsection = getSubsection_(full_name);
      if (!subsection.empty() && (subsections_TOPP_.count(subsection) == 0))
      {
        subsections_TOPP_[subsection] = param.getSectionDescription(subsection);
      }
    }
  }

  void TOPPBase::registerFullParam_(const Param& param)
  {
    // register subsections
    registerParamSubsectionsAsTOPPSubsections_(param);

    // add the actual parameters
    std::vector<ParameterInformation> parameter_information = paramToParameterInformation_(param);
    parameters_.insert(parameters_.end(), parameter_information.begin(), parameter_information.end());
  }

  void TOPPBase::registerStringOption_(const String& name, const String& argument, const String& default_value, const String& description, bool required, bool advanced)
  {
    if (required && default_value != "")
      throw InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Registering a required StringOption param (" + name + ") with a non-empty default is forbidden!", default_value);
    parameters_.push_back(ParameterInformation(name, ParameterInformation::STRING, argument, default_value, description, required, advanced));
  }

  ParameterInformation& TOPPBase::getParameterByName_(const String& name)
  {
    typedef std::vector<ParameterInformation>::iterator TParamInfoIterator;
    //search the right parameter
    for (TParamInfoIterator it = parameters_.begin(); it != parameters_.end(); ++it)
    {
      if (it->name == name)
        return *it;
    }

    //parameter not found
    throw UnregisteredParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
  }

  void TOPPBase::setValidStrings_(const String& name, const std::string vstrings[], int count)
  {
    std::vector<String> vec;
    vec.assign(vstrings, vstrings + count);
    setValidStrings_(name, vec);
  }

  void TOPPBase::setValidStrings_(const String& name, const std::vector<String>& strings)
  {
    //check for commas
    for (Size i = 0; i < strings.size(); ++i)
    {
      if (strings[i].has(','))
      {
        throw InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Comma characters in Param string restrictions are not allowed!");
      }
    }

    // get the matching parameter
    ParameterInformation& p = getParameterByName_(name);

    //check if the type matches
    if (p.type != ParameterInformation::STRING && p.type != ParameterInformation::STRINGLIST)
    {
      throw ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }

    StringList valids = strings;
    StringList defaults;

    if (p.type == ParameterInformation::STRING)
      defaults.push_back(String(p.default_value));
    else
      defaults = p.default_value;

    for (Size j = 0; j < defaults.size(); ++j) // allow the empty string even if not in restrictions
    {
      if (defaults[j].size() > 0 && !ListUtils::contains(valids, defaults[j]))
      {
        throw InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "TO THE DEVELOPER: The TOPP/UTILS tool option '" + name + "' with default value " + String(p.default_value) + " does not meet restrictions!");
      }
    }

    p.valid_strings = strings;
  }

  void TOPPBase::setValidFormats_(const String& name, const std::vector<String>& formats, const bool force_OpenMS_format)
  {
    //check if formats are known
    if (force_OpenMS_format)
    {
      for (Size i = 0; i < formats.size(); ++i)
      {
        if (formats[i] != "fid")
        {
          if (FileHandler::getTypeByFileName(String(".") + formats[i]) == FileTypes::UNKNOWN)
          {
            throw InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The file format '" + formats[i] + "' is invalid!");
          }
        }
      }
    }

    ParameterInformation& p = getParameterByName_(name);

    //check if the type matches
    if (p.type != ParameterInformation::INPUT_FILE
       && p.type != ParameterInformation::OUTPUT_FILE
       && p.type != ParameterInformation::INPUT_FILE_LIST
       && p.type != ParameterInformation::OUTPUT_FILE_LIST)
    {
      throw ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }

    if (p.valid_strings.size() > 0)
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Internal error: Valid formats are already set for '" + name + "'. Please check for typos!");
    }
    p.valid_strings = formats;
  }

  void TOPPBase::setMinInt_(const String& name, Int min)
  {
    ParameterInformation& p = getParameterByName_(name);

    //check if the type matches
    if (p.type != ParameterInformation::INT && p.type != ParameterInformation::INTLIST)
    {
      throw ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }

    IntList defaults;
    if (p.type == ParameterInformation::INT)
      defaults.push_back(Int(p.default_value));
    else
      defaults = p.default_value;
    for (Size j = 0; j < defaults.size(); ++j)
    {
      if (defaults[j] < min)
      {
        throw InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "TO THE DEVELOPER: The TOPP/UTILS tool option '" + name + "' with default value " + String(p.default_value) + " does not meet restrictions!");
      }
    }
    p.min_int = min;
  }

  void TOPPBase::setMaxInt_(const String& name, Int max)
  {
    ParameterInformation& p = getParameterByName_(name);

    //check if the type matches
    if (p.type != ParameterInformation::INT && p.type != ParameterInformation::INTLIST)
    {
      throw ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }
    IntList defaults;
    if (p.type == ParameterInformation::INT)
      defaults.push_back(Int(p.default_value));
    else
      defaults = p.default_value;
    for (Size j = 0; j < defaults.size(); ++j)
    {
      if (defaults[j] > max)
      {
        throw InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "TO THE DEVELOPER: The TOPP/UTILS tool option '" + name + "' with default value " + String(p.default_value) + " does not meet restrictions!");
      }
    }
    p.max_int = max;
  }

  void TOPPBase::setMinFloat_(const String& name, double min)
  {
    ParameterInformation& p = getParameterByName_(name);

    //check if the type matches
    if (p.type != ParameterInformation::DOUBLE && p.type != ParameterInformation::DOUBLELIST)
    {
      throw ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }
    DoubleList defaults;
    if (p.type == ParameterInformation::DOUBLE)
      defaults.push_back(double(p.default_value));
    else
      defaults = p.default_value;
    for (Size j = 0; j < defaults.size(); ++j)
    {
      if (defaults[j] < min)
      {
        throw InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "TO THE DEVELOPER: The TOPP/UTILS tool option '" + name + "' with default value " + String(p.default_value) + " does not meet restrictions!");
      }
    }
    p.min_float = min;
  }

  void TOPPBase::setMaxFloat_(const String& name, double max)
  {
    ParameterInformation& p = getParameterByName_(name);

    //check if the type matches
    if (p.type != ParameterInformation::DOUBLE && p.type != ParameterInformation::DOUBLELIST)
    {
      throw ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }
    DoubleList defaults;
    if (p.type == ParameterInformation::DOUBLE)
      defaults.push_back(double(p.default_value));
    else
      defaults = p.default_value;
    for (Size j = 0; j < defaults.size(); ++j)
    {
      if (defaults[j] > max)
      {
        throw InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "TO THE DEVELOPER: The TOPP/UTILS tool option '" + name + "' with default value " + String(p.default_value) + " does not meet restrictions!");
      }
    }
    p.max_float = max;
  }

  void TOPPBase::registerInputFile_(const String& name, const String& argument, const String& default_value, const String& description, bool required, bool advanced, const StringList& tags)
  {
    if (required && default_value != "" && !ListUtils::contains(tags, "skipexists"))
      throw InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Registering a required InputFile param (" + name + ") with a non-empty default is forbidden!", default_value);
    parameters_.push_back(ParameterInformation(name, ParameterInformation::INPUT_FILE, argument, default_value, description, required, advanced, tags));
  }

  void TOPPBase::registerOutputFile_(const String& name, const String& argument, const String& default_value, const String& description, bool required, bool advanced)
  {
    if (required && default_value != "")
      throw InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Registering a required OutputFile param (" + name + ") with a non-empty default is forbidden!", default_value);
    parameters_.push_back(ParameterInformation(name, ParameterInformation::OUTPUT_FILE, argument, default_value, description, required, advanced));
  }

  void TOPPBase::registerDoubleOption_(const String& name, const String& argument, double default_value, const String& description, bool required, bool advanced)
  {
    if (required)
    {
      throw InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Registering a double param (" + name + ") as 'required' is forbidden (there is no value to indicate it is missing)!", String(default_value));
    }
    parameters_.push_back(ParameterInformation(name, ParameterInformation::DOUBLE, argument, default_value, description, required, advanced));
  }

  void TOPPBase::registerIntOption_(const String& name, const String& argument, Int default_value, const String& description, bool required, bool advanced)
  {
    if (required)
    {
      throw InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Registering an Int param (" + name + ") as 'required' is forbidden (there is no value to indicate it is missing)!", String(default_value));
    }
    parameters_.push_back(ParameterInformation(name, ParameterInformation::INT, argument, default_value, description, required, advanced));
  }

  void TOPPBase::registerOutputFileList_(const String& name, const String& argument, StringList default_value, const String& description, bool required, bool advanced)
  {
    if (required && default_value.size() > 0)
      throw InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Registering a required OutputFileList param (" + name + ") with a non-empty default is forbidden!", ListUtils::concatenate(default_value, ","));
    parameters_.push_back(ParameterInformation(name, ParameterInformation::OUTPUT_FILE_LIST, argument, default_value, description, required, advanced));
  }

  void TOPPBase::registerInputFileList_(const String& name, const String& argument, StringList default_value, const String& description, bool required, bool advanced, const StringList& tags)
  {
    if (required && default_value.size() > 0 && !ListUtils::contains(tags, "skipexists"))
      throw InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Registering a required InputFileList param (" + name + ") with a non-empty default is forbidden!", ListUtils::concatenate(default_value, ","));
    parameters_.push_back(ParameterInformation(name, ParameterInformation::INPUT_FILE_LIST, argument, default_value, description, required, advanced, tags));
  }

  void TOPPBase::registerStringList_(const String& name, const String& argument, StringList default_value, const String& description, bool required, bool advanced)
  {
    if (required && default_value.size() > 0)
      throw InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Registering a required StringList param (" + name + ") with a non-empty default is forbidden!", ListUtils::concatenate(default_value, ","));
    parameters_.push_back(ParameterInformation(name, ParameterInformation::STRINGLIST, argument, default_value, description, required, advanced));
  }

  void TOPPBase::registerIntList_(const String& name, const String& argument, IntList default_value, const String& description, bool required, bool advanced)
  {
    stringstream ss;
    ss << default_value;
    if (required && default_value.size() > 0)
      throw InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Registering a required IntList param (" + name + ") with a non-empty default is forbidden!", String(ss.str()));
    parameters_.push_back(ParameterInformation(name, ParameterInformation::INTLIST, argument, default_value, description, required, advanced));
  }

  void TOPPBase::registerDoubleList_(const String& name, const String& argument, DoubleList default_value, const String& description, bool required, bool advanced)
  {
    stringstream ss;
    ss << default_value;
    if (required && default_value.size() > 0)
      throw InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Registering a required DoubleList param (" + name + ") with a non-empty default is forbidden!", String(ss.str()));
    parameters_.push_back(ParameterInformation(name, ParameterInformation::DOUBLELIST, argument, default_value, description, required, advanced));
  }

  void TOPPBase::registerFlag_(const String& name, const String& description, bool advanced)
  {
    parameters_.push_back(ParameterInformation(name, ParameterInformation::FLAG, "", "", description, false, advanced));
  }

  void TOPPBase::addEmptyLine_()
  {
    parameters_.push_back(ParameterInformation("", ParameterInformation::NEWLINE, "", "", "", false, false));
  }

  void TOPPBase::addText_(const String& text)
  {
    parameters_.push_back(ParameterInformation("", ParameterInformation::TEXT, "", "", text, false, false));
  }

  const ParameterInformation& TOPPBase::findEntry_(const String& name) const
  {
    vector<ParameterInformation>::const_iterator it = parameters_.begin();
    while (it != parameters_.end() && it->name != name)
    {
      ++it;
    }
    if (it == parameters_.end())
    {
      throw UnregisteredParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }
    return *it;
  }

  String TOPPBase::getStringOption_(const String& name) const
  {
    const ParameterInformation& p = findEntry_(name);
    if (p.type != ParameterInformation::STRING && p.type != ParameterInformation::INPUT_FILE && p.type != ParameterInformation::OUTPUT_FILE)
    {
      throw WrongParameterType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }
    if (p.required && (getParam_(name).isEmpty() || getParam_(name) == ""))
    {
      String message = "'" + name + "'";
      if (p.valid_strings.size() > 0)
      {
        message += " [valid: " + ListUtils::concatenate(p.valid_strings, ", ") + "]";
      }
      throw RequiredParameterNotGiven(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message);
    }
    String tmp = getParamAsString_(name, p.default_value);
    writeDebug_(String("Value of string option '") + name + "': " + tmp, 1);

    // if required or set by user, do some validity checks
    if (p.required || (!getParam_(name).isEmpty() && (tmp != p.default_value) &&
                       !tmp.empty()))
    {
      // check if files are readable/writable
      if (p.type == ParameterInformation::INPUT_FILE)
      {
        if (!ListUtils::contains(p.tags, "skipexists"))
        {
          inputFileReadable_(tmp, name);
        }
      }
      else if (p.type == ParameterInformation::OUTPUT_FILE)
      {
        outputFileWritable_(tmp, name);
      }

      // check restrictions
      if (p.valid_strings.size() != 0)
      {
        if (p.type == ParameterInformation::STRING)
        {
          if (find(p.valid_strings.begin(), p.valid_strings.end(), tmp) == p.valid_strings.end())
          {
            String valid_strings = ListUtils::concatenate(p.valid_strings, "', '");
            throw InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Invalid value '") + tmp + "' for string parameter '" + name + "' given. Valid strings are: '" + valid_strings + "'.");
          }
        }
        else if (p.type == ParameterInformation::INPUT_FILE)
        {
          //create upper case list of valid formats
          StringList formats = p.valid_strings;
          StringListUtils::toUpper(formats);
          //determine file type as string
          String format = FileTypes::typeToName(FileHandler::getTypeByFileName(tmp)).toUpper();
          bool invalid = false;
          //Wrong or unknown ending
          if (!ListUtils::contains(formats, format))
          {
            if (format == "UNKNOWN") //Unknown ending => check content
            {
              format = FileTypes::typeToName(FileHandler::getTypeByContent(tmp)).toUpper();
              if (!ListUtils::contains(formats, format))
              {
                if (format == "UNKNOWN") //Unknown format => warning as this might by the wrong format
                {
                  writeLog_("Warning: Could not determine format of input file '" + tmp + "'!");
                }
                else //Wrong ending => invalid
                {
                  invalid = true;
                }
              }
            }
            else //Wrong ending => invalid
            {
              invalid = true;
            }
          }
          if (invalid)
          {
            String valid_formats = "";
            valid_formats.concatenate(p.valid_strings.begin(), p.valid_strings.end(), "','");
            throw InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Input file '" + tmp + "' has invalid format '") + format + "'. Valid formats are: '" + valid_formats + "'.");
          }
        }
        else if (p.type == ParameterInformation::OUTPUT_FILE)
        {
          outputFileWritable_(tmp, name);

          //create upper case list of valid formats
          StringList formats = p.valid_strings;
          StringListUtils::toUpper(formats);
          //determine file type as string
          String format = FileTypes::typeToName(FileHandler::getTypeByFileName(tmp)).toUpper();
          //Wrong or unknown ending
          if (!ListUtils::contains(formats, format) && format != "UNKNOWN")
          {
            String valid_formats = "";
            valid_formats.concatenate(p.valid_strings.begin(), p.valid_strings.end(), "','");
            throw InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Invalid output file extension '") + tmp + "'. Valid file extensions are: '" + valid_formats + "'.");
          }
        }
      }
    }

    return tmp;
  }

  double TOPPBase::getDoubleOption_(const String& name) const
  {
    const ParameterInformation& p = findEntry_(name);
    if (p.type != ParameterInformation::DOUBLE)
    {
      throw WrongParameterType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }
    if (p.required && getParam_(name).isEmpty())
    {
      throw RequiredParameterNotGiven(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }
    double tmp = getParamAsDouble_(name, (double)p.default_value);
    if (p.required && boost::math::isnan(tmp))
    {
      throw RequiredParameterNotGiven(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }
    writeDebug_(String("Value of double option '") + name + "': " + String(tmp), 1);

    //check if in valid range
    if (p.required || (!getParam_(name).isEmpty() && tmp != (double)p.default_value))
    {
      if (tmp < p.min_float || tmp > p.max_float)
      {
        throw InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Invalid value '") + tmp + "' for float parameter '" + name + "' given. Out of valid range: '" + p.min_float + "'-'" + p.max_float + "'.");
      }
    }

    return tmp;
  }

  Int TOPPBase::getIntOption_(const String& name) const
  {
    const ParameterInformation& p = findEntry_(name);
    if (p.type != ParameterInformation::INT)
    {
      throw WrongParameterType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }
    if (p.required && getParam_(name).isEmpty())
    {
      throw RequiredParameterNotGiven(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }
    Int tmp = getParamAsInt_(name, (Int)p.default_value);
    // not checking if NAN here (as done with double, as NAN is not supported for Int)
    writeDebug_(String("Value of int option '") + name + "': " + String(tmp), 1);

    //check if in valid range
    if (p.required || (!getParam_(name).isEmpty() && tmp != (Int)p.default_value))
    {
      if (tmp < p.min_int || tmp > p.max_int)
      {
        throw InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Invalid value '") + tmp + "' for integer parameter '" + name + "' given. Out of valid range: '" + p.min_int + "'-'" + p.max_int + "'.");
      }
    }

    return tmp;
  }

  StringList TOPPBase::getStringList_(const String& name) const
  {
    const ParameterInformation& p = findEntry_(name);
    if (p.type != ParameterInformation::STRINGLIST && p.type != ParameterInformation::INPUT_FILE_LIST && p.type != ParameterInformation::OUTPUT_FILE_LIST)
    {
      throw WrongParameterType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }
    if (p.required && getParam_(name).isEmpty())
    {
      throw RequiredParameterNotGiven(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }
    StringList tmp_list = getParamAsStringList_(name, p.default_value);
    if (p.required && tmp_list.size() == 0)
    {
      throw RequiredParameterNotGiven(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }

    for (StringList::iterator it = tmp_list.begin(); it < tmp_list.end(); ++it)
    {
      String tmp(*it);
      writeDebug_(String("Value of string option '") + name + "': " + tmp, 1);

      // if required or set by user, do some validity checks
      if (p.required || (!getParam_(name).isEmpty() && tmp_list != p.default_value))
      {
        //check if files are readable/writable
        if (p.type == ParameterInformation::INPUT_FILE_LIST)
        {
          if (!ListUtils::contains(p.tags, "skipexists")) inputFileReadable_(tmp, name);
        }
        else if (p.type == ParameterInformation::OUTPUT_FILE_LIST)
        {
          outputFileWritable_(tmp, name);
        }

        //check restrictions
        if (p.valid_strings.size() != 0)
        {
          if (p.type == ParameterInformation::STRINGLIST)
          {
            if (find(p.valid_strings.begin(), p.valid_strings.end(), tmp) == p.valid_strings.end())
            {
              String valid_strings = "";
              valid_strings.concatenate(p.valid_strings.begin(), p.valid_strings.end(), "','");
              throw InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Invalid value '") + tmp + "' for string parameter '" + name + "' given. Valid strings are: '" + valid_strings + "'.");
            }
          }
          else if (p.type == ParameterInformation::INPUT_FILE_LIST)
          {
            //create upper case list of valid formats
            StringList formats = p.valid_strings;
            StringListUtils::toUpper(formats);
            //determine file type as string
            String format = FileTypes::typeToName(FileHandler::getTypeByFileName(tmp)).toUpper();
            bool invalid = false;
            //Wrong or unknown ending
            if (!ListUtils::contains(formats, format))
            {
              if (format == "UNKNOWN") //Unknown ending => check content
              {
                format = FileTypes::typeToName(FileHandler::getTypeByContent(tmp)).toUpper();
                if (!ListUtils::contains(formats, format))
                {
                  if (format == "UNKNOWN") //Unknown format => warning as this might by the wrong format
                  {
                    writeLog_("Warning: Could not determine format of input file '" + tmp + "'!");
                  }
                  else //Wrong ending => invalid
                  {
                    invalid = true;
                  }
                }
              }
              else //Wrong ending => invalid
              {
                invalid = true;
              }
            }
            if (invalid)
            {
              String valid_formats = "";
              valid_formats.concatenate(p.valid_strings.begin(), p.valid_strings.end(), "','");
              throw InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Input file '" + tmp + "' has invalid format '") + format + "'. Valid formats are: '" + valid_formats + "'.");
            }
          }
          else if (p.type == ParameterInformation::OUTPUT_FILE_LIST)
          {
            outputFileWritable_(tmp, name);

            //create upper case list of valid formats
            StringList formats = p.valid_strings;
            StringListUtils::toUpper(formats);
            //determine file type as string
            String format = FileTypes::typeToName(FileHandler::getTypeByFileName(tmp)).toUpper();
            //Wrong or unknown ending
            if (!ListUtils::contains(formats, format) && format != "UNKNOWN")
            {
              String valid_formats = "";
              valid_formats.concatenate(p.valid_strings.begin(), p.valid_strings.end(), "','");
              throw InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Invalid output file extension '") + tmp + "'. Valid file extensions are: '" + valid_formats + "'.");
            }
          }
        }
      }
    }
    return tmp_list;
  }

  DoubleList TOPPBase::getDoubleList_(const String& name) const
  {
    const ParameterInformation& p = findEntry_(name);
    if (p.type != ParameterInformation::DOUBLELIST)
    {
      throw WrongParameterType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }
    if (p.required && getParam_(name).isEmpty())
    {
      throw RequiredParameterNotGiven(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }
    DoubleList tmp_list = getParamAsDoubleList_(name, p.default_value);
    if (p.required && tmp_list.size() == 0)
    {
      throw RequiredParameterNotGiven(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }

    for (DoubleList::iterator it = tmp_list.begin(); it < tmp_list.end(); ++it)
    {
      double tmp = *it;
      writeDebug_(String("Value of string option '") + name + "': " + String(tmp), 1);

      //check if in valid range
      if (p.required || (!getParam_(name).isEmpty() && tmp_list != p.default_value))
      {
        if (tmp < p.min_float || tmp > p.max_float)
        {
          throw InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Invalid value '") + tmp + "' for float parameter '" + name + "' given. Out of valid range: '" + p.min_float + "'-'" + p.max_float + "'.");
        }
      }
    }
    return tmp_list;
  }

  IntList TOPPBase::getIntList_(const String& name) const
  {
    const ParameterInformation& p = findEntry_(name);
    if (p.type != ParameterInformation::INTLIST)
    {
      throw WrongParameterType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }
    if (p.required && getParam_(name).isEmpty())
    {
      throw RequiredParameterNotGiven(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }
    IntList tmp_list = getParamAsIntList_(name, p.default_value);
    if (p.required && tmp_list.size() == 0)
    {
      throw RequiredParameterNotGiven(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }

    Int tmp;
    for (IntList::iterator it = tmp_list.begin(); it < tmp_list.end(); ++it)
    {
      tmp = *it;
      writeDebug_(String("Value of string option '") + name + "': " + String(tmp), 1);

      //check if in valid range
      if (p.required || (!getParam_(name).isEmpty() && tmp_list != p.default_value))
      {
        if (tmp < p.min_int || tmp > p.max_int)
        {
          throw InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Invalid value '") + tmp + "' for integer parameter '" + name + "' given. Out of valid range: '" + p.min_int + "'-'" + p.max_int + "'.");
        }
      }
    }
    return tmp_list;
  }

  bool TOPPBase::getFlag_(const String& name) const
  {
    const ParameterInformation& p = findEntry_(name);
    if (p.type != ParameterInformation::FLAG)
    {
      throw WrongParameterType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }
    bool tmp = getParamAsBool_(name);
    writeDebug_(String("Value of string option '") + name + "': " + String(tmp), 1);
    return tmp;
  }

  void TOPPBase::writeLog_(const String& text) const
  {
    LOG_INFO << text << endl;
    enableLogging_();
    log_ << QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss").toStdString() << ' ' << getIniLocation_() << ": " << text << endl;
  }

  void TOPPBase::writeDebug_(const String& text, UInt min_level) const
  {
    if (debug_level_ >= (Int)min_level)
    {
      writeLog_(text);
    }
  }

  void TOPPBase::writeDebug_(const String& text, const Param& param, UInt min_level) const
  {
    if (debug_level_ >= (Int)min_level)
    {
      LOG_DEBUG << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl
                << QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss").toStdString() << ' ' << getIniLocation_() << " " << text << endl
                << param
                << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl;
      enableLogging_();
      log_ << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl
           << QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss").toStdString() << ' ' << getIniLocation_() << " " << text << endl
           << param
           << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl;
    }
  }

  String TOPPBase::makeTempDirectory_() const
  {
    String temp_dir = QDir::toNativeSeparators((File::getTempDirectory() + "/" + File::getUniqueName() + "/").toQString());
    writeDebug_("Creating temporary directory '" + temp_dir + "'", 1);
    QDir d;
    d.mkpath(temp_dir.toQString());
    return temp_dir;
  }

  void TOPPBase::removeTempDirectory_(const String& temp_dir, Int keep_debug) const
  {
    if (temp_dir.empty()) return; // no temp. dir. created

    if ((keep_debug > 0) && (debug_level_ >= keep_debug))
    {
      writeDebug_("Keeping temporary files in directory '" + temp_dir + "'. Set debug level to " + String(keep_debug) + " or lower to remove them.", keep_debug);
    }
    else
    {
      if ((keep_debug > 0) && (debug_level_ > 0) && (debug_level_ < keep_debug))
      {
        writeDebug_("Deleting temporary directory '" + temp_dir + "'. Set debug level to " + String(keep_debug) + " or higher to keep it.", debug_level_);
      }
      File::removeDirRecursively(temp_dir);
    }
  }

  String TOPPBase::getParamAsString_(const String& key, const String& default_value) const
  {
    const DataValue& tmp = getParam_(key);
    if (!tmp.isEmpty())
    {
      return tmp.toString();
    }
    else
    {
      return default_value;
    }
  }

  Int TOPPBase::getParamAsInt_(const String& key, Int default_value) const
  {
    const DataValue& tmp = getParam_(key);
    if (!tmp.isEmpty())
    {
      if (tmp.valueType() == DataValue::INT_VALUE)
      {
        return (Int)tmp;
      }
      throw WrongParameterType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, key);
    }
    else
    {
      return default_value;
    }
  }

  double TOPPBase::getParamAsDouble_(const String& key, double default_value) const
  {
    const DataValue& tmp = getParam_(key);
    if (!tmp.isEmpty())
    {
      if (tmp.valueType() == DataValue::DOUBLE_VALUE)
      {
        return (double)tmp;
      }
      throw WrongParameterType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, key);
    }
    else
    {
      return default_value;
    }
  }

  StringList TOPPBase::getParamAsStringList_(const String& key, const StringList& default_value) const
  {
    const DataValue& tmp = getParam_(key);
    if (!tmp.isEmpty())
    {
      return tmp;
    }
    else
    {
      return default_value;
    }
  }

  IntList TOPPBase::getParamAsIntList_(const String& key, const IntList& default_value) const
  {
    const DataValue& tmp = getParam_(key);
    if (!tmp.isEmpty())
    {
      if (tmp.valueType() == DataValue::INT_LIST)
      {
        return tmp;
      }
      throw WrongParameterType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, key);
    }
    else
    {
      return default_value;
    }
  }

  DoubleList TOPPBase::getParamAsDoubleList_(const String& key, const DoubleList& default_value) const
  {
    const DataValue& tmp = getParam_(key);
    if (!tmp.isEmpty())
    {
      if (tmp.valueType() == DataValue::DOUBLE_LIST)
      {
        return tmp;
      }
      throw WrongParameterType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, key);
    }
    else
    {
      return default_value;
    }
  }

  bool TOPPBase::getParamAsBool_(const String& key) const
  {
    DataValue tmp = getParam_(key);
    if (tmp.valueType() == DataValue::EMPTY_VALUE)
    {
      return false;
    }
    else if (tmp.valueType() == DataValue::STRING_VALUE)
    {
      if ((String)tmp == "false")
      {
        return false;
      }
      else if ((String)tmp == "true")
      {
        return true;
      }
    }
    throw InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Invalid value '") + tmp.toString() + "' for flag parameter '" + key + "'. Valid values are 'true' and 'false' only.");
  }

  DataValue const& TOPPBase::getParam_(const String& key) const
  {
    if (param_.exists(key))
    {
      return param_.getValue(key);
    }
    else
    {
      // if look up fails everywhere, return EMPTY
      writeDebug_(String("Parameter '") + key + String("' not found."), 1);
      return DataValue::EMPTY;
    }
  }

  Param const& TOPPBase::getParam_() const
  {
    return param_;
  }

  String TOPPBase::getSubsection_(const String& name) const
  {
    size_t pos = name.find_last_of(':');
    if (pos == std::string::npos)
      return ""; // delimiter not found

    return name.substr(0, pos);
  }

  void TOPPBase::enableLogging_() const
  {
    if (!log_.is_open())
    {
      String log_destination = "";
      if (param_cmdline_.exists("log"))
        log_destination = param_cmdline_.getValue("log");
      if (log_destination != "")
      {
        log_.open(log_destination.c_str(), ofstream::out | ofstream::app);
        if (debug_level_ >= 1)
        {
          cout << "Writing to '" << log_destination << '\'' << "\n";
          log_ << QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss").toStdString() << ' ' << getIniLocation_() << ": " << "Writing to '" << log_destination << '\'' <<  "\n";
        }
      }
    }
    return;
  }

  void TOPPBase::checkParam_(const Param& param, const String& filename, const String& location) const
  {
    //cout << endl << "--"<< location<< "--" << endl << param << endl << endl;
    for (Param::ParamIterator it = param.begin(); it != param.end(); ++it)
    {
      // subsections (do not check content, but warn if not registered)
      String subsection = getSubsection_(it.getName());
      if (!subsection.empty() && subsections_TOPP_.count(subsection) == 0) // not found in TOPP subsections
      {
        // for multi-level subsections, check only the first level:
        if (subsections_.count(subsection.substr(0, subsection.find(':'))) == 0) // not found in normal subsections
        {
          if (!(location == "common::" && subsection == tool_name_))
          {
            writeLog_("Warning: Unknown subsection '" + subsection + "' in '" + filename + "' (location '" + location + "')!");
          }
        }
        continue;
      }
      // normal parameter: check its value type
      // if no such parameter is registered an exception is thrown
      try
      {
        //check type
        switch (findEntry_(it.getName()).type)
        {
        case ParameterInformation::STRING:
        case ParameterInformation::INPUT_FILE:
        case ParameterInformation::OUTPUT_FILE:
        case ParameterInformation::FLAG:
          if (it->value.valueType() != DataValue::STRING_VALUE)
          {
            writeLog_("Warning: Wrong parameter type of '" + location + it.getName() + "' in '" + filename + "'. Type should be 'string'!");
          }
          break;

        case ParameterInformation::DOUBLE:
          if (it->value.valueType() != DataValue::DOUBLE_VALUE)
          {
            writeLog_("Warning: Wrong  parameter type of '" + location + it.getName() + "' in '" + filename + "'. Type should be 'double'!");
          }
          break;

        case ParameterInformation::INT:
          if (it->value.valueType() != DataValue::INT_VALUE)
          {
            writeLog_("Warning: Wrong parameter type of '" + location + it.getName() + "' in '" + filename + "'. Type should be 'int'!");
          }
          break;

        case ParameterInformation::STRINGLIST:
        case ParameterInformation::INPUT_FILE_LIST:
        case ParameterInformation::OUTPUT_FILE_LIST:
          if (it->value.valueType() != DataValue::STRING_LIST)
          {
            writeLog_("Warning: Wrong parameter type of '" + location + it.getName() + "' in '" + filename + "'. Type should be 'string list'!");
          }
          break;

        case ParameterInformation::INTLIST:
          if (it->value.valueType() != DataValue::INT_LIST)
          {
            writeLog_("Warning: Wrong parameter type of '" + location + it.getName() + "' in '" + filename + "'. Type should be 'int list'!");
          }
          break;

        case ParameterInformation::DOUBLELIST:
          if (it->value.valueType() != DataValue::DOUBLE_LIST)
          {
            writeLog_("Warning: Wrong parameter type of '" + location + it.getName() + "' in '" + filename + "'. Type should be 'double list'!");
          }
          break;

        default:
          break;
        }
      }
      catch (UnregisteredParameter)
      {
        writeLog_("Warning: Unknown parameter '" + location + it.getName() + "' in '" + filename + "'!");
      }
    }
  }

  void TOPPBase::checkIfIniParametersAreApplicable_(const Param& ini_params)
  {
    Param tool_params = ini_params.copy(getIniLocation_());
    if (tool_params.empty())
    {
      // the ini file does not contain a section for our tool -> warn the user
      writeLog_(String("Warning: The provided INI file does not contain any parameters specific for this tool (expected in '") + getIniLocation_() + "'). Please check your .ini file. The default parameters for this tool will be applied.");
    }
  }

  void TOPPBase::inputFileReadable_(const String& filename, const String& param_name) const
  {
    writeDebug_("Checking input file '" + filename + "'", 2);

    // prepare error message
    String message;
    if (param_name == "")
      message = "Cannot read input file!\n";
    else
      message = "Cannot read input file given from parameter '-" + param_name + "'!\n";

    // check file
    if (!File::exists(filename))
    {
      LOG_ERROR << message;
      throw FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
    if (!File::readable(filename))
    {
      LOG_ERROR << message;
      throw FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
    if (!File::isDirectory(filename) && File::empty(filename))
    {
      LOG_ERROR << message;
      throw FileEmpty(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
  }

  void TOPPBase::outputFileWritable_(const String& filename, const String& param_name) const
  {
    writeDebug_("Checking output file '" + filename + "'", 2);

    // prepare error message
    String message;
    if (param_name == "")
      message = "Cannot write output file!\n";
    else
      message = "Cannot write output file given from parameter '-" + param_name + "'!\n";

    if (!File::writable(filename))
    {
      LOG_ERROR << message;
      throw UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
  }

  void TOPPBase::registerSubsection_(const String& name, const String& description)
  {
    subsections_[name] = description;
  }

  void TOPPBase::registerTOPPSubsection_(const String& name, const String& description)
  {
    subsections_TOPP_[name] = description;
  }

  bool TOPPBase::parseRange_(const String& text, double& low, double& high) const
  {
    bool any_set = false;
    try
    {
      String tmp = text.prefix(':');
      if (!tmp.empty())
      {
        low = tmp.toDouble();
        any_set = true;
      }

      tmp = text.suffix(':');
      if (!tmp.empty())
      {
        high = tmp.toDouble();
        any_set = true;
      }
    }
    catch (Exception::ConversionError&)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "Could not convert string '" + text +
                                       "' to a range of floating point values");
    }
    return any_set;
  }

  bool TOPPBase::parseRange_(const String& text, Int& low, Int& high) const
  {
    bool any_set = false;
    try
    {
      String tmp = text.prefix(':');
      if (!tmp.empty())
      {
        low = tmp.toInt();
        any_set = true;
      }

      tmp = text.suffix(':');
      if (!tmp.empty())
      {
        high = tmp.toInt();
        any_set = true;
      }
    }
    catch (Exception::ConversionError&)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "Could not convert string '" + text +
                                       "' to a range of integer values");
    }
    return any_set;
  }

  Param TOPPBase::getSubsectionDefaults_(const String& /*section*/) const
  {
    throw NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
  }

  Param TOPPBase::getDefaultParameters_() const
  {
    Param tmp;
    String loc = tool_name_ + ":" + String(instance_number_) + ":";
    //parameters
    for (vector<ParameterInformation>::const_iterator it = parameters_.begin(); it != parameters_.end(); ++it)
    {
      if (it->name == "ini" || it->name == "-help" || it->name == "-helphelp" || it->name == "instance" || it->name == "write_ini" || it->name == "write_ctd") // do not store those params in ini file
      {
        continue;
      }
      String name = loc + it->name;
      StringList tags;
      if (it->advanced)
        tags.push_back("advanced");
      if (it->required)
        tags.push_back("required");
      if (it->type == ParameterInformation::INPUT_FILE || it->type == ParameterInformation::INPUT_FILE_LIST)
        tags.push_back("input file");
      if (it->type == ParameterInformation::OUTPUT_FILE || it->type == ParameterInformation::OUTPUT_FILE_LIST)
        tags.push_back("output file");
      switch (it->type)
      {
      case ParameterInformation::STRING:
        tmp.setValue(name, (String)it->default_value, it->description, tags);
        if (it->valid_strings.size() != 0)
        {
          tmp.setValidStrings(name, it->valid_strings);
        }
        break;

      case ParameterInformation::INPUT_FILE:
      case ParameterInformation::OUTPUT_FILE:
        tmp.setValue(name, (String)it->default_value, it->description, tags);
        if (it->valid_strings.size() != 0)
        {
          StringList vss_tmp = it->valid_strings;
          StringList vss;
          foreach(String vs, vss_tmp)
          {
            vss.push_back("*." + vs);
          }
          tmp.setValidStrings(name, vss);
        }
        break;

      case ParameterInformation::DOUBLE:
        tmp.setValue(name, it->default_value, it->description, tags);
        if (it->min_float != -std::numeric_limits<double>::max())
        {
          tmp.setMinFloat(name, it->min_float);
        }
        if (it->max_float != std::numeric_limits<double>::max())
        {
          tmp.setMaxFloat(name, it->max_float);
        }
        break;

      case ParameterInformation::INT:
        tmp.setValue(name, (Int)it->default_value, it->description, tags);
        if (it->min_int != -std::numeric_limits<Int>::max())
        {
          tmp.setMinInt(name, it->min_int);
        }
        if (it->max_int != std::numeric_limits<Int>::max())
        {
          tmp.setMaxInt(name, it->max_int);
        }
        break;

      case ParameterInformation::FLAG:
        tmp.setValue(name, "false", it->description, tags);
        tmp.setValidStrings(name, ListUtils::create<String>("true,false"));
        break;

      case ParameterInformation::INPUT_FILE_LIST:
      case ParameterInformation::OUTPUT_FILE_LIST:
        tmp.setValue(name, it->default_value, it->description, tags);
        if (it->valid_strings.size() != 0)
        {
          StringList vss_tmp = it->valid_strings;
          StringList vss;
          foreach(String vs, vss_tmp)
          {
            vss.push_back("*." + vs);
          }
          tmp.setValidStrings(name, vss);
        }
        break;

      case ParameterInformation::STRINGLIST:
        tmp.setValue(name, it->default_value, it->description, tags);
        if (it->valid_strings.size() != 0)
        {
          tmp.setValidStrings(name, it->valid_strings);
        }
        break;

      case ParameterInformation::INTLIST:
        tmp.setValue(name, it->default_value, it->description, tags);
        if (it->min_int != -std::numeric_limits<Int>::max())
        {
          tmp.setMinInt(name, it->min_int);
        }
        if (it->max_int != std::numeric_limits<Int>::max())
        {
          tmp.setMaxInt(name, it->max_int);
        }
        break;

      case ParameterInformation::DOUBLELIST:
        tmp.setValue(name, it->default_value, it->description, tags);
        if (it->min_float != -std::numeric_limits<double>::max())
        {
          tmp.setMinFloat(name, it->min_float);
        }
        if (it->max_float != std::numeric_limits<double>::max())
        {
          tmp.setMaxFloat(name, it->max_float);
        }
        break;

      default:
        break;
      }
    }

    //subsections intrinsic to TOPP tool (i.e. a command line param with a ':')
    for (map<String, String>::const_iterator it = subsections_TOPP_.begin(); it != subsections_TOPP_.end(); ++it)
    {
      tmp.setSectionDescription(loc + it->first, it->second);
    }

    // set tool version
    tmp.setValue(tool_name_ + ":version", version_, "Version of the tool that generated this parameters file.", ListUtils::create<String>("advanced"));

    // Descriptions
    tmp.setSectionDescription(tool_name_, tool_description_);
    tmp.setSectionDescription(tool_name_ + ":" + String(instance_number_), String("Instance '") + String(instance_number_) + "' section for '" + tool_name_ + "'");

    // add type (as default type is "", but .ini file should have it)
    if (param_cmdline_.exists("type"))
      tmp.setValue(loc + "type", param_cmdline_.getValue("type"));

    // Subsections
    Param sub_sections = getSubsectionDefaults_();
    if (!sub_sections.empty())
    {
      tmp.insert(loc, sub_sections);
    }

    // 2nd stage, use TOPP tool defaults from home (if existing)
    Param tool_user_defaults(getToolUserDefaults_(tool_name_));
    tmp.update(tool_user_defaults);

    // 3rd stage, use OpenMS.ini from library to override settings
    // -> currently disabled as we cannot write back those values to the params

    return tmp;
  }

  Param TOPPBase::getSubsectionDefaults_() const
  {
    Param tmp;

    // Subsections
    for (map<String, String>::const_iterator it = subsections_.begin(); it != subsections_.end(); ++it)
    {
      Param tmp2 = getSubsectionDefaults_(it->first);
      if (!tmp2.empty())
      {
        tmp.insert(it->first + ":", tmp2);
        tmp.setSectionDescription(it->first, it->second);
      }
    }

    return tmp;
  }

  Param TOPPBase::getToolUserDefaults_(const String& tool_name) const
  {
    Param p;
    String ini_name(File::getUserDirectory() + "/" + tool_name + ".ini");
    if (File::readable(ini_name))
    {
      ParamXMLFile paramFile;
      paramFile.load(ini_name, p);
    }
    return p;
  }

  const String& TOPPBase::toolName_() const
  {
    return tool_name_;
  }

  DataProcessing TOPPBase::getProcessingInfo_(DataProcessing::ProcessingAction action) const
  {
    std::set<DataProcessing::ProcessingAction> actions;
    actions.insert(action);

    return getProcessingInfo_(actions);
  }

  DataProcessing TOPPBase::getProcessingInfo_(const std::set<DataProcessing::ProcessingAction>& actions) const
  {
    DataProcessing p;
    //actions
    p.setProcessingActions(actions);
    //software
    p.getSoftware().setName(tool_name_);

    if (test_mode_)
    {
      //version
      p.getSoftware().setVersion("version_string");

      //time
      DateTime date_time;
      date_time.set("1999-12-31 23:59:59");
      p.setCompletionTime(date_time);

      //parameters
      p.setMetaValue("parameter: mode", "test_mode");
    }
    else
    {
      //version
      p.getSoftware().setVersion(version_);
      //time
      p.setCompletionTime(DateTime::now());
      //parameters
      const Param& param = getParam_();
      for (Param::ParamIterator it = param.begin(); it != param.end(); ++it)
      {
        p.setMetaValue(String("parameter: ") + it.getName(), it->value);
      }
    }

    return p;
  }

  void TOPPBase::addDataProcessing_(ConsensusMap& map, const DataProcessing& dp) const
  {
    map.getDataProcessing().push_back(dp);

    //remove absolute map paths
    if (test_mode_)
    {
      for (Size d = 0; d < map.getFileDescriptions().size(); ++d)
      {
        map.getFileDescriptions()[d].filename = File::basename(map.getFileDescriptions()[d].filename);
      }
    }
  }
  
  void TOPPBase::addDataProcessing_(FeatureMap& map, const DataProcessing& dp) const
  {
    map.getDataProcessing().push_back(dp);
  }

  bool TOPPBase::writeCTD_()
  {
    //store ini-file content in ini_file_str
    QString out_dir_str = String(param_cmdline_.getValue("write_ctd")).toQString();
    if (out_dir_str == "")
    {
      out_dir_str = QDir::currentPath();
    }
    StringList type_list = ToolHandler::getTypes(tool_name_);
    if (type_list.size() == 0)
      type_list.push_back(""); // no type for most tools (except GenericWrapper)

    for (Size i = 0; i < type_list.size(); ++i)
    {
      QString write_ctd_file = out_dir_str + QDir::separator() + tool_name_.toQString() + type_list[i].toQString() + ".ctd";
      outputFileWritable_(write_ctd_file, "write_ctd");

      // set type on command line, so that getDefaultParameters_() does not fail (as it calls getSubSectionDefaults() of tool)
      if (type_list[i] != "")
        param_cmdline_.setValue("type", type_list[i]);
      Param default_params = getDefaultParameters_();

      // add type to ini file
      if (type_list[i] != "")
        default_params.setValue(this->ini_location_ + "type", type_list[i]);

      std::stringstream* ss = new std::stringstream();
      ParamXMLFile paramFile;
      paramFile.writeXMLToStream(ss, default_params);
      String ini_file_str(ss->str());

      //
      QString docurl = "", category = "";
      if (official_) // we can only get the docurl/category from registered/official tools
      {
        docurl = "http://ftp.mi.fu-berlin.de/OpenMS/release-documentation/html/TOPP_" + tool_name_.toQString() + ".html";
        category = ToolHandler::getCategory(tool_name_).toQString();
      }
      else if (ToolHandler::getUtilList().count(tool_name_))
      {
        docurl = "http://ftp.mi.fu-berlin.de/OpenMS/release-documentation/html/UTILS_" + tool_name_.toQString() + ".html";
        category = ToolHandler::getCategory(tool_name_).toQString();
      }

      // morph to ctd format
      QStringList lines = ini_file_str.toQString().split("\n", QString::SkipEmptyParts);
      lines.replace(0, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
      lines.insert(1, QString("<tool ctdVersion=\"1.7\" version=\"%1\" name=\"%2\" docurl=\"%3\" category=\"%4\" >").arg(version_.toQString(), tool_name_.toQString(), docurl, category));
      lines.insert(2, QString("<description><![CDATA[") + tool_description_.toQString() + "]]></description>");
      QString html_doc = tool_description_.toQString();
      lines.insert(3, QString("<manual><![CDATA[") + html_doc + "]]></manual>");
      lines.insert(4, QString("<citations>"));
      lines.insert(5, QString("  <citation doi=\"") + QString::fromStdString(cite_openms_.doi) + "\" url=\"\">");
      int l = 5;
      if (!citations_.empty())
      {
        for (Citation c : citations_) 
        {
          lines.insert(++l, QString("  <citation doi=\"") + QString::fromStdString(c.doi) + "\" url=\"\">");
        }
      }
      lines.insert(++l, QString("</citations>"));

      lines.insert(lines.size(), "</tool>");
      String ctd_str = String(lines.join("\n")) + "\n";

      //write to file
      QFile file(write_ctd_file);
      if (!file.open(QIODevice::WriteOnly))
      {
        return false;
      }
      file.write(ctd_str.c_str());
      file.close();
    }

    return true;
  }

  Param TOPPBase::parseCommandLine_(const int argc, const char** argv, const String& misc, const String& unknown)
  {
    Param cmd_params;

    // current state:
    // 'parameters_' contains all commandline params which were registered using 'registerOptionsAndFlags_()' + the common ones (-write_ini etc)
    // .. they are empty/default at this point
    // We now fetch the (so-far unknown) subsection parameters (since they can be addressed on command line as well)

    // special case of GenericWrapper: since we need the subSectionDefaults before pushing the cmd arguments in there
    //                                 but the 'type' is empty currently,
    //                                 we extract and set it beforehand
    StringList sl_args = StringList(argv, argv + argc);
    StringList::iterator it_type = std::find(sl_args.begin(), sl_args.end(), "-type");
    if (it_type != sl_args.end())
    { // found it
      ++it_type; // advance to next argument -- this should be the value of -type
      if (it_type != sl_args.end()) param_.setValue("type", *it_type);
    }

    // prepare map of parameters:
    typedef map<String, vector<ParameterInformation>::const_iterator> ParamMap;
    ParamMap param_map;
    for (vector<ParameterInformation>::const_iterator it = parameters_.begin(); it != parameters_.end(); ++it)
    {
      param_map["-" + it->name] = it;
    }

    vector<ParameterInformation> subsection_param;
    try
    {
      // the parameters from the subsections
      subsection_param = paramToParameterInformation_(getSubsectionDefaults_());
      for (vector<ParameterInformation>::const_iterator it = subsection_param.begin(); it != subsection_param.end(); ++it)
      {
        param_map["-" + it->name] = it;
      }
    }
    catch (BaseException& e)
    { // this only happens for GenericWrapper, if 'type' is not given or invalid (then we do not have subsection params) -- enough to issue a warning
      writeLog_(String("Warning: Unable to fetch subsection parameters! Addressing subsection parameters will not work for this tool (did you forget to specify '-type'?)."));
      writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ")!", 1);
    }

    // list to store "misc"/"unknown" items:
    map<String, StringList> misc_unknown;

    list<String> queue; // queue for arguments
                        // we parse the arguments in reverse order, so that we have arguments already when we encounter the option that uses them!
    for (int i = argc - 1; i > 0; --i)
    {
      String arg = argv[i];
      // options start with "-" or "--" followed by a letter:
      bool is_option = (arg.size() >= 2) && (arg[0] == '-') && (isalpha(arg[1]) || ((arg[1] == '-') && (arg.size() >= 3) &&  isalpha(arg[2])));
      if (is_option) // process content of the queue
      {
        ParamMap::iterator pos = param_map.find(arg);
        if (pos != param_map.end()) // parameter is defined
        {
          DataValue value;
          if (pos->second->type == ParameterInformation::FLAG) // flag
          {
            value = String("true");
          }
          else // option with argument(s)
          {
            switch (pos->second->type)
            {
            case ParameterInformation::STRING:
            case ParameterInformation::INPUT_FILE:
            case ParameterInformation::OUTPUT_FILE:
              if (queue.empty())
                value = String();
              else
                value = queue.front();
              break;

            case ParameterInformation::INT:
              if (!queue.empty())
                value = queue.front().toInt();
              break;

            case ParameterInformation::DOUBLE:
              if (!queue.empty())
                value = queue.front().toDouble();
              break;

            case ParameterInformation::INPUT_FILE_LIST:
            case ParameterInformation::OUTPUT_FILE_LIST:
            case ParameterInformation::STRINGLIST:
            {
              vector<String> arg_list(queue.begin(), queue.end());
              value = StringList(arg_list);
              queue.clear();
              break;
            }

            case ParameterInformation::INTLIST:
            {
              IntList arg_list;
              for (list<String>::iterator it = queue.begin(); it != queue.end(); ++it)
              {
                arg_list.push_back(it->toInt());
              }
              value = arg_list;
              queue.clear();
              break;
            }

            case ParameterInformation::DOUBLELIST:
            {
              DoubleList arg_list;
              for (list<String>::iterator it = queue.begin(); it != queue.end(); ++it)
              {
                arg_list.push_back(it->toDouble());
              }
              value = arg_list;
              queue.clear();
              break;
            }

            default:
              break;
            }
            if (!queue.empty())
              queue.pop_front(); // argument was already used
          }
          LOG_DEBUG << "Command line: setting parameter value: '" << pos->second->name << "' to '" << value << "'" << std::endl;
          cmd_params.setValue(pos->second->name, value);
        }
        else // unknown argument -> append to "unknown" list
        {
          misc_unknown[unknown].push_back(arg);
        }
        // rest of the queue is just text -> insert into "misc" list:
        StringList& misc_list = misc_unknown[misc];
        misc_list.insert(misc_list.begin(), queue.begin(), queue.end());
        queue.clear();
      }
      else // more arguments
      {
        queue.push_front(arg); // order in the queue is not reversed!
      }
    }
    // remaining items in the queue are leading text arguments:
    StringList& misc_list = misc_unknown[misc];
    misc_list.insert(misc_list.begin(), queue.begin(), queue.end());

    // store "misc"/"unknown" items, if there were any:
    for (map<String, StringList>::iterator it = misc_unknown.begin();
         it != misc_unknown.end(); ++it)
    {
      if (it->second.empty())
        continue;

      if (!cmd_params.exists(it->first))
      {
        cmd_params.setValue(it->first, it->second);
      }
      else
      {
        StringList new_value = cmd_params.getValue(it->first);
        new_value.insert(new_value.end(), it->second.begin(), it->second.end());
        cmd_params.setValue(it->first, new_value);
      }
    }

    return cmd_params;
  }

} // namespace OpenMS
