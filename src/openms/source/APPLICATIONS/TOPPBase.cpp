// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl, Johannes Junker, Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/APPLICATIONS/ConsoleUtils.h>
#include <OpenMS/APPLICATIONS/ParameterInformation.h>
#include <OpenMS/APPLICATIONS/ToolHandler.h>

#include <OpenMS/CONCEPT/Colorizer.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/CONCEPT/VersionInfo.h>

#include <OpenMS/DATASTRUCTURES/Date.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/IndentedStream.h>
#include <OpenMS/FORMAT/ParamCTDFile.h>
#include <OpenMS/FORMAT/ParamCWLFile.h>
#include <OpenMS/FORMAT/ParamJSONFile.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/SYSTEM/ExternalProcess.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/SYSTEM/SysInfo.h>
#include <OpenMS/SYSTEM/UpdateCheck.h>

#include <QtCore/QDir>

#include <iostream>

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
  const Citation TOPPBase::cite_openms
    = {"Pfeuffer, J., Bielow, C., Wein, S. et al.", "OpenMS 3 enables reproducible analysis of large-scale mass spectrometry data",
       "Nat Methods (2024)", "10.1038/s41592-024-02197-7"};


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

  String TOPPBase::getToolPrefix() const
  {
    return tool_name_ + ":" + instance_number_ + ":";
  }

  TOPPBase::TOPPBase(const String& tool_name, const String& tool_description, bool official, const std::vector<Citation>& citations, bool toolhandler_test) :
    tool_name_(tool_name),
    tool_description_(tool_description),
    instance_number_(-1),
    official_(official),
    citations_(citations),
    toolhandler_test_(toolhandler_test),
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

    // can be disabled to allow unit tests
    if (toolhandler_test_)
    {
      // check if tool is in official tools list
      if (official_ && tool_name_ != "GenericWrapper" && !ToolHandler::getTOPPToolList().count(tool_name_))
      {
        throw Exception::InvalidValue(__FILE__,
                                      __LINE__,
                                      OPENMS_PRETTY_FUNCTION,
                                      String("If '" + tool_name_ + "' is an official TOPP tool, add it to the tools list in ToolHandler. If it is not, set the 'official' flag of the TOPPBase constructor to false."),
                                      tool_name_);
      }
    }
  }

  TOPPBase::~TOPPBase()
  {
    // delete log file if empty
    const std::string& topplog = getParam_("log").toString();
    if (!topplog.empty() && File::empty(topplog))
    {
      File::remove(topplog);
    }
  }

  TOPPBase::ExitCodes TOPPBase::main(int argc, const char** argv)
  {
    //----------------------------------------------------------
    //parse command line
    //----------------------------------------------------------

    // register values from derived TOPP tool
    registerOptionsAndFlags_();
    addEmptyLine_();
    // common section for all tools
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
    registerStringOption_("write_nested_cwl", "<out_dir>", "", "Writes the Common Workflow Language file(s) (Toolname(s).cwl) to <out_dir>", false, true);
    registerStringOption_("write_cwl", "<out_dir>", "", "Writes the Common Workflow Language file(s) (Toolname(s).cwl) to <out_dir>, but enforce a flat parameter hierarchy", false, true);
    registerStringOption_("write_nested_json", "<out_dir>", "", "Writes the default configuration file", false, true);
    registerStringOption_("write_json", "<out_dir>", "", "Writes the default configuration file, but compatible to the flat hierarchy", false, true);
    registerFlag_("no_progress", "Disables progress logging to command line", true);
    registerFlag_("force", "Overrides tool-specific checks", true);
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
      writeLogError_("Invalid parameter values (" + String(e.getName()) + "): " + String(e.what()) + ". Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    // for now command line is all we have, final assembly will follow below
    param_ = param_cmdline_;

    // assign instance number
    *const_cast<int*>(&instance_number_) = getParamAsInt_("instance", 1);
    writeDebug_(String("Instance: ") + String(instance_number_), 1);

    // assign ini location
    *const_cast<String*>(&ini_location_) = this->getToolPrefix();
    writeDebug_(String("Ini_location: ") + getIniLocation_(), 1);

    // set debug level
    debug_level_ = getParamAsInt_("debug", 0);
    writeDebug_(String("Debug level: ") + String(debug_level_), 1);

    // print command line to console
    StringList args;
    for (int i = 0; i < argc; ++i)
    {
      if (String(argv[i]).has(' '))
      {
        args.push_back(String(argv[i]).quote()); // surround with quotes if argument contains a space
      }
      else
      {
        args.push_back(argv[i]);
      }
    }
    writeDebug_(String(" >> ") + ListUtils::concatenate(args, " "), 1);


    // test if no options were given
    if (argc == 1)
    {
      printUsage_();
      writeLogError_("No options given. Aborting!");
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
      writeLogError_(String("Unknown option(s) '") + getParamAsString_("unknown") + "' given. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    // test if unknown text argument were given (we do not use them)
    if (param_cmdline_.exists("misc"))
    {
      writeLogError_(String("Trailing text argument(s) '") + getParamAsString_("misc") + "' given. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    ExitCodes result;
    try
    {
      //-------------------------------------------------------------
      // store configuration or tool description files
      //-------------------------------------------------------------

      // '-write_ini' given
      if (param_cmdline_.exists("write_ini"))
      {
        String write_ini_file = param_cmdline_.getValue("write_ini").toString();
        outputFileWritable_(write_ini_file, "write_ini");
        Param default_params = getDefaultParameters_();

        // check if augmentation with -ini param is needed
        ParamValue in_ini;
        if (param_cmdline_.exists("ini"))
        {
          in_ini = param_cmdline_.getValue("ini");
          Param ini_params;
          const std::string in_ini_path = in_ini.toString();
          if (FileHandler::getTypeByFileName(in_ini_path) == FileTypes::Type::JSON)
          {
            // The JSON file doesn't carry any information about the parameter tree structure.
            // We hand an additional parameter object with the default values, so we have information
            // about the tree when parsing the JSON file.
            ini_params = getDefaultParameters_();
            if (!ParamJSONFile::load(in_ini_path, ini_params))
            {
              return ILLEGAL_PARAMETERS;
            }
          } else {
            ParamXMLFile().load(in_ini_path, ini_params);
          }

          // check if ini parameters are applicable to this tool
          checkIfIniParametersAreApplicable_(ini_params);
          // update default params with outdated params given in -ini and be verbose
          default_params.update(ini_params, false);
        }
        ParamXMLFile paramFile{};
        paramFile.store(write_ini_file, default_params);
        return EXECUTION_OK;
      }

      // '-write_ctd' given
      if (param_cmdline_.exists("write_ctd"))
      {
        ParamCTDFile paramFile{};
        writeToolDescription_(paramFile, "write_ctd", ".ctd");
        return EXECUTION_OK;
      }

      // '-write_cwl' given
      if (param_cmdline_.exists("write_nested_cwl"))
      {
        ParamCWLFile paramFile{};
        writeToolDescription_(paramFile, "write_nested_cwl", ".cwl");
        return EXECUTION_OK;
      }

      // '-write_flat_cwl' given
      if (param_cmdline_.exists("write_cwl"))
      {
        ParamCWLFile paramFile{};
        paramFile.flatHierarchy = true;
        writeToolDescription_(paramFile, "write_cwl", ".cwl");
        return EXECUTION_OK;
      }

      // '-write_json' given
      if (param_cmdline_.exists("write_nested_json"))
      {
        ParamJSONFile paramFile{};
        writeToolDescription_(paramFile, "write_nested_json", ".json");
        return EXECUTION_OK;
      }

      // '-write_flat_json' given
      if (param_cmdline_.exists("write_json"))
      {
        ParamJSONFile paramFile{};
        paramFile.flatHierarchy = true;
        writeToolDescription_(paramFile, "write_json", ".json");
        return EXECUTION_OK;
      }


      //-------------------------------------------------------------
      // load INI file
      //-------------------------------------------------------------
      {
        String value_ini;

        if (param_cmdline_.exists("ini"))
        {
          value_ini = param_cmdline_.getValue("ini").toString();
          writeDebug_("INI file: " + value_ini, 1);
          writeDebug_("INI location: " + getIniLocation_(), 1);

          if (FileHandler::getTypeByFileName(value_ini) == FileTypes::Type::JSON)
          {
            writeDebug_("Assuming INI is a cwl file", 1);
            // The JSON file doesn't carry any information about the parameter tree structure.
            // We prepopulate the param object with the default values, so we have information
            // about the tree when parsing the JSON file.
            param_inifile_ = getDefaultParameters_();
            if (!ParamJSONFile::load(value_ini, param_inifile_))
            {
              return ILLEGAL_PARAMETERS;
            }

          } else {
            ParamXMLFile().load(value_ini, param_inifile_);
          }
          checkIfIniParametersAreApplicable_(param_inifile_);

          // dissect loaded INI parameters
          param_instance_ = param_inifile_.copy(getIniLocation_(), true);
          writeDebug_("Parameters from instance section:", param_instance_, 2);
          param_common_tool_ = param_inifile_.copy("common:" + tool_name_ + ":", true);
          writeDebug_("Parameters from common section with tool name:", param_common_tool_, 2);
          param_common_ = param_inifile_.copy("common:", true);
          writeDebug_("Parameters from common section without tool name:", param_common_, 2);

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
        if (!param_.update(finalParam, false, false, true, true, OpenMS_Log_warn))
        {
          OPENMS_LOG_ERROR << "Parameters passed to '" << this->tool_name_ << "' are invalid. To prevent usage of wrong defaults, please update/fix the parameters!" << std::endl;
          return ILLEGAL_PARAMETERS;
        }

        if (finalParam.exists("type"))
        {
          param_.setValue("type", finalParam.getValue("type"));
        }

        // check if all parameters are registered and have the correct type
        checkParam_(param_instance_, value_ini, getIniLocation_());
        checkParam_(param_common_tool_, value_ini, "common:" + tool_name_ + "::");
        checkParam_(param_common_, value_ini, "common:");

        // check if the version of the parameters file matches the version of this tool
        // the parameters and values are all ok, but there might be more valid values now or new parameters which are currently not visible in the outdated INI
        String file_version = "";
        if (param_inifile_.exists(tool_name_ + ":version"))
        {
          file_version = param_inifile_.getValue(tool_name_ + ":version").toString();
          if (file_version != version_)
          {
            writeLogInfo_(String("Warning: Parameters file version (") + file_version + ") does not match the version of this tool (" + version_ + ").\n"
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
      if (!test_mode_ && (disable_usage == nullptr || strcmp(disable_usage, "OFF") == 0))
      {
        UpdateCheck::run(tool_name_, version_, debug_level_);
      }
  #endif

      //-------------------------------------------------------------
      // debug level
      //-------------------------------------------------------------
      debug_level_ = getParamAsInt_("debug", 0);
      writeDebug_(String("Debug level (after ini file): ") + String(debug_level_), 1);
      if (debug_level_ > 0) OpenMS_Log_debug.insert(cout); // allows to use OPENMS_LOG_DEBUG << "something" << std::endl;

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
      TOPPBase::setMaxNumberOfThreads(getParamAsInt_("threads", 1));

      //----------------------------------------------------------
      //main
      //----------------------------------------------------------
      StopWatch sw;
      sw.start();
      result = main_(argc, argv);
      sw.stop();
      // useful for benchmarking and for execution on clusters with schedulers
      String mem_usage;
      {
        size_t mem_virtual(0);
        SysInfo::getProcessPeakMemoryConsumption(mem_virtual);
        if (mem_virtual != 0) mem_usage = String("; Peak Memory Usage: ") + (mem_virtual / 1024) + " MB";
      }
      OPENMS_LOG_INFO << this->tool_name_ << " took " << sw.toString() << mem_usage << "." << std::endl;
    } // end try{}
    //----------------------------------------------------------
    //error handling
    //----------------------------------------------------------
    // Errors caused by the user
    catch (UnableToCreateFile& e)
    {
      writeLogError_(String("Error: Unable to write file (") + e.what() + ")");
      writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ")!", 1);
      return CANNOT_WRITE_OUTPUT_FILE;
    }
    catch (FileNotFound& e)
    {
      writeLogError_(String("Error: File not found (") + e.what() + ")");
      writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !", 1);
      return INPUT_FILE_NOT_FOUND;
    }
    catch (ExternalExecutableNotFound& e)
    {
      writeLogError_(String("Error: Executable not found (") + e.what() + ")");
      writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !", 1);
      return EXTERNAL_PROGRAM_NOTFOUND;
    }
    catch (FileNotReadable& e)
    {
      writeLogError_(String("Error: File not readable (") + e.what() + ")");
      writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !", 1);
      return INPUT_FILE_NOT_READABLE;
    }
    catch (FileEmpty& e)
    {
      writeLogError_(String("Error: File empty (") + e.what() + ")");
      writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !", 1);
      return INPUT_FILE_EMPTY;
    }
    catch (ParseError& e)
    {
      writeLogError_(String("Error: Unable to read file (") + e.what() + ")");
      writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !", 1);
      return INPUT_FILE_CORRUPT;
    }
    catch (RequiredParameterNotGiven& e)
    {
      String what = e.what();
      if (!what.hasPrefix("'"))
        what = "'" + what + "'";
      writeLogError_(String("Error: The required parameter ") + what + " was not given or is empty!");
      writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !", 1);
      return MISSING_PARAMETERS;
    }
    catch (InvalidParameter& e)
    {
      writeLogError_(String("Invalid parameter: ") + e.what());
      writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !", 1);
      return ILLEGAL_PARAMETERS;
    }
    // Internal errors because of wrong use of this class
    catch (UnregisteredParameter& e)
    {
      writeLogError_(String("Internal error: Request for unregistered parameter '") + e.what() + "'");
      writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !", 1);
      return INTERNAL_ERROR;
    }
    catch (WrongParameterType& e)
    {
      writeLogError_(String("Internal error: Request for parameter with wrong type '") + e.what() + "'");
      writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !", 1);
      return INTERNAL_ERROR;
    }
    // All other errors
    catch (BaseException& e)
    {
      writeLogError_(String("Error: Unexpected internal error (") + e.what() + ")");
      writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !", 1);
      return UNKNOWN_ERROR;
    }
    log_.close();

    return result;
  }

  void TOPPBase::printUsage_()
  {
    // show advanced options?
    bool verbose = getFlag_("-helphelp");
    String docurl = getDocumentationURL();

    IndentedStream is(cerr, 0, 10);
    // common output
    is << "\n"
       << invert(tool_name_) << " -- " << tool_description_ << "\n"
       << bright("Full documentation: ") << underline(docurl)  // the space is needed, otherwise the remaining line will be underlined on Windows..
       << "\n"
       << bright("Version: ") << verboseVersion_ << "\n"
       << bright("To cite OpenMS:\n") << " + " << is.indent(3) << cite_openms.toString() 
       << is.indent(0) << "\n";
    if (!citations_.empty())
    {
      is << bright() << "To cite " << tool_name_ << ':' << bright().undo() << is.indent(0) << "\n";
      for (const Citation& c : citations_)
        is << " + " << is.indent(3) << c.toString() << is.indent(0) << "\n";
    }
    is << is.indent(0) << "\n";
    is << invert("Usage:") << "\n" // line break needs to be separate, to avoid colored trailing whitespaces
       << "  " << bright(tool_name_) << " <options>" << "\n"
       << "\n";

    // print warning regarding not shown parameters
    if (!subsections_.empty() && !verbose)
    {
      is << "This tool has algorithm parameters that are not shown here! Please check the ini file for a detailed description or use the --helphelp option\n\n";
    }
    

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

    is << bright("Options") << " (" << green("mandatory options marked with '*'") << "):\n";

    // determine max length of parameters (including argument) for indentation
    UInt max_size = 0;
    for (const auto& par : parameters_)
    {
      if (!par.advanced || verbose)
      {
        max_size = max((UInt)max_size, (UInt)(par.name.size() + par.argument.size() + par.required));
      }
    }

    //offset of the descriptions
    UInt offset = 6 + max_size;
    //keep track of the current subsection we are in, to display the subsection help when a new section starts
    String current_TOPP_subsection("");

    // PRINT parameters && description, restrictions and default
    for (vector<ParameterInformation>::const_iterator it = parameters_.begin(); it != parameters_.end(); ++it)
    {
      if (it->advanced && !verbose)
      {
        continue;
      }

      // new subsection?
      String subsection = getSubsection_(it->name);
      if (!subsection.empty() && current_TOPP_subsection != subsection)
      {
        current_TOPP_subsection = subsection;
        map<String, String>::const_iterator subsec_it = subsections_TOPP_.find(current_TOPP_subsection);
        if (subsec_it == subsections_TOPP_.end())
        {
          throw ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "'" + current_TOPP_subsection + "' (TOPP subsection not registered)");
        }
        is << "\n"; // print newline for new subsection

        String subsection_description = subsec_it->second;
        if (subsection_description.length() == 0)
        {
          subsection_description = current_TOPP_subsection;
        }

        is << subsection_description << ":\n"; // print subsection description
      }
      else if (subsection.empty() && !current_TOPP_subsection.empty()) // subsection ended and normal parameters start again
      {
        current_TOPP_subsection = "";
        is << "\n"; // print newline to separate ending subsection
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
      String desc_tmp = String(it->description).firstToUpper();
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
        String tmp_s = ((String)it->default_value.toString()).substitute(", ", " ");
        if (!tmp_s.empty() && tmp_s != "[]")
        {
          addons.push_back(String("default: '") + tmp_s + "'");
        }
      }
      break;

      default:
        break;
      }

      //RESTRICTIONS
      StringList restrictions;
      switch (it->type)
      {
      case ParameterInformation::STRING:
      case ParameterInformation::INPUT_FILE:
      case ParameterInformation::OUTPUT_FILE:
      case ParameterInformation::OUTPUT_PREFIX:
      case ParameterInformation::OUTPUT_DIR:
      case ParameterInformation::STRINGLIST:
      case ParameterInformation::INPUT_FILE_LIST:
      case ParameterInformation::OUTPUT_FILE_LIST:
        if (!it->valid_strings.empty())
        {
          StringList copy = it->valid_strings;
          for (auto& str : copy)
          {
            str.quote('\'');
          }

          String add = "";
          if (it->type == ParameterInformation::INPUT_FILE
            || it->type == ParameterInformation::OUTPUT_FILE
            || it->type == ParameterInformation::OUTPUT_PREFIX
            || it->type == ParameterInformation::OUTPUT_DIR
            || it->type == ParameterInformation::INPUT_FILE_LIST
            || it->type == ParameterInformation::OUTPUT_FILE_LIST)
            add = " formats";

          restrictions.push_back("valid" + add + ": " + ListUtils::concatenate(copy, ", ")); // concatenate restrictions by comma
        }
        break;

      case ParameterInformation::INT:
      case ParameterInformation::INTLIST:
        if (it->min_int != -std::numeric_limits<Int>::max())
        {
          restrictions.push_back(String("min: '") + it->min_int + "'");
        }
        if (it->max_int != std::numeric_limits<Int>::max())
        {
          restrictions.push_back(String("max: '") + it->max_int + "'");
        }
        break;

      case ParameterInformation::DOUBLE:
      case ParameterInformation::DOUBLELIST:
        if (it->min_float != -std::numeric_limits<double>::max())
        {
          restrictions.push_back(String("min: '") + it->min_float + "'");
        }
        if (it->max_float != std::numeric_limits<double>::max())
        {
          restrictions.push_back(String("max: '") + it->max_float + "'");
        }
        break;

      default:
        break;
      }

      string addon_concat;
      //add DEFAULT and RESTRICTIONS
      if (!addons.empty())
      {
        addon_concat = String(" (") + ListUtils::concatenate(addons, " ") + ")";
      }
      string restrict_concat;
      // add DEFAULT and RESTRICTIONS
      if (!restrictions.empty())
      {
        restrict_concat = String(" (") + ListUtils::concatenate(restrictions, " ") + ")";
      }

      if (it->type == ParameterInformation::TEXT)
      {
        is << str_tmp << desc_tmp; // no indentation for text
      }
      else
      {
        is << is.indent(offset);
        if (it->required)
          is << green(str_tmp);
        else
          is << str_tmp;
        is << desc_tmp << cyan(addon_concat) << magenta(restrict_concat);
        is << is.indent(0);
      }
      
      is << "\n";
    }


    // SUBSECTION's at the end
    if (!subsections_.empty() && !verbose)
    {
      //determine indentation of description
      UInt indent = 0;
      for (map<String, String>::const_iterator it = subsections_.begin(); it != subsections_.end(); ++it)
      {
        indent = max((UInt)it->first.size(), indent);
      }
      indent += 6;

      //output
      is << "\n"
         << "The following configuration subsections are valid:\n";
      for (map<String, String>::const_iterator it = subsections_.begin(); it != subsections_.end(); ++it)
      {
        String tmp = String(" - ") + it->first;
        tmp.fillRight(' ', indent);
        is << ConsoleUtils::breakString(tmp + it->second, indent, 10);
        is << "\n";
      }
      is << "\n"
         << "You can write an example INI file using the '-write_ini' option.\n"
         << "Documentation of subsection parameters can be found in the doxygen documentation or the INIFileEditor.\n"
         << "For more information, please consult the online documentation for this tool:\n"
         << "  - " << underline(docurl) << "\n";
    }
    is << endl;
  }

  ParameterInformation TOPPBase::paramEntryToParameterInformation_(const Param::ParamEntry& entry, const String& argument, const String& full_name) const
  {
    String name = full_name.empty() ? entry.name : full_name;
    bool advanced = entry.tags.count("advanced");
    // special case for flags:
    if ((entry.value.valueType() == ParamValue::STRING_VALUE) &&
        /*entry.tags.count("flag") && */ // This would avoid autoconversion from true/false String Params when they default to false
        (entry.value == "false") && // This is the current default
        (entry.valid_strings.size() == 2) &&
        (entry.valid_strings[0] == "true") && (entry.valid_strings[1] == "false"))
    {
      return ParameterInformation(name, ParameterInformation::FLAG, "", "", entry.description, false, advanced);
    }

    const bool input_file = entry.tags.count(TAG_INPUT_FILE);
    const bool output_file = entry.tags.count(TAG_OUTPUT_FILE);
    const bool output_prefix = entry.tags.count(TAG_OUTPUT_PREFIX);
    const bool output_dir = entry.tags.count(TAG_OUTPUT_DIR);
    assert(input_file + output_file + output_prefix + output_dir <= 1); // at most one of these should be true (or none)
    enum ParameterInformation::ParameterTypes type = ParameterInformation::NONE;
    switch (entry.value.valueType())
    {
    case ParamValue::STRING_VALUE:
      if (input_file)
        type = ParameterInformation::INPUT_FILE;
      else if (output_file)
        type = ParameterInformation::OUTPUT_FILE;
      else if (output_prefix)
        type = ParameterInformation::OUTPUT_PREFIX;
      else if (output_dir)
        type = ParameterInformation::OUTPUT_DIR;
      else
        type = ParameterInformation::STRING;
      break;

    case ParamValue::INT_VALUE:
      type = ParameterInformation::INT;
      break;

    case ParamValue::DOUBLE_VALUE:
      type = ParameterInformation::DOUBLE;
      break;

    case ParamValue::STRING_LIST:
      if (input_file)
        type = ParameterInformation::INPUT_FILE_LIST;
      else if (output_file)
        type = ParameterInformation::OUTPUT_FILE_LIST;
      else
        type = ParameterInformation::STRINGLIST;
      break;

    case ParamValue::INT_LIST:
      type = ParameterInformation::INTLIST;
      break;

    case ParamValue::DOUBLE_LIST:
      type = ParameterInformation::DOUBLELIST;
      break;

    case ParamValue::EMPTY_VALUE:
      type = ParameterInformation::NONE;
      break;
    }
    bool required = entry.tags.count("required");
    ParameterInformation param(name, type, argument, entry.value, entry.description, required, advanced);
    param.valid_strings = ListUtils::toStringList<std::string>(entry.valid_strings);
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
    case ParamValue::STRING_VALUE:
      if (entry.valid_strings.empty())
        argument = "<text>"; // name?
      else
        argument = "<choice>";
      break;

    case ParamValue::INT_VALUE:
      argument = "<number>"; // integer?
      break;

    case ParamValue::DOUBLE_VALUE:
      argument = "<value>"; // float?
      break;

    case ParamValue::STRING_LIST:
      argument = "<list>";
      break;

    case ParamValue::INT_LIST:
      argument = "<numbers>";
      break;

    case ParamValue::DOUBLE_LIST:
      argument = "<values>";
      break;

    case ParamValue::EMPTY_VALUE:
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
    if (required && !default_value.empty())
      throw InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Registering a required StringOption param (" + name + ") with a non-empty default is forbidden!", default_value);
    parameters_.emplace_back(name, ParameterInformation::STRING, argument, default_value, description, required, advanced);
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

    const StringList& valids = strings;
    StringList defaults;

    if (p.type == ParameterInformation::STRING)
      defaults.push_back(String(p.default_value.toString()));
    else
      defaults = ListUtils::toStringList<std::string>(p.default_value);

    for (Size j = 0; j < defaults.size(); ++j) // allow the empty string even if not in restrictions
    {
      if (!defaults[j].empty() && !ListUtils::contains(valids, defaults[j]))
      {
        throw InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "TO THE DEVELOPER: The TOPP tool option '" + name + "' with default value " + std::string(p.default_value) + " does not meet restrictions!");
      }
    }
    p.valid_strings = strings;
  }

  void TOPPBase::setValidFormats_(const String& name, const std::vector<String>& formats, const bool force_OpenMS_format)
  {
    // check if formats are known
    if (force_OpenMS_format)
    {
      for (const auto& f : formats)
      {
        if (f != "fid")
        {
          auto ft = FileHandler::getTypeByFileName(String(".") + f);
          if (ft == FileTypes::UNKNOWN)
          {
            throw InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The file format '" + f + "' is invalid!");
          }
        }
      }
    }

    ParameterInformation& p = getParameterByName_(name);

    // check if the type matches
    if (p.type != ParameterInformation::INPUT_FILE
       && p.type != ParameterInformation::OUTPUT_FILE
       && p.type != ParameterInformation::INPUT_FILE_LIST
       && p.type != ParameterInformation::OUTPUT_FILE_LIST
       && p.type != ParameterInformation::OUTPUT_PREFIX)
      // && p.type != ParameterInformation::OUTPUT_DIR ) // output dir is not a file format, hence does not support restricting the format
    {
      throw Exception::WrongParameterType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }

    if (!p.valid_strings.empty())
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Internal error: Valid formats are already set for '" + name + "'. Please check for typos!");
    }
    p.valid_strings = formats;
  }

  void TOPPBase::setMinInt_(const String& name, Int min)
  {
    ParameterInformation& p = getParameterByName_(name);

    // check if the type matches
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
        throw InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "TO THE DEVELOPER: The TOPPS tool option '" + name + "' with default value " + std::string(p.default_value) + " does not meet restrictions!");
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
        throw InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "TO THE DEVELOPER: The TOPPS tool option '" + name + "' with default value " + std::string(p.default_value) + " does not meet restrictions!");
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
        throw InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "TO THE DEVELOPER: The TOPPS tool option '" + name + "' with default value " + std::string(p.default_value) + " does not meet restrictions!");
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
        throw InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "TO THE DEVELOPER: The TOPPS tool option '" + name + "' with default value " + std::string(p.default_value) + " does not meet restrictions!");
      }
    }
    p.max_float = max;
  }

  void TOPPBase::registerInputFile_(const String& name, const String& argument, const String& default_value, const String& description, bool required, bool advanced, const StringList& tags)
  {
    int count_conflicting_tags = (ListUtils::contains(tags, "skipexists") + ListUtils::contains(tags, "is_executable"));
    if (count_conflicting_tags >= 2)
    {
      throw Exception::WrongParameterType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "'skipexists' and 'is_executable' cannot be combined");
    }
    if (required && !default_value.empty() && count_conflicting_tags == 0)
      throw InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Registering a required InputFile param (" + name + ") with a non-empty default is forbidden!", default_value);
    parameters_.emplace_back(name, ParameterInformation::INPUT_FILE, argument, default_value, description, required, advanced, tags);
  }

  void TOPPBase::registerOutputFile_(const String& name, const String& argument, const String& default_value, const String& description, bool required, bool advanced)
  {
    if (required && !default_value.empty())
      throw InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Registering a required OutputFile param (" + name + ") with a non-empty default is forbidden!", default_value);
    parameters_.emplace_back(name, ParameterInformation::OUTPUT_FILE, argument, default_value, description, required, advanced);
  }

  void TOPPBase::registerOutputPrefix_(const String& name, const String& argument, const String& default_value, const String& description, bool required, bool advanced)
  {
    if (required && !default_value.empty())
      throw InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Registering a required OutputPrefix param (" + name + ") with a non-empty default is forbidden!", default_value);
    parameters_.emplace_back(name, ParameterInformation::OUTPUT_PREFIX, argument, default_value, description, required, advanced);
  }
  
  void TOPPBase::registerOutputDir_(const String& name, const String& argument, const String& default_value, const String& description, bool required, bool advanced)
  {
    if (required && !default_value.empty())
      throw InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Registering a required OutputDir param (" + name + ") with a non-empty default is forbidden!", default_value);
    parameters_.emplace_back(name, ParameterInformation::OUTPUT_DIR, argument, default_value, description, required, advanced);
  }

  void TOPPBase::registerDoubleOption_(const String& name, const String& argument, double default_value, const String& description, bool required, bool advanced)
  {
    if (required)
    {
      throw InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Registering a double param (" + name + ") as 'required' is forbidden (there is no value to indicate it is missing)!", String(default_value));
    }
    parameters_.emplace_back(name, ParameterInformation::DOUBLE, argument, default_value, description, required, advanced);
  }

  void TOPPBase::registerIntOption_(const String& name, const String& argument, Int default_value, const String& description, bool required, bool advanced)
  {
    if (required)
    {
      throw InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Registering an Int param (" + name + ") as 'required' is forbidden (there is no value to indicate it is missing)!", String(default_value));
    }
    parameters_.emplace_back(name, ParameterInformation::INT, argument, default_value, description, required, advanced);
  }

  void TOPPBase::registerOutputFileList_(const String& name, const String& argument, const StringList& default_value, const String& description, bool required, bool advanced)
  {
    if (required && !default_value.empty())
      throw InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Registering a required OutputFileList param (" + name + ") with a non-empty default is forbidden!", ListUtils::concatenate(default_value, ","));
    parameters_.emplace_back(name, ParameterInformation::OUTPUT_FILE_LIST, argument, ListUtils::create<std::string>(default_value), description, required, advanced);
  }

  void TOPPBase::registerInputFileList_(const String& name, const String& argument, const StringList& default_value, const String& description, bool required, bool advanced, const StringList& tags)
  {
    int count_conflicting_tags = (ListUtils::contains(tags, "skipexists") + ListUtils::contains(tags, "is_executable"));
    if (count_conflicting_tags >= 2)
    {
      throw Exception::WrongParameterType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "'skipexists' and 'is_executable' cannot be combined");
    }
    if (required && !default_value.empty() && count_conflicting_tags == 0)
      throw InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Registering a required InputFileList param (" + name + ") with a non-empty default is forbidden!", ListUtils::concatenate(default_value, ","));
    parameters_.emplace_back(name, ParameterInformation::INPUT_FILE_LIST, argument, ListUtils::create<std::string>(default_value), description, required, advanced, tags);
  }

  void TOPPBase::registerStringList_(const String& name, const String& argument, const StringList& default_value, const String& description, bool required, bool advanced)
  {
    if (required && !default_value.empty())
      throw InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Registering a required StringList param (" + name + ") with a non-empty default is forbidden!", ListUtils::concatenate(default_value, ","));
    parameters_.emplace_back(name, ParameterInformation::STRINGLIST, argument, ListUtils::create<std::string>(default_value), description, required, advanced);
  }

  void TOPPBase::registerIntList_(const String& name, const String& argument, const IntList& default_value, const String& description, bool required, bool advanced)
  {
    stringstream ss;
    ss << default_value;
    if (required && !default_value.empty())
      throw InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Registering a required IntList param (" + name + ") with a non-empty default is forbidden!", String(ss.str()));
    parameters_.emplace_back(name, ParameterInformation::INTLIST, argument, default_value, description, required, advanced);
  }

  void TOPPBase::registerDoubleList_(const String& name, const String& argument, const DoubleList& default_value, const String& description, bool required, bool advanced)
  {
    stringstream ss;
    ss << default_value;
    if (required && !default_value.empty())
      throw InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Registering a required DoubleList param (" + name + ") with a non-empty default is forbidden!", String(ss.str()));
    parameters_.emplace_back(name, ParameterInformation::DOUBLELIST, argument, default_value, description, required, advanced);
  }

  void TOPPBase::registerFlag_(const String& name, const String& description, bool advanced)
  {
    parameters_.emplace_back(name, ParameterInformation::FLAG, "", "", description, false, advanced);
  }

  void TOPPBase::addEmptyLine_()
  {
    parameters_.emplace_back("", ParameterInformation::NEWLINE, "", "", "", false, false);
  }

  void TOPPBase::addText_(const String& text)
  {
    parameters_.emplace_back("", ParameterInformation::TEXT, "", "", text, false, false);
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
    if (p.type != ParameterInformation::STRING
      && p.type != ParameterInformation::INPUT_FILE
      && p.type != ParameterInformation::OUTPUT_FILE
      && p.type != ParameterInformation::OUTPUT_PREFIX)
    {
      throw WrongParameterType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }
    if (p.required && (getParam_(name).isEmpty() || getParam_(name) == ""))
    {
      String message = "'" + name + "'";
      if (!p.valid_strings.empty())
      {
        message += " [valid: " + ListUtils::concatenate(p.valid_strings, ", ") + "]";
      }
      throw RequiredParameterNotGiven(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message);
    }
    String tmp = getParamAsString_(name, p.default_value.toString());
    writeDebug_(String("Value of string option '") + name + "': " + tmp, 1);

    // if required or set by user, do some validity checks
    if (p.required || !tmp.empty())
    {
      fileParamValidityCheck_(tmp, name, p);
    }

    return tmp;
  }

  String TOPPBase::getOutputDirOption(const String& name) const
  {
    const ParameterInformation& p = findEntry_(name);
    if (p.type != ParameterInformation::OUTPUT_DIR)
    {
      throw WrongParameterType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }
    if (p.required && (getParam_(name).isEmpty() || getParam_(name) == ""))
    {
      String message = "'" + name + "'";
      if (! p.valid_strings.empty()) { message += " [valid: " + ListUtils::concatenate(p.valid_strings, ", ") + "]"; }
      throw RequiredParameterNotGiven(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message);
    }
    String tmp = getParamAsString_(name, p.default_value.toString());
    writeDebug_(String("Value of string(outdir) option '") + name + "': " + tmp, 1);

    // create directory if it does not exist
    File::makeDir(tmp);

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
    if (p.required && std::isnan(tmp))
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

  void TOPPBase::fileParamValidityCheck_(const StringList& param_value, const String& param_name, const ParameterInformation& p) const
  {
    // check if all input files are readable
    if (p.type == ParameterInformation::INPUT_FILE_LIST)
    {
      for (String t : param_value)
      {
        if (!ListUtils::contains(p.tags, "skipexists")) inputFileReadable_(t, param_name);

        // check restrictions
        if (p.valid_strings.empty()) continue;

        // determine file type as string
        FileTypes::Type f_type = FileHandler::getType(t);
        // unknown ending is 'ok'
        if (f_type == FileTypes::UNKNOWN)
        {
          writeLogWarn_("Warning: Could not determine format of input file '" + t + "'!");
        }
        else if (!ListUtils::contains(p.valid_strings, FileTypes::typeToName(f_type).toUpper(), ListUtils::CASE::INSENSITIVE))
        {
            throw InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                   String("Input file '" + t + "' has invalid format '") +
                                   FileTypes::typeToName(f_type) +
                                   "'. Valid formats are: '" + ListUtils::concatenate(p.valid_strings, "','") +
                                   "'.");
        }
      }
    }
  }

  void TOPPBase::fileParamValidityCheck_(String& param_value, const String& param_name, const ParameterInformation& p) const
  {
    // check if files are readable/writable
    if (p.type == ParameterInformation::INPUT_FILE)
    {
      if (ListUtils::contains(p.tags, "is_executable"))
      { // will update to absolute path
        if (File::findExecutable(param_value))
        {
          writeDebug_("Input file resolved to '" + param_value + "'", 2);
        }
        else
        {
          writeLogWarn_("Input file '" + param_value + "' could not be found (by searching on PATH). "
                        "Either provide a full filepath via the '-" +
                          param_name + "' option or fix your PATH environment !" +
                    (p.required ? "" : " Since this file is not strictly required, you might also pass the empty string \"\" as "
                    "argument to prevent its usage (this might limit the usability of the tool)."));
          throw ExternalExecutableNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, param_value);
        }
      }
      if (!ListUtils::contains(p.tags, "skipexists")) inputFileReadable_(param_value, param_name);
    }
    else if (p.type == ParameterInformation::OUTPUT_FILE)
    {
      outputFileWritable_(param_value, param_name);
    }
    else if (p.type == ParameterInformation::OUTPUT_PREFIX)
    {
      outputFileWritable_(param_value + "_0", param_name); // only test one file
    }

    // check restrictions
    if (p.valid_strings.empty()) return;

    switch (p.type)
    {
      case ParameterInformation::STRING:
        if (!ListUtils::contains(p.valid_strings, param_value))
        {
          throw InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            String("Invalid value '") + param_value + "' for string parameter '" + param_name + "' given. Valid strings are: '" +
            ListUtils::concatenate(p.valid_strings, "', '") + "'.");
        }
        break;

      case ParameterInformation::INPUT_FILE:
      {
        // determine file type as string
        FileTypes::Type f_type = FileHandler::getType(param_value);
        // unknown ending is 'ok'
        if (f_type == FileTypes::UNKNOWN)
        {
          writeLogWarn_("Warning: Could not determine format of input file '" + param_value + "'!");
        }
        else if (!ListUtils::contains(p.valid_strings, FileTypes::typeToName(f_type).toUpper(), ListUtils::CASE::INSENSITIVE))
        {
            throw InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                   String("Input file '" + param_value + "' has invalid format '") +
                                   FileTypes::typeToName(f_type) +
                                   "'. Valid formats are: '" + ListUtils::concatenate(p.valid_strings, "','") +
                                   "'.");
        }
        break;
      }

      case ParameterInformation::OUTPUT_FILE:
      {
        // determine file type as string
        FileTypes::Type f_type = FileHandler::getTypeByFileName(param_value);
        // Wrong ending, unknown is is ok.
        if (f_type != FileTypes::UNKNOWN
          && !ListUtils::contains(p.valid_strings, FileTypes::typeToName(f_type).toUpper(), ListUtils::CASE::INSENSITIVE))
        {
          throw InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            String("Invalid output file extension for file '") + param_value + "'. Valid file extensions are: '" +
            ListUtils::concatenate(p.valid_strings, "','") + "'.");
        }
        break;
      }
      case ParameterInformation::OUTPUT_PREFIX: /* no file extension check for out prefixes */
        break;
      default: /*nothing */
        break;
    }
  }

  StringList TOPPBase::getStringList_(const String& name) const
  {
    const ParameterInformation& p = findEntry_(name);
    if (p.type != ParameterInformation::STRINGLIST
      && p.type != ParameterInformation::INPUT_FILE_LIST
      && p.type != ParameterInformation::OUTPUT_FILE_LIST)
    {
      throw WrongParameterType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }
    if (p.required && getParam_(name).isEmpty())
    {
      throw RequiredParameterNotGiven(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }
    StringList tmp_list = getParamAsStringList_(name, ListUtils::toStringList<std::string>(p.default_value));
    if (p.required && tmp_list.empty())
    {
      throw RequiredParameterNotGiven(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }

    for (String& tmp : tmp_list)
    {
      writeDebug_(String("Value of string option '") + name + "': " + tmp, 1);
    }

    // if required or set by user, do some validity checks
    if (p.required || (!getParam_(name).isEmpty() && tmp_list != ListUtils::toStringList<std::string>(p.default_value)))
    {
      fileParamValidityCheck_(tmp_list, name, p);
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
    if (p.required && tmp_list.empty())
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
    if (p.required && tmp_list.empty())
    {
      throw RequiredParameterNotGiven(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }

    for (const Int tmp : tmp_list)
    {
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

  void TOPPBase::writeLogInfo_(const String& text) const
  {
    OPENMS_LOG_INFO << text << endl;
    enableLogging_();
    log_ << QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss").toStdString() << ' ' << getIniLocation_() << ": " << text << endl;
  }

  void TOPPBase::writeLogWarn_(const String& text) const
  {
    OPENMS_LOG_WARN << text << endl;
    enableLogging_();
    log_ << QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss").toStdString() << ' ' << getIniLocation_() << ": " << text << endl;
  }

  void TOPPBase::writeLogError_(const String& text) const
  {
    OPENMS_LOG_ERROR << text << endl;
    enableLogging_();
    log_ << QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss").toStdString() << ' ' << getIniLocation_() << ": " << text << endl;
  }

  void TOPPBase::writeDebug_(const String& text, UInt min_level) const
  {
    if (debug_level_ >= (Int)min_level)
    {
      OPENMS_LOG_DEBUG << text << endl;
      enableLogging_();
      log_ << QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss").toStdString() << ' ' << getIniLocation_() << ": " << text << endl;
    }
  }

  void TOPPBase::writeDebug_(const String& text, const Param& param, UInt min_level) const
  {
    if (debug_level_ >= (Int)min_level)
    {
      OPENMS_LOG_DEBUG << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl
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

  TOPPBase::ExitCodes TOPPBase::runExternalProcess_(const QString& executable, const QStringList& arguments, const QString& workdir) const
  {
    String proc_stdout, proc_stderr; // collect all output (might be useful if program crashes, see below)
    return runExternalProcess_(executable, arguments, proc_stdout, proc_stderr, workdir);
  }

  TOPPBase::ExitCodes TOPPBase::runExternalProcess_(const QString& executable, const QStringList& arguments, String& proc_stdout, String& proc_stderr, const QString& workdir) const
  {
    proc_stdout.clear();
    proc_stderr.clear();

    // callbacks: invoked whenever output is available.
    auto lam_out = [&](const String& out) { proc_stdout += out; if (debug_level_ >= 4) OPENMS_LOG_INFO << out; };
    auto lam_err = [&](const String& out) { proc_stderr += out; if (debug_level_ >= 4) OPENMS_LOG_INFO << out; };
    ExternalProcess ep(lam_out, lam_err);

    const auto& rt = ep.run(executable, arguments, workdir, true); // does automatic escaping etc... start
    if (debug_level_ < 4 && rt != ExternalProcess::RETURNSTATE::SUCCESS)
    { // error occurred: if not written already in callback, do it now
      writeLogError_("Standard output: " + proc_stdout);
      writeLogError_("Standard error: " + proc_stderr);
    }
    switch (rt)
    {
      case ExternalProcess::RETURNSTATE::SUCCESS:
        return EXECUTION_OK;
      case ExternalProcess::RETURNSTATE::NONZERO_EXIT:
      case ExternalProcess::RETURNSTATE::CRASH:
        return EXTERNAL_PROGRAM_ERROR;
      case ExternalProcess::RETURNSTATE::FAILED_TO_START:
        return EXTERNAL_PROGRAM_NOTFOUND;
      default:
        throw Exception::InternalToolError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown return state of external process.");
    }
  }

  String TOPPBase::getParamAsString_(const String& key, const String& default_value) const
  {
    const ParamValue& tmp = getParam_(key);
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
    const ParamValue& tmp = getParam_(key);
    if (!tmp.isEmpty())
    {
      if (tmp.valueType() == ParamValue::INT_VALUE)
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
    const ParamValue& tmp = getParam_(key);
    if (!tmp.isEmpty())
    {
      if (tmp.valueType() == ParamValue::DOUBLE_VALUE)
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
    const ParamValue& tmp = getParam_(key);
    if (!tmp.isEmpty())
    {
      return ListUtils::toStringList<std::string>(tmp);
    }
    else
    {
      return default_value;
    }
  }

  IntList TOPPBase::getParamAsIntList_(const String& key, const IntList& default_value) const
  {
    const ParamValue& tmp = getParam_(key);
    if (!tmp.isEmpty())
    {
      if (tmp.valueType() == ParamValue::INT_LIST)
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
    const ParamValue& tmp = getParam_(key);
    if (!tmp.isEmpty())
    {
      if (tmp.valueType() == ParamValue::DOUBLE_LIST)
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
    ParamValue tmp = getParam_(key);
    if (tmp.valueType() == ParamValue::EMPTY_VALUE)
    {
      return false;
    }
    else if (tmp.valueType() == ParamValue::STRING_VALUE)
    {
      if ((std::string)tmp == "false")
      {
        return false;
      }
      else if ((std::string)tmp == "true")
      {
        return true;
      }
    }
    throw InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, std::string("Invalid value '") + (std::string)tmp + "' for flag parameter '" + key + "'. Valid values are 'true' and 'false' only.");
  }

  ParamValue const& TOPPBase::getParam_(const String& key) const
  {
    if (param_.exists(key))
    {
      return param_.getValue(key);
    }
    else
    {
      // if look up fails everywhere, return EMPTY
      writeDebug_(String("Parameter '") + key + String("' not found."), 1);
      return ParamValue::EMPTY;
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
    if (log_.is_open() || !param_.exists("log")) return;

    std::string log_destination = param_.getValue("log");
    if (log_destination.empty()) return;
    log_.open(log_destination, ofstream::out | ofstream::app);
    if (debug_level_ >= 1)
    {
      cout << "Writing to '" << log_destination << '\'' << "\n";
      log_ << QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss").toStdString() << ' ' << getIniLocation_() << ": " << "Writing to '" << log_destination << '\'' <<  "\n";
    }
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
            writeLogWarn_("Warning: Unknown subsection '" + subsection + "' in '" + filename + "' (location '" + location + "')!");
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
        case ParameterInformation::OUTPUT_PREFIX:
        case ParameterInformation::FLAG:
          if (it->value.valueType() != ParamValue::STRING_VALUE)
          {
            writeLogWarn_("Warning: Wrong parameter type of '" + location + it.getName() + "' in '" + filename + "'. Type should be 'string'!");
          }
          break;

        case ParameterInformation::DOUBLE:
          if (it->value.valueType() != ParamValue::DOUBLE_VALUE)
          {
            writeLogWarn_("Warning: Wrong  parameter type of '" + location + it.getName() + "' in '" + filename + "'. Type should be 'double'!");
          }
          break;

        case ParameterInformation::INT:
          if (it->value.valueType() != ParamValue::INT_VALUE)
          {
            writeLogWarn_("Warning: Wrong parameter type of '" + location + it.getName() + "' in '" + filename + "'. Type should be 'int'!");
          }
          break;

        case ParameterInformation::STRINGLIST:
        case ParameterInformation::INPUT_FILE_LIST:
        case ParameterInformation::OUTPUT_FILE_LIST:
          if (it->value.valueType() != ParamValue::STRING_LIST)
          {
            writeLogWarn_("Warning: Wrong parameter type of '" + location + it.getName() + "' in '" + filename + "'. Type should be 'string list'!");
          }
          break;

        case ParameterInformation::INTLIST:
          if (it->value.valueType() != ParamValue::INT_LIST)
          {
            writeLogWarn_("Warning: Wrong parameter type of '" + location + it.getName() + "' in '" + filename + "'. Type should be 'int list'!");
          }
          break;

        case ParameterInformation::DOUBLELIST:
          if (it->value.valueType() != ParamValue::DOUBLE_LIST)
          {
            writeLogWarn_("Warning: Wrong parameter type of '" + location + it.getName() + "' in '" + filename + "'. Type should be 'double list'!");
          }
          break;

        default:
          break;
        }
      }
      catch (UnregisteredParameter&)
      {
        writeLogWarn_("Warning: Unknown parameter '" + location + it.getName() + "' in '" + filename + "'!");
      }
    }
  }

  void TOPPBase::checkIfIniParametersAreApplicable_(const Param& ini_params)
  {
    Param tool_params = ini_params.copy(getIniLocation_());
    if (tool_params.empty())
    {
      // the ini file does not contain a section for our tool -> warn the user
      writeLogWarn_(String("Warning: The provided INI file does not contain any parameters specific for this tool (expected in '") + getIniLocation_() +
                             "'). Please check your .ini file. The default parameters for this tool will be applied.");
    }
  }

  void TOPPBase::inputFileReadable_(const String& filename, const String& param_name) const
  {
    writeDebug_("Checking input file '" + filename + "'", 2);

    // prepare error message
    String message;
    if (param_name.empty())
      message = "Cannot read input file!\n";
    else
      message = "Cannot read input file given from parameter '-" + param_name + "'!\n";

    // check file existence
    if (!File::exists(filename))
    {
      OPENMS_LOG_ERROR << message;
      throw FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
    if (!File::readable(filename))
    {
      OPENMS_LOG_ERROR << message;
      throw FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
    if (!File::isDirectory(filename) && File::empty(filename))
    {
      OPENMS_LOG_ERROR << message;
      throw FileEmpty(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
  }

  void TOPPBase::outputFileWritable_(const String& filename, const String& param_name) const
  {
    writeDebug_("Checking output file '" + filename + "'", 2);

    // prepare error message
    String message;
    if (param_name.empty())
      message = "Cannot write output file!\n";
    else
      message = "Cannot write output file given from parameter '-" + param_name + "'!\n";

    if (!File::writable(filename))
    {
      OPENMS_LOG_ERROR << message;
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
    String loc = this->getToolPrefix();
    //parameters
    for (vector<ParameterInformation>::const_iterator it = parameters_.begin(); it != parameters_.end(); ++it)
    {
      if (std::unordered_set<std::string>{"ini", "-help", "-helphelp", "instance", "write_ini", "write_ctd", "write_cwl", "write_nested_cwl", "write_json", "write_nested_json"}.count(it->name) > 0) // do not store these params in ini file
      {
        continue;
      }
      String name = loc + it->name;
      std::vector<std::string> tags;
      if (it->advanced)
      {
        tags.emplace_back("advanced");
      }
      if (it->required)
      {
        tags.emplace_back("required");
      }

      if (it->type == ParameterInformation::INPUT_FILE || it->type == ParameterInformation::INPUT_FILE_LIST)
      {
        tags.emplace_back(TAG_INPUT_FILE);
      }

      if (it->type == ParameterInformation::INPUT_FILE && std::find(it->tags.begin(), it->tags.end(), "is_executable") != it->tags.end())
      {
        tags.emplace_back("is_executable");
      }

      if (it->type == ParameterInformation::OUTPUT_FILE || it->type == ParameterInformation::OUTPUT_FILE_LIST) { tags.emplace_back(TAG_OUTPUT_FILE); }
      if (it->type == ParameterInformation::OUTPUT_PREFIX) { tags.emplace_back(TAG_OUTPUT_PREFIX); }
      if (it->type == ParameterInformation::OUTPUT_DIR) { tags.emplace_back(TAG_OUTPUT_DIR); }

      switch (it->type)
      {
      case ParameterInformation::STRING:
        tmp.setValue(name, (String)it->default_value.toString(), it->description, tags);
        if (!it->valid_strings.empty())
        {
          tmp.setValidStrings(name, ListUtils::create<std::string>(it->valid_strings));
        }
        break;

      case ParameterInformation::INPUT_FILE:
      case ParameterInformation::OUTPUT_FILE:
      case ParameterInformation::OUTPUT_PREFIX:
      case ParameterInformation::OUTPUT_DIR:
        tmp.setValue(name, (String)it->default_value.toString(), it->description, tags);
        if (!it->valid_strings.empty())
        {
          StringList vss_tmp = it->valid_strings;
          for (auto& vs : vss_tmp)
          {
            vs = "*." + vs;
          }
          tmp.setValidStrings(name, ListUtils::create<std::string>(vss_tmp));
        }
        break;

      case ParameterInformation::DOUBLE:
        tmp.setValue(name, it->default_value, it->description, tags);
        tmp.setMinFloat(name, it->min_float);
        tmp.setMaxFloat(name, it->max_float);
        break;

      case ParameterInformation::INT:
        tmp.setValue(name, (Int)it->default_value, it->description, tags);
        tmp.setMinInt(name, it->min_int);
        tmp.setMaxInt(name, it->max_int);
        break;

      case ParameterInformation::FLAG:
        tmp.setValue(name, "false", it->description, tags);
        tmp.setValidStrings(name, {"true","false"});
        break;

      case ParameterInformation::INPUT_FILE_LIST:
      case ParameterInformation::OUTPUT_FILE_LIST:
        tmp.setValue(name, it->default_value, it->description, tags);
        if (!it->valid_strings.empty())
        {
          std::vector<std::string> vss = ListUtils::create<std::string>(it->valid_strings);
          std::transform(vss.begin(), vss.end(), vss.begin(), [](const std::string& s) {return "*." + s;});
          tmp.setValidStrings(name, vss);
        }
        break;

      case ParameterInformation::STRINGLIST:
        tmp.setValue(name, it->default_value, it->description, tags);
        if (!it->valid_strings.empty())
        {
          tmp.setValidStrings(name, ListUtils::create<std::string>(it->valid_strings));
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
    tmp.setValue(tool_name_ + ":version", version_, "Version of the tool that generated this parameters file.", {"advanced"});

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
        p.setMetaValue(String("parameter: " + it.getName()), it->value);
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
      for (Size d = 0; d < map.getColumnHeaders().size(); ++d)
      {
        map.getColumnHeaders()[d].filename = File::basename(map.getColumnHeaders()[d].filename);
      }
    }
  }

  void TOPPBase::addDataProcessing_(FeatureMap& map, const DataProcessing& dp) const
  {
    map.getDataProcessing().push_back(dp);
  }

  ///Data processing setter for peak maps

  void TOPPBase::addDataProcessing_(PeakMap& map, const DataProcessing& dp) const
  {
    boost::shared_ptr< DataProcessing > dp_(new DataProcessing(dp));
    for (Size i = 0; i < map.size(); ++i)
    {
      map[i].getDataProcessing().push_back(dp_);
    }
    for (Size i = 0; i < map.getNrChromatograms(); ++i)
    {
      map.getChromatogram(i).getDataProcessing().push_back(dp_);
    }
  }

  String TOPPBase::getDocumentationURL() const
  {
    VersionInfo::VersionDetails ver = VersionInfo::getVersionStruct();
    String tool_prefix = "TOPP_";
    // it is only empty if the GIT_BRANCH inferred or set during CMake config was release/* or master
    // see https://github.com/OpenMS/OpenMS/blob/develop/CMakeLists.txt#L122
    if (ver.pre_release_identifier.empty())
    {
      String release_version = String(ver.version_major) + "." + String(ver.version_minor) + "." + String(ver.version_patch);
      return String("http://www.openms.de/doxygen/release/") + release_version + "/html/"+ tool_prefix + tool_name_ + ".html";
    }
    else
    {
      return String("http://www.openms.de/doxygen/nightly/html/") + tool_prefix + tool_name_ + ".html";
    }
  }

  template <typename Writer>
  void TOPPBase::writeToolDescription_(Writer& writer, std::string write_type, std::string fileExtension)
  {
    //store ini-file content in ini_file_str
    QString out_dir_str = String(param_cmdline_.getValue(write_type).toString()).toQString();
    if (out_dir_str == "")
    {
      out_dir_str = QDir::currentPath();
    }
    StringList type_list = ToolHandler::getTypes(tool_name_);
    if (type_list.empty())
      type_list.push_back(""); // no type for most tools (except GenericWrapper)

    for (Size i = 0; i < type_list.size(); ++i)
    {
      // check file is writable
      QString write_file = out_dir_str + QDir::separator() + tool_name_.toQString() + type_list[i].toQString() + fileExtension.c_str();
      outputFileWritable_(write_file, write_type);

      // set type on command line, so that getDefaultParameters_() does not fail (as it calls getSubSectionDefaults() of tool)
      if (!type_list[i].empty())
        param_cmdline_.setValue("type", type_list[i]);
      Param default_params = getDefaultParameters_();

      // add type to ini file
      if (!type_list[i].empty())
        default_params.setValue(this->ini_location_ + "type", type_list[i]);

      std::stringstream ss;

      // fill program category and docurl
      std::string docurl = getDocumentationURL();
      std::string category;
      if (official_)
      { // we can only get the docurl/category from registered/official tools
        category = ToolHandler::getCategory(tool_name_);
      }

      // collect citation information
      std::vector<std::string> citation_dois;
      citation_dois.reserve(citations_.size() + 1);
      citation_dois.push_back(cite_openms.doi);
      for (auto& citation : citations_)
      {
        citation_dois.push_back(citation.doi);
      }

      // fill tool information
      ToolInfo toolInfo{};
      toolInfo.version_     = version_;
      toolInfo.name_        = tool_name_;
      toolInfo.docurl_      = docurl;
      toolInfo.category_    = category;
      toolInfo.description_ = tool_description_;
      toolInfo.citations_   = citation_dois;

      // this will write the actual data to disk
      writer.store(write_file.toStdString(), default_params, toolInfo);
    }
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
      writeLogWarn_(String("Warning: Unable to fetch subsection parameters! Addressing subsection parameters will not work for this tool (did you forget to specify '-type'?)."));
      writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ")!", 1);
    }

    // list to store "misc"/"unknown" items:
    map<std::string, std::vector<std::string> > misc_unknown;

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
          ParamValue value;
          if (pos->second->type == ParameterInformation::FLAG) // flag
          {
            value = "true";
          }
          else // option with argument(s)
          {
            switch (pos->second->type)
            {
            case ParameterInformation::STRING:
            case ParameterInformation::INPUT_FILE:
            case ParameterInformation::OUTPUT_FILE:
            case ParameterInformation::OUTPUT_PREFIX:
            case ParameterInformation::OUTPUT_DIR:
              if (queue.empty())
                value = std::string();
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
              vector<std::string> arg_list(queue.begin(), queue.end());
              value = arg_list;
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
          OPENMS_LOG_DEBUG << "Command line: setting parameter value: '" << pos->second->name << "' to '" << value << "'" << std::endl;
          cmd_params.setValue(pos->second->name, value);
        }
        else // unknown argument -> append to "unknown" list
        {
          misc_unknown[unknown].push_back(arg);
        }
        // rest of the queue is just text -> insert into "misc" list:
        std::vector<std::string>& misc_list = misc_unknown[misc];
        misc_list.insert(misc_list.begin(), queue.begin(), queue.end());
        queue.clear();
      }
      else // more arguments
      {
        queue.push_front(arg); // order in the queue is not reversed!
      }
    }
    // remaining items in the queue are leading text arguments:
    std::vector<std::string>& misc_list = misc_unknown[misc];
    misc_list.insert(misc_list.begin(), queue.begin(), queue.end());

    // store "misc"/"unknown" items, if there were any:
    for (map<std::string, std::vector<std::string> >::iterator it = misc_unknown.begin();
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
        std::vector<std::string> new_value = cmd_params.getValue(it->first);
        new_value.insert(new_value.end(), it->second.begin(), it->second.end());
        cmd_params.setValue(it->first, new_value);
      }
    }

    return cmd_params;
  }

} // namespace OpenMS
