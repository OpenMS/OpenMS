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
// $Maintainer: Clemens Groepl $
// $Authors: Marc Sturm, Clemens Groepl, Johannes Junker $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/Date.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/VALIDATORS/XMLValidator.h>

#include <iostream>

#include <QDir>
#include <QFile>
#include <QTextStream>

#include <boost/math/special_functions/fpclassify.hpp>

// OpenMP support
#ifdef _OPENMP
	#include <omp.h>
#endif

#include <QtCore/QProcess>

#include <cmath>

#ifdef OPENMS_WINDOWSPLATFORM
#include <windows.h> // for GetConsoleScreenBufferInfo()
#undef min
#undef max
#endif
using namespace std;

namespace OpenMS
{
	using namespace Exception;

	String TOPPBase::topp_ini_file_ = String(QDir::homePath()) + "/.TOPP.ini";

  void TOPPBase::setMaxNumberOfThreads(int
#ifdef _OPENMP
                                       num_threads  // to avoid the unused warning we enable this
                                                    // argument only if openmp is available
#endif
                                       )
  {
#ifdef _OPENMP
    omp_set_num_threads(num_threads);
#endif
  }

  TOPPBase::TOPPBase(const String& tool_name, const String& tool_description, bool official, bool id_tag_support, const String& version)
  	: tool_name_(tool_name),
  		tool_description_(tool_description),
			id_tag_support_(id_tag_support),
			id_tagger_(tool_name),
			instance_number_(-1),
			version_(version),
			log_type_(ProgressLogger::NONE),
			test_mode_(false),
			debug_level_(-1)
	{
		// if version is empty, use the OpenMS/TOPP version and date/time
		if (version_=="")
		{
			version_ = VersionInfo::getVersion() + " " + VersionInfo::getTime();
		}
		// if the revision info is meaningful, show it as well
		if ( !VersionInfo::getRevision().empty() && VersionInfo::getRevision() != "exported" )
		{
			version_ += String(", Revision: ") + VersionInfo::getRevision() + "";
		}

		//check if tool is in official tools list
    if (official && tool_name_!="GenericWrapper" && !ToolHandler::getTOPPToolList().count(tool_name_))
		{
			writeLog_(String("Warning: Message to maintainer - If '") + tool_name_ + "' is an official TOPP tool, add it to the tools list in ToolHandler. If it is not, set the 'official' flag of the TOPPBase constructor to false.");
		}
	}

	TOPPBase::~TOPPBase()
  {
  	//delete log file if empty
  	StringList log_files;
  	if (!getParam_("log").isEmpty()) log_files << (String)(getParam_("log"));
		for (Size i=0; i< log_files.size(); ++i)
		{
  		if (File::empty(log_files[i]))
  		{
  			File::remove(log_files[i]);
  		}
  	}
	}


	TOPPBase::ExitCodes TOPPBase::main(int argc , const char** argv)
	{
		//----------------------------------------------------------
		//parse command line
		//----------------------------------------------------------

		//register values from derived TOPP tool
		registerOptionsAndFlags_();
		addEmptyLine_();
    //common section for all tools
    if (ToolHandler::getTOPPToolList().count(tool_name_)) addText_("Common TOPP options:");
    else addText_("Common UTIL options:");
		registerStringOption_("ini","<file>","","Use the given TOPP INI file",false);
		registerStringOption_("log","<file>","","Name of log file (created only when specified)",false,true);
		registerIntOption_("instance","<n>",1,"Instance number for the TOPP INI file",false,true);
		registerIntOption_("debug","<n>",0,"Sets the debug level",false, true);
		registerIntOption_("threads", "<n>", 1, "Sets the number of threads allowed to be used by the TOPP tool", false);
		registerStringOption_("write_ini","<file>","","Writes the default configuration file",false);
    registerStringOption_("write_ctd","<out_dir>","","Writes the common tool description file(s) (Toolname(s).ctd) to <out_dir>", false, true);
		registerStringOption_("write_wsdl","<file>","","Writes the default WSDL file",false,true);
		registerFlag_("no_progress","Disables progress logging to command line",true);
		if (id_tag_support_)
		{
			registerStringOption_("id_pool","<file>", "", String("ID pool file to DocumentID's for all generated output files. Disabled by default. (Set to 'main' to use ") + String() + id_tagger_.getPoolFile() + ")", false);
		}
		registerFlag_("test","Enables the test mode (needed for internal use only)", true);
		registerFlag_("-help","Shows options");
    registerFlag_("-helphelp","Shows all options (including advanced)",false);

		// parse command line parameters:
		param_cmdline_.parseCommandLine(argc, argv, parameters_);

		// assign instance number
		*const_cast<int*>(&instance_number_) = getParamAsInt_("instance",1);
		writeDebug_(String("Instance: ")+String(instance_number_),1);

		// assign ini location
		*const_cast<String*>(&ini_location_) = tool_name_+':'+String(instance_number_)+':';
		writeDebug_(String("Ini_location: ")+String(ini_location_),1);

		//set debug level
		debug_level_ = getParamAsInt_("debug",0);
		writeDebug_(String("Debug level: ")+String(debug_level_),1);

		// test if no options were given
		if (argc==1)
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

		// '-test' given
		if (param_cmdline_.exists("test"))
		{
			test_mode_ = true;

			// initialize the random generator as early as possible!
      DateTime date_time;
      date_time.set("1999-12-31 23:59:59");
			UniqueIdGenerator::setSeed(date_time);
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
				  ini_params.load( (String)in_ini );

          // check if ini parameters are applicable to this tool
          checkIfIniParametersAreApplicable_(ini_params);
          // update default params with old params given in -ini and be verbose
          default_params.update(ini_params, true);
        }
			  default_params.store(write_ini_file);
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

			// '-write_wsdl' given
			String wsdl_file("");
			if (param_cmdline_.exists("write_wsdl")) wsdl_file = param_cmdline_.getValue("write_wsdl");
      if (wsdl_file != "")
      {
        return writeWSDL_(wsdl_file);
      }

			//-------------------------------------------------------------
			// load INI file
			//-------------------------------------------------------------
			{
				DataValue value_ini;

        if (param_cmdline_.exists("ini")) value_ini = param_cmdline_.getValue("ini");
				if (!value_ini.isEmpty())
				{
					writeDebug_( "INI file: " + (String)value_ini, 1 );
					writeDebug_( "INI location: " + getIniLocation_(), 1);
          Param p_tmp;
          p_tmp.load( (String)value_ini );

          checkIfIniParametersAreApplicable_(p_tmp);
          // set type on commandline if given in .ini file
          if (p_tmp.exists(this->ini_location_ + "type") && !param_cmdline_.exists("type")) param_cmdline_.setValue("type", p_tmp.getValue(this->ini_location_ + "type"));
          // otherwise the following line cannot fetch subsectionDefaults() from tool
          param_inifile_ = this->getDefaultParameters_();
          Logger::LogStream null_stream;
          param_inifile_.update(p_tmp, false, true, null_stream); // silently update (do not trust INI file), but leave unknown params
        }
        else
        { // fill param with default values:
          param_inifile_ = this->getDefaultParameters_();
        }

        // dissect INI file params:
				// (keep in mind: we currently cannot write all parts of an INI file
				// in the same way we read it in, i.e. the "common" section
				// cannot the written, but only be read
				param_instance_ = param_inifile_.copy( getIniLocation_(), true);
				writeDebug_("Parameters from instance section:",param_instance_,2);
				param_common_tool_ = param_inifile_.copy( "common:"+tool_name_+':', true );
				writeDebug_("Parameters from common section with tool name:",param_common_tool_,2);
				param_common_ = param_inifile_.copy( "common:", true );
				writeDebug_("Parameters from common section without tool name:",param_common_,2);

        param_ = param_cmdline_;
				writeDebug_("Applying defaults to instance section:",param_common_,2);
				param_.setDefaults( param_instance_ );
				writeDebug_("Applying defaults to common section with tool name:",param_common_,2);
				param_.setDefaults( param_common_tool_ );
				writeDebug_("Applying defaults to common section without tool name:",param_common_,2);
				param_.setDefaults( param_common_ );

				// check if all parameters are registered and have the correct type
				checkParam_(param_instance_, (String)value_ini, getIniLocation_());
				checkParam_(param_common_tool_, (String)value_ini, "common:" + tool_name_ + "::");
				checkParam_(param_common_, (String)value_ini, "common:" );

				//check if the version of the parameters file matches the version of this tool
				String file_version = "";
				if (param_inifile_.exists(tool_name_ + ":version"))
				{
					file_version = param_inifile_.getValue(tool_name_ + ":version");
					if (file_version!=VersionInfo::getVersion())
					{
						writeLog_(String("Warning: Parameters file version (") + file_version + ") does not match the version of this tool (" + VersionInfo::getVersion() + ").");
					}
				}
			}

			//-------------------------------------------------------------
			// determine and open the real log file
			//-------------------------------------------------------------

			{
				DataValue const & value_log = getParam_("log");
				if (!value_log.isEmpty())
				{
					writeDebug_( "Log file: " + (String)value_log, 1 );
					log_.close();
					log_.open( ((String)value_log) .c_str(), ofstream::out | ofstream::app);
					writeDebug_("Writing to '"+(String)value_log+'\'',1);
				}
			}

			//-------------------------------------------------------------
			// debug level
			//-------------------------------------------------------------
			debug_level_ = getParamAsInt_("debug",0);
			writeDebug_(String("Debug level (after ini file): ")+String(debug_level_),1);

			//-------------------------------------------------------------
			//progress logging
			//-------------------------------------------------------------
			if(!getFlag_("no_progress"))
			{
				log_type_ = ProgressLogger::CMD;
			}

			//-------------------------------------------------------------
			//document ID tagging
			//-------------------------------------------------------------
			if (id_tag_support_ && getStringOption_("id_pool").length()>0)
			{
				// set custom pool file if given
				if (!(getStringOption_("id_pool")==String("main"))) id_tagger_.setPoolFile(getStringOption_("id_pool"));

				//check if there are enough IDs in the pool (we require at least one and warn below 5)
				Int id_count(0);
				if (!id_tagger_.countFreeIDs(id_count))
				{
					writeLog_("Error: Unable to query ID pool! Ending programm (no computation was performed)!");
					return INTERNAL_ERROR;
				}
				if (id_count == 0)
				{
					writeLog_("Error: No Document IDs in the ID pool. Please restock now! Ending programm (no computation was performed)!");
					return INTERNAL_ERROR;
				}
				else if (id_count <= 5)
				{
					writeLog_("Warning: Less than five(!) Document IDs in the ID pool. Please restock soon!");
				}
			}

			//----------------------------------------------------------
			//threads
			//----------------------------------------------------------
			Int threads = getParamAsInt_("threads", 1);
      TOPPBase::setMaxNumberOfThreads(threads);

			//----------------------------------------------------------
			//main
			//----------------------------------------------------------

			result = main_(argc, argv);

#ifndef DEBUG_TOPP
		}

		//----------------------------------------------------------
		//error handling
		//----------------------------------------------------------

		// Errors caused by the user
		catch(UnableToCreateFile& e)
		{
			writeLog_(String("Error: Unable to write file (") + e.what() + ")");
			writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ")!",1);
			return CANNOT_WRITE_OUTPUT_FILE;
		}
		catch(FileNotFound& e)
		{
			writeLog_(String("Error: File not found (") + e.what() + ")");
			writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !",1);
			return INPUT_FILE_NOT_FOUND;
		}
		catch(FileNotReadable& e)
		{
			writeLog_(String("Error: File not readable (") + e.what() + ")");
			writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !",1);
			return INPUT_FILE_NOT_READABLE;
		}
		catch(FileEmpty& e)
		{
			writeLog_(String("Error: File empty (") + e.what() + ")");
			writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !",1);
			return INPUT_FILE_EMPTY;
		}
		catch(ParseError& e)
		{
			writeLog_(String("Error: Unable to read file (") + e.what() + ")");
			writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !",1);
			return INPUT_FILE_CORRUPT;
		}
		catch(RequiredParameterNotGiven& e)
		{
      String what = e.what();
      if (!what.hasPrefix("'")) what = "'"+what+"'";
      writeLog_(String("Error: The required parameter ") + what + " was not given or is empty!");
			writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !",1);
			return MISSING_PARAMETERS;
		}
		catch(InvalidParameter& e)
		{
			writeLog_(String("Invalid parameter: ") + e.what());
			writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !",1);
			return ILLEGAL_PARAMETERS;
		}
		// Internal errors because of wrong use of this class
		catch(UnregisteredParameter& e)
		{
			writeLog_(String("Internal error: Request for unregistered parameter '") + e.what() + "'");
			writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !",1);
			return INTERNAL_ERROR;
		}
		catch(WrongParameterType& e)
		{
			writeLog_(String("Internal error: Request for parameter with wrong type '") + e.what() + "'");
			writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !",1);
			return INTERNAL_ERROR;
		}
		// All other errors
		catch(BaseException& e)
		{
			writeLog_(String("Error: Unexpected internal error (") + e.what() + ")");
			writeDebug_(String("Error occurred in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !",1);
			return UNKNOWN_ERROR;
		}
#endif

		log_.close();

		return result;
	}

  String TOPPBase::breakString_(const String& input, const Size line_len, const Size indentation, const Size max_lines) const
  {
    StringList result;
    Size short_line_len=line_len-indentation;
    if (short_line_len < 1)
    {
      std::cerr << "INTERNAL ERROR: cannot split lines into empty strings! see breakString_()";
      return input;
    }
    for (Size i=0;i<input.size();)
    {
      String line = input.substr(i, result.size()==0 ? line_len : short_line_len); // first line has full length
      Size advance_size = line.size();
      if (line.hasSubstring("\n"))
      {
        advance_size=0;
        while (line.hasPrefix("\n")){ line=line.substr(1);advance_size++;} // advance by # of \n's
        if (line.hasSubstring("\n")) line = line.prefix('\n');
        advance_size += line.size();  // + actual chars
      }
      i+=advance_size;
      String s_intend = (result.size()==0 ? "" : String(indentation,' ')); // first line no intendation
      String r = s_intend + (result.size()==0 ? line : line.trim()); // intended lines get trimmed
      result.push_back(r);//(r.fillRight(' ', (UInt) line_len));
    }
    if (result.size() > max_lines) // remove lines from end if we get too many (but leave the last one)...
    {
      String last = result.back();
      result.erase(result.begin()+max_lines-2, result.end());
      result.push_back((String(indentation,' ') + String("...")));//.fillRight(' ',(UInt) line_len));
      result.push_back(last);
    }
    // remove last " " from last line to prevent automatic linebreak
    //if (result.size()>0 && result[result.size()-1].hasSuffix(" ")) result[result.size()-1] = result[result.size()-1].substr(0,result[result.size()-1].size()-1);
    return result.concatenate("\n");
  }

  void TOPPBase::readConsoleSize_()
  {
    // avoid calling this function more than once
    static bool been_here = false;
    if (been_here) return;
    been_here = true;

    // determine column width of current console (TODO: put that into a dedicated class)
    try
    {
      console_width_ = -1;
      char * p_env;
      p_env = getenv ("COLUMNS");
      if (p_env!=NULL)
      {
				 console_width_ = String(p_env).toInt();
			}
			else
      {
        writeDebug_("output shaping: COLUMNS env does not exist!", 2);
#ifdef OPENMS_WINDOWSPLATFORM
        HANDLE hOut;
        CONSOLE_SCREEN_BUFFER_INFO SBInfo;
        hOut = GetStdHandle(STD_OUTPUT_HANDLE);
        GetConsoleScreenBufferInfo(hOut, &SBInfo);
        console_width_ = SBInfo.dwSize.X;
#else // Linux / MacOS
			  // try "stty size" command
				// don't use QProcess, as stty will not work there
				FILE *fp = popen("stty size", "r" );
				if (fp != NULL)
				{
					char buff[100];
					fgets( buff, sizeof (buff), fp );
					pclose(fp);					
					String output(buff);
					StringList components;
					output.split(' ', components);
					if (components.size()==2) console_width_=components[1].toInt();
				}
				else
				{
					writeDebug_("output shaping: stty size command failed." , 2);
				}
#endif
			}
      --console_width_; // to add the \n at the end of each line without forcing another line break on windows
    }
    catch (...) {}
    // if console_width_ is still -1, we do not use command line reshaping
    if (console_width_<10)
    {
      writeDebug_("Console width could not be determined or is smaller than 10. Not using output shaping!",2);
      console_width_ = std::numeric_limits<int>::max();
    }
  }

	void TOPPBase::printUsage_()
	{
    readConsoleSize_();
		//common output
		cerr << "\n"
	       << tool_name_ << " -- " << tool_description_ << "\n"
	       << "Version: " << version_ << "\n" << "\n"
	       << "Usage:" << "\n"
				 << "  " << tool_name_ << " <options>" << "\n"
				 << "\n"
				 << "Options (mandatory options marked with '*'):" << "\n";

    // show advanced options?
    bool verbose = getFlag_("-helphelp");

		//determine max length of parameters (including argument) for intendation
		UInt max_size = 0;
		for( vector<ParameterInformation>::const_iterator it = parameters_.begin(); it != parameters_.end(); ++it)
		{
      if ((!it->advanced) || (it->advanced && verbose))
      {
			  max_size = max((UInt)max_size,(UInt)(it->name.size()+it->argument.size()+it->required));
      }
		}

		//offset of the descriptions
		UInt offset = 6 + max_size;
    //keep track of the current subsection we are in, to display the subsection help when a new section starts
    String current_TOPP_subsection("");

    // PRINT parameters && description, restrictions and default
		for( vector<ParameterInformation>::const_iterator it = parameters_.begin(); it != parameters_.end(); ++it)
		{
      if (!((!it->advanced) || (it->advanced && verbose))) continue;

      //new subsection?
			String subsection = getSubsection_(it->name);
			if (!subsection.empty() && current_TOPP_subsection != subsection)
      {
        current_TOPP_subsection = subsection;
        map<String,String>::const_iterator it = subsections_TOPP_.find(current_TOPP_subsection);
        if (it==subsections_TOPP_.end()) throw ElementNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,"'"+current_TOPP_subsection+"' (TOPP subsection not registered)");
        cerr << "\n"; // print newline for new subsection
        cerr << breakString_(it->second, console_width_, 0, 10) << ":\n"; // print subsection description
      }
      else if (subsection.empty() && !current_TOPP_subsection.empty())
      { // subsection ended and normal parameters start again
        current_TOPP_subsection="";
        cerr << "\n"; // print newline to separate ending subsection
      }

			//NAME + ARGUMENT
			String tmp = "  -";
			tmp += it->name + " " + it->argument;
			if (it->required) tmp += '*';
			if (it->type == ParameterInformation::NEWLINE) tmp = "";

			//OFFSET
			tmp.fillRight(' ',offset);
			if (it->type == ParameterInformation::TEXT) tmp = "";

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
						String tmp = it->default_value.toString().substitute(", ", " ");
						if (tmp!="" && tmp!="[]")
						{
							addons.push_back(String("default: '") + tmp + "'");
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
					if (it->valid_strings.size()!=0)
					{
						StringList copy = it->valid_strings;
						for (StringList::iterator str_it = copy.begin();
								 str_it != copy.end(); ++str_it)
						{
							str_it->quote('\'');
						}

						String add = "";
						if (it->type == ParameterInformation::INPUT_FILE || it->type == ParameterInformation::OUTPUT_FILE ||
								it->type == ParameterInformation::INPUT_FILE_LIST || it->type == ParameterInformation::OUTPUT_FILE_LIST) add = " formats";

						addons.push_back(String("valid") + add + ": " + copy.concatenate(", ")); // concat restrictions by comma
					}
					break;
				case ParameterInformation::INT:
				case ParameterInformation::INTLIST:
						if (it->min_int!=-std::numeric_limits<Int>::max())
					{
						addons.push_back(String("min: '") + it->min_int + "'");
					}
					if (it->max_int!=std::numeric_limits<Int>::max())
					{
						addons.push_back(String("max: '") + it->max_int + "'");
					}
					break;
				case ParameterInformation::DOUBLE:
				case ParameterInformation::DOUBLELIST:
					if (it->min_float!=-std::numeric_limits<DoubleReal>::max())
					{
						addons.push_back(String("min: '") + it->min_float + "'");
					}
					if (it->max_float!=std::numeric_limits<DoubleReal>::max())
					{
						addons.push_back(String("max: '") + it->max_float + "'");
					}
					break;

				default:
					break;
			}

			//add DEFAULT and RESTRICTIONS
			if (addons.size()!=0)
			{
        desc_tmp += String(" (") + addons.concatenate(" ") + ")";
			}

      if (it->type == ParameterInformation::TEXT) cerr << breakString_(tmp + desc_tmp, console_width_, 0, 10); // no intendation for text
      else  cerr << breakString_(tmp + desc_tmp, console_width_, offset, 10);
			cerr << "\n";
		}

    // SUBSECTION's at the end
		if (subsections_.size()!=0)
		{
			//determine intendation of description
			UInt indent = 0;
			for(map<String,String>::const_iterator it = subsections_.begin(); it!=subsections_.end(); ++it)
			{
				indent = max((UInt)it->first.size(),indent);
			}
			indent += 6;

			//output
			cerr << "\n"
					 << "The following configuration subsections are valid:" << "\n";
			for(map<String,String>::const_iterator it = subsections_.begin(); it!=subsections_.end(); ++it)
			{
				String tmp = String(" - ") + it->first;
				tmp.fillRight(' ',indent);
        cerr << breakString_(tmp  + it->second, console_width_, indent, 10);
				cerr << "\n";
			}
      cerr << "\n"
           << "You can write an example INI file using the '-write_ini' option." << "\n"
					 << "Documentation of subsection parameters can be found in the" << "\n"
					 << "doxygen documentation or the INIFileEditor." << "\n"
					 << "Have a look at OpenMS documentation for more information." << "\n";
		}
		cerr << endl;
	}

	void TOPPBase::registerParamEntry_(const Param::ParamEntry& entry, const String& argument, const String& full_name)
	{
		String name = full_name.empty() ? entry.name : full_name;
		bool advanced = entry.tags.count("advanced");
		// special case for flags:
		if ((entry.value.valueType() == DataValue::STRING_VALUE) && 
				(entry.value == "false") && (entry.valid_strings.size() == 2) && 
				(entry.valid_strings[0] == "true") && (entry.valid_strings[1] == "false"))
		{
			parameters_.push_back(ParameterInformation(name, ParameterInformation::FLAG, "", "", entry.description, false, advanced));
			return;
		}

		bool input_file = entry.tags.count("input file");
		bool output_file = entry.tags.count("output file");
		if (input_file && output_file)
		{
			throw InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Parameter '" + full_name + "' marked as both input and output file");
		}
		enum ParameterInformation::ParameterTypes type;
		switch (entry.value.valueType())
		{
		case DataValue::STRING_VALUE:
			if (input_file) type = ParameterInformation::INPUT_FILE;
			else if (output_file) type = ParameterInformation::OUTPUT_FILE;
			else type = ParameterInformation::STRING;
			break;
		case DataValue::INT_VALUE:
			type = ParameterInformation::INT;
			break;
		case DataValue::DOUBLE_VALUE:
			type = ParameterInformation::DOUBLE;
			break;
		case DataValue::STRING_LIST:
			if (input_file) type = ParameterInformation::INPUT_FILE_LIST;
      else if (output_file) type = ParameterInformation::OUTPUT_FILE_LIST;
			else type = ParameterInformation::STRINGLIST;
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
		parameters_.push_back(param);
	}


	void TOPPBase::registerFullParam_(const Param& param)
	{
		for (Param::ParamIterator it = param.begin(); it != param.end(); ++it)
		{
			String full_name = it.getName();
			String subsection = getSubsection_(full_name);
			// subsection handling:
			if (!subsection.empty() && (subsections_TOPP_.count(subsection) == 0))
			{
				subsections_TOPP_[subsection] = param.getSectionDescription(subsection);
			}
			// make up a value for "argument":
			String argument;
			switch (it->value.valueType())
			{
			case DataValue::STRING_VALUE:
				if (it->valid_strings.empty()) argument = "<text>"; // name?
				else argument = "<choice>";
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
			registerParamEntry_(*it, argument, full_name);
		}
	}


	void TOPPBase::registerStringOption_(const String& name, const String& argument, const String& default_value,const String& description, bool required, bool advanced)
	{
    if (required && default_value!="") throw InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Registering a required StringOption param ("+name+") with a non-empty default is forbidden!", default_value);
    parameters_.push_back(ParameterInformation(name, ParameterInformation::STRING, argument, default_value, description, required, advanced));
	}

	void TOPPBase::setValidStrings_(const String& name, const std::vector<String>& strings)
	{
		//check for commas
		for (Size i=0; i<strings.size(); ++i)
		{
			if (strings[i].has(','))
			{
				throw InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__,"Comma characters in Param string restrictions are not allowed!");
			}
		}
		//search the right parameter
		for (Size i=0; i<parameters_.size(); ++i)
		{
			if (parameters_[i].name==name)
			{
				//check if the type matches
				if (parameters_[i].type!=ParameterInformation::STRING && parameters_[i].type!=ParameterInformation::STRINGLIST)
				{
					throw ElementNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,name);
				}
        StringList valids = strings;
        StringList defaults;
        if (parameters_[i].type==ParameterInformation::STRING) defaults.push_back(String(parameters_[i].default_value));
        else defaults = parameters_[i].default_value;
        for (Size j=0;j<defaults.size();++j)
        {// allow the empty string even if not in restrictions
          if (defaults[j].size()>0 && !valids.contains(defaults[j]))
          {
            throw InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__,"TO THE DEVELOPER: The TOPP/UTILS tool option '" + name + "' with default value "+String(parameters_[i].default_value) + " does not meet restrictions!");
          }
        }

				parameters_[i].valid_strings = strings;
				return;
			}
		}
		//parameter not found
		throw ElementNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,name);
	}

	void TOPPBase::setValidFormats_(const String& name, const std::vector<String>& formats, const bool force_OpenMS_format)
	{
		//check if formats are known
		if (force_OpenMS_format)
		{
			for (Size i=0; i<formats.size(); ++i)
			{
				if (formats[i] != "fid")
				{
					if (FileHandler::getTypeByFileName(String(".")+formats[i])==FileTypes::UNKNOWN)
					{
						throw InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__,"The file format '" + formats[i] + "' is invalid!");
					}
				}
			}
		}
		//search the right parameter
		for (Size i=0; i<parameters_.size(); ++i)
		{
			if (parameters_[i].name==name)
			{
				//check if the type matches
				if (parameters_[i].type!=ParameterInformation::INPUT_FILE && parameters_[i].type!=ParameterInformation::OUTPUT_FILE && parameters_[i].type!=ParameterInformation::INPUT_FILE_LIST && parameters_[i].type!=ParameterInformation::OUTPUT_FILE_LIST)
				{
					throw ElementNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,name);
        }
        parameters_[i].valid_strings = formats;
				return;
			}
		}
		//parameter not found
		throw ElementNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,name);
	}

	void TOPPBase::setMinInt_(const String& name, Int min)
	{
		//search the right parameter
		for (Size i=0; i<parameters_.size(); ++i)
		{
			if (parameters_[i].name==name)
			{
				//check if the type matches
				if (parameters_[i].type!=ParameterInformation::INT && parameters_[i].type!=ParameterInformation::INTLIST)
				{
					throw ElementNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,name);
				}
        IntList defaults;
        if (parameters_[i].type==ParameterInformation::INT) defaults.push_back(Int(parameters_[i].default_value));
        else defaults = parameters_[i].default_value;
        for (Size j=0;j<defaults.size();++j)
        {
          if (defaults[j] < min)
          {
            throw InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__,"TO THE DEVELOPER: The TOPP/UTILS tool option '" + name + "' with default value "+String(parameters_[i].default_value) + " does not meet restrictions!");
          }
        }
				parameters_[i].min_int = min;
				return;
			}
		}
		//parameter not found
		throw ElementNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,name);
	}

	void TOPPBase::setMaxInt_(const String& name, Int max)
	{
		//search the right parameter
		for (Size i=0; i<parameters_.size(); ++i)
		{
			if (parameters_[i].name==name)
			{
				//check if the type matches
				if (parameters_[i].type!=ParameterInformation::INT && parameters_[i].type!=ParameterInformation::INTLIST)
				{
					throw ElementNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,name);
				}
        IntList defaults;
        if (parameters_[i].type==ParameterInformation::INT) defaults.push_back(Int(parameters_[i].default_value));
        else defaults = parameters_[i].default_value;
        for (Size j=0;j<defaults.size();++j)
        {
          if (defaults[j] > max)
          {
            throw InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__,"TO THE DEVELOPER: The TOPP/UTILS tool option '" + name + "' with default value "+String(parameters_[i].default_value) + " does not meet restrictions!");
          }
        }
				parameters_[i].max_int = max;
				return;
			}
		}
		//parameter not found
		throw ElementNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,name);
	}

	void TOPPBase::setMinFloat_(const String& name, DoubleReal min)
	{
		//search the right parameter
		for (Size i=0; i<parameters_.size(); ++i)
		{
			if (parameters_[i].name==name)
			{
				//check if the type matches
				if (parameters_[i].type!=ParameterInformation::DOUBLE && parameters_[i].type!=ParameterInformation::DOUBLELIST)
				{
					throw ElementNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,name);
				}
        DoubleList defaults;
        if (parameters_[i].type==ParameterInformation::DOUBLE) defaults.push_back(DoubleReal(parameters_[i].default_value));
        else defaults = parameters_[i].default_value;
        for (Size j=0;j<defaults.size();++j)
        {
          if (defaults[j] < min)
          {
            throw InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__,"TO THE DEVELOPER: The TOPP/UTILS tool option '" + name + "' with default value "+String(parameters_[i].default_value) + " does not meet restrictions!");
          }
        }
				parameters_[i].min_float = min;
				return;
			}
		}
		//parameter not found
		throw ElementNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,name);
	}

	void TOPPBase::setMaxFloat_(const String& name, DoubleReal max)
	{
		//search the right parameter
		for (Size i=0; i<parameters_.size(); ++i)
		{
			if (parameters_[i].name==name)
			{
				//check if the type matches
				if (parameters_[i].type!=ParameterInformation::DOUBLE && parameters_[i].type!=ParameterInformation::DOUBLELIST)
				{
					throw ElementNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,name);
				}
        DoubleList defaults;
        if (parameters_[i].type==ParameterInformation::DOUBLE) defaults.push_back(DoubleReal(parameters_[i].default_value));
        else defaults = parameters_[i].default_value;
        for (Size j=0;j<defaults.size();++j)
        {
          if (defaults[j] > max)
          {
            throw InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__,"TO THE DEVELOPER: The TOPP/UTILS tool option '" + name + "' with default value "+String(parameters_[i].default_value) + " does not meet restrictions!");
          }
        }
				parameters_[i].max_float = max;
				return;
			}
		}
		//parameter not found
		throw ElementNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,name);
	}


	void TOPPBase::registerInputFile_(const String& name, const String& argument, const String& default_value,const String& description, bool required, bool advanced, const StringList& tags)
	{
    if (required && default_value!="" && !tags.contains("skipexists")) throw InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Registering a required InputFile param ("+name+") with a non-empty default is forbidden!", default_value);
		parameters_.push_back(ParameterInformation(name, ParameterInformation::INPUT_FILE, argument, default_value, description, required, advanced, tags));
	}

	void TOPPBase::registerOutputFile_(const String& name, const String& argument, const String& default_value,const String& description, bool required, bool advanced)
	{
    if (required && default_value!="") throw InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Registering a required OutputFile param ("+name+") with a non-empty default is forbidden!", default_value);
		parameters_.push_back(ParameterInformation(name, ParameterInformation::OUTPUT_FILE, argument, default_value, description, required, advanced));
	}

	void TOPPBase::registerDoubleOption_(const String& name, const String& argument, double default_value, const String& description, bool required, bool advanced)
	{
    if (required && !boost::math::isnan(default_value))
    {
      throw InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Registering a required Double param (" + name + ") with a non-quiet_NaN default is forbidden! Use 'std::numeric_limits<double>::quiet_NaN()' as default value to fix!", String(default_value));
    }
		parameters_.push_back(ParameterInformation(name, ParameterInformation::DOUBLE, argument, default_value, description, required, advanced));
	}

	void TOPPBase::registerIntOption_(const String& name, const String& argument, Int default_value, const String& description, bool required, bool advanced)
	{
    if (required)
    {
      throw InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Registering an Int param (" + name + ") as 'required' is forbidden! (there is no value to indicate it is missing)!", String(default_value));
    }
		parameters_.push_back(ParameterInformation(name, ParameterInformation::INT, argument, default_value, description, required, advanced));
	}

	void TOPPBase::registerOutputFileList_( const String& name, const String& argument, StringList default_value, const String& description, bool required, bool advanced )
	{
    if (required && default_value.size()>0) throw InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Registering a required OutputFileList param ("+name+") with a non-empty default is forbidden!", default_value.concatenate(","));
		parameters_.push_back(ParameterInformation(name,ParameterInformation::OUTPUT_FILE_LIST,argument,default_value,description,required,advanced));
	}

	void TOPPBase::registerInputFileList_( const String& name, const String& argument, StringList default_value, const String& description, bool required, bool advanced )
	{
    if (required && default_value.size()>0) throw InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Registering a required InputFileList param ("+name+") with a non-empty default is forbidden!", default_value.concatenate(","));
		parameters_.push_back(ParameterInformation(name,ParameterInformation::INPUT_FILE_LIST,argument,default_value,description,required,advanced));
	}

	void TOPPBase::registerStringList_( const String& name, const String& argument, StringList default_value, const String& description, bool required, bool advanced)
	{
    if (required && default_value.size()>0) throw InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Registering a required StringList param ("+name+") with a non-empty default is forbidden!", default_value.concatenate(","));
		parameters_.push_back(ParameterInformation(name,ParameterInformation::STRINGLIST,argument,default_value,description,required,advanced));
	}

	void TOPPBase::registerIntList_( const String& name, const String& argument, IntList default_value, const String& description, bool required, bool advanced )
	{
    stringstream ss;
    ss << default_value;
    if (required && default_value.size()>0) throw InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Registering a required IntList param ("+name+") with a non-empty default is forbidden!",  String(ss.str()) );
		parameters_.push_back(ParameterInformation(name,ParameterInformation::INTLIST,argument,default_value,description,required,advanced));
	}

	void TOPPBase::registerDoubleList_( const String& name, const String& argument, DoubleList default_value, const String& description, bool required, bool advanced)
	{
    stringstream ss;
    ss << default_value;
    if (required && default_value.size()>0) throw InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Registering a required DoubleList param ("+name+") with a non-empty default is forbidden!", String(ss.str()));
		parameters_.push_back(ParameterInformation(name,ParameterInformation::DOUBLELIST,argument,default_value,description,required,advanced));
	}

	void TOPPBase::registerFlag_(const String& name, const String& description, bool advanced)
	{
		parameters_.push_back(ParameterInformation(name, ParameterInformation::FLAG, "", "", description, false, advanced));
	}

	void TOPPBase::addEmptyLine_()
	{
		parameters_.push_back(ParameterInformation("",ParameterInformation::NEWLINE, "", "", "", false, false));
	}

	void TOPPBase::addText_(const String& text)
	{
		parameters_.push_back(ParameterInformation("",ParameterInformation::TEXT, "", "", text, false, false));
	}

	const ParameterInformation& TOPPBase::findEntry_(const String& name) const
	{
		vector<ParameterInformation>::const_iterator it = parameters_.begin();
		while(it != parameters_.end() && it->name!=name)
		{
			++it;
		}
		if (it == parameters_.end())
		{
			throw UnregisteredParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
		return *it;
	}

	String TOPPBase::getStringOption_(const String& name) const
	{
		const ParameterInformation& p = findEntry_(name);
		if (p.type!=ParameterInformation::STRING && p.type!=ParameterInformation::INPUT_FILE && p.type!=ParameterInformation::OUTPUT_FILE)
		{
			throw WrongParameterType(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
    String message = "'"+name+"'";
    if (p.valid_strings.size()>0)
    {
      message += " [valid: " + StringList(p.valid_strings).concatenate(", ") + "]";
    }
		if (p.required && (getParam_(name).isEmpty() || getParam_(name)==""))
		{
			throw RequiredParameterNotGiven(__FILE__,__LINE__,__PRETTY_FUNCTION__, message);
		}
		String tmp = getParamAsString_(name, p.default_value);
		writeDebug_(String("Value of string option '") + name + "': " + tmp, 1);

		// if required or set by user, do some validity checks
		if (p.required || ( !getParam_(name).isEmpty() && tmp!=p.default_value))
		{
			//check if files are readable/writable
			if (p.type==ParameterInformation::INPUT_FILE)
			{
				if (!p.tags.contains("skipexists")) inputFileReadable_(tmp, name);
			}
			else if (p.type==ParameterInformation::OUTPUT_FILE)
			{
				outputFileWritable_(tmp, name);
			}

			//check restrictions
			if (p.valid_strings.size()!=0)
			{
				if (p.type==ParameterInformation::STRING)
				{
					if (find(p.valid_strings.begin(),p.valid_strings.end(),tmp)==p.valid_strings.end())
					{
            String valid_strings = StringList(p.valid_strings).concatenate("', '");
						throw InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, String("Invalid value '") + tmp + "' for string parameter '" + name + "' given. Valid strings are: '" + valid_strings + "'.");
					}
				}
				else if (p.type==ParameterInformation::INPUT_FILE)
				{
					if (!p.tags.contains("skipexists")) inputFileReadable_(tmp, name);

					//create upper case list of valid formats
					StringList formats = p.valid_strings;
					formats.toUpper();
					//determine file type as string
					FileHandler fh;
					String format = FileHandler::typeToName(FileHandler::getTypeByFileName(tmp)).toUpper();
					bool invalid = false;
					//Wrong or unknown ending
					if (!formats.contains(format))
					{
						if (format=="UNKNOWN") //Unknown ending => check content
						{
							format = FileHandler::typeToName(FileHandler::getTypeByContent(tmp)).toUpper();
							if (!formats.contains(format))
							{
								if (format=="UNKNOWN") //Unknown format => warning as this might by the wrong format
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
						valid_formats.concatenate(p.valid_strings.begin(),p.valid_strings.end(),"','");
						throw InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, String("Input file '" + tmp + "' has invalid format '") + format + "'. Valid formats are: '" + valid_formats + "'.");
					}
				}
				else if (p.type==ParameterInformation::OUTPUT_FILE)
				{
					outputFileWritable_(tmp, name);

					//create upper case list of valid formats
					StringList formats = p.valid_strings;
					formats.toUpper();
					//determine file type as string
					FileHandler fh;
					String format = FileHandler::typeToName(FileHandler::getTypeByFileName(tmp)).toUpper();
					//Wrong or unknown ending
					if (!formats.contains(format) && format!="UNKNOWN")
					{
						String valid_formats = "";
						valid_formats.concatenate(p.valid_strings.begin(),p.valid_strings.end(),"','");
						throw InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, String("Invalid output file extension '") + tmp + "'. Valid file extensions are: '" + valid_formats + "'.");
					}
				}
			}
		}

		return tmp;
	}

	DoubleReal TOPPBase::getDoubleOption_(const String& name) const
	{
		const ParameterInformation& p = findEntry_(name);
		if (p.type != ParameterInformation::DOUBLE)
		{
			throw WrongParameterType(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
    if (p.required && getParam_(name).isEmpty() )
		{
			throw RequiredParameterNotGiven(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
		double tmp = getParamAsDouble_(name, (DoubleReal)p.default_value);
    if (p.required && boost::math::isnan(tmp ))
		{
			throw RequiredParameterNotGiven(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
		writeDebug_(String("Value of double option '") + name + "': " + String(tmp), 1);

		//check if in valid range
		if (p.required || ( !getParam_(name).isEmpty() && tmp!=(DoubleReal)p.default_value))
		{
			if (tmp<p.min_float || tmp>p.max_float)
			{
				throw InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, String("Invalid value '") + tmp + "' for float parameter '" + name + "' given. Out of valid range: '" + p.min_float + "'-'" + p.max_float + "'.");
			}
		}

		return tmp;
	}

	Int TOPPBase::getIntOption_(const String& name) const
	{
		const ParameterInformation& p = findEntry_(name);
		if (p.type != ParameterInformation::INT)
		{
			throw WrongParameterType(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
		if (p.required && getParam_(name).isEmpty() )
		{
			throw RequiredParameterNotGiven(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
		Int tmp = getParamAsInt_(name, (Int)p.default_value);
    // not checking if NAN here (as done with DoubleReal, as NAN is not supported for Int)
		writeDebug_(String("Value of int option '") + name + "': " + String(tmp), 1);

		//check if in valid range
		if (p.required || ( !getParam_(name).isEmpty() && tmp!=(Int)p.default_value))
		{
			if (tmp<p.min_int || tmp>p.max_int)
			{
				throw InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, String("Invalid value '") + tmp + "' for integer parameter '" + name + "' given. Out of valid range: '" + p.min_int + "'-'" + p.max_int + "'.");
			}
		}

		return tmp;
	}

	StringList TOPPBase::getStringList_(const String& name) const
	{
		const ParameterInformation& p = findEntry_(name);
		if (p.type!=ParameterInformation::STRINGLIST && p.type!=ParameterInformation::INPUT_FILE_LIST && p.type!=ParameterInformation::OUTPUT_FILE_LIST)
		{
			throw WrongParameterType(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
		if (p.required && getParam_(name).isEmpty() )
		{
			throw RequiredParameterNotGiven(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
		StringList tmp_list = getParamAsStringList_(name, (StringList)p.default_value);
		if (p.required && tmp_list.size()==0)
		{
			throw RequiredParameterNotGiven(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}

    for (StringList::iterator it = tmp_list.begin(); it < tmp_list.end(); ++it)
		{
			String tmp(*it);
			writeDebug_(String("Value of string option '") + name + "': " + tmp, 1);

			// if required or set by user, do some validity checks
			if (p.required || ( !getParam_(name).isEmpty()&&tmp_list!=p.default_value))
			{
				//check if files are readable/writeable
				if (p.type==ParameterInformation::INPUT_FILE_LIST)
				{
					inputFileReadable_(tmp, name);
				}
				else if (p.type==ParameterInformation::OUTPUT_FILE_LIST)
				{
					outputFileWritable_(tmp, name);
				}

				//check restrictions
				if (p.valid_strings.size()!=0)
				{
					if (p.type==ParameterInformation::STRINGLIST)
					{
						if (find(p.valid_strings.begin(),p.valid_strings.end(),tmp)==p.valid_strings.end())
						{
							String valid_strings = "";
							valid_strings.concatenate(p.valid_strings.begin(),p.valid_strings.end(),"','");
							throw InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, String("Invalid value '") + tmp + "' for string parameter '" + name + "' given. Valid strings are: '" + valid_strings + "'.");
						}
					}
					else if (p.type==ParameterInformation::INPUT_FILE_LIST)
					{
						inputFileReadable_(tmp, name);

						//create upper case list of valid formats
						StringList formats = p.valid_strings;
						formats.toUpper();
						//determine file type as string
						FileHandler fh;
						String format = FileHandler::typeToName(FileHandler::getTypeByFileName(tmp)).toUpper();
						bool invalid = false;
						//Wrong or unknown ending
						if (!formats.contains(format))
						{
							if (format=="UNKNOWN") //Unknown ending => check content
							{
								format = FileHandler::typeToName(FileHandler::getTypeByContent(tmp)).toUpper();
								if (!formats.contains(format))
								{
									if (format=="UNKNOWN") //Unknown format => warning as this might by the wrong format
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
							valid_formats.concatenate(p.valid_strings.begin(),p.valid_strings.end(),"','");
							throw InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, String("Input file '" + tmp + "' has invalid format '") + format + "'. Valid formats are: '" + valid_formats + "'.");
						}
					}
					else if (p.type==ParameterInformation::OUTPUT_FILE_LIST)
					{
						outputFileWritable_(tmp, name);

						//create upper case list of valid formats
						StringList formats = p.valid_strings;
						formats.toUpper();
						//determine file type as string
						FileHandler fh;
						String format = FileHandler::typeToName(FileHandler::getTypeByFileName(tmp)).toUpper();
						//Wrong or unknown ending
						if (!formats.contains(format) && format!="UNKNOWN")
						{
							String valid_formats = "";
							valid_formats.concatenate(p.valid_strings.begin(),p.valid_strings.end(),"','");
							throw InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, String("Invalid output file extension '") + tmp + "'. Valid file extensions are: '" + valid_formats + "'.");
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
			throw WrongParameterType(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
		if (p.required && getParam_(name).isEmpty() )
		{
			throw RequiredParameterNotGiven(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
		DoubleList tmp_list = getParamAsDoubleList_(name, (DoubleList)p.default_value);
		if (p.required && tmp_list.size()==0)
		{
			throw RequiredParameterNotGiven(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
		DoubleReal tmp;
		for(DoubleList::iterator it = tmp_list.begin(); it < tmp_list.end(); ++it)
		{
			tmp = *it;
			writeDebug_(String("Value of string option '") + name + "': " + String(tmp), 1);

			//check if in valid range
			if (p.required || ( !getParam_(name).isEmpty() && tmp_list!=(DoubleList)p.default_value))
			{
				if (tmp<p.min_float || tmp>p.max_float)
				{
					throw InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, String("Invalid value '") + tmp + "' for float parameter '" + name + "' given. Out of valid range: '" + p.min_float + "'-'" + p.max_float + "'.");
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
			throw WrongParameterType(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
		if (p.required && getParam_(name).isEmpty() )
		{
			throw RequiredParameterNotGiven(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
		IntList tmp_list = getParamAsIntList_(name, (IntList)p.default_value);
		if (p.required && tmp_list.size()==0)
		{
			throw RequiredParameterNotGiven(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}

		Int tmp;
		for(IntList::iterator it = tmp_list.begin(); it < tmp_list.end();++it)
		{
			tmp = *it;
			writeDebug_(String("Value of string option '") + name + "': " + String(tmp), 1);

			//check if in valid range
			if (p.required || ( !getParam_(name).isEmpty() && tmp_list!=(IntList)p.default_value))
			{
				if (tmp<p.min_int || tmp>p.max_int)
				{
					throw InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, String("Invalid value '") + tmp + "' for integer parameter '" + name + "' given. Out of valid range: '" + p.min_int + "'-'" + p.max_int + "'.");
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
			throw WrongParameterType(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
		bool tmp = getParamAsBool_(name);
		writeDebug_(String("Value of string option '") + name + "': " + String(tmp), 1);
		return tmp;
	}

	void TOPPBase::writeLog_(const String& text) const
	{
		LOG_INFO << text << endl;
		enableLogging_();
		log_ << QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss").toStdString() << ' ' << getIniLocation_() << ": " << text<< endl;
	}

	void TOPPBase::writeDebug_(const String& text, UInt min_level) const
	{
		if (debug_level_>=(Int)min_level)
		{
		  writeLog_(text);
		}
	}

	void TOPPBase::writeDebug_(const String& text, const Param& param, UInt min_level) const
	{
		if (debug_level_>=(Int)min_level)
		{
			LOG_DEBUG << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl
						    << QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss").toStdString() << ' ' << getIniLocation_() << " " << text<< endl
						    << param
						    << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl;
			enableLogging_();
			log_  << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl
						<< QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss").toStdString() << ' ' << getIniLocation_() << " " << text<< endl
						<< param
						<< " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl;
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
			if (tmp.valueType()==DataValue::INT_VALUE)
			{
				return (Int)tmp;
			}
			throw WrongParameterType(__FILE__, __LINE__, __PRETTY_FUNCTION__, key);
		}
		else
		{
			return default_value;
		}
	}

	DoubleReal TOPPBase::getParamAsDouble_(const String& key, DoubleReal default_value) const
	{
		const DataValue& tmp = getParam_(key);
		if (!tmp.isEmpty())
		{
			if (tmp.valueType()==DataValue::DOUBLE_VALUE)
			{
				return (DoubleReal)tmp;
			}
			throw WrongParameterType(__FILE__, __LINE__, __PRETTY_FUNCTION__, key);
		}
		else
		{
			return default_value;
		}
	}

	StringList TOPPBase::getParamAsStringList_(const String& key, const StringList& default_value) const
	{
		const DataValue& tmp = getParam_(key);
		if(!tmp.isEmpty())
		{
			return (StringList)tmp;
		}
		else
		{
			return default_value;
		}
	}

	IntList TOPPBase::getParamAsIntList_(const String& key,const IntList& default_value) const
	{
		const DataValue& tmp = getParam_(key);
		if(!tmp.isEmpty())
		{
			if (tmp.valueType()==DataValue::INT_LIST)
			{
				return (IntList)tmp;
			}
			throw WrongParameterType(__FILE__, __LINE__, __PRETTY_FUNCTION__, key);
		}
		else
		{
			return default_value;
		}
	}

	DoubleList TOPPBase::getParamAsDoubleList_(const String& key,const DoubleList& default_value) const
	{
		const DataValue& tmp = getParam_(key);
		if(!tmp.isEmpty())
		{
			if (tmp.valueType()==DataValue::DOUBLE_LIST)
			{
				return (DoubleList)tmp;
			}
			throw WrongParameterType(__FILE__, __LINE__, __PRETTY_FUNCTION__, key);
		}
		else
		{
			return default_value;
		}
	}

	bool TOPPBase::getParamAsBool_(const String& key) const
	{
		DataValue tmp = getParam_(key);
		if (tmp.valueType()==DataValue::EMPTY_VALUE)
		{
			return false;
		}
		else if (tmp.valueType()==DataValue::STRING_VALUE)
		{
			if ((String)tmp=="false")
			{
				return false;
			}
			else if ((String)tmp=="true")
			{
				return true;
			}
		}
		throw InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__,String("Invalid value '") + tmp.toString() + "' for flag parameter '" + key +"'. Valid values are 'true' and 'false' only.");
	}

	DataValue const& TOPPBase::getParam_(const String& key) const
	{
		// look up in command line
		{
			if (param_cmdline_.exists(key))
			{
				writeDebug_(String("Parameter '")+key+String("' from COMMAND LINE: ")+String(param_cmdline_.getValue(key)),3);
				return param_cmdline_.getValue(key);
			}
		}

		// look up in instance section
		{
			if (param_instance_.exists(key))
			{
				writeDebug_(String("Parameter '")+key+String("' from INSTANCE SECTION: ")+String(param_instance_.getValue(key)),3);
				return param_instance_.getValue(key);
			}
		}

		// look up in common section with tool name
		{
			if (param_common_tool_.exists(key))
			{
				writeDebug_(String("Parameter '")+key+String("' from COMMON SECTION (TOOL SPECIFIC): ")+String(param_common_tool_.getValue(key)),3);
				return param_common_tool_.getValue( key );
			}
		}

		// look up in common section without tool name
		{
			if (param_common_.exists(key))
			{
				writeDebug_(String("Parameter '")+key+String("' from COMMON SECTION: ")+String(param_common_.getValue(key)),3);
				return param_common_.getValue( key );
			}
		}

		// if look up fails everywhere, return EMPTY
		writeDebug_( String("Parameter '")+key+String("' not found."), 1 );
		return DataValue::EMPTY;
	}

	Param const& TOPPBase::getParam_() const
	{
		return param_;
	}

	String TOPPBase::getSubsection_(const String& name) const
	{
		size_t pos = name.find_last_of(':');
		if (pos == string::npos) return ""; // delimiter not found
		return name.substr(0, pos);
	}

	void TOPPBase::enableLogging_() const
	{
		if ( !log_.is_open() )
		{
			String log_destination = "";
			if(param_cmdline_.exists("log")) log_destination = param_cmdline_.getValue("log");
			if ( log_destination!="" )
			{
				log_.open( log_destination.c_str(), ofstream::out | ofstream::app);
				if (debug_level_>=1)
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
		for (Param::ParamIterator it = param.begin(); it!=param.end(); ++it)
		{
			// subsections (do not check content, but warn if not registered)
			String subsection = getSubsection_(it.getName());
			if (!subsection.empty() && subsections_TOPP_.count(subsection) == 0) // not found in TOPP subsections
			{
				// for multi-level subsections, check only the first level:
				if (subsections_.count(subsection.substr(0, subsection.find(':'))) == 0) // not found in normal subsections
				{
					if (!(location == "common::" && subsection==tool_name_) )
					{
						writeLog_("Warning: Unknown subsection '" + subsection + "' in '" + filename + "' (location '"+location+"')!");
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
						if (it->value.valueType()!=DataValue::STRING_VALUE)
						{
							writeLog_("Warning: Wrong parameter type of '" + location + it.getName() + "' in '" + filename + "'. Type should be 'string'!");
						}
						break;
					case ParameterInformation::DOUBLE:
						if (it->value.valueType()!=DataValue::DOUBLE_VALUE)
						{
							writeLog_("Warning: Wrong  parameter type of '" + location + it.getName() + "' in '" + filename + "'. Type should be 'double'!");
						}
						break;
					case ParameterInformation::INT:
						if (it->value.valueType()!=DataValue::INT_VALUE)
						{
							writeLog_("Warning: Wrong parameter type of '" + location + it.getName() + "' in '" + filename + "'. Type should be 'int'!");
						}
						break;
					case ParameterInformation::STRINGLIST:
					case ParameterInformation::INPUT_FILE_LIST:
					case ParameterInformation::OUTPUT_FILE_LIST:
						if (it->value.valueType()!=DataValue::STRING_LIST)
						{
							writeLog_("Warning: Wrong parameter type of '" + location + it.getName() + "' in '" + filename + "'. Type should be 'string list'!");
						}
						break;
					case ParameterInformation::INTLIST:
						if (it->value.valueType()!=DataValue::INT_LIST)
						{
							writeLog_("Warning: Wrong parameter type of '" + location + it.getName() + "' in '" + filename + "'. Type should be 'int list'!");
						}
						break;
					case ParameterInformation::DOUBLELIST:
						if (it->value.valueType()!=DataValue::DOUBLE_LIST)
						{
							writeLog_("Warning: Wrong parameter type of '" + location + it.getName() + "' in '" + filename + "'. Type should be 'double list'!");
						}
						break;
					default:
						break;
				};
			}
			catch (UnregisteredParameter)
			{
				writeLog_("Warning: Unknown parameter '" + location + it.getName() + "' in '" + filename + "'!");
			}
		}
	}

  void TOPPBase::checkIfIniParametersAreApplicable_(const Param &ini_params)
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
		writeDebug_( "Checking input file '" + filename + "'", 2 );

    // prepare error message
    String message;
    if (param_name == "") message = "Cannot read input file!\n";
    else message = "Cannot read input file given from parameter '-" + param_name + "'!\n";

    // check file
    if (!File::exists(filename))
	  {
      LOG_ERROR << message;
		  throw FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
	  }
	  if (!File::readable(filename))
	  {
      LOG_ERROR << message;
		  throw FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
	  }
    if (!File::isDirectory(filename) && File::empty(filename))
    {
      LOG_ERROR << message;
      throw FileEmpty(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }
	}

	void TOPPBase::outputFileWritable_(const String& filename, const String& param_name) const
	{
		writeDebug_( "Checking output file '" + filename + "'", 2 );

    // prepare error message
    String message;
    if (param_name == "") message = "Cannot write output file!\n";
    else message = "Cannot write output file given from parameter '-" + param_name + "'!\n";

    if (!File::writable(filename))
		{
      LOG_ERROR << message;
			throw UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
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

	void TOPPBase::parseRange_(const String& text, double& low, double& high) const
	{
		String tmp = text.prefix(':');
		if (tmp!="")
		{
			low = tmp.toDouble();
		}
		tmp = "";
		tmp = text.suffix(':');
		if (tmp!="")
		{
			high = tmp.toDouble();
		}
	}

	Param TOPPBase::getSubsectionDefaults_(const String& /*section*/) const
	{
		throw NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);

		return Param();
	}

	Param TOPPBase::getDefaultParameters_() const
	{
		Param tmp;
		String loc = tool_name_ + ":" + String(instance_number_) + ":";
		//parameters
		for( vector<ParameterInformation>::const_iterator it = parameters_.begin(); it != parameters_.end(); ++it)
		{
			if (it->name=="ini" || it->name=="-help" || it->name=="-helphelp" || it->name=="instance" || it->name=="write_ini" || it->name=="write_wsdl" || it->name=="write_ctd")
			{ // do not store those params in ini file
        continue;
      }
			String name = loc + it->name;
			StringList tags;
			if (it->advanced) tags.push_back("advanced");
      if (it->required) tags.push_back("required");
			if (it->type == ParameterInformation::INPUT_FILE || it->type == ParameterInformation::INPUT_FILE_LIST) tags.push_back("input file");
			if (it->type == ParameterInformation::OUTPUT_FILE || it->type == ParameterInformation::OUTPUT_FILE_LIST) tags.push_back("output file");
			switch (it->type)
			{
				case ParameterInformation::STRING:
          tmp.setValue(name,(String)it->default_value, it->description, tags);
          if (it->valid_strings.size() != 0)
          {
            tmp.setValidStrings(name,it->valid_strings);
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
					if (it->min_float!=-std::numeric_limits<DoubleReal>::max())
					{
						tmp.setMinFloat(name, it->min_float);
					}
					if (it->max_float!=std::numeric_limits<DoubleReal>::max())
					{
						tmp.setMaxFloat(name, it->max_float);
					}
					break;
				case ParameterInformation::INT:
					tmp.setValue(name,(Int)it->default_value, it->description, tags);
					if (it->min_int!=-std::numeric_limits<Int>::max())
					{
						tmp.setMinInt(name, it->min_int);
					}
					if (it->max_int!=std::numeric_limits<Int>::max())
					{
						tmp.setMaxInt(name, it->max_int);
					}
					break;
				case ParameterInformation::FLAG:
					tmp.setValue(name,"false", it->description, tags);
					tmp.setValidStrings(name,StringList::create("true,false"));
					break;
				case ParameterInformation::INPUT_FILE_LIST:
				case ParameterInformation::OUTPUT_FILE_LIST:
          tmp.setValue(name,(StringList)it->default_value, it->description, tags);
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
					tmp.setValue(name,(StringList)it->default_value, it->description, tags);
					if (it->valid_strings.size()!=0)
					{
						tmp.setValidStrings(name,it->valid_strings);
					}
					break;
				case ParameterInformation::INTLIST:
					tmp.setValue(name,(IntList)it->default_value, it->description, tags);
					if (it->min_int!=-std::numeric_limits<Int>::max())
					{
						tmp.setMinInt(name, it->min_int);
					}
					if (it->max_int!=std::numeric_limits<Int>::max())
					{
						tmp.setMaxInt(name, it->max_int);
					}
					break;
				case ParameterInformation::DOUBLELIST:
					tmp.setValue(name,(DoubleList)it->default_value, it->description, tags);
					if (it->min_float!=-std::numeric_limits<DoubleReal>::max())
					{
						tmp.setMinFloat(name, it->min_float);
					}
					if (it->max_float!=std::numeric_limits<DoubleReal>::max())
					{
						tmp.setMaxFloat(name, it->max_float);
					}
					break;
				default:
					break;
			}
		}

		//subsections intrinsic to TOPP tool (i.e. a commandline param with a ':')
		for(map<String,String>::const_iterator it = subsections_TOPP_.begin(); it!=subsections_TOPP_.end(); ++it)
		{
			tmp.setSectionDescription(loc + it->first, it->second);
		}

		// set tool version
		tmp.setValue(tool_name_ + ":version", VersionInfo::getVersion(), "Version of the tool that generated this parameters file.", StringList::create("advanced"));

		// Descriptions
		tmp.setSectionDescription(tool_name_, tool_description_);
		tmp.setSectionDescription(tool_name_ + ":" + String(instance_number_), String("Instance '") + String(instance_number_) + "' section for '" + tool_name_ + "'");

    // add type (as default type is "", but .ini file should have it)
    if (param_cmdline_.exists("type")) tmp.setValue(loc + "type", param_cmdline_.getValue("type"));

		// Subsections
		for(map<String,String>::const_iterator it = subsections_.begin(); it!=subsections_.end(); ++it)
		{
			Param tmp2 = getSubsectionDefaults_(it->first);
			if (!tmp2.empty())
			{
				tmp.insert(loc + it->first + ":",tmp2);
				tmp.setSectionDescription(loc + it->first, it->second);
			}
		}

		// 2nd stage, use TOPP tool defaults from home (if existing)
		Param tool_user_defaults(getToolUserDefaults_(tool_name_));
		tmp.update(tool_user_defaults);

		// 3rd stage, use OpenMS.ini from library to override settings
		Param system_defaults(File::getSystemParameters());
    // this currently writes to the wrong part of the ini-file (revise) or remove altogether:
    //   there should be no section which already contains these params (-> thus a lot of warnings are emitted)
    //   furthermore, entering those params will not allow us to change settings in OpenMS.ini and expect them to be effective, as they will be overridden by the tools' ini file
		//tmp.update(system_defaults);

		return tmp;
	}

	Param TOPPBase::getToolUserDefaults_(const String& tool_name) const
	{
		Param p;
		String ini_name(File::getUserDirectory() + "/" + tool_name + ".ini");
		if (File::readable(ini_name))
		{
			p.load(ini_name);
		}
		return p;
	}

	const DocumentIDTagger& TOPPBase::getDocumentIDTagger_() const
	{
		if (!id_tag_support_)
		{
			writeLog_(String("Error: Message to maintainer - You created your TOPP tool without id_tag_support and query the ID Pool class! Decide what you want!"));
			exit(INTERNAL_ERROR);
		}
		else if (id_tag_support_ && getStringOption_("id_pool").length()==0)
		{
			writeLog_(String("Error: Message to maintainer - You created your TOPP tool with id_tag_support and query the ID Pool class without the user actually requesting it (-id_pool is not set)!"));
			exit(INTERNAL_ERROR);
		}
		return id_tagger_;
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
			p.setMetaValue("parameter: mode" , "test_mode");
		}
		else
		{
			//version
			p.getSoftware().setVersion(VersionInfo::getVersion());
			//time
			p.setCompletionTime(DateTime::now());
			//parameters
			const Param& param = getParam_();
			for (Param::ParamIterator it=param.begin(); it!=param.end(); ++it)
			{
			   p.setMetaValue(String("parameter: ") + it.getName() , it->value);
			}
		}

		return p;
	}

  void TOPPBase::addDataProcessing_(ConsensusMap& map, const DataProcessing& dp) const
  {
  	map.getDataProcessing().push_back(dp);

  	//remove abolute map paths
  	if (test_mode_)
		{
			for (Size d=0; d<map.getFileDescriptions().size(); ++d)
			{
				map.getFileDescriptions()[d].filename = File::basename(map.getFileDescriptions()[d].filename);
			}
		}
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
    if (type_list.size()==0) type_list.push_back(""); // no type for most tools (except GenericWrapper)

    for (Size i=0; i<type_list.size(); ++i)
    {
      QString write_ctd_file = out_dir_str + QDir::separator() + tool_name_.toQString() + type_list[i].toQString() + ".ctd";
      outputFileWritable_(write_ctd_file, "write_ctd");

      // set type on commandline, so that getDefaultParameters_() does not fail (as it calls getSubSectionDefaults() of tool)
      if (type_list[i]!="") param_cmdline_.setValue("type", type_list[i]);
      Param default_params = getDefaultParameters_();
      // add type to ini file
      if (type_list[i]!="") default_params.setValue(this->ini_location_ + "type", type_list[i]);

      std::stringstream* ss = new std::stringstream();
      default_params.writeXMLToStream(ss);
      String ini_file_str(ss->str());

      //morph to ctd format
      QStringList lines = ini_file_str.toQString().split("\n", QString::SkipEmptyParts);
      lines.replace(0, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
      lines.insert(1, "<tool status=\"internal\">");
      lines.insert(2, QString("<name>")+tool_name_.toQString()+"</name>");
      lines.insert(3, QString("<version>")+VersionInfo::getVersion().toQString()+"</version>");
      lines.insert(4, QString("<description><![CDATA[")+tool_description_.toQString()+"]]></description>");
      QString html_doc = tool_description_.toQString();
      lines.insert(5, QString("<manual><![CDATA[")+html_doc+"]]></manual>");
      lines.insert(6, "<docurl>http://www-bs2.informatik.uni-tuebingen.de/services/OpenMS/OpenMS-release/html/TOPP__"+tool_name_.toQString()+".html</docurl>");
      lines.insert(7, "<category>"+ToolHandler::getCategory(tool_name_).toQString()+"</category>");
      lines.insert(8, "<type></type>");
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

  TOPPBase::ExitCodes TOPPBase::writeWSDL_(const String& filename)
  {
    outputFileWritable_(filename, "write_wsdl");
    ofstream os(filename.c_str());

    //write header
    os << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
    os << "<wsdl:definitions targetNamespace=\"http://www-bs.informatik.uni-tuebingen.de/compas\" xmlns:ns1=\"http://org.apache.axis2/xsd\" xmlns:plnk=\"http://schemas.xmlsoap.org/ws/2003/05/partner-link/\" xmlns:soap=\"http://schemas.xmlsoap.org/wsdl/soap/\" xmlns:tns=\"http://www-bs.informatik.uni-tuebingen.de/compas\" xmlns:wsdl=\"http://schemas.xmlsoap.org/wsdl/\" xmlns:xs=\"http://www.w3.org/2001/XMLSchema\">" << endl;
    os << "  <wsdl:types>" << endl;
    os << "    <xs:schema attributeFormDefault=\"unqualified\" elementFormDefault=\"qualified\" targetNamespace=\"http://org.apache.axis2/xsd\" xmlns:ns1=\"http://org.apache.axis2/xsd\" xmlns:xs=\"http://www.w3.org/2001/XMLSchema\">" << endl;
    os << "      <xs:element name=\"" << tool_name_ << "Request\">" << endl;
    os << "        <xs:complexType>" << endl;
    os << "          <xs:sequence>" << endl;

    //write types (forward declaration for readability only. Could be defined in the message as well.
    Param param = getDefaultParameters_();
    param = param.copy(tool_name_ + ":" + String(instance_number_) + ":",true);
    for (Param::ParamIterator it=param.begin(); it!=param.end(); ++it)
    {
      //find out if the value is restricted
      bool restricted = false;
      if (it->value.valueType()==DataValue::STRING_VALUE  && !it->valid_strings.empty())
      {
        restricted = true;
      }
      else if (it->value.valueType()==DataValue::STRING_LIST || it->value.valueType()==DataValue::INT_LIST || it->value.valueType()==DataValue::DOUBLE_LIST)
      {
        restricted = true;
      }
      else if (it->value.valueType()==DataValue::INT_VALUE && (it->min_int!=-std::numeric_limits<Int>::max() || it->max_int!=std::numeric_limits<Int>::max()))
      {
        restricted = true;
      }
      else if (it->value.valueType()==DataValue::DOUBLE_VALUE && (it->min_float!=-std::numeric_limits<DoubleReal>::max() || it->max_float!=std::numeric_limits<DoubleReal>::max()))
      {
        restricted = true;
      }

      //name, default (and type if not restricted)
      os << "            <xs:element name=\"" << it.getName() << "\"";
      if (!restricted)
      {
        if (it->value.valueType()==DataValue::STRING_VALUE) os << " type=\"xs:string\"";
        if (it->value.valueType()==DataValue::DOUBLE_VALUE) os << " type=\"xs:double\"";
        if (it->value.valueType()==DataValue::INT_VALUE) os << " type=\"xs:integer\"";
      }
      os << " default=\"" << it->value.toString() << "\">" << endl;
      //docu
      if (it->description!="")
      {
        String description = it->description;
        description.substitute("<","&lt;");
        description.substitute(">","&gt;");
        os << "              <xs:annotation>" << endl;
        os << "                <xs:documentation>" << description << "</xs:documentation>" << endl;
        os << "              </xs:annotation>" << endl;
      }
      //restrictions
      if (restricted)
      {
        os << "              <xs:simpleType>" << endl;
        if (it->value.valueType()==DataValue::STRING_LIST)
        {
          os << "                <xs:restriction base=\"xs:stringlist\">" << endl;
          os << "                  <xs:pattern value=\"^$|[^,](,[^,]+)*\"/>" << endl;
        }
        else if (it->value.valueType()==DataValue::INT_LIST)
        {
          os << "                <xs:restriction base=\"xs:intlist\">" << endl;
          os << "                  <xs:pattern value=\"^$|[^,](,[^,]+)*\"/>" << endl;
        }
        else if (it->value.valueType()==DataValue::DOUBLE_LIST)
        {
          os << "                <xs:restriction base=\"xs:doublelist\">" << endl;
          os << "                  <xs:pattern value=\"^$|[^,](,[^,]+)*\"/>" << endl;
        }
        else if (it->value.valueType()==DataValue::STRING_VALUE)
        {
          os << "                <xs:restriction base=\"xs:string\">" << endl;
          for (Size i=0; i<it->valid_strings.size(); ++i)
          {
            os << "                  <xs:enumeration value=\"" << it->valid_strings[i] << "\"/>" << endl;
          }
        }
        else if (it->value.valueType()==DataValue::DOUBLE_VALUE)
        {
          os << "                <xs:restriction base=\"xs:double\">" << endl;
          if (it->min_float!=-std::numeric_limits<DoubleReal>::max())
          {
            os << "                  <xs:minInclusive value=\"" << it->min_float << "\"/>" << endl;
          }
          if (it->max_float!=std::numeric_limits<DoubleReal>::max())
          {
            os << "                  <xs:maxInclusive value=\"" << it->max_float << "\"/>" << endl;
          }
        }
        else if (it->value.valueType()==DataValue::INT_VALUE)
        {
          os << "                <xs:restriction base=\"xs:integer\">" << endl;
          if (it->min_int!=-std::numeric_limits<Int>::max())
          {
            os << "                  <xs:minInclusive value=\"" << it->min_int << "\"/>" << endl;
          }
          if (it->max_int!=std::numeric_limits<Int>::max())
          {
            os << "                  <xs:maxInclusive value=\"" << it->max_int << "\"/>" << endl;
          }
        }
        os << "                </xs:restriction>" << endl;
        os << "              </xs:simpleType>" << endl;
      }
      os << "            </xs:element>" << endl;
    }
    os << "          </xs:sequence>" << endl;
    os << "        </xs:complexType>" << endl;
    os << "      </xs:element>" << endl;
    os << "    </xs:schema>" << endl;
    os << "  </wsdl:types>" << endl;
    //message
    os << "  <wsdl:message name=\"" << tool_name_ << "RequestMessage\">" << endl;
    os << "    <wsdl:part element=\"ns1:" << tool_name_ << "Request\" name=\"part1\"/>" << endl;
    os << "  </wsdl:message>" << endl;
    //port
    os << "  <wsdl:portType name=\"SVMHCProcessPortType\">" << endl;
    os << "    <wsdl:operation name=\"request\">" << endl;
    os << "      <wsdl:input message=\"tns:" << tool_name_ << "RequestMessage\"/>" << endl;
    os << "    </wsdl:operation>" << endl;
    os << "  </wsdl:portType>" << endl;
    //binding
    os << "  <wsdl:binding name=\"" << tool_name_ << "ProviderServiceBinding\" type=\"tns:" << tool_name_ << "PortType\">" << endl;
    os << "    <soap:binding style=\"rpc\" transport=\"http://schemas.xmlsoap.org/soap/http\" xmlns:soap=\"http://schemas.xmlsoap.org/wsdl/soap/\"/>" << endl;
    os << "    <wsdl:operation name=\"request\">" << endl;
    os << "      <soap:operation soapAction=\"\" style=\"rpc\" xmlns:soap=\"http://schemas.xmlsoap.org/wsdl/soap/\"/>" << endl;
    os << "      <wsdl:input>" << endl;
    os << "        <soap:body encodingStyle=\"http://schemas.xmlsoap.org/soap/encoding/\" use=\"encoded\" xmlns:soap=\"http://schemas.xmlsoap.org/wsdl/soap/\"/>" << endl;
    os << "      </wsdl:input>" << endl;
    os << "      <wsdl:output>" << endl;
    os << "        <soap:body encodingStyle=\"http://schemas.xmlsoap.org/soap/encoding/\" use=\"encoded\" xmlns:soap=\"http://schemas.xmlsoap.org/wsdl/soap/\"/>" << endl;
    os << "      </wsdl:output>" << endl;
    os << "    </wsdl:operation>" << endl;
    os << "  </wsdl:binding>" << endl;
    //service
    os << "  <wsdl:service name=\"" << tool_name_ << "ProviderService\">" << endl;
    os << "    <wsdl:port binding=\"tns:" << tool_name_ << "ProviderServiceBinding\" name=\"" << tool_name_ << "ProviderServicePort\">" << endl;
    os << "     <soap:address location=\"http://trypsin.informatik.uni-tuebingen.de:30090/active-bpel/services/" << tool_name_ << "ProviderService\" xmlns:soap=\"http://schemas.xmlsoap.org/wsdl/soap/\"/>" << endl;
    os << "    </wsdl:port>" << endl;
    os << "  </wsdl:service>" << endl;
    //end
    os << "</wsdl:definitions>" << endl;
    os.close();

    //validate written file
    XMLValidator validator;
    if (!validator.isValid(filename, File::find("SCHEMAS/WSDL_20030211.xsd")))
    {
      writeLog_("Error: The written WSDL file does not validate against the XML schema. Please report this bug!");
      return INTERNAL_ERROR;
    }

    return EXECUTION_OK;
  }

} // namespace OpenMS

