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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include "TOPPBase.h"

using namespace std;

namespace OpenMS
{
  TOPPBase::TOPPBase(const String& tool_name)
  	: tool_name_(tool_name),
		debug_level_(-1),
		param_(),
		log_(),
		options_(), 
		flags_(),
		instance_number_(-1)
	{
		
	}

	TOPPBase::~TOPPBase()
  {
  	
  }

	TOPPBase::ExitCodes TOPPBase::main(int argc , char** argv)
	{
		//----------------------------------------------------------
		//parse command line
		//----------------------------------------------------------
		setOptionsAndFlags_();
		options_["-ini"] = "ini";
		options_["-log"] = "log";
		options_["-n"] = "instance";
		options_["-d"] = "debug";

		flags_["--help"] = "help";
		flags_["--help-opt"] = "helpopt";		
		
		parseCommandLine_(argc,argv);
		
		//start logging to default location
		log_.open(getParamAsString_("log","TOPP.log").c_str(), std::ofstream::out | std::ofstream::app);
		log_ << "-----------------------------------------------------------" << std::endl;
		//set debug level
		debug_level_ = getParamAsInt_("debug",0);
		writeDebug_(String("Debug level: ")+String(debug_level_),1);
		
		//set instance number
		instance_number_ = getParamAsInt_("instance",1);
		writeDebug_(String("Instance: ")+String(instance_number_),1);

		// test if no options were given
		if (argc==1)
		{
			writeLog_("No options given. Aborting!");
			printUsage_();
			return ILLEGAL_PARAMETERS;
		}

		// '--help' given
		if (!(param_.getValue("help").isEmpty()))
		{
			printUsage_();
			return OK;
		}
	
		// '--help-opt' given
		if (!(param_.getValue("helpopt").isEmpty()))
		{
			printHelpOpt_();
			return OK;
		}
			
		// test if unknown options were given
		if (!param_.getValue("unknown").isEmpty())
		{

			writeLog_(String("Unknown option(s) '") + getParamAsString_("unknown") + "' given. Aborting!");
			printUsage_();
			return ILLEGAL_PARAMETERS;
		}
		
		// test if unknown text argument were given (we do not use them)
		if (!param_.getValue("misc").isEmpty())
		{
			writeLog_(String("Trailing text argument(s) '") + getParamAsString_("misc") + "' given. Aborting!");
			printUsage_();
			return ILLEGAL_PARAMETERS;
		}
		
		ExitCodes result;
		try
		{
			//-------------------------------------------------------------
			// loading INI file
			//-------------------------------------------------------------
			if (!param_.getValue("ini").isEmpty())
			{
				writeDebug_(String("INI file: ") + getParamAsString_("ini"),1);
				param_.load((String)(param_.getValue("ini")));
			}

			
			//-------------------------------------------------------------
			// determine and open the real log file
			//-------------------------------------------------------------
			if (!getParam_("log").isEmpty())
			{
				writeDebug_(String("Log file: ")+getParamAsString_("log"),1);
				log_.close();
				log_.open(getParamAsString_("log").c_str(), std::ofstream::out | std::ofstream::app);
			}
			

			
			//----------------------------------------------------------
			//main
			//----------------------------------------------------------
			
			result = main_(argc, argv);
		}

		//----------------------------------------------------------
		//error handling
		//----------------------------------------------------------

		catch(Exception::UnableToCreateFile& e)
		{
			writeLog_(String("Error: Unable to write file (") + e.what() + ")");
			return CANNOT_WRITE_OUTPUT_FILE;
		}	
		catch(Exception::FileNotFound& e)
		{
			writeLog_(String("Error: File not found (") + e.what() + ")");
			return INPUT_FILE_NOT_FOUND;
		}
		catch(Exception::FileNotReadable& e)
		{
			writeLog_(String("Error: File not readable (") + e.what() + ")");
			return INPUT_FILE_NOT_READABLE;
		}
		catch(Exception::FileEmpty& e)
		{
			writeLog_(String("Error: File empty (") + e.what() + ")");
			return INPUT_FILE_EMPTY;
		}
		catch(Exception::ParseError& e)
		{
			writeLog_(String("Error: Unable to read file (") + e.what() + ")");
			return INPUT_FILE_CORRUPT;
		}
		catch(Exception::Base& e)
		{
			writeLog_(String("Error: Unexpected error (") + e.what() + ")");
			return UNKNOWN_ERROR;
		}
	  
		log_.close();
		
		return result;
	}

	void TOPPBase::printUsage_()
	{
		printToolUsage_();
		
		std::cerr << std::endl
							<< "Common TOPP options are:" << std::endl
							<< "  -ini <file>       Use the given TOPP INI file" << std::endl
							<< "  -log <file>       log file (default: TOPP.log)" << std::endl
							<< "  -n <int>          instance number (default: 1)" << std::endl
							<< "  -d <level>        set debug level (default: 0)" << std::endl
							<< "  --help            show this help" << std::endl
							<< "  --help-opt        show help on the INI options accepted" << std::endl
							<< std::endl ;	
	}


	void TOPPBase::printHelpOpt_()
	{
		printToolHelpOpt_();
		
		std::cerr << std::endl
							<< "Common TOPP INI options are:" << std::endl
							<< "  log       log file (default: TOPP.log)" << std::endl
							<< std::endl ;	
	}

	void TOPPBase::parseCommandLine_(int argc , char** argv)
	{
		param_.parseCommandLine(argc,argv,options_,flags_,"misc","unknown");
	}
	
	void TOPPBase::writeLog_(const String& text)
	{
		std::cout << text << std::endl;
		log_ << Date::now() << " " << tool_name_ << ":" << instance_number_ << ": " << text<< std::endl;
	}
	
	void TOPPBase::writeDebug_(const String& text, UnsignedInt min_level)
	{
		if (debug_level_>=(SignedInt)min_level)
		{
			log_ << Date::now() << " " << tool_name_ << ":" << instance_number_ << ": " << text<< std::endl;
		}
	}
	
	String TOPPBase::getParamAsString_(const String& key, const String& default_value)
	{
		const DataValue& tmp = getParam_(key);
		if (!tmp.isEmpty())
		{
			return (String)(tmp);
		}
		else
		{
			return default_value;
		}
	}

	SignedInt TOPPBase::getParamAsInt_(const String& key, SignedInt default_value)
	{
		const DataValue& tmp = getParam_(key);
		if (!tmp.isEmpty())
		{
			return (SignedInt)(tmp);
		}
		else
		{
			return default_value;
		}
	}

	double TOPPBase::getParamAsDouble_(const String& key, double default_value)
	{
		const DataValue& tmp = getParam_(key);
		if (!tmp.isEmpty())
		{
			return (double)(tmp);
		}
		else
		{
			return default_value;
		}
	}

	DataValue TOPPBase::getParam_(const String& key)
	{
		//command line
		String key_string = key;
		if (!param_.getValue(key_string).isEmpty())
		{
			writeDebug_(String("Parameter '")+key+String("' from COMMAND LINE: ")+(String)(param_.getValue(key_string)),3);
			return param_.getValue(key_string);
		}
		//instance section
		key_string = String(tool_name_) + ":" + String(instance_number_) + ":" + key;
		if (!param_.getValue(key_string).isEmpty())
		{
			writeDebug_(String("Parameter '")+key+String("' from INSTANCE SECTION: ")+(String)(param_.getValue(key_string)),3);
			return param_.getValue(key_string);
		}
		//common secion
		key_string = String("common:")+String(tool_name_)+":"+key;
		if (!param_.getValue(key_string).isEmpty())
		{
			writeDebug_(String("Parameter '")+key+String("' from COMMON SECTION: ")+(String)(param_.getValue(key_string)),3);
			return param_.getValue(key_string);
		}
		writeDebug_(String("Parameter '")+key+String("' NOT FOUND!"),3);
		return DataValue::EMPTY;
	}
	
	Param TOPPBase::getParamCopy_(const std::string& prefix, bool remove_prefix, const std::string& new_prefix)
	{
		return param_.copy(prefix,remove_prefix,new_prefix);
	}


} // namespace OpenMS
