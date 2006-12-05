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
// $Maintainer: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/Date.h>

using namespace std;

namespace OpenMS
{

	namespace
	{
		char const log_separator_[] = "================================================================================";
	}

  TOPPBase::TOPPBase(const String& tool_name)
  	: tool_name_(tool_name),
			instance_number_(-1),
			debug_level_(-1)
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
		
		param_cmdline_.parseCommandLine(argc,argv,options_,flags_,"misc","unknown");

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
		if (!(param_cmdline_.getValue("help").isEmpty()))
		{
			printUsage_();
			return EXECUTION_OK;
		}
	
		// '--help-opt' given
		if (!(param_cmdline_.getValue("helpopt").isEmpty()))
		{
			printHelpOpt_();
			return EXECUTION_OK;
		}
			
		// test if unknown options were given
		if (!param_cmdline_.getValue("unknown").isEmpty())
		{

			writeLog_(String("Unknown option(s) '") + getParamAsString_("unknown") + "' given. Aborting!");
			printUsage_();
			return ILLEGAL_PARAMETERS;
		}
		
		// test if unknown text argument were given (we do not use them)
		if (!param_cmdline_.getValue("misc").isEmpty())
		{
			writeLog_(String("Trailing text argument(s) '") + getParamAsString_("misc") + "' given. Aborting!");
			printUsage_();
			return ILLEGAL_PARAMETERS;
		}
		
		ExitCodes result;
		try
		{
			//-------------------------------------------------------------
			// load INI file
			//-------------------------------------------------------------
			{
				DataValue const & value_ini = param_cmdline_.getValue("ini");
				if (!value_ini.isEmpty())
				{
					writeDebug_( "INI file: " + (String)value_ini, 1 );
					param_inifile_.load( (String)value_ini );
					param_instance_ = param_inifile_.copy( getIniLocation(), true, "" );
					writeDebug_("Parameters from instance section:",param_instance_,2);
					param_instance_inherited_ = param_inifile_.copyWithInherit( getIniLocation(), true, "" );
					writeDebug_("Parameters from instance section, including inherited ones:",param_instance_inherited_,2);
					param_common_tool_ = param_inifile_.copy( "common:"+tool_name_+':', true, "" );
					writeDebug_("Parameters from common section with tool name:",param_common_tool_,2);
					param_common_ = param_inifile_.copy( "common:", true, "" );
					writeDebug_("Parameters from common section without tool name:",param_common_,2);
				}
				param_ = param_cmdline_;
				param_.setDefaults( param_instance_inherited_ );
				param_.setDefaults( param_common_tool_ );
				param_.setDefaults( param_common_ );
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
					log_ << log_separator_ << endl;
					writeDebug_("Writing to '"+(String)value_log+'\'',1);
				}
			}

			//-------------------------------------------------------------
			// debug level
			//-------------------------------------------------------------
			debug_level_ = getParamAsInt_("debug",0);
			writeDebug_(String("Debug level (after ini file): ")+String(debug_level_),1);	

			
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

	void TOPPBase::printUsage_() const
	{
		printToolUsage_();
		
		cerr << endl
							<< "Common TOPP options are:" << endl
							<< "  -ini <file>       Use the given TOPP INI file" << endl
							<< "  -log <file>       log file (default: TOPP.log)" << endl
							<< "  -n <int>          instance number (default: 1)" << endl
							<< "  -d <level>        set debug level (default: 0)" << endl
							<< "  --help            show this help" << endl
							<< "  --help-opt        show help on the INI options accepted" << endl
							<< endl ;	
	}


	void TOPPBase::printHelpOpt_() const
	{
		printToolHelpOpt_();
		
		cerr << endl
							<< "Common TOPP INI options are:" << endl
							<< "  log       log file (default: TOPP.log)" << endl
							<< endl ;	
	}

	void TOPPBase::writeLog_(const String& text) const
	{
		cout << text << endl;
		enableLogging_();
		log_ << Date::now() << ' ' << getIniLocation() << ": " << text<< endl;
	}
	
	void TOPPBase::writeDebug_(const String& text, UnsignedInt min_level) const
	{
		if (debug_level_>=(SignedInt)min_level)
		{
			cout << text << endl;
			enableLogging_();
			log_ << Date::now() << ' ' << getIniLocation() << ": " << text<< endl;
		}
	}

	void TOPPBase::writeDebug_(const String& text, const Param& param, UnsignedInt min_level) const
	{
		if (debug_level_>=(SignedInt)min_level)
		{
			cout << text << endl << param;
			enableLogging_();
			log_ 
				<< " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl
				<< Date::now() << ' ' << getIniLocation() << ": " << text<< endl
				<< param 
				<< " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl;
		}
	}
	
	String TOPPBase::getParamAsString_(const String& key, const String& default_value) const
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

	SignedInt TOPPBase::getParamAsInt_(const String& key, SignedInt default_value) const
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

	double TOPPBase::getParamAsDouble_(const String& key, double default_value) const
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

	bool TOPPBase::getParamAsBool_(const String& key, bool default_value) const
	{
		const DataValue& tmp = getParam_(key);
		switch (tmp.valueType())
		{
			case DataValue::STRVALUE:
				{
					String tmp2 = (string)(tmp);
					tmp2.toLower();
					if (tmp2=="off" || tmp2=="false")
					{
						return false;
					}
					else if (tmp2=="on" || tmp2=="true")
					{
						return true;
					}
				}
				break;
			case DataValue::INTVALUE:
			case DataValue::SHOVALUE:
			case DataValue::LONVALUE:
				{
					SignedInt tmp2 = (SignedInt)(tmp);
					if (tmp2==0)
					{
						return false;
					}
					else if (tmp2==1)
					{
						return true;
					}
				}
				break;
			case DataValue::FLOVALUE:
			case DataValue::DOUVALUE:
				{
					float tmp2 = (float)(tmp);
					if (tmp2==0.0)
					{
						return false;
					}
					else if (tmp2==1.0)
					{
						return true;
					}
				}
				break;
			case DataValue::EMPTYVALUE:
			  break;
		}
		return default_value; 
	}

	DataValue const& TOPPBase::getParam_(const String& key) const
	{
		// look up in command line
		{
			DataValue const & value = param_cmdline_.getValue( key );
			if (!value.isEmpty())
			{
				writeDebug_(String("Parameter '")+key+String("' from COMMAND LINE: ")+(String)(value),3);
				return value;
			}
		}

		// look up in instance section
		{
			DataValue const & value = param_instance_.getValue( key );
			if (!value.isEmpty())
			{
				writeDebug_(String("Parameter '")+key+String("' from INSTANCE SECTION: ")+(String)(value),3);
				return value;
			}
		}

		// look up in instance section, with inheritance
		{
			DataValue const & value = param_instance_inherited_.getValue( key );
			if (!value.isEmpty())
			{
				writeDebug_(String("Parameter '")+key+String("' from INSTANCE SECTION (INHERITED): ")+(String)(value),3);
				return value;
			}
		}

		// look up in common secion with tool name
		{
			DataValue const & value = param_common_tool_.getValue( key );
			if (!value.isEmpty())
			{
				writeDebug_(String("Parameter '")+key+String("' from COMMON SECTION (TOOL SPECIFIC): ")+(String)(value),3);
				return value;
			}
		}

		// look up in common secion without tool name
		{
			DataValue const & value = param_common_.getValue( key );
			if (!value.isEmpty())
			{
				writeDebug_(String("Parameter '")+key+String("' from COMMON SECTION: ")+(String)(value),3);
				return value;
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

	Param TOPPBase::getParamCopy_(const string& prefix) const
	{
		return param_inifile_.copyWithInherit(prefix,true,"");
	}

	void TOPPBase::enableLogging_() const
	{
		if ( !log_ )
		{
			DataValue const & log_destination = param_cmdline_.getValue("log");
			if ( log_destination.isEmpty() )
			{
				log_.open("TOPP.log", ofstream::out | ofstream::app);
				log_ << log_separator_ << endl;
				writeDebug_("Writing to 'TOPP.log'",1);
			}
			else
			{
				log_.open( ((String)log_destination) .c_str(), ofstream::out | ofstream::app);
				log_ << log_separator_ << endl;
				writeDebug_("Writing to '"+(String)log_destination+'\'',1);
			}
		}
		return;
	}

} // namespace OpenMS
