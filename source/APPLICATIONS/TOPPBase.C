// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/Date.h>
#include <OpenMS/CONCEPT/VersionInfo.h>

#include <math.h>

using namespace std;

namespace OpenMS
{

	namespace
	{
		char const log_separator_[] = "================================================================================";
	}

  TOPPBase::TOPPBase(const String& tool_name, const String& tool_description)
  	: tool_name_(tool_name),
  		tool_description_(tool_description),
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
		
		//register values
		registerOptionsAndFlags_();
		addEmptyLine_();
		addText_("Common TOPP options:");
		registerStringOption_("ini","<file>","","Use the given TOPP INI file",false);
		registerStringOption_("log","<file>","TOPP.log","Location of the log file",false);
		registerIntOption_("instance","<n>",1,"Instance number for the TOPP INI file",false);
		registerIntOption_("debug","<n>",0,"Sets the debug level",false);
		registerStringOption_("write_ini","<file>","","Writes an example INI file",false);
		registerFlag_("-help","Shows this help");
		
		// prepare options and flags for command line parsing
		map<string,string> options;
		map<string,string> flags;
		for( vector<ParameterInformation>::const_iterator it = parameters_.begin(); it != parameters_.end(); ++it)
		{
			switch (it->type)
			{		
				case ParameterInformation::TEXT:
				case ParameterInformation::NEWLINE:
					break;
				case ParameterInformation::FLAG:
					flags[string("-") + it->name] = it->name;
					break;
				case ParameterInformation::NONE:
					writeLog_("Warning: undefined type of parameter '" + it->name + "'");
					break;
				default:
					options[string("-") + it->name] = it->name;
					break;
			}
		}
		param_cmdline_.parseCommandLine(argc,argv,options,flags,"misc","unknown");
		
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
		if (!(param_cmdline_.getValue("-help").isEmpty()))
		{
			printUsage_();
			return EXECUTION_OK;
		}
		
		// '-write_ini' given
		String ini_file = (String)param_cmdline_.getValue("write_ini");
		if (ini_file != "")
		{
			outputFileWritable_(ini_file);
			Param tmp;
			String loc = tool_name_ + ":1:";
			//parameters
			for( vector<ParameterInformation>::const_iterator it = parameters_.begin(); it != parameters_.end(); ++it)
			{
				if (it->name!="ini" && it->name!="-help" && it->name!="instance" && it->name!="write_ini")
				{
					switch(it->type)
					{
						case ParameterInformation::STRING:
							tmp.setValue(loc + it->name,it->default_value);
							break;
						case ParameterInformation::DOUBLE:
							tmp.setValue(loc + it->name,String(it->default_value).toDouble());
							break;
						case ParameterInformation::INT:
							tmp.setValue(loc + it->name,String(it->default_value).toInt());
							break;
						case ParameterInformation::FLAG:
							tmp.setValue(loc + it->name,"off");
							break;
						default:
							break;
					}
				}
			}
			//subsections
			for(vector<String>::const_iterator it = subsections_.begin(); it!=subsections_.end(); ++it)
			{
				tmp.setValue(loc + *it + ":dummy" , "Replace this ITEM with your seettings!");
			}
			tmp.store(ini_file);
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
#ifndef DEBUG_TOPP
		try
		{
#endif
			//-------------------------------------------------------------
			// load INI file
			//-------------------------------------------------------------
			{
				DataValue const & value_ini = param_cmdline_.getValue("ini");
				if (!value_ini.isEmpty())
				{
					writeDebug_( "INI file: " + (String)value_ini, 1 );
					writeDebug_( "INI location: " + getIniLocation_(), 1);
					param_inifile_.load( (String)value_ini );
					param_instance_ = param_inifile_.copy( getIniLocation_(), true, "" );
					writeDebug_("Parameters from instance section:",param_instance_,2);
					param_instance_inherited_ = param_inifile_.copyWithInherit( getIniLocation_(), "" );
					writeDebug_("Parameters from instance section, including inherited ones:",param_instance_inherited_,2);
					param_common_tool_ = param_inifile_.copy( "common:"+tool_name_+':', true, "" );
					writeDebug_("Parameters from common section with tool name:",param_common_tool_,2);
					param_common_ = param_inifile_.copy( "common:", true, "" );
					writeDebug_("Parameters from common section without tool name:",param_common_,2);
				}
				param_ = param_cmdline_;
				writeDebug_("Applying defaults to instance section, including inherited ones:",param_common_,2);
				param_.setDefaults( param_instance_inherited_ );
				writeDebug_("Applying defaults to common section with tool name:",param_common_,2);
				param_.setDefaults( param_common_tool_ );
				writeDebug_("Applying defaults to common section without tool name:",param_common_,2);
				param_.setDefaults( param_common_ );

				// check if all parameters are registered and have the correct type
				checkParam_(param_instance_inherited_, (String)value_ini, getIniLocation_());
				checkParam_(param_common_tool_, (String)value_ini, "common:" + tool_name_ + "::");
				checkParam_(param_common_, (String)value_ini, "common::" );
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

#ifndef DEBUG_TOPP
		}

		//----------------------------------------------------------
		//error handling
		//----------------------------------------------------------
		
		// Errors caused by the user
		catch(Exception::UnableToCreateFile& e)
		{
			writeLog_(String("Error: Unable to write file (") + e.what() + ")");
			writeDebug_(String("Error occured in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ")!",1);
			return CANNOT_WRITE_OUTPUT_FILE;
		}	
		catch(Exception::FileNotFound& e)
		{
			writeLog_(String("Error: File not found (") + e.what() + ")");
			writeDebug_(String("Error occured in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !",1);
			return INPUT_FILE_NOT_FOUND;
		}
		catch(Exception::FileNotReadable& e)
		{
			writeLog_(String("Error: File not readable (") + e.what() + ")");
			writeDebug_(String("Error occured in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !",1);
			return INPUT_FILE_NOT_READABLE;
		}
		catch(Exception::FileEmpty& e)
		{
			writeLog_(String("Error: File empty (") + e.what() + ")");
			writeDebug_(String("Error occured in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !",1);
			return INPUT_FILE_EMPTY;
		}
		catch(Exception::ParseError& e)
		{
			writeLog_(String("Error: Unable to read file (") + e.what() + ")");
			writeDebug_(String("Error occured in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !",1);
			return INPUT_FILE_CORRUPT;
		}
		catch(Exception::RequiredParameterNotGiven& e)
		{
			writeLog_(String("Error: The required parameter '") + e.what() + "' was not given!");
			writeDebug_(String("Error occured in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !",1);
			return MISSING_PARAMETERS;
		}
		// Internal errors because of wrong use of this class
		catch(Exception::UnregisteredParameter& e)
		{
			writeLog_(String("Internal error: Request for unregistered parameter '") + e.what() + ")");
			writeDebug_(String("Error occured in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !",1);
			return INTERNAL_ERROR;
		}
		catch(Exception::WrongParameterType& e)
		{
			writeLog_(String("Internal error: Request for parameter with wrong type '") + e.what() + ")");
			writeDebug_(String("Error occured in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !",1);
			return INTERNAL_ERROR;
		}
		// All other errors
		catch(Exception::Base& e)
		{
			writeLog_(String("Error: Unexpected internal error (") + e.what() + ")");
			writeDebug_(String("Error occured in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !",1);
			return UNKNOWN_ERROR;
		}
#endif
	  
		log_.close();
		
		return result;
	}

	void TOPPBase::printUsage_() const
	{
		//common output
		cerr << endl
	       << tool_name_ << " -- " << tool_description_ << endl
	       << "Version: " << VersionInfo::getVersion() << endl
	       << endl
	       << "Usage:" << endl
				 << "  " << tool_name_ << " <options>" << endl
				 << endl
				 << "Options (mandatory options marked with '*'):" << endl;
		
		//determine max length of parameters (including argument) for indentation
		string::size_type max_size = 0;
		for( vector<ParameterInformation>::const_iterator it = parameters_.begin(); it != parameters_.end(); ++it)
		{
			max_size = max(max_size,it->name.size()+it->argument.size()+it->required);
		}
		
		//offset of the descriptions
		UnsignedInt offset = 6 + max_size;
		
		for( vector<ParameterInformation>::const_iterator it = parameters_.begin(); it != parameters_.end(); ++it)
		{
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
			
			//handle newlines in description
			vector<String> parts;
			if (!desc_tmp.split('\n',parts))
			{
				cerr << tmp << desc_tmp;
			}
			else
			{
				vector<String>::iterator it2 = parts.begin();
				cerr << tmp << *it2 << endl;
				it2++;
				for (; (it2+1)!=parts.end(); ++it2)
				{
					it2->firstToUpper();
					if (it->type != ParameterInformation::TEXT) cerr << String(offset,' ');
					cerr << *it2 << endl;
				}
				if (it->type != ParameterInformation::TEXT) cerr << String(offset,' ');
				cerr << *it2;
			}
			

			//DEFAULT
			switch (it->type)
			{
				case ParameterInformation:: STRING:
					if (it->default_value!="") 
					{
						cerr << " (default: '" << it->default_value << "')";
					}
					break;
				case ParameterInformation::DOUBLE:
					cerr << " (default: '" << it->default_value << "')";
					break;
				case ParameterInformation::INT:
					cerr << " (default: '" << it->default_value << "')";
					break;				
				default:
					break;
			}
			
			cerr << endl;

		} 
	}


	void TOPPBase::registerStringOption_(const String& name, const String& argument, const String& default_value,const String& description, bool required)
	{
		parameters_.push_back(ParameterInformation(name, ParameterInformation::STRING, argument, default_value, description, required));
	}
	

	void TOPPBase::registerDoubleOption_(const String& name, const String& argument, double default_value, const String& description, bool required)
	{
		parameters_.push_back(ParameterInformation(name, ParameterInformation::DOUBLE, argument, String(default_value), description, required));
	}
	

	void TOPPBase::registerIntOption_(const String& name, const String& argument, SignedInt default_value, const String& description, bool required)
	{
		parameters_.push_back(ParameterInformation(name, ParameterInformation::INT, argument, String(default_value), description, required));
	}

	void TOPPBase::registerFlag_(const String& name, const String& description)
	{
		parameters_.push_back(ParameterInformation(name, ParameterInformation::FLAG, "", "", description, false));
	}

	void TOPPBase::addEmptyLine_()
	{
		parameters_.push_back(ParameterInformation("",ParameterInformation::NEWLINE, "", "", "", false));
	}

	void TOPPBase::addText_(const String& text)
	{
		parameters_.push_back(ParameterInformation("",ParameterInformation::TEXT, "", "", text, false));
	}

	const TOPPBase::ParameterInformation& TOPPBase::findEntry_(const String& name) const throw (Exception::UnregisteredParameter)
	{
		vector<ParameterInformation>::const_iterator it = parameters_.begin();
		while(it != parameters_.end() && it->name!=name)
		{ 
			++it;
		}
		if (it == parameters_.end())
		{
			throw Exception::UnregisteredParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
		return *it;
	}
	
	String TOPPBase::getStringOption_(const String& name) const throw (Exception::UnregisteredParameter, Exception::RequiredParameterNotGiven, Exception::WrongParameterType )
	{
		const ParameterInformation& p = findEntry_(name);
		if (p.type != ParameterInformation::STRING)
		{
			throw Exception::WrongParameterType(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
		String tmp = getParamAsString_(name, p.default_value);
		writeDebug_(String("Value of string option '") + name + "': " + tmp, 1);
		
		if (p.required && (tmp==p.default_value) )
		{
			throw Exception::RequiredParameterNotGiven(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
		
		return tmp;
	}
	
	double TOPPBase::getDoubleOption_(const String& name) const throw (Exception::UnregisteredParameter, Exception::RequiredParameterNotGiven, Exception::WrongParameterType )
	{
		const ParameterInformation& p = findEntry_(name);
		if (p.type != ParameterInformation::DOUBLE)
		{
			throw Exception::WrongParameterType(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
		double tmp = getParamAsDouble_(name, String(p.default_value).toFloat());
		writeDebug_(String("Value of string option '") + name + "': " + String(tmp), 1);

		if (p.required && fabs(tmp-String(p.default_value).toDouble()< 0.0001) )
		{
			throw Exception::RequiredParameterNotGiven(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
		
		return tmp;
	}
	
	SignedInt TOPPBase::getIntOption_(const String& name) const throw (Exception::UnregisteredParameter, Exception::RequiredParameterNotGiven, Exception::WrongParameterType )
	{
		const ParameterInformation& p = findEntry_(name);
		if (p.type != ParameterInformation::INT)
		{
			throw Exception::WrongParameterType(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
		SignedInt tmp = getParamAsInt_(name, String(p.default_value).toInt());
		writeDebug_(String("Value of string option '") + name + "': " + String(tmp), 1);

		if (p.required && (tmp==String(p.default_value).toInt()) )
		{
			throw Exception::RequiredParameterNotGiven(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}

		return tmp;
	}
	
	bool TOPPBase::getFlag_(const String& name) const throw (Exception::UnregisteredParameter, Exception::WrongParameterType )
	{
		const ParameterInformation& p = findEntry_(name);
		if (p.type != ParameterInformation::FLAG)
		{
			throw Exception::WrongParameterType(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
		bool tmp = getParamAsBool_(name);
		writeDebug_(String("Value of string option '") + name + "': " + String(tmp), 1);
		return tmp;
	}

	void TOPPBase::writeLog_(const String& text) const
	{
		cout << text << endl;
		enableLogging_();
		log_ << Date::now() << ' ' << getIniLocation_() << ": " << text<< endl;
	}
	
	void TOPPBase::writeDebug_(const String& text, UnsignedInt min_level) const
	{
		if (debug_level_>=(SignedInt)min_level)
		{
			cout << text << endl;
			enableLogging_();
			log_ << Date::now() << ' ' << getIniLocation_() << ": " << text<< endl;
		}
	}

	void TOPPBase::writeDebug_(const String& text, const Param& param, UnsignedInt min_level) const
	{
		if (debug_level_>=(SignedInt)min_level)
		{
			cout 	<< " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl
						<< Date::now() << ' ' << getIniLocation_() << " " << text<< endl
						<< param 
						<< " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl;
			enableLogging_();
			log_  << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl
						<< Date::now() << ' ' << getIniLocation_() << " " << text<< endl
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

	void TOPPBase::enableLogging_() const
	{
		if ( !log_.is_open() )
		{
			DataValue const & log_destination = param_cmdline_.getValue("log");
			if ( log_destination.isEmpty() )
			{
				log_.open("TOPP.log", ofstream::out | ofstream::app);
				log_ << log_separator_ << endl;
				if (debug_level_>=1)
				{
					cout << "Writing to 'TOPP.log'" << endl;
					log_ << Date::now() << ' ' << getIniLocation_() << ": " << "Writing to 'TOPP.log'"<< endl;
				}
			}
			else
			{
				log_.open( ((String)log_destination) .c_str(), ofstream::out | ofstream::app);
				log_ << log_separator_ << endl;
				if (debug_level_>=1)
				{
					cout << "Writing to '" << (String)log_destination << '\'' << endl;
					log_ << Date::now() << ' ' << getIniLocation_() << ": " << "Writing to '" << (String)log_destination << '\'' <<  endl;
				}
			}
		}
		return;
	}

	void TOPPBase::checkParam_(const Param& param, const String& filename, const String& location) const
	{
		//cout << endl << "--"<< location<< "--" << endl << param << endl << endl;
		for (Param::ConstIterator it = param.begin(); it!=param.end(); ++it)
		{
			// subsections
			if (String(it->first).has(':'))
			{
				String sec = String(it->first).prefix(':');
				if (find(subsections_.begin(), subsections_.end(), sec)==subsections_.end())
				{
					if (!(location == "common::" && sec==tool_name_) )
					{
						writeLog_("Warning: Unknown subsection '" + sec + "' in '" + filename + "' (location '"+location+"')!");
					}
				}
				continue;
			}
			// if no such parameter is registered an exception is thrown
			try
			{
				//check type
				switch (findEntry_(it->first).type)
				{
					case ParameterInformation::STRING:
						if (it->second.valueType()!=DataValue::STRVALUE) 
						{
							writeLog_("Warning: Wrong parameter type of '" + location + it->first + "' in '" + filename + "'. Type should be 'string'!");
						}
						break;
					case ParameterInformation::DOUBLE:
						if (it->second.valueType()!=DataValue::DOUVALUE && it->second.valueType()!=DataValue::FLOVALUE)
						{
							writeLog_("Warning: Wrong  parameter type of '" + location + it->first + "' in '" + filename + "'. Type should be 'double'!");
						}
						break;						
					case ParameterInformation::INT:
						if (it->second.valueType()!=DataValue::INTVALUE)
						{
							writeLog_("Warning: Wrong parameter type of '" + location + it->first + "' in '" + filename + "'. Type should be 'int'!");
						}
						break;
					case ParameterInformation::FLAG:
						switch (it->second.valueType())
						{
							case DataValue::STRVALUE:
								{
									String tmp = it->second;
									if (tmp!="on" && tmp!="off" && tmp!="true" && tmp!="false") 
									{
										writeLog_("Warning: Unrecognized value for '" + location + it->first + "' in '" + filename + "'. It should be 'on' or 'off'!");
									}
								}
								break;
							case DataValue::SHOVALUE:
							case DataValue::LONVALUE:
							case DataValue::INTVALUE:
								{
									SignedInt tmp = it->second;
									if (tmp!=1 && tmp!=0) 
									{
										writeLog_("Warning: Unrecognized value for '" + location + it->first + "' in '" + filename + "'. It should be '0' or '1'!");
									}
								}
								break;
							default:
								writeLog_("Warning: Wrong parameter type of '" + location + it->first + "' in '" + filename + "'. Type should be 'string' or 'int'!");
								break;
						};
						break;
					default:
						break;
				};
			}
			catch (Exception::UnregisteredParameter)
			{
				writeLog_("Warning: Unknown parameter '" + location + it->first + "' in '" + filename + "'!");
			}
		}
	}

	void TOPPBase::inputFileReadable_(const String& filename) const throw (Exception::FileNotFound, Exception::FileNotReadable, Exception::FileEmpty)
	{
		if (!File::exists(filename))
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
		}
		if (!File::readable(filename))
		{
			throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);			
		}
    if (File::empty(filename))
    {
      throw Exception::FileEmpty(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }
	}
	
	void TOPPBase::outputFileWritable_(const String& filename) const throw (Exception::UnableToCreateFile)
	{
		if (!File::writable(filename))
		{
			throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
		}
	}

	void TOPPBase::registerSubsection_(const String& name)
	{
		subsections_.push_back(name);
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

} // namespace OpenMS

