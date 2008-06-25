// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/FORMAT/FileHandler.h>

#include <math.h>

using namespace std;

namespace OpenMS
{
	using namespace Exception;
	
  TOPPBase::TOPPBase(const String& tool_name, const String& tool_description)
  	: tool_name_(tool_name),
  		tool_description_(tool_description),
			instance_number_(-1),
			debug_level_(-1),
			log_type_(ProgressLogger::NONE)
	{
	}

	TOPPBase::~TOPPBase()
  {
  	//delete TOPP.log and custom log file if they are empty
  	StringList log_files;
  	log_files << "TOPP.log";
  	if (!getParam_("log").isEmpty()) log_files << (String)(getParam_("log"));
		for(UInt i=0; i< log_files.size(); ++i)
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

		//register values
		registerOptionsAndFlags_();
		addEmptyLine_();
		addText_("Common TOPP options:");
		registerStringOption_("ini","<file>","","Use the given TOPP INI file",false);
		registerStringOption_("log","<file>","TOPP.log","Location of the log file",false);
		registerIntOption_("instance","<n>",1,"Instance number for the TOPP INI file",false);
		registerIntOption_("debug","<n>",0,"Sets the debug level",false);
		registerStringOption_("write_ini","<file>","","Writes an example INI file",false);
		registerFlag_("no_progress","Disables progress logging to command line");
		registerFlag_("-help","Shows this help");

		// prepare options and flags for command line parsing
		map<String,String> options;
		map<String,String> flags;
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
		if (param_cmdline_.exists("-help"))
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
			String ini_file("");
			if (param_cmdline_.exists("write_ini")) ini_file = param_cmdline_.getValue("write_ini");
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
						String name = loc + it->name;
						switch(it->type)
						{
							case ParameterInformation::STRING:
								tmp.setValue(name,it->default_value, it->description);
								if (it->valid_strings.size()!=0)
								{
									tmp.setValidStrings(name,it->valid_strings);
								}
								break;
							case ParameterInformation::INPUT_FILE:
							case ParameterInformation::OUTPUT_FILE:
								{
									String formats;
									if (it->valid_strings.size()!=0)
									{
										formats.implode(it->valid_strings.begin(),it->valid_strings.end(),",");
										formats = String("(valid formats: '") + formats + "')";
									}
									tmp.setValue(name,it->default_value, it->description + formats);
								}
								break;
							case ParameterInformation::DOUBLE:
								tmp.setValue(name,String(it->default_value).toDouble(), it->description);
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
								tmp.setValue(name,String(it->default_value).toInt(), it->description);
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
								tmp.setValue(name,"false", it->description);
								tmp.setValidStrings(name,StringList::create("true,false"));
								break;
							default:
								break;
						}
					}
				}
				//subsections
				for(map<String,String>::const_iterator it = subsections_.begin(); it!=subsections_.end(); ++it)
				{
					Param tmp2 = getSubsectionDefaults_(it->first);
					if (!tmp2.empty())
					{
						tmp.insert(loc + it->first + ":",tmp2);
						tmp.setSectionDescription(loc + it->first, it->second);
					}
				}
				tmp.setSectionDescription(tool_name_ + ":1", String("Instance '1' section for '") + tool_name_ + "'");
				
				// store "type" in INI-File (if given)
				if (param_cmdline_.exists("type")) tmp.setValue(loc + "type", (String) param_cmdline_.getValue("type"));
							
				tmp.store(ini_file);
				return EXECUTION_OK;
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
					param_inifile_.load( (String)value_ini );
					param_instance_ = param_inifile_.copy( getIniLocation_(), true);
					writeDebug_("Parameters from instance section:",param_instance_,2);
					param_common_tool_ = param_inifile_.copy( "common:"+tool_name_+':', true );
					writeDebug_("Parameters from common section with tool name:",param_common_tool_,2);
					param_common_ = param_inifile_.copy( "common:", true );
					writeDebug_("Parameters from common section without tool name:",param_common_,2);
				}
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
		catch(Exception::InvalidParameter& e)
		{
			writeLog_(String("Invalid parameter: ") + e.what());
			writeDebug_(String("Error occured in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !",1);
			return ILLEGAL_PARAMETERS;
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
		catch(Exception::BaseException& e)
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
	       << "Version: " << VersionInfo::getVersionAndTime() << endl
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
		UInt offset = 6 + max_size;

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

			//DEFAULT
			StringList addons;
			switch (it->type)
			{
				case ParameterInformation::STRING:
					if (it->default_value!="")
					{
						addons.push_back(String("default: '") + it->default_value + "'");
					}
					break;
				case ParameterInformation::DOUBLE:
					addons.push_back(String("default: '") + it->default_value + "'");
					break;
				case ParameterInformation::INT:
					addons.push_back(String("default: '") + it->default_value + "'");
					break;
				default:
					break;
			}
			
			//RESTRICTIONS
			if (it->type == ParameterInformation::STRING || it->type == ParameterInformation::INPUT_FILE || it->type == ParameterInformation::OUTPUT_FILE)
			{
				if (it->valid_strings.size()!=0)
				{
					String tmp;
					tmp.implode(it->valid_strings.begin(),it->valid_strings.end(),",");

					String add = "";
					if (it->type == ParameterInformation::INPUT_FILE || it->type == ParameterInformation::OUTPUT_FILE) add = " formats";

					addons.push_back(String("valid") + add + ": '" + tmp + "'");
				}
			}
			else if (it->type == ParameterInformation::INT)
			{
				if (it->min_int!=-std::numeric_limits<Int>::max())
				{
					addons.push_back(String("min: '") + it->min_int + "'");
				}
				if (it->max_int!=std::numeric_limits<Int>::max())
				{
					addons.push_back(String("max: '") + it->max_int + "'");
				}
			}
			else if (it->type == ParameterInformation::DOUBLE)
			{
				if (it->min_float!=-std::numeric_limits<DoubleReal>::max())
				{
					addons.push_back(String("min: '") + it->min_float + "'");
				}
				if (it->max_float!=std::numeric_limits<DoubleReal>::max())
				{
					addons.push_back(String("max: '") + it->max_float + "'");
				}
			}
			//add DEFAULT and RESTRICTIONS
			if (addons.size()!=0)
			{
				String output;
				output.implode(addons.begin(),addons.end()," ");
				if (desc_tmp[desc_tmp.size()-1]!='\n') desc_tmp += " ";
				desc_tmp += String("(") + output + ")";
			}

			//handle newlines in description
			vector<String> parts;
			if (!desc_tmp.split('\n',parts))
			{
				cerr << tmp << desc_tmp;
			}
			else
			{
				vector<String>::iterator it2 = parts.begin();
				it2->firstToUpper();
				cerr << tmp << *it2 << endl;
				it2++;
				for (; (it2+1)!=parts.end(); ++it2)
				{
					if (it->type != ParameterInformation::TEXT) cerr << String(offset,' ');
					cerr << *it2 << endl;
				}
				if (it->type != ParameterInformation::TEXT)
				{
					// Note: one space less here, if default will appear at beginning of line
					cerr << String(offset-it2->empty(),' ');
				}
				cerr << *it2;
			}

			cerr << endl;
		}

		if (subsections_.size()!=0)
		{
			//determin indentation of description
			UInt indent = 0;
			for(map<String,String>::const_iterator it = subsections_.begin(); it!=subsections_.end(); ++it)
			{
				indent = max((UInt)it->first.size(),indent);
			}
			indent += 6;

			//output
			cerr << endl
					 << "The following configuration subsections are valid:" << endl;
			for(map<String,String>::const_iterator it = subsections_.begin(); it!=subsections_.end(); ++it)
			{
				String tmp = String(" - ") + it->first;
				tmp.fillRight(' ',indent);
				cerr << tmp << it->second << endl;
			}
			cerr << endl
					 << "You can write an example INI file using the '-write_ini' option." << endl
					 << "Documentation of subsection parameters can be found in the" << endl
					 << "doxygen documentation or the INIFileEditor." << endl
					 << "Have a look at OpenMS/doc/index.html for more information." << endl;
		}
		cerr << endl;
	}


	void TOPPBase::registerStringOption_(const String& name, const String& argument, const String& default_value,const String& description, bool required)
	{
		parameters_.push_back(ParameterInformation(name, ParameterInformation::STRING, argument, default_value, description, required));
	}

	void TOPPBase::setValidStrings_(const String& name, const std::vector<String>& strings)
	{
		//check for commas
		for (UInt i=0; i<strings.size(); ++i)
		{
			if (strings[i].has(','))
			{
				throw InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__,"Comma characters in Param string restrictions are not allowed!");
			}
		}
		//search the right parameter
		for (UInt i=0; i<parameters_.size(); ++i)
		{
			if (parameters_[i].name==name)
			{
				//check if the type matches
				if (parameters_[i].type!=ParameterInformation::STRING)
				{
					throw ElementNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,name);
				}
				parameters_[i].valid_strings = strings;
				return;
			}
		}
		//parameter not found
		throw ElementNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,name);		
	}

	void TOPPBase::setValidFormats_(const String& name, const std::vector<String>& formats)
	{
		//check for commas
		for (UInt i=0; i<formats.size(); ++i)
		{
			if (FileHandler::getTypeByFileName(String(".")+formats[i])==FileHandler::UNKNOWN)
			{
				throw InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__,"The file format '" + formats[i] + "' is invalid!");
			}
		}
		//search the right parameter
		for (UInt i=0; i<parameters_.size(); ++i)
		{
			if (parameters_[i].name==name)
			{
				//check if the type matches
				if (parameters_[i].type!=ParameterInformation::INPUT_FILE && parameters_[i].type!=ParameterInformation::OUTPUT_FILE)
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
		for (UInt i=0; i<parameters_.size(); ++i)
		{
			if (parameters_[i].name==name)
			{
				//check if the type matches
				if (parameters_[i].type!=ParameterInformation::INT)
				{
					throw ElementNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,name);
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
		for (UInt i=0; i<parameters_.size(); ++i)
		{
			if (parameters_[i].name==name)
			{
				//check if the type matches
				if (parameters_[i].type!=ParameterInformation::INT)
				{
					throw ElementNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,name);
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
		for (UInt i=0; i<parameters_.size(); ++i)
		{
			if (parameters_[i].name==name)
			{
				//check if the type matches
				if (parameters_[i].type!=ParameterInformation::DOUBLE)
				{
					throw ElementNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,name);
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
		for (UInt i=0; i<parameters_.size(); ++i)
		{
			if (parameters_[i].name==name)
			{
				//check if the type matches
				if (parameters_[i].type!=ParameterInformation::DOUBLE)
				{
					throw ElementNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,name);
				}
				parameters_[i].max_float = max;
				return;
			}
		}
		//parameter not found
		throw ElementNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,name);	
	}
	

	void TOPPBase::registerInputFile_(const String& name, const String& argument, const String& default_value,const String& description, bool required)
	{
		parameters_.push_back(ParameterInformation(name, ParameterInformation::INPUT_FILE, argument, default_value, description, required));
	}
	
	void TOPPBase::registerOutputFile_(const String& name, const String& argument, const String& default_value,const String& description, bool required)
	{
		parameters_.push_back(ParameterInformation(name, ParameterInformation::OUTPUT_FILE, argument, default_value, description, required));
	}	
	
	void TOPPBase::registerDoubleOption_(const String& name, const String& argument, double default_value, const String& description, bool required)
	{
		parameters_.push_back(ParameterInformation(name, ParameterInformation::DOUBLE, argument, String(default_value), description, required));
	}

	void TOPPBase::registerIntOption_(const String& name, const String& argument, Int default_value, const String& description, bool required)
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

	const TOPPBase::ParameterInformation& TOPPBase::findEntry_(const String& name) const
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

	String TOPPBase::getStringOption_(const String& name) const
	{
		const ParameterInformation& p = findEntry_(name);
		if (p.type!=ParameterInformation::STRING && p.type!=ParameterInformation::INPUT_FILE && p.type!=ParameterInformation::OUTPUT_FILE)
		{
			throw Exception::WrongParameterType(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
		if (p.required && !setByUser_(name) )
		{
			throw Exception::RequiredParameterNotGiven(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
		String tmp = getParamAsString_(name, p.default_value);
		writeDebug_(String("Value of string option '") + name + "': " + tmp, 1);
		
		// if required or set by user, do some validity checks
		if (p.required || ( setByUser_(name) && tmp!=p.default_value))
		{
			//check if files are readable/writeable
			if (p.type==ParameterInformation::INPUT_FILE)
			{
				writeDebug_( "Checking input file '" + name + "': '" + tmp + "'", 2 );
				inputFileReadable_(tmp);
			}
			else if (p.type==ParameterInformation::OUTPUT_FILE)
			{
				writeDebug_( "Checking output file '" + name + "': '" + tmp + "'", 2 );
				outputFileWritable_(tmp);
			}
			
			//check restrictions
			if (p.valid_strings.size()!=0)
			{
				if (p.type==ParameterInformation::STRING)
				{
					if (find(p.valid_strings.begin(),p.valid_strings.end(),tmp)==p.valid_strings.end())
					{
						String valid_strings = "";
						valid_strings.implode(p.valid_strings.begin(),p.valid_strings.end(),"','");
						throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, String("Invalid value '") + tmp + "' for string parameter '" + name + "' given. Valid strings are: '" + valid_strings + "'.");
					}
				}
				else if (p.type==ParameterInformation::INPUT_FILE)
				{
					writeDebug_( "Checking input file '" + name + "': '" + tmp + "'", 2 );
					inputFileReadable_(tmp);
					
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
						valid_formats.implode(p.valid_strings.begin(),p.valid_strings.end(),"','");
						throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, String("Input file '" + tmp + "' has invalid format '") + format + "'. Valid formats are: '" + valid_formats + "'.");					
					}
				}
				else if (p.type==ParameterInformation::OUTPUT_FILE)
				{
					writeDebug_( "Checking output file '" + name + "': '" + tmp + "'", 2 );
					outputFileWritable_(tmp);
	
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
						valid_formats.implode(p.valid_strings.begin(),p.valid_strings.end(),"','");
						throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, String("Invalid output file extension '") + tmp + "'. Valid file extensions are: '" + valid_formats + "'.");					
					}
				}
			}
		}
		
		return tmp;
	}

	bool TOPPBase::setByUser_(const String& name) const
	{
		//look up because of possible exception only
		findEntry_(name);

		if (param_cmdline_.exists(name))
		{
			return true;
		}

		if (param_instance_.exists(name))
		{
			return true;
		}

		if (param_common_tool_.exists(name))
		{
			return true;
		}

		if (param_common_.exists(name))
		{
			return true;
		}

		return false;
	}

	double TOPPBase::getDoubleOption_(const String& name) const
	{
		const ParameterInformation& p = findEntry_(name);
		if (p.type != ParameterInformation::DOUBLE)
		{
			throw Exception::WrongParameterType(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
		if (p.required && !setByUser_(name) )
		{
			throw Exception::RequiredParameterNotGiven(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
		double tmp = getParamAsDouble_(name, String(p.default_value).toDouble());
		writeDebug_(String("Value of string option '") + name + "': " + String(tmp), 1);

		//check if in valid range
		if (p.required || ( setByUser_(name) && tmp!=p.default_value.toDouble()))
		{
			if (tmp<p.min_float || tmp>p.max_float)
			{
				throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, String("Invalid value '") + tmp + "' for float parameter '" + name + "' given. Out of valid range: '" + p.min_float + "'-'" + p.max_float + "'.");
			}
		}

		return tmp;
	}

	Int TOPPBase::getIntOption_(const String& name) const
	{
		const ParameterInformation& p = findEntry_(name);
		if (p.type != ParameterInformation::INT)
		{
			throw Exception::WrongParameterType(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
		if (p.required && !setByUser_(name) )
		{
			throw Exception::RequiredParameterNotGiven(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
		Int tmp = getParamAsInt_(name, String(p.default_value).toInt());
		writeDebug_(String("Value of string option '") + name + "': " + String(tmp), 1);

		//check if in valid range
		if (p.required || ( setByUser_(name) && tmp!=p.default_value.toInt()))
		{
			if (tmp<p.min_int || tmp>p.max_int)
			{
				throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, String("Invalid value '") + tmp + "' for integer parameter '" + name + "' given. Out of valid range: '" + p.min_int + "'-'" + p.max_int + "'.");
			}
		}

		return tmp;
	}

	bool TOPPBase::getFlag_(const String& name) const
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
		log_ << QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss").toStdString() << ' ' << getIniLocation_() << ": " << text<< endl;
	}

	void TOPPBase::writeDebug_(const String& text, UInt min_level) const
	{
		if (debug_level_>=(Int)min_level)
		{
			cout << text << endl;
			enableLogging_();
			log_ << QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss").toStdString() << ' ' << getIniLocation_() << ": " << text<< endl;
		}
	}

	void TOPPBase::writeDebug_(const String& text, const Param& param, UInt min_level) const
	{
		if (debug_level_>=(Int)min_level)
		{
			cout 	<< " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl
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
			return tmp.toString().toInt();
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
			return tmp.toString().toDouble();
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

		// look up in common secion with tool name
		{
			if (param_common_tool_.exists(key))
			{
				writeDebug_(String("Parameter '")+key+String("' from COMMON SECTION (TOOL SPECIFIC): ")+String(param_common_tool_.getValue(key)),3);
				return param_common_tool_.getValue( key );
			}
		}

		// look up in common secion without tool name
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

	void TOPPBase::enableLogging_() const
	{
		if ( !log_.is_open() )
		{
			String log_destination = "";
			if(param_cmdline_.exists("log")) log_destination = param_cmdline_.getValue("log");
			if ( log_destination!="" )
			{
				log_.open("TOPP.log", ofstream::out | ofstream::app);
				if (debug_level_>=1)
				{
					cout << "Writing to 'TOPP.log'" << endl;
					log_ << QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss").toStdString() << ' ' << getIniLocation_() << ": " << "Writing to 'TOPP.log'"<< endl;
				}
			}
			else
			{
				log_.open( log_destination.c_str(), ofstream::out | ofstream::app);
				if (debug_level_>=1)
				{
					cout << "Writing to '" << log_destination << '\'' << endl;
					log_ << QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss").toStdString() << ' ' << getIniLocation_() << ": " << "Writing to '" << log_destination << '\'' <<  endl;
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
			// subsections
			if (it.getName().has(':'))
			{
				String sec = it.getName().prefix(':');
				if (subsections_.find(sec)==subsections_.end())
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
				switch (findEntry_(it.getName()).type)
				{
					case ParameterInformation::STRING:
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
					case ParameterInformation::FLAG:
						if (it->value.valueType()!=DataValue::STRING_VALUE)
						{
							writeLog_("Warning: Wrong parameter type of '" + location + it.getName() + "' in '" + filename + "'. Type should be 'string'!");
						}
						break;
					default:
						break;
				};
			}
			catch (Exception::UnregisteredParameter)
			{
				writeLog_("Warning: Unknown parameter '" + location + it.getName() + "' in '" + filename + "'!");
			}
		}
	}

	void TOPPBase::inputFileReadable_(const String& filename) const
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

	void TOPPBase::outputFileWritable_(const String& filename) const
	{
		if (!File::writable(filename))
		{
			throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
		}
	}

	void TOPPBase::registerSubsection_(const String& name, const String& description)
	{
		subsections_[name] = description;
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
		throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);

		return Param();
	}


} // namespace OpenMS

