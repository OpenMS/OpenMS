// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/Date.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/VALIDATORS/XMLValidator.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder_impl.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h>
#include <OpenMS/FILTERING/TRANSFORMERS/PreprocessingFunctor.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithm.h>

// OpenMP support
#ifdef _OPENMP
	#include <omp.h>
#endif

#include <cmath>

using namespace std;

namespace OpenMS
{
	using namespace Exception;

  TOPPBase::TOPPBase(const String& tool_name, const String& tool_description, bool official, bool id_tag_support, const String& version)
  	: tool_name_(tool_name),
  		tool_description_(tool_description),
			id_tag_support_(id_tag_support),
			id_tagger_(tool_name),
			instance_number_(-1),
			debug_level_(-1),
			version_(version),
			log_type_(ProgressLogger::NONE),
			test_mode_(false)
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
		if (official && !getToolList().has(tool_name_))
		{
			writeLog_(String("Error: Message to maintainer - If '") + tool_name_ + "' is an official TOPP tool, add it to the TOPPBase tools list. If it is not, set the 'official' bool of the TOPPBase constructor to false.");
		}
	}

	TOPPBase::~TOPPBase()
  {
  	//delete TOPP.log and custom log file if they are empty
  	StringList log_files;
  	log_files << "TOPP.log";
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

		//register values
		registerOptionsAndFlags_();
		addEmptyLine_();
		addText_("Common TOPP options:");
		registerStringOption_("ini","<file>","","Use the given TOPP INI file",false);
		registerStringOption_("log","<file>","TOPP.log","Location of the log file",false, true);
		registerIntOption_("instance","<n>",1,"Instance number for the TOPP INI file",false);
		registerIntOption_("debug","<n>",0,"Sets the debug level",false, true);
		registerIntOption_("threads", "<n>", 1, "Sets the number of threads allowed to be used by the TOPP tool", false);
		registerStringOption_("write_ini","<file>","","Writes the default configuration file",false);
		registerStringOption_("write_wsdl","<file>","","Writes the default WSDL file",false);
		registerFlag_("no_progress","Disables progress logging to command line");
		if (id_tag_support_)
		{
			registerStringOption_("id_pool","<file>",
														"",
														String("ID pool file to DocumentID's for all generated output files. Disabled by default. (Set to 'main' to use ") + String() + id_tagger_.getPoolFile() + ")"
														,false);
		}
		registerFlag_("test","Enables the test mode (needed for software testing only)", true);
		registerFlag_("-help","Shows this help");

		// prepare options and flags for command line parsing
		Map<String,String> options;
		Map<String,String> flags;
		Map<String,String> multi_options;
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
				case ParameterInformation::INPUT_FILE_LIST:
				case ParameterInformation::OUTPUT_FILE_LIST:
				case ParameterInformation::INTLIST:
				case ParameterInformation::DOUBLELIST:
				case ParameterInformation::STRINGLIST:
					multi_options[string("-") + it->name] = it->name;
					break;
				default:
					options[string("-") + it->name] = it->name;
					break;
			}
		}
		param_cmdline_.parseCommandLine(argc,argv,options,flags,multi_options);

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
			String write_ini_file("");
			if (param_cmdline_.exists("write_ini")) write_ini_file = param_cmdline_.getValue("write_ini");
			if (write_ini_file != "")
			{
				Param default_params = getDefaultParameters_();
				// check if augmentation with -ini param is needed
				DataValue in_ini;
				if (param_cmdline_.exists("ini")) in_ini = param_cmdline_.getValue("ini");
				if (!in_ini.isEmpty())
				{
					Param ini_params;
					ini_params.load( (String)in_ini );
					// augment
					for(Param::ParamIterator it = ini_params.begin(); it != ini_params.end();++it)
					{
						if (default_params.exists(it.getName()))
						{
							// param 'type': do not override!
							if (it.getName().hasSuffix(":type"))
							{	
								if (default_params.getValue(it.getName()) != it->value)
								{
									writeLog_("Warning: Augmented and Default Ini-File differ in value. Default value will not be altered!");
									continue;
								}
							}
							
							// all other parameters:
							Param::ParamEntry entry = default_params.getEntry (it.getName());
							if (entry.value.valueType() == it->value.valueType())
							{
								if (entry.value != it->value)
								{
									// check entry for consistency (in case restrictions have changed)							
									entry.value = it->value;
									String s;
									if (entry.isValid(s))
									{
										// overwrite default value
										writeDebug_(String("Overriding Default-Parameter ") + it.getName() + " with new value "+String(it->value) ,10); 
										default_params.setValue(it.getName(),it->value, entry.description, default_params.getTags(it.getName()));
									}
									else
									{
										writeDebug_(String("Parameter ") + it.getName() + " does not fit into new restriction settings! Ignoring...",10); 
									}
								}
								else
								{
									// value stayed the same .. nothing to be done
								}
							}
							else
							{
								writeDebug_(String("Parameter ") + it.getName() + " has changed value type! Ignoring...",10); 
							}
						}
						else
						{
							writeDebug_(String("Deprecated Parameter ") + it.getName() + " given in -ini argument! Ignoring...",10); 
						}
					}
				
				}
				outputFileWritable_(write_ini_file);
				
				default_params.store(write_ini_file);
				return EXECUTION_OK;
			}

			// '-write_wsdl' given
			String wsdl_file("");
			if (param_cmdline_.exists("write_wsdl")) wsdl_file = param_cmdline_.getValue("write_wsdl");
			if (wsdl_file != "")
			{
				outputFileWritable_(wsdl_file);
				ofstream os(wsdl_file.c_str());

				//write header
				os << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
				os << "<wsdl:definitions targetNamespace=\"http://www-bs.informatik.uni-tuebingen.de/compas\" xmlns:ns1=\"http://org.apache.axis2/xsd\" xmlns:plnk=\"http://schemas.xmlsoap.org/ws/2003/05/partner-link/\" xmlns:soap=\"http://schemas.xmlsoap.org/wsdl/soap/\" xmlns:tns=\"http://www-bs.informatik.uni-tuebingen.de/compas\" xmlns:wsdl=\"http://schemas.xmlsoap.org/wsdl/\" xmlns:xs=\"http://www.w3.org/2001/XMLSchema\">" << endl;
				os << "  <wsdl:types>" << endl;
				os << "    <xs:schema attributeFormDefault=\"unqualified\" elementFormDefault=\"qualified\" targetNamespace=\"http://org.apache.axis2/xsd\" xmlns:ns1=\"http://org.apache.axis2/xsd\" xmlns:xs=\"http://www.w3.org/2001/XMLSchema\">" << endl;
				os << "      <xs:element name=\"" << tool_name_ << "Request\">" << endl;
				os << "        <xs:complexType>" << endl;
				os << "          <xs:sequence>" << endl;

				//write types (forward declaration for readablility only. Could be defined in the message as well.
				Param param = getDefaultParameters_();
				param = param.copy(tool_name_ + ":1:",true);
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

				//validate written file
				XMLValidator validator;
				if (!validator.isValid(wsdl_file,File::find("SCHEMAS/WSDL_20030211.xsd")))
				{
					writeLog_("Error: The written WSDL file does not validate against the XML schema. Please report this bug!");
					return INTERNAL_ERROR;
				}

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
			#ifdef _OPENMP
			Int threads = getParamAsInt_("threads", 1);
			omp_set_num_threads(threads);
			#endif

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
			writeLog_(String("Internal error: Request for unregistered parameter '") + e.what() + "'");
			writeDebug_(String("Error occured in line ") + e.getLine() + " of file " + e.getFile() + " (in function: " + e.getFunction() + ") !",1);
			return INTERNAL_ERROR;
		}
		catch(Exception::WrongParameterType& e)
		{
			writeLog_(String("Internal error: Request for parameter with wrong type '") + e.what() + "'");
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
	       << "Version: " << version_ << endl << endl
	       << "Usage:" << endl
				 << "  " << tool_name_ << " <options>" << endl
				 << endl
				 << "Options (mandatory options marked with '*'):" << endl;

		//determine max length of parameters (including argument) for indentation
		UInt max_size = 0;
		for( vector<ParameterInformation>::const_iterator it = parameters_.begin(); it != parameters_.end(); ++it)
		{
			max_size = max((UInt)max_size,(UInt)(it->name.size()+it->argument.size()+it->required));
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
				case ParameterInformation::DOUBLE:
				case ParameterInformation::INT:
				case ParameterInformation::STRINGLIST:
				case ParameterInformation::INTLIST:
				case ParameterInformation::DOUBLELIST:
					{
						String tmp = it->default_value.toString();
						if (tmp!="" && tmp!="[]")
						{
							addons.push_back(String("default: \"") + tmp + "\"");
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
						String tmp;
						tmp.concatenate(it->valid_strings.begin(),it->valid_strings.end(),",");

						String add = "";
						if (it->type == ParameterInformation::INPUT_FILE || it->type == ParameterInformation::OUTPUT_FILE ||
								it->type == ParameterInformation::INPUT_FILE_LIST || it->type == ParameterInformation::OUTPUT_FILE_LIST) add = " formats";

						addons.push_back(String("valid") + add + ": \"" + tmp + "\"");
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
				String output;
				output.concatenate(addons.begin(),addons.end()," ");
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
			//determine indentation of description
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


	void TOPPBase::registerStringOption_(const String& name, const String& argument, const String& default_value,const String& description, bool required, bool advanced)
	{
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
		for (Size i=0; i<formats.size(); ++i)
		{
			if (FileHandler::getTypeByFileName(String(".")+formats[i])==FileTypes::UNKNOWN)
			{
				throw InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__,"The file format '" + formats[i] + "' is invalid!");
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
				parameters_[i].max_float = max;
				return;
			}
		}
		//parameter not found
		throw ElementNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,name);
	}


	void TOPPBase::registerInputFile_(const String& name, const String& argument, const String& default_value,const String& description, bool required, bool advanced)
	{
		parameters_.push_back(ParameterInformation(name, ParameterInformation::INPUT_FILE, argument, default_value, description, required, advanced));
	}

	void TOPPBase::registerOutputFile_(const String& name, const String& argument, const String& default_value,const String& description, bool required, bool advanced)
	{
		parameters_.push_back(ParameterInformation(name, ParameterInformation::OUTPUT_FILE, argument, default_value, description, required, advanced));
	}

	void TOPPBase::registerDoubleOption_(const String& name, const String& argument, double default_value, const String& description, bool required, bool advanced)
	{
		parameters_.push_back(ParameterInformation(name, ParameterInformation::DOUBLE, argument, default_value, description, required, advanced));
	}

	void TOPPBase::registerIntOption_(const String& name, const String& argument, Int default_value, const String& description, bool required, bool advanced)
	{
		parameters_.push_back(ParameterInformation(name, ParameterInformation::INT, argument, default_value, description, required, advanced));
	}

	void TOPPBase::registerOutputFileList_( const String& name, const String& argument, StringList default_value, const String& description, bool required, bool advanced )
	{
		parameters_.push_back(ParameterInformation(name,ParameterInformation::OUTPUT_FILE_LIST,argument,default_value,description,required,advanced));
	}

	void TOPPBase::registerInputFileList_( const String& name, const String& argument, StringList default_value, const String& description, bool required, bool advanced )
	{
		parameters_.push_back(ParameterInformation(name,ParameterInformation::INPUT_FILE_LIST,argument,default_value,description,required,advanced));
	}

	void TOPPBase::registerStringList_( const String& name, const String& argument, StringList default_value, const String& description, bool required, bool advanced)
	{
		parameters_.push_back(ParameterInformation(name,ParameterInformation::STRINGLIST,argument,default_value,description,required,advanced));
	}

	void TOPPBase::registerIntList_( const String& name, const String& argument, IntList default_value, const String& description, bool required, bool advanced )
	{
		parameters_.push_back(ParameterInformation(name,ParameterInformation::INTLIST,argument,default_value,description,required,advanced));
	}

	void TOPPBase::registerDoubleList_( const String& name, const String& argument, DoubleList default_value, const String& description, bool required, bool advanced)
	{
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
						valid_strings.concatenate(p.valid_strings.begin(),p.valid_strings.end(),"','");
						throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, String("Invalid value '") + tmp + "' for string parameter \"" + name + "\" given. Valid strings are: \"" + valid_strings + "\".");
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
						valid_formats.concatenate(p.valid_strings.begin(),p.valid_strings.end(),"','");
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
						valid_formats.concatenate(p.valid_strings.begin(),p.valid_strings.end(),"','");
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

	DoubleReal TOPPBase::getDoubleOption_(const String& name) const
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
		double tmp = getParamAsDouble_(name, (DoubleReal)p.default_value);
		writeDebug_(String("Value of string option '") + name + "': " + String(tmp), 1);

		//check if in valid range
		if (p.required || ( setByUser_(name) && tmp!=(DoubleReal)p.default_value))
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
		Int tmp = getParamAsInt_(name, (Int)p.default_value);
		writeDebug_(String("Value of string option '") + name + "': " + String(tmp), 1);

		//check if in valid range
		if (p.required || ( setByUser_(name) && tmp!=(Int)p.default_value))
		{
			if (tmp<p.min_int || tmp>p.max_int)
			{
				throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, String("Invalid value '") + tmp + "' for integer parameter '" + name + "' given. Out of valid range: '" + p.min_int + "'-'" + p.max_int + "'.");
			}
		}

		return tmp;
	}

	StringList TOPPBase::getStringList_(const String& name) const
	{
		const ParameterInformation& p = findEntry_(name);
		if (p.type!=ParameterInformation::STRINGLIST && p.type!=ParameterInformation::INPUT_FILE_LIST && p.type!=ParameterInformation::OUTPUT_FILE_LIST)
		{
			throw Exception::WrongParameterType(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
		if (p.required && !setByUser_(name) )
		{
			throw Exception::RequiredParameterNotGiven(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
		StringList tmp_list = getParamAsStringList_(name, (StringList)p.default_value);
		String tmp;
		for(StringList::iterator it = tmp_list.begin(); it < tmp_list.end(); ++it)
		{
			tmp = *it;
			writeDebug_(String("Value of string option '") + name + "': " + tmp, 1);

			// if required or set by user, do some validity checks
			if (p.required || ( setByUser_(name) && tmp_list!=p.default_value))
			{
				//check if files are readable/writeable
				if (p.type==ParameterInformation::INPUT_FILE_LIST)
				{
					writeDebug_( "Checking input file '" + name + "': '" + tmp + "'", 2 );
					inputFileReadable_(tmp);
				}
				else if (p.type==ParameterInformation::OUTPUT_FILE_LIST)
				{
					writeDebug_( "Checking output file '" + name + "': '" + tmp + "'", 2 );
					outputFileWritable_(tmp);
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
							throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, String("Invalid value '") + tmp + "' for string parameter '" + name + "' given. Valid strings are: '" + valid_strings + "'.");
						}
					}
					else if (p.type==ParameterInformation::INPUT_FILE_LIST)
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
							valid_formats.concatenate(p.valid_strings.begin(),p.valid_strings.end(),"','");
							throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, String("Input file '" + tmp + "' has invalid format '") + format + "'. Valid formats are: '" + valid_formats + "'.");
						}
					}
					else if (p.type==ParameterInformation::OUTPUT_FILE_LIST)
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
							valid_formats.concatenate(p.valid_strings.begin(),p.valid_strings.end(),"','");
							throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, String("Invalid output file extension '") + tmp + "'. Valid file extensions are: '" + valid_formats + "'.");
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
			throw Exception::WrongParameterType(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
		if (p.required && !setByUser_(name) )
		{
			throw Exception::RequiredParameterNotGiven(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
		DoubleList tmp_list = getParamAsDoubleList_(name, (DoubleList)p.default_value);
		DoubleReal tmp;
		for(DoubleList::iterator it = tmp_list.begin(); it < tmp_list.end(); ++it)
		{
			tmp = *it;
			writeDebug_(String("Value of string option '") + name + "': " + String(tmp), 1);

			//check if in valid range
			if (p.required || ( setByUser_(name) && tmp_list!=(DoubleList)p.default_value))
			{
				if (tmp<p.min_float || tmp>p.max_float)
				{
					throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, String("Invalid value '") + tmp + "' for float parameter '" + name + "' given. Out of valid range: '" + p.min_float + "'-'" + p.max_float + "'.");
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
			throw Exception::WrongParameterType(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
		if (p.required && !setByUser_(name) )
		{
			throw Exception::RequiredParameterNotGiven(__FILE__,__LINE__,__PRETTY_FUNCTION__, name);
		}
		IntList tmp_list = getParamAsIntList_(name, (IntList)p.default_value);

		Int tmp;
		for(IntList::iterator it = tmp_list.begin(); it < tmp_list.end();++it)
		{
			tmp = *it;
			writeDebug_(String("Value of string option '") + name + "': " + String(tmp), 1);

			//check if in valid range
			if (p.required || ( setByUser_(name) && tmp_list!=(IntList)p.default_value))
			{
				if (tmp<p.min_int || tmp>p.max_int)
				{
					throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, String("Invalid value '") + tmp + "' for integer parameter '" + name + "' given. Out of valid range: '" + p.min_int + "'-'" + p.max_int + "'.");
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
			//if the parameter has the correct type (ini file) no conversion is necessary
			if (tmp.valueType()==DataValue::INT_VALUE)
			{
				return (Int)tmp;
			}
			//for the command line a conversion is necessary
			try
			{
				return tmp.toString().toInt();
			}
			catch(Exception::ConversionError&)
			{
				throw InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__,String("Invalid value '") + tmp.toString() + "' for integer parameter '" + key + "' given.");
			}
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
			//if the parameter has the correct type (ini file) no conversion is necessary
			if (tmp.valueType()==DataValue::DOUBLE_VALUE)
			{
				return (DoubleReal)tmp;
			}
			//for the command line a conversion is necessary
			try
			{
				return tmp.toString().toDouble();
			}
			catch(Exception::ConversionError&)
			{
				throw InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__,String("Invalid value '") + tmp.toString() + "' for double parameter '" + key + "' given.");
			}
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
			//if the parameter has the correct type (ini file) no conversion is necessary
			if (tmp.valueType()==DataValue::INT_LIST)
			{
				return (IntList)tmp;
			}
			//for the command line a conversion is necessary
			StringList sl = (StringList)tmp;
			IntList il;
			il.resize(sl.size());
		 	for (Size i = 0; i < sl.size(); ++i)
			{
				try
				{
					il[i] = sl[i].toInt();
				}
				catch(Exception::ConversionError&)
				{
					throw InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__,String("Invalid value '") + sl[i] + "' in integer list parameter '" + key + "' given.");
				}
			}

			return il;
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
			//if the parameter has the correct type (ini file) no conversion is necessary
			if (tmp.valueType()==DataValue::DOUBLE_LIST)
			{
				return (DoubleList)tmp;
			}
			//for the command line a conversion is necessary
			StringList sl = (StringList)tmp;
			DoubleList dl;
			dl.resize(sl.size());
		 	for (Size i = 0; i < sl.size(); ++i)
			{
				try
				{
					dl[i] = sl[i].toDouble();
				}
				catch(Exception::ConversionError&)
				{
					throw InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__,String("Invalid value '") + sl[i] + "' in double list parameter '" + key + "' given.");
				}
			}

			return dl;
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
							writeLog_("Warning: Wrong parameter type of '" + location + it.getName() + "' in '" + filename + "'. Type should be 'string'!");
						}
						break;
					case ParameterInformation::INTLIST:
						if (it->value.valueType()!=DataValue::INT_LIST)
						{
							writeLog_("Warning: Wrong parameter type of '" + location + it.getName() + "' in '" + filename + "'. Type should be 'string'!");
						}
						break;
					case ParameterInformation::DOUBLELIST:
						if (it->value.valueType()!=DataValue::DOUBLE_LIST)
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
    if (!File::isDirectory(filename) && File::empty(filename))
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

	Param TOPPBase::getDefaultParameters_() const
	{
		Param tmp;
		String loc = tool_name_ + ":1:";
		//parameters
		for( vector<ParameterInformation>::const_iterator it = parameters_.begin(); it != parameters_.end(); ++it)
		{
			if (it->name!="ini" && it->name!="-help" && it->name!="instance" && it->name!="write_ini" && it->name!="write_wsdl")
			{
				String name = loc + it->name;
				StringList tags;
				if (it->advanced) tags.push_back("advanced");
				if(it->type == ParameterInformation::INPUT_FILE || it->type == ParameterInformation::INPUT_FILE_LIST) tags.push_back("input file");
				if(it->type == ParameterInformation::OUTPUT_FILE || it->type == ParameterInformation::OUTPUT_FILE_LIST) tags.push_back("output file");
				switch(it->type)
				{
					case ParameterInformation::STRING:
						tmp.setValue(name,(String)it->default_value, it->description, tags);
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
								formats.concatenate(it->valid_strings.begin(),it->valid_strings.end(),",");
								formats = String("(valid formats: '") + formats + "')";
							}
							tmp.setValue(name,(String)it->default_value, it->description + formats, tags);
						}
						break;
					case ParameterInformation::DOUBLE:
						tmp.setValue(name,(DoubleReal)it->default_value, it->description, tags);
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
						{
							String formats;
							if (it->valid_strings.size()!=0)
							{
								formats.concatenate(it->valid_strings.begin(),it->valid_strings.end(),",");
								formats = String("(valid formats: '") + formats + "')";
							}
							tmp.setValue(name,(StringList)it->default_value, it->description + formats, tags);
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
		
		//set tool version
		tmp.setValue(tool_name_ + ":version", VersionInfo::getVersion(), "Version of the tool that generated this parameters file.", StringList::create("advanced"));
		
		//descriptions
		tmp.setSectionDescription(tool_name_, tool_description_);
		tmp.setSectionDescription(tool_name_ + ":1", String("Instance '1' section for '") + tool_name_ + "'");

		// store "type" in INI-File (if given)
		if (param_cmdline_.exists("type")) tmp.setValue(loc + "type", (String) param_cmdline_.getValue("type"));

		return tmp;
	}


	const IDTagger& TOPPBase::getIDTagger_() const
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

	Map<String,StringList> TOPPBase::getToolList()
	{
		Map<String,StringList> tools_map;

		tools_map["AdditiveSeries"] = StringList::create("");
		tools_map["BaselineFilter"] = StringList::create("");
		tools_map["ConsensusID"] = StringList::create("");
		tools_map["DBExporter"] = StringList::create("");
		tools_map["DBImporter"] = StringList::create("");
		tools_map["DTAExtractor"] = StringList::create("");
		tools_map["Decharger"] = StringList::create("");
		tools_map["FalseDiscoveryRate"] = StringList::create("");
		tools_map["FeatureFinder"] = Factory<FeatureFinderAlgorithm<Peak1D,Feature> >::registeredProducts();
		tools_map["FeatureLinker"] = Factory<FeatureGroupingAlgorithm>::registeredProducts();
		tools_map["FileConverter"] = StringList::create("");
		tools_map["FileFilter"] = StringList::create("");
		tools_map["FileInfo"] = StringList::create("");
		tools_map["FileMerger"] = StringList::create("");
		tools_map["IDDecoyProbability"] = StringList::create("");
		tools_map["IDFileConverter"] = StringList::create("");
		tools_map["IDFilter"] = StringList::create("");
		tools_map["IDMapper"] = StringList::create("");
		tools_map["IDMerger"] = StringList::create("");
		tools_map["IDRTCalibration"] = StringList::create("");
		tools_map["ITRAQAnalyzer"] = StringList::create("4plex,8plex");
		tools_map["InspectAdapter"] = StringList::create("");
		tools_map["InternalCalibration"] = StringList::create("");
		tools_map["MapAligner"] = Factory<MapAlignmentAlgorithm>::registeredProducts();
		tools_map["MapNormalizer"] = StringList::create("");
		tools_map["MascotAdapter"] = StringList::create("");
		tools_map["MascotAdapterOnline"] = StringList::create("");
		tools_map["NoiseFilter"] = StringList::create("sgolay,gaussian");
		tools_map["OMSSAAdapter"] = StringList::create("");
		tools_map["PILISModel"] = StringList::create("");
		tools_map["PILISIdentification"] = StringList::create("");
		tools_map["PTModel"] = StringList::create("");
		tools_map["PTPredict"] = StringList::create("");
		tools_map["PeakPicker"] = StringList::create("wavelet,high_res");
		tools_map["PepNovoAdapter"] = StringList::create("");
		tools_map["RTModel"] = StringList::create("");
		tools_map["RTPredict"] = StringList::create("");
		tools_map["Resampler"] = StringList::create("");
		tools_map["SILACAnalyzer"] = StringList::create("");
		tools_map["SequestAdapter"] = StringList::create("");
		tools_map["SpectraFilter"] = Factory<PreprocessingFunctor>::registeredProducts();
		tools_map["SpecLibSearcher"] = StringList::create("");
		tools_map["TOFCalibration"] = StringList::create("");
		tools_map["TextExporter"] = StringList::create("");
		tools_map["TextImporter"] = StringList::create("");
		tools_map["XTandemAdapter"] = StringList::create("");
		tools_map["PrecursorIonSelector"] = StringList::create("");
		tools_map["PeptideIndexer"] = StringList::create("");
		tools_map["CompNovo"] = StringList::create("CompNovo,CompNovoCID");
    tools_map["SILACAnalyzer"] = StringList::create("double,triple");

		return tools_map;
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

} // namespace OpenMS

