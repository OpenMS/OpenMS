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

#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Date.h>
#include <OpenMS/CONCEPT/VersionInfo.h>

#include <iostream>
#include <fstream>

namespace OpenMS
{
	/**
		 @brief Base class for TOPP Applications.
		
		 You have to implement the virtual methods @ref printToolUsage_, @ref printToolHelpOpt_, @ref setOptionsAndFlags_ 
		 and @ref main_ only.
		 
		 @todo complete the tests (Clemens)
	*/
  class TOPPBase
  {
	 public:
			
		/// Exit codes
		enum ExitCodes
			{
				EXECUTION_OK,
				INPUT_FILE_NOT_FOUND,
				INPUT_FILE_NOT_READABLE,
				INPUT_FILE_CORRUPT,
				INPUT_FILE_EMPTY,
				CANNOT_WRITE_OUTPUT_FILE,
				ILLEGAL_PARAMETERS,
				UNKNOWN_ERROR,
				EXTERNAL_PROGRAM_ERROR,
				PARSE_ERROR
			};
			
		/// Construtor
		TOPPBase(const String& tool_name);

		/// Destructor
		virtual ~TOPPBase();
		
		/// Main routine of all TOPP applications
		ExitCodes main(int argc , char** argv);

	 public:
		
		/**@brief Returns the name of the tool, e.g. "TOPPTool".

			 This is assigned once and for all in the constructor.
		*/
		String const& getToolName() const { return tool_name_; }

		/**@brief Returns the location of the ini file where parameters are taken
			 from.  E.g. if the command line was <code>TOPPTool -n 17</code>, then
			 this will be <code>"TOPPTool:17:"</code>.  Note the ':' at the end.

			 This is assigned during tool startup, depending on the command line but
			 (of course) not depending on ini files.
		*/
		String const& getIniLocation() const { return ini_location_; }

	 private:

		/**@brief Tool name.  This is assigned once and for all in the constructor.
		 */
		String const tool_name_;

		/**@brief Instance number
		 */
		SignedInt const instance_number_;

		/**@brief Location in the ini file where to look for parameters.
		*/
		String const ini_location_;
		
		/// No default constructor.  It is "declared away".
		TOPPBase();
		
		/// No default copy constructor.  It is "declared away".
		TOPPBase(const TOPPBase&);
			
		/// Debug level
		SignedInt debug_level_;

		/// All parameters relevant to this invocation of the program.
		Param param_;

		/// All parameters specified in the ini file
		Param param_inifile_;

		/// Parameters from command line
		Param param_cmdline_;

		/**@brief Parameters from instance section

		   @note We might as well skip this (since the instance section with
			 inheritance contains a superset of it) but I leave it so that we can
			 trace where the parameters were obtained from. -- (cg) 2006-10-05
		*/
		Param param_instance_;

		/// Parameters from instance section, including inherited ones
		Param param_instance_inherited_;

		/// Parameters from common section with tool name.
		Param param_common_tool_;

		/// Parameters from common section without tool name.
		Param param_common_;

		/** @brief Log file stream.  Use the writeLog_() and writeDebug_() methods
				to access it.
		*/
		mutable std::ofstream log_;
		
		/** @brief Ensures that at least some default logging destination is
				opened for writing in append mode.

				@note This might be invoked at various places early in the startup
				process of the TOPP tool.  Thus we cannot consider the ini file here.
				The final logging destination is determined in main().
		*/
		void enableLogging_() const;

	 protected:
		
		/// Command line options with argument (options)
		std::map<std::string,std::string> options_;

		/// Command line options without argument (flags)
		std::map<std::string,std::string> flags_;

		/// Prints the help for the command line options and usage.  Do not list the common options here.
		virtual void printToolUsage_() const =0 ; 
			
		/// Prints the tool-specific command line options <b>and</b> appends the common options.
		void printUsage_() const;
			
		/**
			 @brief Prints the help for the INI-File options and a sample entry.
			
			 Be careful the types of the sample entries. Do not list the common options.
		*/
		virtual void printToolHelpOpt_() const =0;			

		/// Prints the tool-specific INI options and flags and appends the common options and flags.
		void printHelpOpt_() const;

		/**
			 @brief Sets the valid commannd line options (with argument) and flags (without argument).
				
			 The following values are automatically set:<BR>		
			 options_["-ini"] = "ini";<BR>
			 options_["-log"] = "log";<BR>
			 options_["-n"] = "instance";<BR>
			 options_["-d"] = "debug";<BR>
			 flags_["--help"] = "help";<BR>
			 flags_["--help-opt"] = "helpopt";<BR>
			 
			 @note Make sure the name command line name and the internal name are the same!
			       Otherwise the lookup from the ini file does not work!
		*/
		virtual void setOptionsAndFlags_()=0;
			
		/// The actual "main" method.  main_() is invoked by main().
		virtual ExitCodes main_(int argc , char** argv)=0;
		
		
		/** @name Debug and Log output
		 */
		//@{			
		/// Writes a string to the log file and to std::cout
		void writeLog_(const String& text) const;
		
		/// Writes a @p text to the log file and to std::cout if the debug level is at least @p min_level
		void writeDebug_(const String& text, UnsignedInt min_level) const;

		/**@brief Writes a String followed by a Param to the log file and to
			 std::cout if the debug level is at least @p min_level
		*/
		void writeDebug_(const String& text, const Param& param, UnsignedInt min_level) const;
		//@}

		/** 
			@name Parameter lookup (command line and ini file)
		 */
		//@{
		/**
			 @brief Return the value of parameter @p key as a string or @p default_value if this value is not set.
				
			 @note See getParam_(const String&) const for the order in which parameters are searched.
		*/
		String getParamAsString_(const String& key, const String& default_value="") const;
		/**
			 @brief Return the value of parameter @p key as an integer or @p default_value if this value is not set.
				
			 @note See getParam_(const String&) const for the order in which parameters are searched.
		*/
		SignedInt getParamAsInt_(const String& key, SignedInt default_value=0) const;
		/**
			 @brief Return the value of parameter @p key as a double or @p default_value if this value is not set.
				
			 @note See getParam_(const String&) const for the order in which parameters are searched.
		*/
		double getParamAsDouble_(const String& key, double default_value=0) const;
		/**
			 @brief Return the value of parameter @p key as a bool or @p default_value if this value is not set.
			 
			 If the DataValue is a string, the values 'off', 'on', 'true' and 'false' are interpreted. 
			 If the DataValue is a numerical value, the values '0' and '1' interpreted.
			 For all other values and when the value of key @p key is not set, the @p default_value is returned.
			 
			 @note See getParam_(const String&) const for the order in which parameters are searched.
		*/
		bool getParamAsBool_(const String& key, bool default_value=false) const;
		/**
			 @brief Return the value @p key of parameters as DataValue.
			 DataValue::EMPTY indicates that a parameter was not found.
				
			 Parameters are searched in this order:
			 -# command line
			 -# instance section, e.g. "TOPPTool:1:some_key", see getIniLocation().
			 -# inherited from instance section
			 -# common section with tool name,  e.g. "common:ToolName:some_key"
			 -# common section without tool name,  e.g. "common:some_key"
			 .
			 where "some_key" == key in the examples.
			 
			 @note The search key for the commandline command parameter is the internal name! See @ref setOptionsAndFlags_.
		*/
		DataValue const& getParam_(const String& key) const;

		/**@brief Return <em>all</em> parameters relevant to this TOPP tool.
			 Returns a Param that contains everything you can get by the
			 getParamAs...() methods.
		 */
		Param const& getParam_() const;

		/**
			 @brief Returns a new Param object containing all entries that start with @p prefix.
  	
			 @p prefix should contain a ':' at the end if you want to extract a
			 subtree.  In this case, "inherit" tags are supported.
			 This agrees with the convention that getIniLocation() ends with a ':'.

			 Inheritance of parameters is supported as follows: If the subtree
			 specified by prefix contains an <code>&lt;ITEM name="inherit"
			 value="other:place" type="string"/&gt;</code>, then everything from
			 <code>other:place</code> is inherited.  This works recursively, but at
			 most 15 steps.  Otherwise an exception is thrown (e.g. to detect
			 cycles).  (BTW, it is not an error if <code>other:place</code> is not an
			 existing <code>&lt;NODE&gt;</code>)

			 Otherwise not only nodes, but as well values with that prefix are
			 copied.  (I am yet to see a useful application for this...)

		*/
		Param getParamCopy_( const std::string& prefix ) const;

  };

} // namespace OpenMS

