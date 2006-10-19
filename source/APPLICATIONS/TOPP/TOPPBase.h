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
	*/
  class TOPPBase
  {
	 public:
			
		/// Exit codes
		enum ExitCodes
			{
				OK,
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
			
		/// Default construtor
		TOPPBase(const String& tool_name);

		/// Destructor
		~TOPPBase();
		
		/// Main routine of all TOPP applications
		ExitCodes main(int argc , char** argv);
			
	 protected:

		/// Tool name
		String tool_name_;

		/// Debug level
		SignedInt debug_level_;

		/// Settings
		Param param_;

		/// Log file stream
		std::ofstream log_;

		/// Command line options with argument (options)
		std::map<std::string,std::string> options_;

		/// Command line options without argument (flags)
		std::map<std::string,std::string> flags_;

		/// Current Instance number
		SignedInt instance_number_;

		/// Prints the help for the command line options and usage. Do not list the common options.
		virtual void printToolUsage_()=0;
			
		/// Prints the tool-specific command line options and appends the common options.
		void printUsage_();
			
		/**
			 @brief Prints the help for the INI-File options and a sample entry.
			
			 Be carefull the types of the sample entries. Do not list the common options.
		*/
		virtual void printToolHelpOpt_()=0;			

		/// Prints the tool-specific INI options and flags and appends the common options and flags.
		void printHelpOpt_();

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
			
		/// Actual main method.
		virtual ExitCodes main_(int argc , char** argv)=0;
			
		/// Parses the command line.
		void parseCommandLine_(int argc , char** argv);
		
		/** @name Debug and Log output
		 */
		//@{
		/// Writes a string to the log file and to std::out
		void writeLog_(const String& text);
		/// Writes a string to the log file if the debug level is at least @p min_level
		void writeDebug_(const String& text, UnsignedInt min_level);
		/// Writes a Param to the log file if the debug level is at least @p min_level
		void writeDebug_(const String& text,const Param& param, UnsignedInt min_level);
		//@}

		/** 
			@name Parameter lookup (command line and ini file)
		 */
		//@{
		/**
			 @brief Return the value @key of param_ as a string or @p default_value when this value is not set.
				
			 @note See @ref getParam_ for the order of the lookup.
		*/
		String getParamAsString_(const String& key, const String& default_value="");
		/**
			 @brief Return the value @key of param_ as an integer or @p default_value when this value is not set.
				
			 @note See @ref getParam_ for the order of the lookup.
		*/
		SignedInt getParamAsInt_(const String& key, SignedInt default_value=0);
		/**
			 @brief Return the value @key of param_ as a double or @p default_value when this value is not set.
				
			 @note See @ref getParam_ for the order of the lookup.
		*/
		double getParamAsDouble_(const String& key, double default_value=0);
		/**
			 @brief Return the value @key of param_ as a bool.
			 
			 If the DataValue is a string, the values 'off', 'on', 'true' and 'false' are interpreted. 
			 If the DataValue is a numerical value, the values '0' and '1' interpreted.
			 For all other values and when the value of key @p key is not set, the @p default_value is returned.
			 
			 @note See @ref getParam_ for the order of the lookup.
		*/
		bool getParamAsBool_(const String& key, bool default_value=false);
		/**
			 @brief Return a value of param_ as DataValue.
	
			 All methods in this section search for the paramters in the following order:
			 <UL>
			 	<LI>bla (command line)
				<LI>&lt;toolname&gt;:&lt;instance&gt;:bla
				<LI>common:&lt;toolname&gt;:bla
				<LI>common:bla
			 </UL>
			 
			 @note The search key for the commandline command parameter is the internal name! See @ref setOptionsAndFlags_.

		*/
		const DataValue& getParam_(const String& key);
		//@}
		
		/**
			 @brief Returns a new Param object containing all entries that start with @p prefix.
  	
			 @param prefix should contain a ':' at the end if you want to extract a subtree.
			 Otherwise not only nodes, but as well values with that prefix are copied.
			 @param remove_prefix indicates if the prefix is removed before adding entries to the new Param
			 @param new_prefix is added to the front of all keys

			 Inheritance of parameters is supported as follows: If the subtree
			 specified by prefix contains an <code>&lt;ITEM name="inherit"
			 value="other:place" type="string"/&gt;</code>, then everything from
			 <code>other:place</code> is inherited.  This works recursively, but at
			 most 15 steps.  Otherwise an exception is thrown (e.g. to detect
			 cycles).  It is not an error if <code>other:place</code> is not an
			 existing <code>&lt;NODE&gt;</code> - we just won't inherit anything
			 from there in this case.

			 @todo Support for new_prefix is not fully tested and should be considered experimental. (Clemens)

		*/
		Param getParamCopy_(const std::string& prefix, bool remove_prefix=false, const std::string& new_prefix="");
  };

} // namespace OpenMS

