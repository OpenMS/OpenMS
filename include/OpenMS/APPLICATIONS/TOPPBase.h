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

#ifndef OPENMS_APPLICATIONS_TOPPBASE_H
#define OPENMS_APPLICATIONS_TOPPBASE_H

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <iostream>
#include <fstream>
class QStringList;

namespace OpenMS
{

  namespace Exception
  {
    /// An unregistered parameter was accessed
    class UnregisteredParameter
          : public Exception::Base
    {
      public:
        UnregisteredParameter( const char* file, int line, const char* function, const String& parameter )
            : Base( file, line, function, "UnregisteredParameter", parameter )
        {
          globalHandler.setMessage( what_ );
        }
    };
    /// A parameter was accessed with the wrong type
    class WrongParameterType
          : public Exception::Base
    {
      public:
        WrongParameterType( const char* file, int line, const char* function, const String& parameter )
            : Base( file, line, function, "WrongParameterType", parameter )
        {
          globalHandler.setMessage( what_ );
        }
    };
    /// A required parameter was not given
    class RequiredParameterNotGiven
          : public Exception::Base
    {
      public:
        RequiredParameterNotGiven( const char* file, int line, const char* function, const String& parameter )
            : Base( file, line, function, "RequiredParameterNotGiven", parameter )
        {
          globalHandler.setMessage( what_ );
        }
    };
  }

  /**
  	 @brief Base class for TOPP applications.

  	This base class implements functionality used in most TOPP tools:
  	- parameter handling
  	- file handling
  	- progress logging

  	If you want to create a new TOPP tool, please take care of the following:
  	- derive a new class from this class
  	- impelment the registerOptionsAndFlags_ and main_ methods
  	- add a doxygen page for the tool and add the page to TOPP.doxygen
  	- hide the derived class in the OpenMS documentation by using doxygen codition macros.

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
        MISSING_PARAMETERS,
        UNKNOWN_ERROR,
        EXTERNAL_PROGRAM_ERROR,
        PARSE_ERROR,
        INCOMPATIBLE_INPUT_DATA,
        INTERNAL_ERROR
    };

      /// Construtor
      TOPPBase( const String& tool_name, const String& tool_description );

      /// Destructor
      virtual ~TOPPBase();

      /// Main routine of all TOPP applications
      ExitCodes main( int argc , char** argv );

      /// Stuct that captures all information of a parameter
      struct ParameterInformation
      {
        /// Parameter types
        enum ParameterTypes
        {
          NONE = 0,         ///< Undefined type
          STRING,         ///< String parameter
          DOUBLE,         ///< Floating point number parameter
          INT,            ///< Integer parameter
          FLAG,           ///< Parameter without argument
          TEXT,           ///< Left aligned text, see addText_
          NEWLINE  ///< An empty line, see addEmptyLine_
        };

        /// name of the parameter (internal and external)
        String name;
        /// type of the parameter
        ParameterTypes type;
        /// default value of the parameter stored as string
        String default_value;
        /// description of the parameter
        String description;
        /// argument in the description
        String argument;
        /// flag that indicates if this parameter is required i.e. it must differ from the default value
        bool required;

        /// Constructor that takes all members in declaration order
        ParameterInformation( const String& n, ParameterTypes t, const String& arg, const String& def, const String& desc, bool req )
        {
          name = n;
          type = t;
          default_value = def;
          description = desc;
          argument = arg;
          required = req;
        }

        ParameterInformation()
            : name(),
            type( NONE ),
            default_value(),
            description(),
            argument(),
            required( true )
        {}

        ParameterInformation& operator=( const ParameterInformation& rhs )
        {
          if ( &rhs == this ) return * this;

          name = rhs.name;
          type = rhs.type;
          default_value = rhs.default_value;
          description = rhs.description;
          argument = rhs.argument;
          required = rhs.required;
          return *this;
        }


      };


    private:

      /// Tool name.  This is assigned once and for all in the constructor.
      String const tool_name_;

      /// Tool description. This is assigned once and for all in the constructor.
      String const tool_description_;

      ///Instance number
      Int const instance_number_;

      ///Location in the ini file where to look for parameters.
      String const ini_location_;

      /// No default constructor.  It is "declared away".
      TOPPBase();

      /// No default copy constructor.  It is "declared away".
      TOPPBase( const TOPPBase& );

      /// Debug level
      Int debug_level_;

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

      ///Log file stream.  Use the writeLog_() and writeDebug_() methods to access it.
      mutable std::ofstream log_;

      /** @brief Ensures that at least some default logging destination is
      		opened for writing in append mode.

      		@note This might be invoked at various places early in the startup
      		process of the TOPP tool.  Thus we cannot consider the ini file here.
      		The final logging destination is determined in main().
      */
      void enableLogging_() const;

      /// Storage location for parameter information
      std::vector<ParameterInformation> parameters_;

      /**
      	@brief This method should return the default parameters for subsections.

      	It is called once for each registered subsection, when writing the an example ini file.

      	Reimplement this method to set the defaults written in the 'write_ini' method.
      */
      virtual Param getSubsectionDefaults_( const String& section ) const;

      /// Storage location and description for allowed subsections
      std::map<String, String> subsections_;

      /**
      	@name Internal parameter handling
       */
      //@{
      /**
      	 @brief Return the value of parameter @p key as a string or @p default_value if this value is not set.

      	 @note See getParam_(const String&) const for the order in which parameters are searched.
      */
      String getParamAsString_( const String& key, const String& default_value = "" ) const;

      /**
      	 @brief Return the value of parameter @p key as an integer or @p default_value if this value is not set.

      	 @note See getParam_(const String&) const for the order in which parameters are searched.
      */
      Int getParamAsInt_( const String& key, Int default_value = 0 ) const;

      /**
      	 @brief Return the value of parameter @p key as a double or @p default_value if this value is not set.

      	 @note See getParam_(const String&) const for the order in which parameters are searched.
      */
      double getParamAsDouble_( const String& key, double default_value = 0 ) const;

      /**
      	 @brief Return the value of parameter @p key as a bool or @p default_value if this value is not set.

      	 If the DataValue is a string, the values 'off', 'on', 'true' and 'false' are interpreted.
      	 If the DataValue is a numerical value, the values '0' and '1' interpreted.
      	 For all other values and when the value of key @p key is not set, the @p default_value is returned.

      	 @note See getParam_(const String&) const for the order in which parameters are searched.
      */
      bool getParamAsBool_( const String& key, bool default_value = false ) const;

      /**
      	 @brief Return the value @p key of parameters as DataValue.
      	 DataValue::EMPTY indicates that a parameter was not found.

      	 Parameters are searched in this order:
      	 -# command line
      	 -# instance section, e.g. "TOPPTool:1:some_key", see getIniLocation_().
      	 -# inherited from instance section
      	 -# common section with tool name,  e.g. "common:ToolName:some_key"
      	 -# common section without tool name,  e.g. "common:some_key"
      	 .
      	 where "some_key" == key in the examples.

      */
      DataValue const& getParam_( const String& key ) const;
      //@}

    protected:

      /**
      	@brief Returns the location of the ini file where parameters are taken
      	from.  E.g. if the command line was <code>TOPPTool -instance 17</code>, then
      	this will be <code>"TOPPTool:17:"</code>.  Note the ':' at the end.

      	This is assigned during tool startup, depending on the command line but (of course) not depending on ini files.
      */
      String const& getIniLocation_() const
      {
        return ini_location_;
      }

      /**
      	@name Parameter handling

      	Use the methods registerStringOption_, registerDoubleOption_, registerIntOption_ and registerFlag_
      	in order to register parameters in registerOptionsAndFlags_.

      	To access the values of registered parameters in the main_ method use methods
      	getStringOption_, getDoubleOption_, getIntOption_ and getFlag_.

      	In order to format the help output the methods addEmptyLine_ and addText_ can be used.
       */
      //@{
      /**
      	 @brief Sets the valid command line options (with argument) and flags (without argument).

      	 The options '-ini' '-log' '-instance' '-debug' and the flag '--help' are automatically registered.
      */
      virtual void registerOptionsAndFlags_() = 0;

      /**
      	@brief Registers a string option. Indentation of newline is done automatically

      	@param name Name of the option in the command line and the INI file
      	@param argument Argument description text for the help output
      	@param default_value Default argument
      	@param description Description of the parameter. Indentation of newline is done automatically.
      	@param required If the user has to provide a value i.e. if the value has to differ from the default (checked in get-method)
      */
      void registerStringOption_( const String& name, const String& argument, const String& default_value, const String& description, bool required = true );

      /**
      	@brief Registers a double option. Indentation of newline is done automatically.

      	@param name Name of the option in the command line and the INI file
      	@param argument Argument description text for the help output
      	@param default_value Default argument
      	@param description Description of the parameter. Indentation of newline is done automatically.
      	@param required If the user has to provide a value i.e. if the value has to differ from the default (checked in get-method)
      */
      void registerDoubleOption_( const String& name, const String& argument, double default_value, const String& description, bool required = true );

      /**
      	@brief Registers an integer option.

      	@param name Name of the option in the command line and the INI file
      	@param argument Argument description text for the help output
      	@param default_value Default argument
      	@param description Description of the parameter. Indentation of newline is done automatically.
      	@param required If the user has to provide a value i.e. if the value has to differ from the default (checked in get-method)
      */
      void registerIntOption_( const String& name, const String& argument, Int default_value, const String& description, bool required = true );

      /// Registers a flag
      void registerFlag_( const String& name, const String& description );

      /**
      	@brief Registers an allowed subsection in the INI file.

      	Use this method to register subsections that are passed to algorithms

      	@see checkParam_
      */
      void registerSubsection_( const String& name, const String& description );

      /// Adds an empty line between registered variables in the documentation.
      void addEmptyLine_();

      /// Adds a left aligned text between registered variables in the documentation e.g. for subdividing the documentation.
      void addText_( const String& text );

      ///Returns the value of a previously registered string option
      String getStringOption_( const String& name ) const throw ( Exception::UnregisteredParameter, Exception::RequiredParameterNotGiven, Exception::WrongParameterType );

      /**
      	@brief Returns the value of a previously registered double option

      	If you want to find out if a value was really set or is a default value, use the setByUser_(String) method.
      */
      double getDoubleOption_( const String& name ) const throw ( Exception::UnregisteredParameter, Exception::RequiredParameterNotGiven, Exception::WrongParameterType );

      /**
      	@brief Returns the value of a previously registered integer option

      	If you want to find out if a value was really set or is a default value, use the setByUser_(String) method.
      */
      Int getIntOption_( const String& name ) const throw ( Exception::UnregisteredParameter, Exception::RequiredParameterNotGiven, Exception::WrongParameterType );

      ///Returns the value of a previously registered flag
      bool getFlag_( const String& name ) const throw ( Exception::UnregisteredParameter, Exception::WrongParameterType );

      ///Returns if an option was set by the user (needed to distinguish between user-set and default value)
      bool setByUser_( const String& name ) const throw ( Exception::UnregisteredParameter );

      /// Finds the entry in the parameters_ array that has the name @p name
      const ParameterInformation& findEntry_( const String& name ) const throw ( Exception::UnregisteredParameter );

      /**
      	@brief Return <em>all</em> parameters relevant to this TOPP tool.

      	Returns a Param that contains everything you can get by the getParamAs...() methods.
      */
      Param const& getParam_() const;

      /**
      	@brief Checks top-level entries of @p param according to the the information during registration

      	Only top-lvel entries and allowed subsections are checked.
      	Checking the content of the subsection is the duty of the algorithm it is passed to.

      	This method does not abort execution of the tool, but will warn the user through stderr!
      	It is called automatically in the main method

      	@param param Parameters to check
      	@param filename The source file name
      	@param location Exact location inside the source file
      */
      void checkParam_( const Param& param, const String& filename, const String& location ) const;
      //@}

      /// Prints the tool-specific command line options and appends the common options.
      void printUsage_() const;

      /// The actual "main" method.  main_() is invoked by main().
      virtual ExitCodes main_( int argc , char** argv ) = 0;

      /** @name Debug and Log output
       */
      //@{
      /// Writes a string to the log file and to std::cout
      void writeLog_( const String& text ) const;

      /// Writes a @p text to the log file and to std::cout if the debug level is at least @p min_level
      void writeDebug_( const String& text, UInt min_level ) const;

      /// Writes a String followed by a Param to the log file and to std::cout if the debug level is at least @p min_level
      void writeDebug_( const String& text, const Param& param, UInt min_level ) const;
      //@}


      /**
      	@name File IO checking methods

      	Methods used to check the validity of input and output files in main_.
      	The exceptions thrown in these methods are catched in the main method of this class.
      	They do not have to be handled in the tool itself!
      */
      //@{
      /// Checks if an input file exists, is readable and is not empty
      void inputFileReadable_( const String& filename ) const throw ( Exception::FileNotFound, Exception::FileNotReadable, Exception::FileEmpty );

      /// Checks if an output file is writable
      void outputFileWritable_( const String& filename ) const throw ( Exception::UnableToCreateFile );
      //@}

      /// Helper function that parses a range string ([a]:[b]) into to variables
      void parseRange_( const String& text, double& low, double& high ) const;

      ///Type of progress logging
      ProgressLogger::LogType log_type_;
  };



} // namespace OpenMS

#endif //OPENMS_APPLICATIONS_TOPPBASE_H
