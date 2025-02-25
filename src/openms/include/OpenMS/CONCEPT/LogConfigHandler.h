// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/StreamHandler.h>

#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/DATASTRUCTURES/Param.h>

namespace OpenMS
{

  /**
    @brief The LogConfigHandler provides the functionality to configure the internal logging of OpenMS algorithms that use the
    global instances of LogStream.

    @ingroup Concept
   */
  class OPENMS_DLLAPI LogConfigHandler
  {
public:

    static String PARAM_NAME; ///< Name of the parameter in which the configuration should be stored

    /**
      @brief Translates the given list of parameter settings into a LogStream configuration

      Translates the given list of parameter settings into a LogStream configuration.
      Usually this list stems from a command line call.

      Each element in the stringlist should follow this naming convention

      &lt;LOG_NAME&gt; &lt;ACTION&gt; &lt;PARAMETER&gt;

      with
      LOG_NAME: DEBUG,INFO,WARNING,ERROR,FATAL_ERROR
      ACTION: add,remove,clear
      PARAMETER: for 'add'/'remove' it is the stream name (cout, cerr or a filename), 'clear' does not require any further parameter

      Example:
      <code>DEBUG add debug.log</code><br>

      This function will <b>not</b> apply to settings to the log handlers. Use configure() for that.

      @param setting StringList containing the configuration options
      @throw Exception::ParseError In case of an invalid configuration.
      @return Param object containing all settings, that can be applied using the LogConfigHandler::configure() method
     */
    Param parse(const StringList & setting);

    /**
      @brief Applies the given parameters (@p param) to the current configuration


      &lt;LOG_NAME&gt; &lt;ACTION&gt; &lt;PARAMETER&gt; &lt;STREAMTYPE&gt;

      LOG_NAME: DEBUG,INFO,WARNING,ERROR,FATAL_ERROR
      ACTION: add,remove,clear
      PARAMETER: for 'add'/'remove' it is the stream name (cout, cerr or a filename), 'clear' does not require any further parameter
      STREAMTYPE: FILE, STRING (for a StringStream, which you can grab by this name using getStream() )

      You cannot specify a file named "cout" or "cerr" even if you specify streamtype 'FILE' - the handler will mistake this for the
      internal streams, but you can use "./cout" to print to a file named cout.

      A classical configuration would contain a list of settings e.g.

      <code>DEBUG add debug.log FILE</code><br>
      <code>INFO remove cout FILE</code> (FILE will be ignored)<br>
      <code>INFO add string_stream1 STRING</code><br>

      @throw Exception::ElementNotFound If the LogStream (first argument) does not exist.
      @throw Exception::FileNotWritable If a file (or stream) should be opened as log file (or stream) that is not accessible.
      @throw Exception::IllegalArgument If a stream should be registered, that was already registered with a different type.
     */
    void configure(const Param & param);

    /**
      @brief Returns a reference to the registered stream with the name @p stream_name.

      @throw Exception::IllegalArgument If no stream with the name @p stream_name was registered before.

      @return Reference to the stream.
     */
    std::ostream & getStream(const String & stream_name);

    /**
      @brief Sets a minimum @p log_level by removing all streams from loggers lower than that level.
      order of log_level: "DEBUG", "INFO", "WARNING", "ERROR", "FATAL_ERROR"
     */
    void setLogLevel(const String & log_level);

    /**
      @brief Returns the instance of LogConfigHandler.
     */
    static LogConfigHandler * getInstance();

    /// Destructor
    virtual ~LogConfigHandler();

protected:

    /**
      @brief Returns the named global instance of the LogStream. (OpenMS::OpenMS_Log_debug, OpenMS::OpenMS_Log_info, OpenMS::OpenMS_Log_warn, OpenMS::OpenMS_Log_error, OpenMS::OpenMS_Log_fatal)

      @param stream_name Name of the stream. Should be DEBUG,INFO,WARNING,ERROR,FATAL_ERROR.

      @throw ElementNotFoundException if the given @p stream_name does not correspond to one of the known LogStream instances

      @return A reference to the named LogStream
     */
    Logger::LogStream & getLogStreamByName_(const String & stream_name);

    /**
      @brief Returns the correct set of registered streams for the given stream type (e.g. DEBUG, INFO, ..)

      @param stream_type String representation of the stream type (DEBUG, INFO, ..)

      @throw ElementNotFoundException if the given @p stream_type does not correspond to one of the known LogStreams
     */
    std::set<String> & getConfigSetByName_(const String & stream_type);

    /**
      @brief Translates the given @p stream_type String into a valid StreamHandler::StreamType

      @param stream_type String representation of the StreamHandler::StreamType

      @throw Exception::IllegalArgument is thrown when the passed @p stream_type does not correspond to an existing StreamHandler::StreamType

      @return The requested StreamHandler::StreamType
     */
    StreamHandler::StreamType getStreamTypeByName_(const String & stream_type);

    std::set<String> debug_streams_; ///< List of all streams that were appended to OpenMS::OpenMS_Log_debug
    std::set<String> info_streams_; ///< List of all streams that were appended to OpenMS::OpenMS_Log_info
    std::set<String> warn_streams_; ///< List of all streams that were appended to OpenMS::OpenMS_Log_warn
    std::set<String> error_streams_; ///< List of all streams that were appended to OpenMS::OpenMS_Log_error
    std::set<String> fatal_streams_; ///< List of all streams that were appended to OpenMS::OpenMS_Log_fatal

    std::map<String, StreamHandler::StreamType> stream_type_map_; ///< Maps the registered streams to a StreamHandler::StreamType

private:
    friend OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, LogConfigHandler const & lch);

    // the pointer to the single instance
    static LogConfigHandler * instance_;

    // it is a singleton so we hide all c' and d'tors

    /// Default constructor
    LogConfigHandler();

    /// Copy constructor
    LogConfigHandler(const LogConfigHandler & source);

    /// Assignment operator
    virtual LogConfigHandler & operator=(const LogConfigHandler & source);
  };

  /// Overload for the \a insertion \a operator (operator<<) to have a formatted output of the LogConfigHandler
  OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, LogConfigHandler const & lch);

} // end namespace OpenMS

