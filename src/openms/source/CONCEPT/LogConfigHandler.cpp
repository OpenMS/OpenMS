// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <iostream>

#include <OpenMS/CONCEPT/LogConfigHandler.h>

using std::cout;
using std::cerr;
using std::endl;

namespace OpenMS
{
  String LogConfigHandler::PARAM_NAME = "log";

  LogConfigHandler::LogConfigHandler()
  {
    // add default configuration
    fatal_streams_.insert("cerr");
    error_streams_.insert("cerr");

    warn_streams_.insert("cout");
    info_streams_.insert("cout");
  }

  LogConfigHandler::LogConfigHandler(const LogConfigHandler & other)
  {
    debug_streams_ = other.debug_streams_;
    info_streams_ = other.info_streams_;
    warn_streams_ = other.warn_streams_;
    error_streams_ = other.error_streams_;
    fatal_streams_ = other.fatal_streams_;

    stream_type_map_ = other.stream_type_map_;
  }

  LogConfigHandler::~LogConfigHandler() = default;

  LogConfigHandler & LogConfigHandler::operator=(const LogConfigHandler & source) = default;

  Param LogConfigHandler::parse(const StringList & settings)
  {
    Param p;
    String suffix = " FILE";
    std::vector<std::string> commands;
    for (StringList::const_iterator iter = settings.begin(); iter != settings.end(); ++iter)
    {
      // split by " " to get all keywords
      StringList l;
      (*iter).split(' ', l, true);

      if (l.size() < 2 || l.size() > 3)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, (*iter), "Error while parsing logger config. Setting can only have 2 or 3 arguments.");
      }

      // we parse a command line here, so we append a FILE to each of the arguments
      // to indicate, that all of these streams are FILE streams
      // for cout/cerr the type parameter is ignored
      String new_command = *iter + suffix;
      commands.push_back(new_command);
    }

    p.setValue(LogConfigHandler::PARAM_NAME, commands, "List of all settings that should be applied to the current Logging Configuration");

    return p;
  }

  void LogConfigHandler::configure(const Param & param)
  {
    StringList configurations = ListUtils::toStringList<std::string>(param.getValue(LogConfigHandler::PARAM_NAME));

    for (StringList::const_iterator iter = configurations.begin(); iter != configurations.end(); ++iter)
    {
      // split by " " to get the commands
      StringList commands;
      iter->split(' ', commands, true);

      Logger::LogStream & log = getLogStreamByName_(commands[0]);

      // convenience variables
      String & command = commands[1];

      // identify action
      if (command == "add")
      {
        // convenience variables
        const String & stream_name = commands[2];

        // add the stream given by the 3rd argument to the defined log
        if (stream_name == "cout")
        {
          log.insert(cout);
        }
        else if (stream_name == "cerr")
        {
          log.insert(cerr);
        }
        else
        {
          if (commands.size() <= 3) // write error to cerr and not a LogStream (because we're just configuring it...)
          {
            std::cerr << "Error during configuring logging: the command '" << (*iter) << "' requires 4 entries but has only " << commands.size() << "\n";
            continue;
          }
          const String & stream_type = commands[3];

          // check if a stream with the same name, but different type was already registered
          if (stream_type_map_.count(stream_name) != 0)
          {
            if (stream_type_map_[stream_name] != getStreamTypeByName_(stream_type))
            {
              throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "A stream with the same name but different type was already registered.");
            }
          }

          StreamHandler::StreamType type = getStreamTypeByName_(stream_type);
          Int status = STREAM_HANDLER.registerStream(type, stream_name);

          if (!status)
          {
            // operation failed
            throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, commands[2]);
          }

          log.insert(STREAM_HANDLER.getStream(type, stream_name));
          log.setPrefix(STREAM_HANDLER.getStream(type, stream_name), "[%S] ");

          stream_type_map_[stream_name] = type;

        }

        // register the stream internally, so that the LogConfigHandler knows
        // which streams were added
        getConfigSetByName_(commands[0]).insert(stream_name);
      }
      else if (command == "remove")
      {
        // convenience variables
        const String & stream_name = commands[2];

        // add the stream given by the 3rd argument to the defined log
        if (stream_name == "cout")
        {
          log.remove(cout);
        }
        else if (stream_name == "cerr")
        {
          log.remove(cerr);
        }
        else
        {
          if (commands.size() <= 3) // write error to cerr and not a LogStream (because we're just configuring it...)
          {
            std::cerr << "Error during configuring logging: the command '" << (*iter) << "' requires 4 entries but has only " << commands.size() << "\n";
            continue;
          }
          const String & stream_type = commands[3];
          StreamHandler::StreamType type = getStreamTypeByName_(stream_type);

          // it is a file, get the ostream from the StreamHandler
          if (STREAM_HANDLER.hasStream(type, stream_name))
          {
            log.remove(STREAM_HANDLER.getStream(type, stream_name));
            STREAM_HANDLER.unregisterStream(type, stream_name);
          }
        }
        // unregister the stream internally, so that the LogConfigHandler knows
        // which streams were added
        getConfigSetByName_(commands[0]).erase(stream_name);

        // remove the type from the stream_type_map if there is no
        // stream referencing it anymore
        if (!STREAM_HANDLER.hasStream(stream_type_map_[stream_name], stream_name))
        {
          stream_type_map_.erase(stream_name);
        }
      }
      else if (command == "clear")
      {
        // remove all streams from the given log
        for (std::set<String>::iterator it = getConfigSetByName_(commands[0]).begin(); it != getConfigSetByName_(commands[0]).end(); ++it)
        {
          if (*it == "cout")
          {
            log.remove(cout);
          }
          else if (*it == "cerr")
          {
            log.remove(cerr);
          }
          else // handle the file streams
          {
            log.remove(STREAM_HANDLER.getStream(stream_type_map_[*it], *it));
            STREAM_HANDLER.unregisterStream(stream_type_map_[*it], *it);

            // remove the type from the stream_type_map if there is no
            // stream referencing it anymore
            if (!STREAM_HANDLER.hasStream(stream_type_map_[*it], *it))
            {
              stream_type_map_.erase(*it);
            }
          }
        }

        // clean the set
        getConfigSetByName_(commands[0]).clear();
      }
    }
  }

  void LogConfigHandler::setLogLevel(const String & log_level)
  {
    std::vector<String> lvls = {"DEBUG", "INFO", "WARNING", "ERROR", "FATAL_ERROR"};
    for (const auto& lvl : lvls)
    {
      if (lvl == log_level) break;
      getLogStreamByName_(lvl).removeAllStreams();
    }
  }

  Logger::LogStream & LogConfigHandler::getLogStreamByName_(const String & stream_name)
  {
    Logger::LogStream * log = &OpenMS_Log_debug; // default

    if (stream_name == "DEBUG")
    {
      log = &OpenMS_Log_debug;
    }
    else if (stream_name == "INFO")
    {
      log = &OpenMS_Log_info;
    }
    else if (stream_name == "WARNING")
    {
      log = &OpenMS_Log_warn;
    }
    else if (stream_name == "ERROR")
    {
      log = &OpenMS_Log_error;
    }
    else if (stream_name == "FATAL_ERROR")
    {
      log = &OpenMS_Log_fatal;
    }
    else
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, stream_name);
    }
    return *log;
  }

  std::set<String> & LogConfigHandler::getConfigSetByName_(const String & stream_type)
  {
    std::set<String> * s = &debug_streams_;
    if (stream_type == "DEBUG")
    {
      s = &debug_streams_;
    }
    else if (stream_type == "INFO")
    {
      s = &info_streams_;
    }
    else if (stream_type == "WARNING")
    {
      s = &warn_streams_;
    }
    else if (stream_type == "ERROR")
    {
      s = &error_streams_;
    }
    else if (stream_type == "FATAL_ERROR")
    {
      s = &fatal_streams_;
    }
    else
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, stream_type);
    }

    return *s;
  }

  std::ostream & LogConfigHandler::getStream(const String & name)
  {
    if (stream_type_map_.count(name) != 0)
    {
      return STREAM_HANDLER.getStream(stream_type_map_[name], name);
    }
    else
    {
      // there is no stream with this name
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "There is no stream with the given name.");
    }
  }

  StreamHandler::StreamType LogConfigHandler::getStreamTypeByName_(const String & stream_type)
  {
    StreamHandler::StreamType type;
    if (stream_type == "FILE")
    {
      type = StreamHandler::FILE;
    }
    else if (stream_type == "STRING")
    {
      type = StreamHandler::STRING;
    }
    else
    {
      // unsupported log type
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The log type " + stream_type + " is not supported");
    }

    return type;
  }

  void printStreamConfig_(std::ostream & os, const String & name, const std::set<String> & stream_names, const std::map<String, StreamHandler::StreamType> & stream_type_map);
  void printStreamConfig_(std::ostream & os, const String & name, const std::set<String> & stream_names, const std::map<String, StreamHandler::StreamType> & stream_type_map)
  {
    os << name << endl;
    for (std::set<String>::const_iterator it = stream_names.begin(); it != stream_names.end(); ++it)
    {
      os << "->" << "\t" << *it;
      // append stream type
      os << " (";

      switch ((stream_type_map.find(*it))->second)
      {
      case StreamHandler::STRING:
        os << "STRINGSTREAM";
        break;

      case StreamHandler::FILE:
      default:
        os << "FILE";
        break;
      }
      os <<  ")";
      os << std::endl;
    }
  }

  std::ostream & operator<<(std::ostream & os, LogConfigHandler const & lch)
  {

    printStreamConfig_(os, "OPENMS_LOG_DEBUG", lch.debug_streams_, lch.stream_type_map_);
    printStreamConfig_(os, "OPENMS_LOG_INFO", lch.info_streams_, lch.stream_type_map_);
    printStreamConfig_(os, "OPENMS_LOG_WARN", lch.warn_streams_, lch.stream_type_map_);
    printStreamConfig_(os, "OPENMS_LOG_ERROR", lch.error_streams_, lch.stream_type_map_);
    printStreamConfig_(os, "OPENMS_LOG_FATAL_ERROR", lch.fatal_streams_, lch.stream_type_map_);

    return os;
  }

  LogConfigHandler * LogConfigHandler::instance_ = nullptr;

  LogConfigHandler * LogConfigHandler::getInstance()
  {
    if (LogConfigHandler::instance_ == nullptr)
    {
      LogConfigHandler::instance_ = new LogConfigHandler();
    }
    return LogConfigHandler::instance_;
  }

} // end namespace OpenMS
